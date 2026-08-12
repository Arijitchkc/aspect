#include <aspect/material_model/reaction_model/padawan.h>

#include <boost/property_tree/json_parser.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>


std::vector<std::string>
padawan::read_string_array(const boost::property_tree::ptree &pt,
                           const std::string &key)
{
  std::vector<std::string> values;
  for (const auto &item : pt.get_child(key))
    values.push_back(item.second.get_value<std::string>());
  return values;
}


bool padawan::load(const std::string &directory)
{
  try
    {
      boost::property_tree::ptree metadata;
      boost::property_tree::read_json(directory + "/metadata.json", metadata);

      minerals = read_string_array(metadata, "phases");
      oxides = read_string_array(metadata, "oxides");
      const std::vector<std::string> input_columns =
        read_string_array(metadata, "input_columns");

      n_total_oxides = static_cast<int>(oxides.size());
      n_inputs = 2 + n_total_oxides;
      if (input_columns.size() != static_cast<std::size_t>(n_inputs) ||
          input_columns[0] != "P_kbar" || input_columns[1] != "T_C" ||
          !std::equal(oxides.begin(), oxides.end(), input_columns.begin() + 2))
        throw std::runtime_error("metadata input_columns do not match [P_kbar, T_C, oxides]");

      module = torch::jit::load(directory + "/" + metadata.get<std::string>("artifact"));
      module.eval();
      return true;
    }
  catch (const std::exception &error)
    {
      std::cerr << "Neural-network loading failed: " << error.what() << std::endl;
      return false;
    }
}


padawan::Prediction
padawan::predict_batch(const std::vector<double> &P_kbar,
                       const std::vector<double> &T_C,
                       const std::vector<std::vector<double>> &bulk_oxides,
                       const bool validate_physical_outputs) const
{
  if (P_kbar.size() != T_C.size() || P_kbar.size() != bulk_oxides.size())
    throw std::invalid_argument("padawan::predict_batch received arrays of different sizes");

  Prediction prediction;
  const long batch_size = static_cast<long>(P_kbar.size());
  if (batch_size == 0)
    return prediction;

  std::vector<float> values(static_cast<std::size_t>(batch_size) * n_inputs, 0.0f);
  std::vector<int> invalid(static_cast<std::size_t>(batch_size), 0);
  for (long batch = 0; batch < batch_size; ++batch)
    {
      const std::size_t offset = static_cast<std::size_t>(batch) * n_inputs;
      invalid[batch] = !std::isfinite(P_kbar[batch]) || !std::isfinite(T_C[batch]) ||
                       bulk_oxides[batch].size() != static_cast<std::size_t>(n_total_oxides);
      values[offset] = std::isfinite(P_kbar[batch])
                       ? static_cast<float>(P_kbar[batch]) : 0.0f;
      values[offset + 1] = std::isfinite(T_C[batch])
                           ? static_cast<float>(T_C[batch]) : 0.0f;

      for (int oxide = 0; oxide < n_total_oxides; ++oxide)
        {
          float value = oxide < static_cast<int>(bulk_oxides[batch].size())
                        ? static_cast<float>(bulk_oxides[batch][oxide]) : 0.0f;
          if (!std::isfinite(value) || value > 2.0f)
            {
              invalid[batch] = 1;
              value = 0.0f;
            }
          values[offset + 2 + oxide] = std::max(value, 0.0f);
        }
    }

  const torch::Tensor inputs = torch::from_blob(
                                 values.data(), {batch_size, n_inputs}, torch::kFloat32).clone();
  torch::NoGradGuard no_grad;
  const torch::jit::IValue output = module.forward({inputs});
  if (!output.isTuple())
    throw std::runtime_error("The neural network did not return a tuple");

  const auto &items = output.toTuple()->elements();
  if (items.size() != 6)
    throw std::runtime_error("The neural network must return six tensors");

  const torch::Tensor proportions =
    items.at(0).toTensor().to(torch::kCPU, torch::kFloat32).contiguous();
  const torch::Tensor mineral_oxides =
    items.at(1).toTensor().to(torch::kCPU, torch::kFloat32).contiguous();
  const torch::Tensor physical =
    items.at(5).toTensor().to(torch::kCPU, torch::kFloat32).contiguous();

  if (proportions.dim() != 2 || mineral_oxides.dim() != 3 || physical.dim() != 2 ||
      proportions.size(0) != batch_size || mineral_oxides.size(0) != batch_size ||
      physical.size(0) != batch_size ||
      proportions.size(1) != static_cast<long>(minerals.size()) ||
      mineral_oxides.size(1) != static_cast<long>(minerals.size()) ||
      mineral_oxides.size(2) != n_total_oxides || physical.size(1) != 4)
    throw std::runtime_error("The neural-network output tensor shapes do not match metadata");

  const auto modal = proportions.accessor<float, 2>();
  const auto chemistry = mineral_oxides.accessor<float, 3>();
  const auto properties = physical.accessor<float, 2>();
  const long n_minerals = proportions.size(1);

  prediction.proportions.assign(batch_size, std::vector<float>(n_minerals));
  prediction.mineral_oxides.assign(
    batch_size,
    std::vector<std::vector<float>>(n_minerals,
                                     std::vector<float>(n_total_oxides)));
  prediction.physical.assign(batch_size, std::vector<float>(4));
  prediction.is_ood = invalid;

  for (long batch = 0; batch < batch_size; ++batch)
    {
      float modal_sum = 0.0f;
      for (long mineral = 0; mineral < n_minerals; ++mineral)
        {
          const float raw_proportion = modal[batch][mineral];
          if (!std::isfinite(raw_proportion))
            {
              prediction.is_ood[batch] = 1;
              continue;
            }

          const float proportion = raw_proportion < 1.0e-5f
                                   ? 0.0f : std::clamp(raw_proportion, 0.0f, 1.0f);
          prediction.proportions[batch][mineral] = proportion;
          modal_sum += proportion;

          for (int oxide = 0; oxide < n_total_oxides; ++oxide)
            {
              const float value = chemistry[batch][mineral][oxide];
              if (!std::isfinite(value))
                prediction.is_ood[batch] = 1;
              prediction.mineral_oxides[batch][mineral][oxide] =
                proportion == 0.0f || !std::isfinite(value)
                ? 0.0f : std::clamp(value, 0.0f, 1.0f);
            }
        }

      if (modal_sum < 0.9f || modal_sum > 1.1f)
        prediction.is_ood[batch] = 1;

      for (long property = 0; property < 4; ++property)
        {
          const float value = properties[batch][property];
          prediction.physical[batch][property] = value;
          if (validate_physical_outputs &&
              (!std::isfinite(value) || value <= 0.0f))
            prediction.is_ood[batch] = 1;
        }
    }

  return prediction;
}
