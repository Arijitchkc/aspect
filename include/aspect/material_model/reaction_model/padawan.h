#ifndef PADAWAN_H
#define PADAWAN_H
#pragma push_macro("Assert")
#pragma push_macro("AssertThrow")
#undef Assert
#undef AssertThrow
#include <torch/script.h>
#pragma pop_macro("AssertThrow")
#pragma pop_macro("Assert")

#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

class padawan
{
  public:
    struct Prediction
    {
      std::vector<std::vector<float>> proportions;
      std::vector<std::vector<std::vector<float>>> mineral_oxides;
      // [rho_solid, rho_liquid, Cp_solid, Cp_liquid]
      std::vector<std::vector<float>> physical;
      std::vector<int> is_ood;
    };

    bool load(const std::string &file_path);
    Prediction predict_batch(const std::vector<double> &P_kbar,
                             const std::vector<double> &T_C,
                             const std::vector<std::vector<double>> &bulk_oxides,
                             const bool validate_physical_outputs = true) const;

    const std::vector<std::string> &mineral_names() const
    {
      return minerals;
    }
    const std::vector<std::string> &oxide_names()   const
    {
      return oxides;
    }

  private:
    static std::vector<std::string> read_string_array(const boost::property_tree::ptree &pt,
                                                      const std::string &key);

    mutable torch::jit::script::Module module;
    int n_inputs = 0;
    int n_total_oxides = 0;
    std::vector<std::string> oxides;
    std::vector<std::string> minerals;
};
#endif
