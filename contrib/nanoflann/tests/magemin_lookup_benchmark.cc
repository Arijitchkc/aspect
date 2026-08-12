#include <MAGEMin_cpp.h>
#include <nanoflann.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <vector>


namespace
{
  constexpr unsigned int dimension = 13;
  constexpr unsigned int n_oxides = 11;

  const std::array<double, dimension> tolerances =
  {{
      0.05, 5.0,
      0.001, 0.0005, 0.0005, 0.001, 0.001, 0.0001,
      0.0002, 0.0001, 0.00005, 0.0001, 0.0001
    }};


  struct LookupPoint
  {
    std::array<double, dimension> coordinates;
    std::uint64_t assemblage_mask = 0;
    stableAssemblage result;
  };


  struct PointCloud
  {
    std::vector<LookupPoint> points;

    std::size_t kdtree_get_point_count() const
    {
      return points.size();
    }

    double kdtree_get_pt(const std::size_t point_index,
                         const std::size_t coordinate) const
    {
      return points[point_index].coordinates[coordinate] / tolerances[coordinate];
    }

    template <class BoundingBox>
    bool kdtree_get_bbox(BoundingBox &) const
    {
      return false;
    }
  };


  using Distance = nanoflann::L2_Simple_Adaptor<double, PointCloud>;
  using IncrementalIndex = nanoflann::KDTreeSingleIndexIncrementalAdaptor<
    Distance, PointCloud, dimension, std::size_t>;
  using StaticIndex = nanoflann::KDTreeSingleIndexAdaptor<
    Distance, PointCloud, dimension, std::size_t>;


  std::uint64_t
  assemblage_mask(const stableAssemblage &result)
  {
    std::uint64_t mask = 0;
    for (unsigned int i = 0;
         i < result.mineral_proportions.size() && i < 64;
         ++i)
      if (result.mineral_proportions[i] > 1e-4)
        mask |= (std::uint64_t(1) << i);
    return mask;
  }


  double
  melt_fraction(const stableAssemblage &result)
  {
    const auto mineral = std::find(result.mineral_names.begin(),
                                   result.mineral_names.end(),
                                   "liq");
    if (mineral == result.mineral_names.end())
      return 0.0;

    const std::size_t index = std::distance(result.mineral_names.begin(), mineral);
    return index < result.mineral_proportions.size()
           ? result.mineral_proportions[index] : 0.0;
  }


  std::array<double, dimension>
  make_coordinates(const double pressure_GPa,
                   const double temperature_C,
                   const std::vector<double> &composition)
  {
    std::array<double, dimension> coordinates = {{0.0}};
    coordinates[0] = pressure_GPa;
    coordinates[1] = temperature_C;
    for (unsigned int i = 0; i < n_oxides; ++i)
      coordinates[i+2] = composition[i];
    return coordinates;
  }


  std::vector<double>
  normalized_composition(const std::array<double, dimension> &coordinates)
  {
    std::vector<double> composition(n_oxides);
    for (unsigned int i = 0; i < n_oxides; ++i)
      composition[i] = std::max(0.0, coordinates[i+2]);

    const double sum = std::accumulate(composition.begin(), composition.end(), 0.0);
    for (double &value : composition)
      value /= sum;
    return composition;
  }


  std::vector<std::size_t>
  box_query(const IncrementalIndex &index,
            const std::array<double, dimension> &coordinates)
  {
    IncrementalIndex::BoundingBox box;
    for (unsigned int d = 0; d < dimension; ++d)
      {
        const double scaled_coordinate = coordinates[d] / tolerances[d];
        box[d].low = scaled_coordinate - 1.0;
        box[d].high = scaled_coordinate + 1.0;
      }

    std::vector<std::size_t> candidates;
    nanoflann::BoxResultSet<std::size_t> result(candidates);
    const std::size_t n_candidates = index.findWithinBox(result, box);
    (void) n_candidates;
    std::sort(candidates.begin(), candidates.end());
    return candidates;
  }


  std::vector<std::size_t>
  brute_force_query(const PointCloud &cloud,
                    const std::array<double, dimension> &coordinates)
  {
    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < cloud.points.size(); ++i)
      {
        bool inside = true;
        for (unsigned int d = 0; d < dimension; ++d)
          if (std::abs(cloud.points[i].coordinates[d] - coordinates[d]) > tolerances[d])
            {
              inside = false;
              break;
            }
        if (inside)
          candidates.push_back(i);
      }
    return candidates;
  }


  std::size_t
  nearest_candidate(const PointCloud &cloud,
                    const std::vector<std::size_t> &candidates,
                    const std::array<double, dimension> &coordinates)
  {
    double nearest_distance = std::numeric_limits<double>::max();
    std::size_t nearest = std::numeric_limits<std::size_t>::max();
    for (const std::size_t index : candidates)
      {
        double distance = 0.0;
        for (unsigned int d = 0; d < dimension; ++d)
          {
            const double difference =
              (cloud.points[index].coordinates[d] - coordinates[d]) / tolerances[d];
            distance += difference * difference;
          }
        if (distance < nearest_distance)
          {
            nearest_distance = distance;
            nearest = index;
          }
      }
    return nearest;
  }


  unsigned int
  neighbourhood_status(const PointCloud &cloud,
                       const std::vector<std::size_t> &candidates)
  {
    if (candidates.size() < 3)
      return 1;

    const std::uint64_t mask = cloud.points[candidates[0]].assemblage_mask;
    double minimum_melt = std::numeric_limits<double>::max();
    double maximum_melt = -std::numeric_limits<double>::max();

    for (const std::size_t index : candidates)
      {
        if (cloud.points[index].assemblage_mask != mask)
          return 2;
        const double melt = melt_fraction(cloud.points[index].result);
        minimum_melt = std::min(minimum_melt, melt);
        maximum_melt = std::max(maximum_melt, melt);
      }

    if (maximum_melt - minimum_melt > 0.001)
      return 3;

    const unsigned int n_minerals = cloud.points[candidates[0]].result.mineral_proportions.size();
    for (unsigned int mineral = 0; mineral < n_minerals; ++mineral)
      {
        double minimum = std::numeric_limits<double>::max();
        double maximum = -std::numeric_limits<double>::max();
        for (const std::size_t index : candidates)
          {
            const double value = cloud.points[index].result.mineral_proportions[mineral];
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
          }
        if (maximum - minimum > 0.10)
          return 4;
      }

    for (unsigned int mineral = 0; mineral < n_minerals; ++mineral)
      {
        if ((mask & (std::uint64_t(1) << mineral)) == 0)
          continue;

        for (unsigned int oxide = 0; oxide < n_oxides; ++oxide)
          {
            double minimum = std::numeric_limits<double>::max();
            double maximum = -std::numeric_limits<double>::max();
            for (const std::size_t index : candidates)
              {
                const double value =
                  cloud.points[index].result.mineral_oxides[mineral][oxide];
                minimum = std::min(minimum, value);
                maximum = std::max(maximum, value);
              }
            if (maximum - minimum > 0.01)
              return 5;
          }
      }

    return 0;
  }
}


int main()
{
  const std::vector<double> base_composition =
  {
    0.4486, 0.03512, 0.0307, 0.39524, 0.08202, 0.000183,
    0.003, 0.00155, 0.000298, 0.00321, 0.0
  };

  std::mt19937 generator(1984);
  std::uniform_real_distribution<double> pressure_distribution(0.5, 5.0);
  std::uniform_real_distribution<double> temperature_distribution(1050.0, 1650.0);
  std::uniform_real_distribution<double> unit_distribution(-1.0, 1.0);

  constexpr unsigned int n_clusters = 12;
  constexpr unsigned int training_points_per_cluster = 10;
  constexpr unsigned int query_points_per_cluster = 5;

  std::vector<std::array<double, dimension>> training_coordinates;
  std::vector<std::array<double, dimension>> query_coordinates;

  for (unsigned int cluster = 0; cluster < n_clusters; ++cluster)
    {
      std::array<double, dimension> center =
        make_coordinates(pressure_distribution(generator),
                         temperature_distribution(generator),
                         base_composition);
      for (unsigned int d = 2; d < dimension; ++d)
        center[d] = std::max(0.0, center[d] + 0.25 * tolerances[d] * unit_distribution(generator));

      for (unsigned int point = 0; point < training_points_per_cluster; ++point)
        {
          std::array<double, dimension> coordinates = center;
          const double perturbation = cluster % 2 == 0 ? 0.02 : 0.60;
          for (unsigned int d = 0; d < dimension; ++d)
            coordinates[d] += perturbation * tolerances[d] * unit_distribution(generator);
          const std::vector<double> composition = normalized_composition(coordinates);
          training_coordinates.push_back(make_coordinates(coordinates[0], coordinates[1], composition));
        }

      for (unsigned int point = 0; point < query_points_per_cluster; ++point)
        {
          std::array<double, dimension> coordinates = center;
          const double perturbation = cluster % 2 == 0 ? 0.01 : 0.35;
          for (unsigned int d = 0; d < dimension; ++d)
            coordinates[d] += perturbation * tolerances[d] * unit_distribution(generator);
          const std::vector<double> composition = normalized_composition(coordinates);
          query_coordinates.push_back(make_coordinates(coordinates[0], coordinates[1], composition));
        }
    }

  std::vector<std::array<double, dimension>> all_coordinates = training_coordinates;
  all_coordinates.insert(all_coordinates.end(), query_coordinates.begin(), query_coordinates.end());

  std::vector<double> temperatures;
  std::vector<double> pressures;
  std::vector<std::vector<double>> compositions;
  for (const auto &coordinates : all_coordinates)
    {
      temperatures.push_back(coordinates[1] + 273.15);
      pressures.push_back(coordinates[0] * 10.0);
      compositions.push_back(normalized_composition(coordinates));
    }

  MAGEMin_wrapper wrapper;
  std::vector<stableAssemblage> results;
  const auto magemin_start = std::chrono::steady_clock::now();
  wrapper.executeMAGEMin(0, nullptr, temperatures, pressures, n_oxides,
                         const_cast<char *>("ig"), compositions, results);
  const auto magemin_end = std::chrono::steady_clock::now();

  PointCloud cloud;
  cloud.points.reserve(training_coordinates.size());
  for (unsigned int i = 0; i < training_coordinates.size(); ++i)
    cloud.points.push_back({training_coordinates[i], assemblage_mask(results[i]), results[i]});

  IncrementalIndex index(dimension, cloud);
  const auto build_start = std::chrono::steady_clock::now();
  index.addPoints(0, cloud.points.size()-1);
  const auto build_end = std::chrono::steady_clock::now();

  unsigned int exact_candidate_sets = 0;
  unsigned int accepted_queries = 0;
  std::array<unsigned int, 6> neighbourhood_counts = {{0}};
  double density_error_sum = 0.0;
  double density_error_max = 0.0;
  double cp_error_sum = 0.0;
  double cp_error_max = 0.0;
  double enthalpy_error_sum = 0.0;
  double enthalpy_error_max = 0.0;
  double melt_error_sum = 0.0;
  double melt_error_max = 0.0;

  const auto query_start = std::chrono::steady_clock::now();
  for (unsigned int i = 0; i < query_coordinates.size(); ++i)
    {
      const std::vector<std::size_t> candidates = box_query(index, query_coordinates[i]);
      const std::vector<std::size_t> brute_force = brute_force_query(cloud, query_coordinates[i]);
      if (candidates == brute_force)
        ++exact_candidate_sets;

      const unsigned int status = neighbourhood_status(cloud, candidates);
      ++neighbourhood_counts[status];
      if (status != 0)
        continue;

      const std::size_t nearest = nearest_candidate(cloud, candidates, query_coordinates[i]);
      const stableAssemblage &approximation = cloud.points[nearest].result;
      const stableAssemblage &truth = results[training_coordinates.size()+i];

      const double density_error = std::abs(approximation.bulk_density - truth.bulk_density);
      const double cp_error = std::abs(approximation.bulk_cp - truth.bulk_cp);
      const double enthalpy_error = std::abs(approximation.bulk_enthalpy - truth.bulk_enthalpy);
      const double melt_error = std::abs(melt_fraction(approximation) - melt_fraction(truth));

      density_error_sum += density_error;
      density_error_max = std::max(density_error_max, density_error);
      cp_error_sum += cp_error;
      cp_error_max = std::max(cp_error_max, cp_error);
      enthalpy_error_sum += enthalpy_error;
      enthalpy_error_max = std::max(enthalpy_error_max, enthalpy_error);
      melt_error_sum += melt_error;
      melt_error_max = std::max(melt_error_max, melt_error);
      ++accepted_queries;
    }
  const auto query_end = std::chrono::steady_clock::now();

  const std::vector<std::size_t> oxide_guard_candidates =
    box_query(index, query_coordinates.front());
  bool oxide_guard_rejected = false;
  if (!oxide_guard_candidates.empty())
    {
      stableAssemblage &changed = cloud.points[oxide_guard_candidates.front()].result;
      for (unsigned int mineral = 0;
           mineral < changed.mineral_proportions.size();
           ++mineral)
        if (changed.mineral_proportions[mineral] > 1e-4)
          {
            changed.mineral_oxides[mineral][0] += 0.02;
            oxide_guard_rejected =
              neighbourhood_status(cloud, oxide_guard_candidates) == 5;
            break;
          }
    }

  const auto milliseconds = [](const auto begin, const auto end)
  {
    return std::chrono::duration<double, std::milli>(end-begin).count();
  };

  std::cout << std::fixed << std::setprecision(6)
            << "MAGEMin points: " << results.size() << '\n'
            << "MAGEMin batch time [ms]: " << milliseconds(magemin_start, magemin_end) << '\n'
            << "nanoflann build time [ms]: " << milliseconds(build_start, build_end) << '\n'
            << "nanoflann query time [ms]: " << milliseconds(query_start, query_end) << '\n'
            << "box-query agreement with brute force: " << exact_candidate_sets
            << '/' << query_coordinates.size() << '\n'
            << "thermodynamically accepted queries: " << accepted_queries
            << '/' << query_coordinates.size() << '\n'
            << "synthetic oxide-gradient rejection: "
            << (oxide_guard_rejected ? "pass" : "fail") << '\n'
            << "safe/too-few/phase/melt/mineral/oxide: "
            << neighbourhood_counts[0] << ' '
            << neighbourhood_counts[1] << ' '
            << neighbourhood_counts[2] << ' '
            << neighbourhood_counts[3] << ' '
            << neighbourhood_counts[4] << ' '
            << neighbourhood_counts[5] << '\n';

  if (accepted_queries > 0)
    std::cout << "density error mean/max: "
              << density_error_sum/accepted_queries << ' ' << density_error_max << '\n'
              << "Cp error mean/max: "
              << cp_error_sum/accepted_queries << ' ' << cp_error_max << '\n'
              << "enthalpy error mean/max: "
              << enthalpy_error_sum/accepted_queries << ' ' << enthalpy_error_max << '\n'
              << "melt-fraction error mean/max: "
              << melt_error_sum/accepted_queries << ' ' << melt_error_max << '\n';

  return exact_candidate_sets == query_coordinates.size() && oxide_guard_rejected
         ? 0 : 1;
}
