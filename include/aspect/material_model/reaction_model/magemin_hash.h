/**
 * Hash-based cache for MAGEMin equilibrium results.
 */

#ifndef _aspect_material_model_reaction_model_magemin_hash_h
#define _aspect_material_model_reaction_model_magemin_hash_h

#include <MAGEMin_cpp.h>

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace MAGEMin_hash
      {
        constexpr std::uint64_t timestep_unavailable =
          std::numeric_limits<std::uint64_t>::max();
        constexpr std::uint64_t timestep_preseed = timestep_unavailable - 1;

        struct Entry
        {
          double pressure;
          double temperature;
          std::vector<double> bulkComposition;
          stableAssemblage result;
          double meltFraction;
          std::uint64_t timestep;
          bool occupied;
          int access_count;
          std::uint64_t cached_hash;

          Entry();
        };


        class MAGEMin_Hash_table
        {
          public:
            MAGEMin_Hash_table();

            void set_oxide_names(const std::vector<std::string> &names);
            void set_context_signature(const std::string &signature)
            {
              context_signature = signature;
            }

            void set_scales(const double pressure_scale,
                            const double temperature_scale,
                            const double composition_scale)
            {
              P_scale = pressure_scale;
              T_scale = temperature_scale;
              X_scale.assign(n_oxides, composition_scale);
            }

            void set_oxide_scales(const std::vector<double> &composition_scales)
            {
              X_scale = composition_scales;
            }

            void snap_query(double &pressure_GPa,
                            double &temperature_C,
                            std::vector<double> &bulk_composition) const;

            const stableAssemblage *lookup(double pressure_GPa,
                                           double temperature_C,
                                           const std::vector<double> &bulk_composition,
                                           double &melt_fraction,
                                           bool update_statistics = true);

            void insert(double pressure_GPa,
                        double temperature_C,
                        double melt_fraction,
                        const std::vector<double> &bulk_composition,
                        const stableAssemblage &result,
                        const std::uint64_t timestep = timestep_unavailable);

            void resize_hash_table(std::size_t new_size);
            void reduce_hash_table_size(double target_load);
            void reset_access_counts();
            std::size_t getHits() const;
            std::size_t getMisses() const;
            std::size_t getMaxProbe() const;
            std::size_t getOccupancy() const;
            const std::vector<Entry> &getEntries() const
            {
              return table;
            }

            bool save_binary(const std::string &filepath) const;
            bool load_binary(const std::string &filepath);
            std::string serialize_entries() const;
            bool merge_serialized_entries(const std::string &data,
                                          const bool clear_first = false);

          private:
            void clear();
            std::uint64_t binValue(double value, double scale) const;
            std::uint64_t computeFullHash(double pressure_GPa,
                                          double temperature_C,
                                          const std::vector<double> &bulk_composition) const;
            bool binnedMatch(const Entry &entry,
                             double pressure_GPa,
                             double temperature_C,
                             const std::vector<double> &bulk_composition) const;
            void place_entry(Entry &&entry);

            static constexpr std::uint64_t FIBONACCI_MULT = 11400714819323198485ULL;

            std::size_t table_size;
            std::size_t n_oxides;
            std::vector<std::string> oxide_names;
            std::string context_signature;
            double P_scale;
            double T_scale;
            std::vector<double> X_scale;

            std::vector<Entry> table;
            std::size_t hits;
            std::size_t misses;
            std::size_t shift;
            mutable std::size_t max_probe_depth;
        };
      }
    }
  }
}

#endif
