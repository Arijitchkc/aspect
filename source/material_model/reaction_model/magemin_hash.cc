/**
 * Hash-based cache for MAGEMin equilibrium results.
 */

#include <aspect/material_model/reaction_model/magemin_hash.h>
#include <aspect/utilities.h>

#include <deal.II/base/exceptions.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace MAGEMin_hash
      {
        Entry::Entry()
          : pressure(0.0),
            temperature(0.0),
            meltFraction(0.0),
            timestep(timestep_unavailable),
            occupied(false),
            access_count(0),
            cached_hash(0)
        {}


        MAGEMin_Hash_table::MAGEMin_Hash_table()
          : table_size(0),
            n_oxides(0),
            P_scale(50.0),
            T_scale(1.0),
            table(),
            hits(0),
            misses(0),
            shift(64),
            max_probe_depth(0)
        {
        }

        void MAGEMin_Hash_table::set_oxide_names(const std::vector<std::string> &names)
        {
          AssertThrow(!names.empty(),
                      ExcMessage("A MAGEMin cache needs at least one oxide."));
          n_oxides = names.size();
          X_scale.assign(n_oxides, 200.0);
          clear();
          oxide_names = names;
        }

        void MAGEMin_Hash_table::snap_query(double &pressure_GPa,
                                            double &temperature_C,
                                            std::vector<double> &bulk_composition) const
        {
          pressure_GPa = std::round(std::max(0.0, pressure_GPa) * P_scale) / P_scale;
          temperature_C = std::round(std::max(0.0, temperature_C) * T_scale) / T_scale;
          AssertDimension(bulk_composition.size(), n_oxides);
          for (std::size_t i = 0; i < n_oxides; ++i)
            bulk_composition[i] = std::round(std::max(0.0, bulk_composition[i]) * X_scale[i]) / X_scale[i];
        }

        void MAGEMin_Hash_table::resize_hash_table(size_t new_size)
        {
          AssertThrow(n_oxides > 0,
                      ExcMessage("Set the MAGEMin cache oxide count before resizing it."));
          AssertThrow( (new_size>0) && (new_size & (new_size - 1)) == 0,
                       ExcMessage("Hash table size MUST be a power of 2"));
          table_size = new_size;
          table.clear();
          table.resize(new_size);
          hits = 0;
          misses = 0;
          shift = 64 - static_cast<int>(std::log2(table_size));
        }

        uint64_t MAGEMin_Hash_table::binValue(double v, double scale) const
        {
          return static_cast<uint64_t>(std::round(std::max(0.0, v) * scale));
        }

        bool MAGEMin_Hash_table::binnedMatch(const Entry &e, double P, double T,
                                             const std::vector<double> &bulk) const
        {
          AssertDimension(bulk.size(), n_oxides);
          AssertDimension(e.bulkComposition.size(), n_oxides);
          if (binValue(e.pressure, P_scale)    != binValue(P, P_scale))    return false;
          if (binValue(e.temperature, T_scale) != binValue(T, T_scale))    return false;
          for (size_t i = 0; i < n_oxides; ++i)
            if (binValue(e.bulkComposition[i], X_scale[i]) != binValue(bulk[i], X_scale[i]))
              return false;
          return true;
        }

        uint64_t MAGEMin_Hash_table::computeFullHash(double P, double T,
                                                     const std::vector<double> &bulk) const
        {
          AssertDimension(bulk.size(), n_oxides);
          uint64_t p_bin = binValue(P, P_scale);
          uint64_t t_bin = binValue(T, T_scale);
          uint64_t h = (t_bin << 32) | (p_bin & 0xFFFFFFFF);

          // Include all oxides, particularly O and H2O, so chemically distinct
          // states do not share the same initial probe chain.
          for (std::size_t i = 0; i < n_oxides; ++i)
            {
              const uint64_t ox_bin = binValue(bulk[i], X_scale[i]);
              h ^= ox_bin;
              h *= 0xc6a4a7935bd1e995ULL;
              h ^= (h >> 47);
            }
          return h * FIBONACCI_MULT;   // pre-shift
        }








        // Lookup entry in cache
        const stableAssemblage *MAGEMin_Hash_table::lookup(double P_GPa,
                                                           double T_celsius,
                                                           const std::vector<double> &bulkComp,
                                                           double &out_meltFraction,
                                                           const bool update_statistics)
        {
          AssertThrow(table_size > 0, ExcMessage("MAGEMin cache has not been initialized."));

          uint64_t full_hash = computeFullHash(P_GPa, T_celsius, bulkComp);
          size_t idx = static_cast<size_t>(full_hash >> shift);
          size_t current_depth = 0;
          // Linear probing for collision resolution
          for (size_t probe = 0; probe < table_size; ++probe)
            {
              current_depth = probe;
              size_t probe_idx = (idx + probe) & (table_size - 1);

              // If we hit an unoccupied slot, the entry definitely isn't here
              if (!table[probe_idx].occupied)
                {
                  if (update_statistics)
                    {
                      misses++;
                      if (current_depth > max_probe_depth)
                        max_probe_depth = current_depth;
                    }
                  return nullptr;
                }


              if (table[probe_idx].cached_hash != full_hash)
                {
                  continue;
                }

              if (binnedMatch(table[probe_idx], P_GPa,T_celsius, bulkComp))
                {
                  if (update_statistics)
                    {
                      hits++;
                      table[probe_idx].access_count++;
                      if (current_depth > max_probe_depth)
                        max_probe_depth = current_depth;
                    }
                  out_meltFraction = table[probe_idx].meltFraction;
                  return &table[probe_idx].result;
                }
            }


          // values is not found after full probe
          if (update_statistics)
            {
              misses++;
              if (current_depth > max_probe_depth) max_probe_depth = current_depth;
            }
          return nullptr;
        }

        void MAGEMin_Hash_table::insert(double P_GPa,
                                        double T_celsius,
                                        double melt_MAGEMin,
                                        const std::vector<double> &bulkComp,
                                        const stableAssemblage &mage_data_structure,
                                        const std::uint64_t solve_timestep)
        {
          AssertDimension(bulkComp.size(), n_oxides);
          // Composition is part of the key, preventing pressure-temperature
          // collisions from concentrating entries in a few probe chains.
          uint64_t full_hash = computeFullHash(P_GPa, T_celsius, bulkComp);
          size_t idx = static_cast<size_t>(full_hash >> shift);


          // Linear probing to find slot
          for (size_t probe = 0; probe < table_size; ++probe)
            {
              size_t probe_idx = (idx + probe) & (table_size - 1);

              // --- Diagnostic Tracking ---
              if (probe > max_probe_depth) max_probe_depth = probe;

              // Case 1: Empty slot - insert here
              if (!table[probe_idx].occupied)
                {
                  table[probe_idx].pressure = P_GPa;
                  table[probe_idx].temperature = T_celsius;

                  table[probe_idx].bulkComposition = bulkComp;
                  table[probe_idx].result = mage_data_structure;
                  table[probe_idx].meltFraction = melt_MAGEMin;
                  table[probe_idx].timestep = solve_timestep;
                  table[probe_idx].cached_hash = full_hash;
                  table[probe_idx].occupied = true;
                  table[probe_idx].access_count = 0;


                  return;
                }

              if (binnedMatch(table[probe_idx], P_GPa, T_celsius, bulkComp))
                {

                  // Update existing entry with new result
                  table[probe_idx].bulkComposition = bulkComp;
                  table[probe_idx].result = mage_data_structure;
                  table[probe_idx].meltFraction = melt_MAGEMin;
                  table[probe_idx].timestep = solve_timestep;
                  table[probe_idx].cached_hash = full_hash;
                  table[probe_idx].occupied = true;
                  return;
                }
            }


          static int full_warn_count = 0;
          if (full_warn_count < 10)
            {
              std::cerr << "WARNING: Hash table neighborhood full at " << idx << "\n";
              full_warn_count++;
            }
        }




        void MAGEMin_Hash_table::reduce_hash_table_size(double target_load)
        {
          const size_t target_count = static_cast<size_t>(table_size * target_load);

          // Pass 1: Welford mean/std of access_count across occupied entries.
          double access_mean = 0.0;
          double access_M2 = 0.0;
          size_t access_n = 0;
          for (const auto &entry : table)
            {
              if (!entry.occupied) continue;
              access_n++;
              double delta = entry.access_count - access_mean;
              access_mean += delta / access_n;
              double delta2 = entry.access_count - access_mean;
              access_M2 += delta * delta2;
            }
          if (access_n == 0) return;

          const double access_std = (access_n > 1) ? std::sqrt(access_M2 / access_n) : 0.0;
          const double threshold  = std::max(access_mean - 0.5 * access_std, 0.5);

          // Pass 2: extract survivors by MOVE (entries hold a full stableAssemblage —
          // copying them here would be the most expensive line in the function).
          std::vector<Entry> survivors;
          survivors.reserve(access_n);
          for (auto &entry : table)
            {
              if (!entry.occupied) continue;
              if (entry.access_count >= threshold)
                survivors.push_back(std::move(entry));
              entry.occupied = false;   // moved-from or evicted either way
            }

          // Pass 3: if the threshold wasn't aggressive enough, drop the least-used
          // survivors until we are at target. nth_element: O(N), no full sort needed.
          if (survivors.size() > target_count)
            {
              std::nth_element(survivors.begin(),
                               survivors.begin() + target_count,
                               survivors.end(),
                               [](const Entry &a, const Entry &b)
              {
                return a.access_count > b.access_count;
              });
              survivors.resize(target_count);
            }

          // Pass 4: rebuild. Every survivor is re-placed from its home slot, so every
          // probe chain is contiguous again. This is the fix for the unreachable-entry
          // bug: previously, punching holes in place broke linear probing and produced
          // false misses + duplicate inserts.
          for (auto &e : survivors)
            place_entry(std::move(e));
        }

        void MAGEMin_Hash_table::place_entry(Entry &&e)
        {
          // cached_hash stores the full pre-shift 64-bit hash, so no recompute needed,
          // and this stays correct even if shift changed via resize_hash_table.
          const size_t idx = static_cast<size_t>(e.cached_hash >> shift);
          for (size_t probe = 0; probe < table_size; ++probe)
            {
              const size_t pi = (idx + probe) & (table_size - 1);
              if (!table[pi].occupied)
                {
                  table[pi] = std::move(e);
                  table[pi].occupied = true;
                  return;
                }
            }
          // Unreachable when target_load < 1.0
        }

        void MAGEMin_Hash_table::reset_access_counts()
        {
          hits = 0;
          misses = 0;
          max_probe_depth = 0;
          for (auto &entry : table)
            {
              entry.access_count = 0;
            }
        }

        void MAGEMin_Hash_table::clear()
        {
          for (auto &entry : table)
            {
              entry.occupied = false;
            }
          hits = 0;
          misses = 0;
          max_probe_depth = 0;
        }

        size_t MAGEMin_Hash_table::getHits() const
        {
          return hits;
        }

        size_t MAGEMin_Hash_table::getMisses() const
        {
          return misses;
        }
        size_t MAGEMin_Hash_table::getMaxProbe() const
        {
          return max_probe_depth;
        }

        size_t MAGEMin_Hash_table::getOccupancy() const
        {
          size_t count = 0;
          for (const auto &entry : table)
            {
              if (entry.occupied)
                {
                  count++;
                }
            }
          return count;
        }


        namespace
        {
          // One writer/reader pair carries every field: save_binary and
          // load_binary list the same fields in the same order, so the two
          // functions stay symmetric by inspection.
          struct Writer
          {
            std::string buf;
            void raw(const void *src, std::size_t n)
            {
              buf.append(static_cast<const char *>(src), n);
            }
            void put(double v)
            {
              raw(&v, sizeof v);
            }
            void put(std::int32_t v)
            {
              raw(&v, sizeof v);
            }
            void put(std::uint32_t v)
            {
              raw(&v, sizeof v);
            }
            void put(std::uint64_t v)
            {
              raw(&v, sizeof v);
            }
            void put(const std::string &s)
            {
              put(static_cast<std::uint32_t>(s.size()));
              raw(s.data(), s.size());
            }
            void put(const std::vector<double> &v)
            {
              put(static_cast<std::uint32_t>(v.size()));
              raw(v.data(), v.size() * sizeof(double));
            }
            void put(const std::vector<std::string> &v)
            {
              put(static_cast<std::uint32_t>(v.size()));
              for (const auto &s : v) put(s);
            }
          };

          struct Reader
          {
            const char *p;
            const char *end;
            bool ok = true;
            void raw(void *dst, std::size_t n)
            {
              if (!ok || static_cast<std::size_t>(end - p) < n)
                {
                  ok = false;
                  return;
                }
              std::memcpy(dst, p, n);
              p += n;
            }
            void get(double &v)
            {
              raw(&v, sizeof v);
            }
            void get(std::int32_t &v)
            {
              raw(&v, sizeof v);
            }
            void get(std::uint32_t &v)
            {
              raw(&v, sizeof v);
            }
            void get(std::uint64_t &v)
            {
              raw(&v, sizeof v);
            }
            void get(std::string &s)
            {
              std::uint32_t n = 0;
              get(n);
              if (!ok || static_cast<std::size_t>(end - p) < n)
                {
                  ok = false;
                  return;
                }
              s.assign(p, n);
              p += n;
            }
            void get(std::vector<double> &v)
            {
              std::uint32_t n = 0;
              get(n);
              if (!ok || static_cast<std::size_t>(end - p) < n * sizeof(double))
                {
                  ok = false;
                  return;
                }
              v.resize(n);
              raw(v.data(), n * sizeof(double));
            }
            void get(std::vector<std::string> &v)
            {
              std::uint32_t n = 0;
              get(n);
              if (!ok) return;
              v.resize(n);
              for (auto &s : v) get(s);
            }
          };

          constexpr std::uint32_t MAGE_MAGIC   = 0x4D414745u;  // 'M','A','G','E'
          // Version 9 also records when each exact state was solved.
          constexpr std::uint32_t MAGE_VERSION = 9u;

          struct DecodedEntry
          {
            double pressure;
            double temperature;
            double melt_fraction;
            std::vector<double> bulk_composition;
            std::uint64_t timestep;
            stableAssemblage result;
          };

          std::uint64_t checksum64(const std::string &data)
          {
            std::uint64_t hash = 1469598103934665603ULL;
            for (const unsigned char byte : data)
              {
                hash ^= byte;
                hash *= 1099511628211ULL;
              }
            return hash;
          }
        }


        std::string MAGEMin_Hash_table::serialize_entries() const
        {
          std::uint32_t n_entries = 0;
          for (const auto &e : table)
            if (e.occupied) ++n_entries;

          Writer header;
          header.put(MAGE_MAGIC);
          header.put(MAGE_VERSION);
          header.put(static_cast<std::uint32_t>(n_oxides));
          header.put(oxide_names);
          header.put(P_scale);
          header.put(T_scale);
          header.put(X_scale);
          header.put(context_signature);
          header.put(n_entries);

          Writer w;
          for (const auto &e : table)
            {
              if (!e.occupied) continue;
              AssertDimension(e.bulkComposition.size(), n_oxides);

              w.put(e.pressure);
              w.put(e.temperature);
              w.put(e.meltFraction);
              w.put(e.bulkComposition);
              w.put(e.timestep);

              const stableAssemblage &a = e.result;
              w.put(static_cast<std::int32_t>(a.len_oxides));
              w.put(a.bulk_composition);
              w.put(a.meltFraction);
              w.put(a.mineral_names);
              w.put(a.mineral_proportions);
              w.put(static_cast<std::uint32_t>(a.mineral_oxides.size()));
              for (const auto &row : a.mineral_oxides) w.put(row);
              w.put(a.stability_of_solution);
              w.put(a.bulk_density);
              w.put(a.bulk_entropy);
              w.put(a.bulk_enthalpy);
              w.put(a.bulk_specific_enthalpy);
              w.put(a.bulk_cp);
              w.put(a.fO2);
              w.put(a.dQFM);
              w.put(a.fixed_density);
              w.put(a.fixed_entropy);
              w.put(a.fixed_enthalpy);
              w.put(a.fixed_cp);
              w.put(a.fixed_volume_molar);
              // Preserve the version-9 record layout without retaining the
              // rejected warm-start state.
              w.put(std::uint32_t(0));
            }

          const std::string compressed = dealii::Utilities::compress(w.buf);
          header.put(static_cast<std::uint64_t>(w.buf.size()));
          header.put(static_cast<std::uint64_t>(compressed.size()));
          header.put(checksum64(w.buf));

          return header.buf + compressed;
        }


        bool MAGEMin_Hash_table::save_binary(const std::string &filepath) const
        {
          std::ofstream out(filepath, std::ios::binary);
          if (!out.is_open()) return false;

          const std::string data = serialize_entries();
          out.write(data.data(), static_cast<std::streamsize>(data.size()));
          return out.good();
        }


        bool MAGEMin_Hash_table::merge_serialized_entries(
          const std::string &file,
          const bool clear_first)
        {
          Reader h {file.data(), file.data() + file.size()};
          std::uint32_t magic = 0, version = 0, saved_oxides = 0, n_entries = 0;
          std::vector<std::string> saved_oxide_names;
          std::vector<double> saved_X_scale;
          std::string saved_context_signature;
          double saved_P_scale = 0.0, saved_T_scale = 0.0;
          std::uint64_t payload_uncompressed = 0, payload_on_disk = 0;
          std::uint64_t saved_checksum = 0;
          h.get(magic);
          h.get(version);
          h.get(saved_oxides);
          h.get(saved_oxide_names);
          h.get(saved_P_scale);
          h.get(saved_T_scale);
          h.get(saved_X_scale);
          h.get(saved_context_signature);
          h.get(n_entries);
          h.get(payload_uncompressed);
          h.get(payload_on_disk);
          h.get(saved_checksum);
          if (!h.ok ||
              magic != MAGE_MAGIC || version != MAGE_VERSION ||
              saved_oxides != n_oxides ||
              (!oxide_names.empty() && saved_oxide_names != oxide_names) ||
              saved_P_scale != P_scale || saved_T_scale != T_scale ||
              saved_X_scale != X_scale ||
              saved_context_signature != context_signature ||
              static_cast<std::size_t>(h.end - h.p) != payload_on_disk)
            return false;

          std::string payload;
          try
            {
              payload = dealii::Utilities::decompress(
                          std::string(h.p, payload_on_disk));
            }
          catch (...)
            {
              return false;
            }
          if (payload.size() != payload_uncompressed) return false;
          if (checksum64(payload) != saved_checksum) return false;

          Reader r {payload.data(), payload.data() + payload.size()};
          std::vector<DecodedEntry> decoded;
          decoded.reserve(n_entries);
          for (std::uint32_t i = 0; i < n_entries; ++i)
            {
              DecodedEntry entry;
              r.get(entry.pressure);
              r.get(entry.temperature);
              r.get(entry.melt_fraction);
              r.get(entry.bulk_composition);
              r.get(entry.timestep);

              stableAssemblage &a = entry.result;
              std::int32_t len_oxides = 0;
              r.get(len_oxides);
              a.len_oxides = len_oxides;
              r.get(a.bulk_composition);
              r.get(a.meltFraction);
              a.oxideNames = saved_oxide_names;
              r.get(a.mineral_names);
              r.get(a.mineral_proportions);
              std::uint32_t n_rows = 0;
              r.get(n_rows);
              if (!r.ok) return false;
              a.mineral_oxides.resize(n_rows);
              for (auto &row : a.mineral_oxides) r.get(row);
              r.get(a.stability_of_solution);
              r.get(a.bulk_density);
              r.get(a.bulk_entropy);
              r.get(a.bulk_enthalpy);
              r.get(a.bulk_specific_enthalpy);
              r.get(a.bulk_cp);
              r.get(a.fO2);
              r.get(a.dQFM);
              r.get(a.fixed_density);
              r.get(a.fixed_entropy);
              r.get(a.fixed_enthalpy);
              r.get(a.fixed_cp);
              r.get(a.fixed_volume_molar);
              std::uint32_t n_guesses = 0;
              r.get(n_guesses);
              if (!r.ok) return false;
              for (std::uint32_t guess = 0; guess < n_guesses; ++guess)
                {
                  std::string phase_name, phase_type;
                  std::int32_t id = -1, em = -1;
                  std::vector<double> composition_variables;
                  r.get(phase_name);
                  r.get(phase_type);
                  r.get(id);
                  r.get(em);
                  r.get(composition_variables);
                }
              if (!r.ok || entry.bulk_composition.size() != n_oxides ||
                  a.len_oxides != static_cast<int>(n_oxides) ||
                  a.bulk_composition.size() != n_oxides)
                return false;
              decoded.emplace_back(std::move(entry));
            }

          // Do not alter a live cache unless the complete payload was valid.
          if (!r.ok || r.p != r.end || decoded.size() > table_size)
            return false;
          if (clear_first)
            clear();
          for (auto &entry : decoded)
            insert(entry.pressure, entry.temperature, entry.melt_fraction,
                   entry.bulk_composition, entry.result, entry.timestep);
          return true;
        }


        bool MAGEMin_Hash_table::load_binary(const std::string &filepath)
        {
          std::ifstream in(filepath, std::ios::binary);
          if (!in.is_open()) return false;

          in.seekg(0, std::ios::end);
          std::string file(static_cast<std::size_t>(in.tellg()), '\0');
          in.seekg(0, std::ios::beg);
          if (!file.empty())
            in.read(&file[0], static_cast<std::streamsize>(file.size()));
          if (!in.good()) return false;

          return merge_serialized_entries(file, true);
        }
      }
    }
  }
}
