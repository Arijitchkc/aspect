/**
 * Implementation of MAGEMin in ASPECT as an external library
 */
#include "aspect/material_model/interface.h"
#include "aspect/material_model/reactive_fluid_transport.h"
#include <aspect/postprocess/visualization/named_additional_outputs.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/reaction_model/meltMagemin.h>
#include <aspect/utilities.h>
#include <cstddef>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe_data.h>
#include <type_traits>
#include <vector>
#include "aspect/simulator_signals.h"

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {

      namespace
      {
        std::uint64_t cache_timestep(const unsigned int timestep,
                                     const bool store_timestep)
        {
          return (!store_timestep || timestep == numbers::invalid_unsigned_int)
                 ? MAGEMin_hash::timestep_unavailable
                 : timestep;
        }
      }


      template <int dim>
      class MageDiagnosticOutputs : public NamedAdditionalMaterialOutputs<dim>
      {
        public:
          // Register the names that will appear in ParaView.
          MageDiagnosticOutputs(const unsigned int n_points)
            : NamedAdditionalMaterialOutputs<dim>({"proportion_ol", "proportion_cpx","proportion_g","proportion_spl","proportion_opx","proportion_fsp", "mage_enthalpy_of_fusion", "mage_bulk_specific_enthalpy", "mage_log10_fO2", "mage_delta_QFM", "output_mechanism"}),
          proportion_ol(n_points, 0.0),
          proportion_cpx(n_points, 0.0),
          proportion_g(n_points, 0.0),
          proportion_spl(n_points, 0.0),
          proportion_opx(n_points, 0.0),
          proportion_fsp(n_points, 0.0),
          mage_enthalpy_of_fusion(n_points, 0.0),
          mage_bulk_specific_enthalpy(n_points, 0.0),
          mage_log10_fO2(n_points, std::numeric_limits<double>::quiet_NaN()),
          mage_delta_QFM(n_points, std::numeric_limits<double>::quiet_NaN()),
          output_mechanism(n_points, 0.0)
          {}

          // Route the output index to its data vector.
          std::vector<double> get_nth_output(const unsigned int idx) const override
          {
            AssertIndexRange(idx, 11);
            if (idx == 0) return proportion_ol;
            if (idx == 1) return proportion_cpx;
            if (idx == 2) return proportion_g;
            if (idx == 3) return proportion_spl;
            if (idx == 4) return proportion_opx;
            if (idx == 5) return proportion_fsp;
            if (idx == 6) return mage_enthalpy_of_fusion;
            if (idx == 7) return mage_bulk_specific_enthalpy;
            if (idx == 8) return mage_log10_fO2;
            if (idx == 9) return mage_delta_QFM;
            if (idx == 10) return output_mechanism;
            return proportion_ol;
          }

          // Store the named output data.
          std::vector<double> proportion_ol;
          std::vector<double> proportion_cpx;
          std::vector<double> proportion_g;
          std::vector<double> proportion_spl;
          std::vector<double> proportion_opx;
          std::vector<double> proportion_fsp;
          std::vector<double> mage_enthalpy_of_fusion;
          std::vector<double> mage_bulk_specific_enthalpy;
          std::vector<double> mage_log10_fO2;
          std::vector<double> mage_delta_QFM;
          std::vector<double> output_mechanism;
      };

      template <int dim>
      void meltMagemin<dim>::create_additional_named_outputs(
        typename aspect::MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        out.additional_outputs.push_back(std::make_shared<MageDiagnosticOutputs<dim>>(out.densities.size()));
      }



      template <int dim> meltMagemin<dim>::meltMagemin() = default;

      template <int dim>
      void meltMagemin<dim>::initializeNewandModern()
      {
        std::vector<std::string> mineral_names;
        if (surrogateMode == 1)
          {
            AssertThrow(nn.load(surrogate_directory),
                        ExcMessage("Failed to load neural-network models from: "
                                   + surrogate_directory));
            oxide_names = nn.oxide_names();
            mineral_names = nn.mineral_names();
          }
        else
          {
            wrap.initializeMAGEMin(const_cast<char *>(database_name.c_str()), &sAssemblage);
            oxide_names = wrap.get_oxide_names(database_name);
            mineral_names = wrap.get_mineral_names(database_name);
          }
        connect_to_signals();

        solid_oxide_field_names = configured_solid_oxide_field_names;
        if (solid_oxide_field_names.empty())
          for (const std::string &oxide_name : oxide_names)
            solid_oxide_field_names.push_back(oxide_name == "FeOt" ? "FeO" : oxide_name);
        liquid_oxide_field_names = configured_liquid_oxide_field_names;
        if (liquid_oxide_field_names.empty())
          {
            liquid_oxide_field_names.reserve(oxide_names.size());
            for (const std::string &oxide_name : oxide_names)
              liquid_oxide_field_names.push_back(
                (oxide_name == "FeOt" ? "FeO" : oxide_name) + "l");
          }

        AssertThrow(solid_oxide_field_names.size() == oxide_names.size(),
                    ExcMessage("'Solid oxide field names' must contain one field for each native "
                               "MAGEMin oxide in database '" + database_name + "'."));
        AssertThrow(liquid_oxide_field_names.size() == oxide_names.size(),
                    ExcMessage("'Liquid oxide field names' must contain one field for each native "
                               "MAGEMin oxide in database '" + database_name + "'."));
        for (unsigned int i = 0; i < oxide_names.size(); ++i)
          {
            AssertThrow(this->introspection().compositional_name_exists(solid_oxide_field_names[i]),
                        ExcMessage("The native MAGEMin oxide '" + oxide_names[i] +
                                   "' is mapped to missing solid compositional field '" +
                                   solid_oxide_field_names[i] + "'."));
            AssertThrow(this->introspection().compositional_name_exists(liquid_oxide_field_names[i]),
                        ExcMessage("The native MAGEMin oxide '" + oxide_names[i] +
                                   "' is mapped to missing liquid compositional field '" +
                                   liquid_oxide_field_names[i] + "'."));
          }

        water_oxide_index = std::numeric_limits<unsigned int>::max();
        for (unsigned int i = 0; i < oxide_names.size(); ++i)
          if (oxide_names[i] == "H2O")
            water_oxide_index = i;

        if (have_per_oxide_scales)
          AssertThrow(hash_X_scales.size() == oxide_names.size(),
                      ExcMessage("'Hash table X scales per oxide' must contain one value for each "
                                 "native MAGEMin oxide in database '" + database_name + "'."));

        if (surrogateMode == 0 && preseed_kdtree_enabled == 1)
          AssertThrow(preseed_bulk_composition.size() == oxide_names.size(),
                      ExcMessage("'Preseed bulk composition' must contain one value for each "
                                 "native MAGEMin oxide in database '" + database_name + "'."));

        setup_hash_tiers();
        if (surrogateMode == 0)
          clear_kdtree();

        this->get_pcout()
            << std::endl
            << "==========================================================" << std::endl
            << "   " << (surrogateMode == 1 ? "Neural-network" : "MAGEMin")
            << " thermodynamics initialized" << std::endl
            << "==========================================================" << std::endl;
        if (surrogateMode == 1)
          this->get_pcout() << "  Model directory : " << surrogate_directory << std::endl;
        else
          this->get_pcout() << "  Database        : " << database_name << std::endl;
        this->get_pcout()
            << "  Fine hash size  : " << size_of_hash_table << " entries" << std::endl
            << "  Coarse hash     : "
            << (surrogateMode == 0 && enable_two_tier == 1 ? "ON" : "OFF");
        if (surrogateMode == 0 && enable_two_tier == 1)
          this->get_pcout() << " (" << size_of_coarse_hash_table << " entries)";
        this->get_pcout()
            << std::endl
            << "  Hash scales     : P=" << hash_table_P_scale
            << ", T=" << hash_table_T_scale
            << ", X=" << hash_table_X_scale << std::endl
            << "  KD-tree         : "
            << (surrogateMode == 0 && useKd_tree == 1 ? "ON" : "OFF") << std::endl
            << "  Oxides (" << oxide_names.size() << ")    : ";
        for (const std::string &oxide : oxide_names)
          this->get_pcout() << oxide << " ";
        this->get_pcout() << std::endl
                          << "  Phases (" << mineral_names.size() << ")    : ";
        for (const std::string &mineral : mineral_names)
          this->get_pcout() << mineral << " ";
        this->get_pcout() << std::endl
                          << "=========================================================="
                          << std::endl;

        kdtree_hits = 0;
        preseed_done = false;

      }


      /**
      * Pre-seed the KD-tree by running MAGEMin on a P-T grid with the
      * initial bulk composition. This maps phase boundaries before the
      * simulation starts, so the KD-tree can detect them from timestep 1.
      */
      template <int dim>
      void meltMagemin<dim>::preseed_kdtree()
      {
        const MPI_Comm comm = this->get_mpi_communicator();
        const unsigned int rank = Utilities::MPI::this_mpi_process(comm);
        const unsigned int n_ranks = Utilities::MPI::n_mpi_processes(comm);

        // Use cutoff depth as max depth (no point seeding below melting range).
        const double max_depth = cutOff_depth;  // meters
        const double min_depth = 0.0;

        // Convert depth to pressure (approximate: rho * g * depth)
        const double min_P_GPa = min_depth * 3.3e4 / 1e9;  // ~0 GPa at surface
        const double max_P_GPa = max_depth * 3.3e4 / 1e9;

        // Build P-T grid
        std::vector<double> P_values, T_values;
        for (double P = min_P_GPa; P <= max_P_GPa + 1e-10; P += preseed_P_spacing)
          P_values.push_back(P);
        for (double T = preseed_T_min; T <= preseed_T_max + 1e-10; T += preseed_T_spacing)
          T_values.push_back(T);

        // Normalize the bulk composition
        std::vector<double> comp = preseed_bulk_composition;
        double comp_sum = 0;
        for (const double value : comp)
          comp_sum += value;
        if (comp_sum > 0.01)
          for (double &value : comp)
            value /= comp_sum;

        struct PreseedPoint
        {
          double pressure_GPa;
          double temperature_C;
          bool near_solidus;
        };
        std::vector<PreseedPoint> grid;
        for (const double P_GPa : P_values)
          for (const double T_C : T_values)
            {
              const double P_Pa = P_GPa * 1e9;
              const double dry_solidus = A1 + 273.15 + A2 * P_Pa
                                         + A3 * P_Pa * P_Pa;
              const double saturation = Xie1 * std::pow(P_Pa, Lambda)
                                        + Xie2 * P_Pa;
              const double bulk_water = water_oxide_index < comp.size()
                                        ? comp[water_oxide_index]*100 : 0.0;
              const double water_shift = std::min(
                                           K * std::pow(std::max(0.0, bulk_water), gamma),
                                           K * std::pow(std::max(0.0, saturation), gamma));
              const double dT = T_C + 273.15 - (dry_solidus - water_shift);
              if (dT >= -fine_band_below_solidus)
                grid.push_back({P_GPa, T_C, dT <= fine_band_above_solidus});
            }

        // Keep the distributed preseed separate from the live caches while it
        // is built. This avoids exchanging unrelated run-time entries.
        MAGEMin_hash::MAGEMin_Hash_table local_fine;
        local_fine.set_oxide_names(oxide_names);
        local_fine.set_context_signature(cache_context_signature);
        local_fine.set_scales(hash_table_P_scale, hash_table_T_scale,
                              hash_table_X_scale);
        if (have_per_oxide_scales)
          local_fine.set_oxide_scales(hash_X_scales);
        local_fine.resize_hash_table(size_of_hash_table);

        MAGEMin_hash::MAGEMin_Hash_table local_coarse;
        if (enable_two_tier == 1)
          {
            local_coarse.set_oxide_names(oxide_names);
            local_coarse.set_context_signature(cache_context_signature);
            local_coarse.set_scales(hash_table_P_scale / coarse_tier_PT_factor,
                                    hash_table_T_scale / coarse_tier_PT_factor,
                                    hash_table_X_scale / coarse_tier_X_factor);
            if (have_per_oxide_scales)
              {
                std::vector<double> coarse_x = hash_X_scales;
                for (double &scale : coarse_x)
                  scale /= coarse_tier_X_factor;
                local_coarse.set_oxide_scales(coarse_x);
              }
            local_coarse.resize_hash_table(size_of_coarse_hash_table);
          }

        // Run each eligible grid point on exactly one rank.
        int dummyArg1 = 0;
        char **dummyArgc1 = nullptr;
        const int len_oxides = oxide_names.size();
        char *database = const_cast<char *>(database_name.c_str());
        std::size_t local_solves = 0;
        std::size_t local_kept = 0;
        std::size_t local_count_with_melt = 0;

        const auto t_start = std::chrono::high_resolution_clock::now();

        for (std::size_t point_index = rank;
             point_index < grid.size(); point_index += n_ranks)
          {
            const PreseedPoint &point = grid[point_index];
            const double P_kbar = point.pressure_GPa * 10.0;
            const double T_K = point.temperature_C + 273.15;

            const std::vector<std::vector<double>> bulkComposition = {comp};
            const std::vector<double> temperatures = {T_K};
            const std::vector<double> pressures = {P_kbar};

            wrap.executeMAGEMin(dummyArg1, dummyArgc1, temperatures,
                                pressures, len_oxides, database,
                                bulkComposition, &sAssemblage);
            ++local_solves;

            double meltFraction = 0.0;
            for (size_t k = 0; k < sAssemblage.mineral_names.size(); ++k)
              if (sAssemblage.mineral_names[k] == "liq")
                {
                  meltFraction = sAssemblage.mineral_proportions[k];
                  break;
                }
            if (meltFraction > 1e-5)
              ++local_count_with_melt;

            if (meltFraction > 1e-5 || point.near_solidus)
              {
                auto &cache = (enable_two_tier == 1 && !point.near_solidus)
                              ? local_coarse : local_fine;
                double Pq = point.pressure_GPa, Tq = point.temperature_C;
                std::vector<double> compq = comp;
                if (snap_queries_to_bins == 1)
                  cache.snap_query(Pq, Tq, compq);
                cache.insert(Pq, Tq, meltFraction, compq, sAssemblage,
                             store_cache_timestep
                             ? MAGEMin_hash::timestep_preseed
                             : MAGEMin_hash::timestep_unavailable);
                ++local_kept;
              }
          }

        // All ranks receive the same exact preseed cache. The serialized form
        // is already compressed and is also used by checkpoint/restart.
        const auto fine_data = Utilities::MPI::all_gather(
                                 comm, local_fine.serialize_entries());
        for (const std::string &data : fine_data)
          AssertThrow(magemin_cache.merge_serialized_entries(data),
                      ExcMessage("Failed to exchange the fine MAGEMin preseed cache."));

        if (enable_two_tier == 1)
          {
            const auto coarse_data = Utilities::MPI::all_gather(
                                       comm, local_coarse.serialize_entries());
            for (const std::string &data : coarse_data)
              AssertThrow(magemin_cache_coarse.merge_serialized_entries(data),
                          ExcMessage("Failed to exchange the coarse MAGEMin preseed cache."));
          }

        // Rebuild identical KD-tree metadata from the merged cache. Repeating
        // cheap lookups is far smaller than repeating the thermodynamic solves.
        std::size_t shared_points = 0;
        for (const PreseedPoint &point : grid)
          {
            auto &cache = (enable_two_tier == 1 && !point.near_solidus)
                          ? magemin_cache_coarse : magemin_cache;
            double Pq = point.pressure_GPa;
            double Tq = point.temperature_C;
            double melt_fraction = 0.0;
            std::vector<double> compq = comp;
            if (snap_queries_to_bins == 1)
              cache.snap_query(Pq, Tq, compq);
            const stableAssemblage *result =
              cache.lookup(Pq, Tq, compq, melt_fraction, false);
            if (result == nullptr) continue;

            std::uint64_t assemblage_mask = 0;
            for (std::size_t m = 0;
                 m < result->mineral_proportions.size() && m < 64; ++m)
              if (result->mineral_proportions[m] > 1e-4)
                assemblage_mask |= (std::uint64_t(1) << m);
            add_kdtree_point(Pq, Tq, compq, assemblage_mask,
                             enable_two_tier != 1 || point.near_solidus);
            ++shared_points;
          }

        const std::size_t count = Utilities::MPI::sum(local_solves, comm);
        const std::size_t kept = Utilities::MPI::sum(local_kept, comm);
        const std::size_t count_with_melt =
          Utilities::MPI::sum(local_count_with_melt, comm);
        const double elapsed_s = Utilities::MPI::max(
                                   std::chrono::duration<double>(
                                     std::chrono::high_resolution_clock::now() - t_start).count(), comm);

        this->get_pcout() << "MAGEMin preseed: solved=" << count
                          << ", retained=" << shared_points << "/" << kept
                          << ", melt=" << count_with_melt
                          << ", seconds=" << elapsed_s << std::endl;
      }

      template <int dim> void meltMagemin<dim>::connect_to_signals()
      {
        // Get ASPECT's signal system
        SimulatorSignals<dim> &signals = this->get_signals();

        // Connect our cache clearing callback to start_timestep signal
        signals.start_timestep.connect(
          std::bind(&meltMagemin<dim>::clear_cache_on_timestep,
                    this,
                    std::placeholders::_1)
        );

        this->get_pcout()
            << (surrogateMode == 1 ? "Neural-network" : "MAGEMin")
            << " cache connected to ASPECT timestep signals" << std::endl;

        // Connect to checkpoint/restart signals for cache persistence
        if (save_hash_table == 1 && surrogateMode == 0)
          {
            signals.pre_checkpoint_store_user_data.connect(
              std::bind(&meltMagemin<dim>::save_caches_on_checkpoint,
                        this,
                        std::placeholders::_1)
            );
            signals.post_resume_load_user_data.connect(
              std::bind(&meltMagemin<dim>::load_caches_on_resume,
                        this,
                        std::placeholders::_1)
            );
          }
      }


      /// Function to save caches in binary
      template <int dim>
      void meltMagemin<dim>::save_caches_on_checkpoint(
        typename parallel::distributed::Triangulation<dim> &/*tria*/)
      {
        const int rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

        const std::string hash_file = hash_storage_location + "magemin_cache_rank_"
                                      + std::to_string(rank) + ".bin";
        const std::string hash_file_c = hash_storage_location + "magemin_cache_coarse_rank_"
                                        + std::to_string(rank) + ".bin";

        const bool hash_ok = magemin_cache.save_binary(hash_file);
        const bool coarse_ok = (enable_two_tier == 1)
                               ? magemin_cache_coarse.save_binary(hash_file_c) : true;

        if (rank == 0)
          this->get_pcout() << "   MAGEMin cache checkpoint: hash="
                            << (hash_ok ? "OK" : "FAILED")
                            << ", coarse=" << (coarse_ok ? "OK" : "FAILED")
                            << " (" << magemin_cache.getOccupancy() << " entries)"
                            << std::endl;
      }

      // Function to load caches in binary
      template <int dim>
      void meltMagemin<dim>::load_caches_on_resume(
        typename parallel::distributed::Triangulation<dim> &/*tria*/)
      {
        if (cache_loaded_before_resume)
          {
            cache_loaded_before_resume = false;
            return;
          }

        const int rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

        const std::string hash_file = hash_storage_location + "magemin_cache_rank_"
                                      + std::to_string(rank) + ".bin";
        const std::string hash_file_c = hash_storage_location + "magemin_cache_coarse_rank_"
                                        + std::to_string(rank) + ".bin";

        const bool hash_ok = magemin_cache.load_binary(hash_file);
        const bool coarse_ok = (enable_two_tier == 1)
                               ? magemin_cache_coarse.load_binary(hash_file_c) : true;

        clear_kdtree();
        const auto add_cache_entries = [this](
                                         const MAGEMin_hash::MAGEMin_Hash_table &cache,
                                         const bool uses_fine_cache)
        {
          for (const MAGEMin_hash::Entry &entry : cache.getEntries())
            {
              if (!entry.occupied)
                continue;
              std::uint64_t assemblage_mask = 0;
              for (std::size_t mineral = 0;
                   mineral < entry.result.mineral_proportions.size() && mineral < 64;
                   ++mineral)
                if (entry.result.mineral_proportions[mineral] > 1e-4)
                  assemblage_mask |= (std::uint64_t(1) << mineral);
              add_kdtree_point(entry.pressure, entry.temperature,
                               entry.bulkComposition, assemblage_mask,
                               uses_fine_cache);
            }
        };
        if (hash_ok)
          add_cache_entries(magemin_cache, true);
        if (enable_two_tier == 1 && coarse_ok)
          add_cache_entries(magemin_cache_coarse, false);

        // A saved cache already contains any preseed grid from the original
        // run. Only skip pre-seeding when every rank loaded its files. This
        // keeps the later collective preseed call consistent when a restart
        // uses a different number of MPI ranks.
        const unsigned int all_caches_loaded = Utilities::MPI::min(
                                                 static_cast<unsigned int>(hash_ok && coarse_ok),
                                                 this->get_mpi_communicator());
        preseed_done = (all_caches_loaded == 1);

        if (rank == 0)
          this->get_pcout() << "   MAGEMin cache resume: hash="
                            << (hash_ok ? "OK" : "not found/failed")
                            << ", coarse=" << (coarse_ok ? "OK" : "not found/failed")
                            << " (" << magemin_cache.getOccupancy() << " entries)"
                            << "; rebuilt KD-tree with "
                            << kdtree_points.points.size() << " points; "
                            << (preseed_done ? "complete on all ranks"
                                : "missing rank files, preseed scheduled")
                            << std::endl;
      }

      template <int dim>
      void meltMagemin<dim>::clear_cache_on_timestep(
        const SimulatorAccess<dim> &simulator_access)
      {
        const MPI_Comm comm = this->get_mpi_communicator();
        const std::size_t local_hits =
          magemin_cache.getHits() + magemin_cache_coarse.getHits();
        const std::size_t local_misses =
          magemin_cache.getMisses() + magemin_cache_coarse.getMisses();
        const std::size_t local_occupancy =
          magemin_cache.getOccupancy() + magemin_cache_coarse.getOccupancy();
        const std::size_t local_kd_size = kdtree_points.points.size();

        const std::size_t total_hits = Utilities::MPI::sum(local_hits, comm);
        const std::size_t total_misses = Utilities::MPI::sum(local_misses, comm);
        const std::size_t total_kd_hits = Utilities::MPI::sum(kdtree_hits, comm);
        const std::size_t total_calls =
          Utilities::MPI::sum(total_actual_magemin_calls, comm);
        const std::size_t total_queries =
          Utilities::MPI::sum(total_points_evaluated, comm);
        const std::size_t total_occupancy = Utilities::MPI::sum(local_occupancy, comm);
        const std::size_t total_kd_size = Utilities::MPI::sum(local_kd_size, comm);
        const std::size_t total_nn = Utilities::MPI::sum(count_nn_calls, comm);
        const std::size_t total_ood = Utilities::MPI::sum(count_nn_ood, comm);

        const std::size_t max_queries = Utilities::MPI::max(total_points_evaluated, comm);
        const std::size_t min_queries = Utilities::MPI::min(total_points_evaluated, comm);
        const std::size_t max_hits = Utilities::MPI::max(local_hits, comm);
        const std::size_t min_hits = Utilities::MPI::min(local_hits, comm);
        const std::size_t max_misses = Utilities::MPI::max(local_misses, comm);
        const std::size_t min_misses = Utilities::MPI::min(local_misses, comm);
        const std::size_t max_occupancy = Utilities::MPI::max(local_occupancy, comm);
        const std::size_t min_occupancy = Utilities::MPI::min(local_occupancy, comm);
        const std::size_t max_probe = Utilities::MPI::max(
                                        std::max(magemin_cache.getMaxProbe(), magemin_cache_coarse.getMaxProbe()), comm);
        const std::size_t max_kd_size = Utilities::MPI::max(local_kd_size, comm);
        const std::size_t min_kd_size = Utilities::MPI::min(local_kd_size, comm);
        const std::size_t max_kd_hits = Utilities::MPI::max(kdtree_hits, comm);
        const std::size_t min_kd_hits = Utilities::MPI::min(kdtree_hits, comm);
        const std::size_t max_calls = Utilities::MPI::max(total_actual_magemin_calls, comm);
        const std::size_t min_calls = Utilities::MPI::min(total_actual_magemin_calls, comm);
        const std::size_t max_nn = Utilities::MPI::max(count_nn_calls, comm);
        const std::size_t min_nn = Utilities::MPI::min(count_nn_calls, comm);
        const std::size_t max_ood = Utilities::MPI::max(count_nn_ood, comm);
        const std::size_t min_ood = Utilities::MPI::min(count_nn_ood, comm);

        const double max_hash_time = Utilities::MPI::max(time_hash_us, comm);
        const double min_hash_time = Utilities::MPI::min(time_hash_us, comm);
        const double max_kd_time = Utilities::MPI::max(time_kdtree_us, comm);
        const double min_kd_time = Utilities::MPI::min(time_kdtree_us, comm);
        const double max_magemin_time = Utilities::MPI::max(time_magemin_us, comm);
        const double min_magemin_time = Utilities::MPI::min(time_magemin_us, comm);
        const double max_nn_time = Utilities::MPI::max(time_nn_us, comm);
        const double min_nn_time = Utilities::MPI::min(time_nn_us, comm);
        const double total_magemin_time = Utilities::MPI::sum(time_magemin_us, comm);
        const double total_nn_time = Utilities::MPI::sum(time_nn_us, comm);

        const auto percentage = [total_queries](const std::size_t value)
        {
          return total_queries > 0 ? 100.0 * value / total_queries : 0.0;
        };
        const double hit_rate = total_hits + total_misses > 0
                                ? 100.0 * total_hits / (total_hits + total_misses)
                                : 0.0;

        this->get_pcout()
            << std::endl
            << "==========================================================" << std::endl
            << "   " << (surrogateMode == 1 ? "Neural-network" : "MAGEMin")
            << " Statistics (Timestep "
            << simulator_access.get_timestep_number() << ")" << std::endl
            << "   [max|min] = per-rank extremes" << std::endl
            << "==========================================================" << std::endl
            << "  Total queries : " << total_queries
            << "  [" << max_queries << "|" << min_queries << "]" << std::endl
            << "  Resolution:" << std::endl
            << "    Hash table  : " << total_hits << " (" << percentage(total_hits) << "%)" << std::endl;

        if (surrogateMode == 1)
          this->get_pcout()
              << "    Neural net  : " << total_nn << " (" << percentage(total_nn) << "%)" << std::endl;
        else
          this->get_pcout()
              << "    KD-tree     : " << total_kd_hits << " (" << percentage(total_kd_hits) << "%)" << std::endl
              << "    MAGEMin     : " << total_calls << " (" << percentage(total_calls) << "%)" << std::endl;

        this->get_pcout()
            << "----------------------------------------------------------" << std::endl
            << "  [1] Hash table:" << std::endl
            << "      Hits       : " << total_hits << "  [" << max_hits << "|" << min_hits
            << "]  " << hit_rate << "%" << std::endl
            << "      Misses     : " << total_misses << "  [" << max_misses << "|" << min_misses << "]" << std::endl
            << "      Occupancy  : " << total_occupancy << "  [" << max_occupancy << "|" << min_occupancy << "]" << std::endl
            << "      Max probe  : " << max_probe << std::endl;

        if (surrogateMode == 1)
          {
            const std::size_t total_inlier = total_nn >= total_ood ? total_nn - total_ood : 0;
            this->get_pcout()
                << "----------------------------------------------------------" << std::endl
                << "  [2] Neural network:" << std::endl
                << "      States      : " << total_nn << "  [" << max_nn << "|" << min_nn << "]" << std::endl
                << "      In-domain   : " << total_inlier << " ("
                << (total_nn > 0 ? 100.0 * total_inlier / total_nn : 0.0) << "% of NN)" << std::endl
                << "      OOD flagged : " << total_ood << "  [" << max_ood << "|" << min_ood << "] ("
                << (total_nn > 0 ? 100.0 * total_ood / total_nn : 0.0) << "% of NN)" << std::endl
                << "      Average     : " << (total_nn > 0 ? total_nn_time / total_nn : 0.0)
                << " us/state" << std::endl
                << "----------------------------------------------------------" << std::endl
                << "  [3] Timing (s) [max|min]:" << std::endl
                << "      Hash        : [" << max_hash_time / 1e6 << "|" << min_hash_time / 1e6 << "]" << std::endl
                << "      Neural net  : [" << max_nn_time / 1e6 << "|" << min_nn_time / 1e6 << "]" << std::endl;
          }
        else
          {
            this->get_pcout()
                << "----------------------------------------------------------" << std::endl
                << "  [2] KD-tree:" << std::endl
                << "      Points      : " << total_kd_size << "  [" << max_kd_size << "|" << min_kd_size << "]" << std::endl
                << "      Hits        : " << total_kd_hits << "  [" << max_kd_hits << "|" << min_kd_hits << "]" << std::endl
                << "----------------------------------------------------------" << std::endl
                << "  [3] MAGEMin solver:" << std::endl
                << "      Calls       : " << total_calls << "  [" << max_calls << "|" << min_calls << "]" << std::endl
                << "      Average     : " << (total_calls > 0 ? total_magemin_time / total_calls : 0.0)
                << " us/call" << std::endl
                << "----------------------------------------------------------" << std::endl
                << "  [4] Timing (s) [max|min]:" << std::endl
                << "      Hash        : [" << max_hash_time / 1e6 << "|" << min_hash_time / 1e6 << "]" << std::endl
                << "      KD-tree     : [" << max_kd_time / 1e6 << "|" << min_kd_time / 1e6 << "]" << std::endl
                << "      MAGEMin     : [" << max_magemin_time / 1e6 << "|" << min_magemin_time / 1e6 << "]" << std::endl;
          }
        this->get_pcout() << "==========================================================" << std::endl;

        const auto evict_if_full = [](MAGEMin_hash::MAGEMin_Hash_table &cache,
                                      const std::size_t capacity)
        {
          if (capacity > 0 &&
              static_cast<double>(cache.getOccupancy()) / capacity > 0.75)
            {
              cache.reduce_hash_table_size(0.5);
              return true;
            }
          return false;
        };

        bool evicted = evict_if_full(magemin_cache, size_of_hash_table);
        if (enable_two_tier == 1)
          evicted = evict_if_full(magemin_cache_coarse,
                                  size_of_coarse_hash_table) || evicted;
        if (evicted)
          this->get_pcout() << "Hash table eviction triggered on this timestep" << std::endl;
        if (evicted && useKd_tree == 1)
          clear_kdtree();

        magemin_cache.reset_access_counts();
        magemin_cache_coarse.reset_access_counts();
        if (useKd_tree == 1 && kdtree_points.points.size() > 1000000)
          {
            clear_kdtree();
            this->get_pcout() << "KD-tree limit reached; tree cleared" << std::endl;
          }

        if (!preseed_done && this->get_timestep_number() >= 1 &&
            preseed_kdtree_enabled == 1 && useKd_tree == 1 &&
            surrogateMode == 0)
          {
            preseed_kdtree();
            preseed_done = true;
          }

        total_actual_magemin_calls = 0;
        total_points_evaluated = 0;
        kdtree_hits = 0;
        time_hash_us = time_kdtree_us = time_magemin_us = time_nn_us = 0.0;
        count_nn_calls = count_nn_ood = 0;
      }


      /**
      * Function to guess melt fraction based on Katz 2003
      */


      template <int dim>
      double meltMagemin<dim>::guess_MeltFraction(const double pressure,
                                                  const double temperature,
                                                  const double bulk_water) const
      {
        // Katz parameters are assumed to be ASPECT-side pressure units:
        // pressure in Pa, temperature in K.
        const double T_solidus_dry =
          A1 + 273.15 + A2 * pressure + A3 * pressure * pressure;

        const double T_lherz_liquidus_dry =
          B1 + 273.15 + B2 * pressure + B3 * pressure * pressure;

        const double T_liquidus_dry =
          C1 + 273.15 + C2 * pressure + C3 * pressure * pressure;

        // The hydrous formulation expects bulk_water in wt%.
        const auto X_H2O_partitioned = [&](const double F) -> double
        {
          return bulk_water / (D_H2O + F * (1.0 - D_H2O));
        };

        const auto X_H2O_saturated = [&]() -> double
        {
          return Xie1 * std::pow(std::max(0.0, pressure), Lambda)
          + Xie2 * std::max(0.0, pressure);
        };

        const auto deltaT_H2O = [&](const double X_H2O) -> double
        {
          return K * std::pow(std::max(0.0, X_H2O), gamma);
        };

        const double deltaT_saturated = deltaT_H2O(X_H2O_saturated());

        const auto melt_fraction_shifted = [&](const double deltaT) -> double
        {
          const double T_solidus = T_solidus_dry - deltaT;
          const double T_lherz_liquidus = T_lherz_liquidus_dry - deltaT;
          const double T_liquidus = T_liquidus_dry - deltaT;

          double peridotite_melt_fraction;

          if (temperature < T_solidus || pressure > 1.3e10)
            peridotite_melt_fraction = 0.0;
          else if (temperature > T_lherz_liquidus)
            peridotite_melt_fraction = 1.0;
          else
            peridotite_melt_fraction =
            std::pow((temperature - T_solidus)
            / (T_lherz_liquidus - T_solidus),
            beta);

          const double R_cpx = r1 + r2 * std::max(0.0, pressure);
          const double F_max = M_cpx / R_cpx;

          if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
            {
              const double T_max =
              std::pow(F_max, 1.0 / beta)
              * (T_lherz_liquidus - T_solidus)
              + T_solidus;

              peridotite_melt_fraction =
              F_max
              + (1.0 - F_max)
              * std::pow((temperature - T_max)
              / (T_liquidus - T_max),
              beta);
            }

          return std::max(0.0, std::min(1.0, peridotite_melt_fraction));
        };

        // The anhydrous branch does not require the implicit partitioning solve.
        if (bulk_water <= 0.0)
          return melt_fraction_shifted(0.0);

        // Solve F = melt_fraction(T, P, deltaT(F)).
        const auto root_equation = [&](const double F) -> double
        {
          const double X_H2O = X_H2O_partitioned(F);

          const double deltaT =
          std::min(deltaT_H2O(X_H2O), deltaT_saturated);

          return melt_fraction_shifted(deltaT) - F;
        };

        double a = 0.0;
        double b = 1.0;

        double fa = root_equation(a);
        double fb = root_equation(b);

        if (fa * fb > 0.0)
          return std::abs(fa) < std::abs(fb) ? a : b;

        constexpr double tolerance = 1e-6;
        constexpr unsigned int max_iterations = 100;

        for (unsigned int iteration = 0; iteration < max_iterations; ++iteration)
          {
            const double m = 0.5 * (a + b);
            const double fm = root_equation(m);

            if (std::abs(fm) < tolerance || 0.5 * (b - a) < tolerance)
              return m;

            if (fa * fm < 0.0)
              {
                b = m;
                fb = fm;
              }
            else
              {
                a = m;
                fa = fm;
              }
          }

        return 0.5 * (a + b);
      }




      // ASPECT's EnthalpyOutputs stores latent enthalpy in J/kg. Convert
      // MAGEMin's molar phase entropies and return T(S_liquid-S_solid).
      // latent_heat_melt.cc applies the sign needed by the energy equation.
      template <int dim>
      double meltMagemin<dim>::compute_enthalpy_of_fusion(const std::vector<std::string> &names,
                                                          const std::vector<double> &props,
                                                          const std::vector<double> &entropies,
                                                          const std::vector<double> &densities,
                                                          const std::vector<double> &volumes,
                                                          const double temperature) const
      {
        double S_liq_specific = 0.0;
        double S_solid_specific = 0.0;
        double total_solid_wt = 0.0;
        bool found_liq = false;

        for (size_t k = 0; k < names.size(); ++k)
          {
            double prop = (k < props.size()) ? props[k] : 0.0;
            if (prop < 1e-5) continue;

            double s   = (k < entropies.size()) ? entropies[k] : 0.0;
            double rho = (k < densities.size())  ? densities[k]  : 0.0;
            double vol = (k < volumes.size())     ? volumes[k]    : 0.0;

            if (rho<1000.0 || rho > 10000.0 || vol <= 0.0 || !std::isfinite(s)) continue;

            const double molar_mass = rho * vol * 1e-3; // g/mol
            const double S_specific = s * 1e6 / molar_mass; // J/(kg K); s is kJ/(mol K)

            if (names[k] == "liq")
              {
                S_liq_specific = S_specific;
                found_liq = true;
              }
            else
              {
                S_solid_specific += prop * S_specific;
                total_solid_wt += prop;
              }
          }

        if (!found_liq || total_solid_wt < 1e-5)
          {
            return 0.0;
          }

        S_solid_specific /= total_solid_wt;
        return std::max(temperature * (S_liq_specific - S_solid_specific), 0.0);
      }


      template <int dim>
      double
      meltMagemin<dim>::
      compute_bulk_density (const double porosity,
                            const double solid_density,
                            const double fluid_density) const
      {
        return (1 - porosity) * solid_density + porosity * fluid_density;
      }

      template <int dim>
      double
      meltMagemin<dim>::
      compute_mass_fraction (const double volume_frac,
                             const double material_density,
                             const double bulk_density) const
      {
        return volume_frac * material_density / bulk_density;
      }


      /**
      * This function fills the data struct with all data initially
      */
      template<int dim> void meltMagemin<dim>::load_variables_from_present_time_step(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in, unsigned int q) const
      {
        qp.oxideLength = oxide_names.size();

        qp.idx_phi_porosity = this->introspection().compositional_index_for_name("porosity");
        if (this->introspection().compositional_name_exists("depletion"))
          qp.idx_depletion = this->introspection().compositional_index_for_name("depletion");
        if (this->introspection().compositional_name_exists("crust"))
          qp.idx_phi_crust = this->introspection().compositional_index_for_name("crust");

        for (std::size_t i=0; i<static_cast<std::size_t>(qp.oxideLength); ++i)
          {
            qp.idx_X_solid[i] = this->introspection().compositional_index_for_name(solid_oxide_field_names[i]);
            qp.idx_X_melt[i] = this->introspection().compositional_index_for_name(liquid_oxide_field_names[i]);
            qp.X_solid[i] = std::clamp(in.composition[q][qp.idx_X_solid[i]], 0.0, 1.0);
            qp.X_melt[i] = std::clamp(in.composition[q][qp.idx_X_melt[i]], 0.0, 1.0);
          }

        double volume_fraction_porosity = std::clamp(in.composition[q][qp.idx_phi_porosity], 0.0, 1.0);
        qp.bulk_density_current = compute_bulk_density(volume_fraction_porosity, qp.density_solid_current, qp.density_liquid_current);

        // COnverting volume fraction porosity to mass fraction as MAGEMin uses mass fractions
        qp.phi_melt = compute_mass_fraction(volume_fraction_porosity, qp.density_liquid_current, qp.bulk_density_current);
        qp.phi_solid = 1.0 - qp.phi_melt;

        qp.depletion = std::clamp(in.composition[q][qp.idx_depletion], 0.0, 1.0);
        qp.phi_crust = std::clamp(in.composition[q][qp.idx_phi_crust], 0.0, 1.0);

        qp.depth_m = this->get_geometry_model().depth(in.position[q]);
        qp.P_Pa = std::max(this->get_adiabatic_conditions().pressure(in.position[q]),0.0);
        qp.T_K = std::max(in.temperature[q], 273.15);
        qp.P_kbar = qp.P_Pa/1e8;
        qp.T_C = qp.T_K - 273.15;

        qp.sum_liquid=0.0;
        qp.sum_solid=0.0;
        qp.amount_of_melt_reacting_with_solid = 0.0;
        qp.amount_of_free_melt = 0.0;
      }

      // template<int dim> void meltMagemin<dim>::normalize_drifted_compositions(quad_point_data &qp) const
      // {
      //   // Advection and diffusion can make phase compositions drift away from a
      //   // unit sum. Correct only substantial drift so small numerical changes are
      //   // not amplified before the equilibrium calculation.
      //   qp.liquid_valid = true;
      //
      //   for (int i = 0; i < qp.oxideLength; ++i)
      //     {
      //       qp.sum_liquid += qp.X_melt[i];
      //       qp.sum_solid += qp.X_solid[i];
      //     }
      //
      //
      //
      //   // Normalize a substantially drifted solid composition.
      //   if (qp.phi_solid > 0.01 && qp.sum_solid > 0.0 &&
      //       (qp.sum_solid<0.9 || qp.sum_solid>1.1))
      //     {
      //       double s_scale = 1.0 / qp.sum_solid;
      //       for (int i = 0; i < qp.oxideLength; i++)
      //         {
      //           double current =qp.X_solid[i];
      //           qp.X_solid[i] = current * s_scale;
      //         }
      //       qp.sum_solid=1.0;
      //     }
      //
      //   // Normalize a substantially drifted liquid composition.
      //   if (qp.phi_melt > 0.00001 && qp.sum_liquid > 0.6 && (qp.sum_liquid<0.9 || qp.sum_liquid>1.1))
      //     {
      //       double l_scale = 1.0 / qp.sum_liquid;
      //       for (int i = 0; i < qp.oxideLength; i++)
      //         {
      //           double current = qp.X_melt[i];
      //           qp.X_melt[i] = current * l_scale;
      //         }
      //       qp.sum_liquid=1.0;
      //     }
      // }

      template<int dim> void meltMagemin<dim>::normalize_compositions(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in,
                                                                      typename Interface<dim>::MaterialModelOutputs &out,
                                                                      const unsigned int q) const
      {
        // Normalize a substantially drifted solid composition.
        qp.sum_solid = 0.0;
        qp.sum_liquid = 0.0;
        for (int i = 0; i < qp.oxideLength; ++i)
          {
            qp.sum_solid += qp.X_solid[i];
            qp.sum_liquid += qp.X_melt[i];
          }
        if (qp.phi_melt <= 1e-5)
          {
            qp.phi_melt = 0.0;
            qp.phi_solid = 1.0;

            for (int i = 0; i < qp.oxideLength; ++i)
              {
                out.reaction_terms[q][qp.idx_X_melt[i]] = -in.composition[q][qp.idx_X_melt[i]];
              }
            std::fill(qp.X_melt.begin(),
                      qp.X_melt.end(),
                      0.0);
            qp.sum_liquid = 0.0;
          }
        if (qp.phi_solid > 1e-5)
          {
            AssertThrow(std::isfinite(qp.sum_solid) &&
                        qp.sum_solid > 1e-5,
                        ExcMessage("Nonzero solid has zero oxide composition."));
            for (int i = 0; i < qp.oxideLength; ++i)
              qp.X_solid[i] /= qp.sum_solid;

            qp.sum_solid = 1.0;
          }

        if (qp.phi_melt > 1e-5)
          {
            AssertThrow(std::isfinite(qp.sum_liquid) &&
                        qp.sum_liquid > 1e-5,
                        ExcMessage("Nonzero melt has zero oxide composition."));

            for (int i = 0; i < qp.oxideLength; ++i)
              qp.X_melt[i] /= qp.sum_liquid;

            qp.sum_liquid = 1.0;
          }


      }


      /**
      * If melt goes below a certain temperature, convert all liquid at that point to solid
      */
      template<int dim> double meltMagemin<dim>::freeze_melt_below_certain_tenperature(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in,typename Interface<dim>::MaterialModelOutputs &out, unsigned q) const
      {
        // Default: no change
        for (int i = 0; i < qp.oxideLength; i++)
          {
            out.reaction_terms[q][qp.idx_X_solid[i]] = 0.0;
            out.reaction_terms[q][qp.idx_X_melt[i]] = 0.0;
          }
        if (qp.phi_melt < 1e-5)
          {
            return qp.phi_melt;
          }
        // const bool valid = qp.liquid_valid && qp.sum_liquid > 0.8 && qp.sum_liquid < 1.2;

        // if (!valid)
        //   {
        //     for (int i = 0; i < qp.oxideLength; ++i)
        //       {
        //         out.reaction_terms[q][qp.idx_X_melt[i]] = -in.composition[q][qp.idx_X_melt[i]];
        //         out.reaction_terms[q][qp.idx_X_solid[i]] = 0.0;
        //       }
        //     return 0.0;
        //   }


        // valid melt: mix into solid, remove melt
        for (int i = 0; i < qp.oxideLength; ++i)
          {
            const double bulk_new = (qp.X_solid[i] * (1-qp.phi_melt) + qp.X_melt[i] * qp.phi_melt);
            out.reaction_terms[q][qp.idx_X_solid[i]] = bulk_new - in.composition[q][qp.idx_X_solid[i]];
            out.reaction_terms[q][qp.idx_X_melt[i]] = -in.composition[q][qp.idx_X_melt[i]];
          }
        qp.phi_crust = qp.phi_crust+ qp.phi_melt; //std::clamp(qp.phi_melt, 0.0, 1.0);
        return 0.0;
      }



      template<int dim> void meltMagemin<dim>::handle_cases_when_there_is_no_solid_left(quad_point_data &qp, typename Interface<dim>::MaterialModelOutputs &out, unsigned q) const
      {

        for (int i=0; i<qp.oxideLength; i++)
          {
            out.reaction_terms[q][qp.idx_X_solid[i]] =0;
            out.reaction_terms[q][qp.idx_X_melt[i]] =0;
          }
      }


      template<int dim> bool meltMagemin<dim>::skip_magemin_evaluation(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out, unsigned int q) const
      {

        const double bulk_water = water_oxide_index < qp.X_solid.size()
                                  ? qp.X_solid[water_oxide_index]*100 : 0.0;
        const double guessedMeltFraction =
          guess_MeltFraction(qp.P_Pa, qp.T_K, bulk_water);

        if (guessedMeltFraction > 1e-5 || qp.phi_melt > 1e-5 ||
            (qp.T_K > qp.T_solidus_K &&
             this->get_adiabatic_conditions().pressure(in.position[q]) > 0.0))
          return false;
        else
          return true;

        // // The solids have already been normalized previously.
        // if (qp.sum_solid<0.9 || qp.sum_solid > 1.1)
        //   {
        //     double s_sum = 0;
        //     for (int i = 0; i < qp.oxideLength; i++)
        //       s_sum += std::max(in.composition[q][qp.idx_X_solid[i]], 0.0);
        //     if (s_sum > 0.5 &&  (s_sum<0.9 || s_sum>1.2))
        //       {
        //         const double s_scale = 1.0 / s_sum;
        //         for (int i = 0; i < qp.oxideLength; i++)
        //           {
        //             const double current = std::max(in.composition[q][qp.idx_X_solid[i]], 0.0);
        //             out.reaction_terms[q][qp.idx_X_solid[i]] = (current * s_scale) - current;
        //           }
        //       }
        //   }
        // else
        //   for (int i = 0; i < qp.oxideLength; i++)
        //     out.reaction_terms[q][qp.idx_X_solid[i]] =
        //       qp.X_solid[i] - std::max(in.composition[q][qp.idx_X_solid[i]], 0.0);
        //
        //
        // Finite-element liquid fields are cleared below the solidus. Liquid
        // particles retain a normalized dormant composition so cell-average
        // interpolation cannot dilute valid liquid chemistry with zeros.

        // for (int i = 0; i < qp.oxideLength; i++)
        //   out.reaction_terms[q][qp.idx_X_melt[i]] = -in.composition[q][qp.idx_X_melt[i]];

        // return true;
      }

      template<int dim> void meltMagemin<dim>::calculate_bulk_composition(quad_point_data &qp) const
      {

        qp.amount_of_melt_reacting_with_solid = std::min(qp.phi_melt, melt_retention_threshold);
        qp.amount_of_free_melt = std::max(qp.phi_melt - qp.amount_of_melt_reacting_with_solid, 0.0);


        // if (qp.sum_liquid<0.8 || qp.sum_liquid>1.2)
        //   {
        //     qp.liquid_valid=false;
        //     qp.amount_of_melt_reacting_with_solid = 0.0;
        //     qp.amount_of_free_melt = qp.phi_melt;
        //   }

        const double participating_volume = qp.phi_solid + qp.amount_of_melt_reacting_with_solid;

        for (size_t i=0; i<static_cast<size_t>(qp.oxideLength); i++)
          {
            // ASPECT transports porosity as volume fractions.
            // Preserve that convention when constructing the bulk
            // composition passed to MAGEMin. With complete retention,
            // participating_volume is one and this reduces to
            // X_solid*(1-porosity) + X_liquid*porosity.
            qp.X_bulk_passed_to_magemin_individual[i] = std::clamp((qp.X_solid[i] * qp.phi_solid
                                                                    + qp.X_melt[i]
                                                                    * qp.amount_of_melt_reacting_with_solid)/participating_volume, 0.0, 1.0);
          }

        double bulk_comp_passed_to_magemin_sum = 0.0;
        for (int i = 0; i < qp.oxideLength; ++i)
          bulk_comp_passed_to_magemin_sum += qp.X_bulk_passed_to_magemin_individual[i];

        for (int i = 0; i <qp.oxideLength; ++i)
          qp.X_bulk_passed_to_magemin_individual[i] /= bulk_comp_passed_to_magemin_sum;


        // These are the final form of input to any surrogate/magemin
        qp.X_bulk_passed_to_magemin= {qp.X_bulk_passed_to_magemin_individual};
      }


      template<int dim>
      void meltMagemin<dim>::fill_stableAssemblage_from_nn(const padawan::Prediction &pred,
                                                           const std::size_t prediction_index,
                                                           quad_point_data &qp) const
      {
        const auto &mins = nn.mineral_names();
        const auto &oxs = nn.oxide_names();

        qp.assemblage.len_oxides = static_cast<int>(oxs.size());
        qp.assemblage.oxideNames = oxs;
        qp.assemblage.meltFraction = 0.0;
        qp.assemblage.bulk_specific_enthalpy =
          std::numeric_limits<double>::quiet_NaN();
        // The current neural-network surrogate does not predict redox outputs.
        qp.assemblage.fO2 = std::numeric_limits<double>::quiet_NaN();
        qp.assemblage.dQFM = std::numeric_limits<double>::quiet_NaN();

        static const std::vector<std::string> ig_mage_mins=
        {
          // Solution phases
          "spl",   "bi",    "cd",    "cpx",   "ep",    "g",     "amp",
          "ilm",   "liq",   "ol",    "opx",   "fsp",   "fl",    "mu",
          "fper",  "chl",   "ne",
          // Pure phases
          "q",     "crst",  "trd",   "coe",   "stv",   "ky",    "sill",
          "and",   "ru",    "sph",
          // Buffers / activities / special
          "O2",    "H2O",   "cor",   "qfm",   "mw",    "qif",   "nno",
          "hm",    "iw",    "cco",   "aH2O",  "aO2",   "aMgO",  "aFeO",
          "aAl2O3","aTiO2"
        };


        const size_t M = ig_mage_mins.size();
        qp.assemblage.mineral_proportions.assign(M, 0.0);
        qp.assemblage.mineral_oxides.assign(M, std::vector<double>(oxs.size(), 0.0));
        qp.assemblage.fixed_density.assign(M, 0.0);
        qp.assemblage.fixed_entropy.assign(M, 0.0);
        qp.assemblage.fixed_enthalpy.assign(M, 0.0);
        qp.assemblage.fixed_cp.assign(M, 0.0);
        qp.assemblage.fixed_volume_molar.assign(M, 0.0);
        qp.assemblage.mineral_names = ig_mage_mins;
        qp.fraction_of_melt = 0;

        const bool has_physical = prediction_index < pred.physical.size() &&
                                  pred.physical[prediction_index].size() == 4;
        const double rho_solid = has_physical ? pred.physical[prediction_index][0] : 0.0;
        const double rho_liquid = has_physical ? pred.physical[prediction_index][1] : 0.0;
        const double cp_solid = has_physical ? pred.physical[prediction_index][2] : 0.0;
        const double cp_liquid = has_physical ? pred.physical[prediction_index][3] : 0.0;

        for (size_t n = 0; n < mins.size(); n++)
          {
            const auto it = std::find(ig_mage_mins.begin(), ig_mage_mins.end(), mins[n]);
            if (it == ig_mage_mins.end()) continue;
            const size_t m = static_cast<size_t>(it - ig_mage_mins.begin());

            const double prop = (prediction_index < pred.proportions.size() &&
                                 n < pred.proportions[prediction_index].size())
                                ? static_cast<double>(pred.proportions[prediction_index][n]) : 0.0;
            qp.assemblage.mineral_proportions[m] = prop;
            qp.assemblage.fixed_density[m] = mins[n] == "liq" ? rho_liquid : rho_solid;
            qp.assemblage.fixed_cp[m] = mins[n] == "liq" ? cp_liquid : cp_solid;

            if (prediction_index < pred.mineral_oxides.size() &&
                n < pred.mineral_oxides[prediction_index].size())
              for (size_t k = 0;
                   k < oxs.size() && k < pred.mineral_oxides[prediction_index][n].size();
                   ++k)
                {
                  qp.assemblage.mineral_oxides[m][k] =
                    static_cast<double>(pred.mineral_oxides[prediction_index][n][k]);
                }

            if (mins[n] == "liq")
              {
                qp.fraction_of_melt = prop;
                qp.assemblage.meltFraction = prop;
              }
          }

      }




      template<int dim>
      void meltMagemin<dim>::run_neural_network_batch(
        std::vector<quad_point_data> &quad_points,
        const std::vector<unsigned int> &point_indices) const
      {
        if (point_indices.empty())
          return;

        std::vector<double> pressures;
        std::vector<double> temperatures;
        std::vector<std::vector<double>> compositions;
        pressures.reserve(point_indices.size());
        temperatures.reserve(point_indices.size());
        compositions.reserve(point_indices.size());

        for (const unsigned int point_index : point_indices)
          {
            quad_point_data &qp = quad_points[point_index];
            std::vector<double> composition = qp.X_bulk_passed_to_magemin_individual;
            double sum = 0.0;
            for (const double value : composition)
              sum += value;
            if (sum > 0.5)
              for (double &value : composition)
                value /= sum;

            pressures.push_back(qp.P_kbar);
            temperatures.push_back(qp.T_C);
            compositions.push_back(std::move(composition));
          }

        const auto start_time = std::chrono::high_resolution_clock::now();
        const padawan::Prediction prediction =
          nn.predict_batch(pressures, temperatures, compositions,
                           use_MAGEMin_density == 1 ||
                           use_MAGEMin_specific_heat == 1);
        time_nn_us += std::chrono::duration<double, std::micro>(
                        std::chrono::high_resolution_clock::now() - start_time).count();
        count_nn_calls += point_indices.size();

        for (std::size_t batch_index = 0; batch_index < point_indices.size(); ++batch_index)
          {
            quad_point_data &qp = quad_points[point_indices[batch_index]];
            if (batch_index >= prediction.is_ood.size() || prediction.is_ood[batch_index] == 1)
              {
                ++count_nn_ood;
                qp.run_type = 20;
                qp.fraction_of_melt = qp.phi_melt;
                continue;
              }

            fill_stableAssemblage_from_nn(prediction, batch_index, qp);
            qp.run_type = 5;
            qp.active = &qp.assemblage;
            extract_stable_assemblages_and_physical_parameters(qp);

            // This cache belongs to the NN execution path and starts empty.
            // Repeated states therefore avoid another TorchScript evaluation.
            magemin_cache.insert(qp.P_kbar / 10.0,
                                 qp.T_C,
                                 qp.fraction_of_melt,
                                 qp.X_bulk_passed_to_magemin_individual,
                                 qp.assemblage,
                                 cache_timestep(this->get_timestep_number(),
                                                store_cache_timestep));
          }
      }

      template<int dim> bool meltMagemin<dim>::try_running_hash_table(quad_point_data &qp) const
      {
        auto t0 = std::chrono::high_resolution_clock::now();
        double cached_melt_fraction = 0.0;
        const stableAssemblage *a = pick_tier(qp).lookup(qp.P_kbar/10.0,
                                                         qp.T_C,
                                                         qp.X_bulk_passed_to_magemin_individual,
                                                         cached_melt_fraction);
        auto t1 = std::chrono::high_resolution_clock::now();
        time_hash_us += std::chrono::duration<double, std::micro>(t1 - t0).count();
        if (a==nullptr)
          {
            return false;
          }
        qp.run_type=10;
        qp.active = a;
        return true;
      }



      template <int dim>
      void meltMagemin<dim>::clear_kdtree() const
      {
        kdtree_points.points.clear();
        kdtree = std::make_unique<KDTree>(oxide_names.size() + 2,
                                          kdtree_points);
      }



      template <int dim>
      void meltMagemin<dim>::add_kdtree_point(
        const double pressure_GPa,
        const double temperature_C,
        const std::vector<double> &composition,
        const std::uint64_t assemblage_mask,
        const bool uses_fine_cache) const
      {
        if (useKd_tree != 1)
          return;

        if (kdtree == nullptr)
          clear_kdtree();

        MAGEMinLookupPoint point;
        point.coordinates.assign(oxide_names.size() + 2, 0.0);
        point.coordinates[0] = pressure_GPa * hash_table_P_scale;
        point.coordinates[1] = temperature_C * hash_table_T_scale;
        AssertDimension(composition.size(), oxide_names.size());
        for (unsigned int i = 0; i < oxide_names.size(); ++i)
          {
            const double scale = have_per_oxide_scales
                                 ? hash_X_scales[i] : hash_table_X_scale;
            point.coordinates[i+2] = composition[i] * scale;
          }
        point.assemblage_mask = assemblage_mask;
        point.uses_fine_cache = uses_fine_cache;

        kdtree_points.points.push_back(point);
        kdtree->addPoint(kdtree_points.points.size()-1);
      }



      template<int dim>
      void meltMagemin<dim>::fallback_to_magemin_batch(
        std::vector<quad_point_data> &quad_points,
        const std::vector<unsigned int> &point_indices) const
      {
        if (point_indices.empty())
          return;

        // Equal snapped inputs need only one MAGEMin calculation. Sorting a vector
        // keeps the grouping deterministic and avoids another owning container.
        std::vector<unsigned int> sorted_indices = point_indices;
        std::sort(sorted_indices.begin(), sorted_indices.end(),
                  [&quad_points, this](const unsigned int left_index,
                                       const unsigned int right_index)
        {
          const quad_point_data &left = quad_points[left_index];
          const quad_point_data &right = quad_points[right_index];
          const bool left_uses_fine_cache = enable_two_tier != 1 || left.near_solidus;
          const bool right_uses_fine_cache = enable_two_tier != 1 || right.near_solidus;

          if (left_uses_fine_cache != right_uses_fine_cache)
            return left_uses_fine_cache < right_uses_fine_cache;
          if (left.P_kbar != right.P_kbar)
            return left.P_kbar < right.P_kbar;
          if (left.T_C != right.T_C)
            return left.T_C < right.T_C;
          return std::lexicographical_compare(
                   left.X_bulk_passed_to_magemin_individual.begin(),
                   left.X_bulk_passed_to_magemin_individual.end(),
                   right.X_bulk_passed_to_magemin_individual.begin(),
                   right.X_bulk_passed_to_magemin_individual.end());
        });

        const auto have_same_key = [&quad_points, this](const unsigned int left_index,
                                                        const unsigned int right_index)
        {
          const quad_point_data &left = quad_points[left_index];
          const quad_point_data &right = quad_points[right_index];
          const bool left_uses_fine_cache = enable_two_tier != 1 || left.near_solidus;
          const bool right_uses_fine_cache = enable_two_tier != 1 || right.near_solidus;
          return left_uses_fine_cache == right_uses_fine_cache &&
                 left.P_kbar == right.P_kbar &&
                 left.T_C == right.T_C &&
                 left.X_bulk_passed_to_magemin_individual ==
                 right.X_bulk_passed_to_magemin_individual;
        };

        std::vector<std::vector<unsigned int>> groups;
        for (const unsigned int point_index : sorted_indices)
          {
            if (groups.empty() || !have_same_key(groups.back().front(), point_index))
              groups.emplace_back();
            groups.back().push_back(point_index);
          }

        std::vector<double> temperatures;
        std::vector<double> pressures;
        std::vector<std::vector<double>> compositions;
        temperatures.reserve(groups.size());
        pressures.reserve(groups.size());
        compositions.reserve(groups.size());
        for (const std::vector<unsigned int> &group : groups)
          {
            const quad_point_data &qp = quad_points[group.front()];
            temperatures.push_back(qp.T_K);
            pressures.push_back(qp.P_kbar);
            compositions.push_back(qp.X_bulk_passed_to_magemin.front());
          }

        int dummy_argument = 0;
        char **dummy_arguments = nullptr;
        std::vector<stableAssemblage> results;
        const quad_point_data &first_point = quad_points[groups.front().front()];
        const auto start_time = std::chrono::high_resolution_clock::now();
        wrap.executeMAGEMin(dummy_argument, dummy_arguments, temperatures,
                            pressures, first_point.oxideLength,
                            first_point.database, compositions, results);
        time_magemin_us += std::chrono::duration<double, std::micro>(
                             std::chrono::high_resolution_clock::now() - start_time).count();
        total_actual_magemin_calls += groups.size();

        AssertThrow(results.size() == groups.size(),
                    ExcMessage("MAGEMin returned a different number of results than requested."));

        for (std::size_t result_index = 0; result_index < groups.size(); ++result_index)
          {
            const unsigned int representative_index = groups[result_index].front();
            quad_point_data &representative = quad_points[representative_index];
            double melt_fraction = 0.0;
            for (std::size_t mineral = 0;
                 mineral < results[result_index].mineral_names.size();
                 ++mineral)
              if (results[result_index].mineral_names[mineral] == "liq")
                {
                  melt_fraction = results[result_index].mineral_proportions[mineral];
                  break;
                }

            MAGEMin_hash::MAGEMin_Hash_table &cache = pick_tier(representative);
            cache.insert(representative.P_kbar / 10.0,
                         representative.T_C,
                         melt_fraction,
                         representative.X_bulk_passed_to_magemin_individual,
                         results[result_index],
                         cache_timestep(this->get_timestep_number(),
                                        store_cache_timestep));

            if (useKd_tree == 1 && this->get_timestep_number() > 0)
              {
                std::uint64_t assemblage_mask = 0;
                for (std::size_t mineral = 0;
                     mineral < results[result_index].mineral_proportions.size() && mineral < 64;
                     ++mineral)
                  if (results[result_index].mineral_proportions[mineral] > 1e-4)
                    assemblage_mask |= (std::uint64_t(1) << mineral);

                add_kdtree_point(representative.P_kbar / 10.0,
                                 representative.T_C,
                                 representative.X_bulk_passed_to_magemin_individual,
                                 assemblage_mask,
                                 enable_two_tier != 1 || representative.near_solidus);
              }

            for (const unsigned int point_index : groups[result_index])
              {
                quad_point_data &qp = quad_points[point_index];
                qp.run_type = 15;
                qp.fraction_of_melt = melt_fraction;
                double cached_melt_fraction = 0.0;
                qp.active = pick_tier(qp).lookup(qp.P_kbar / 10.0,
                                                 qp.T_C,
                                                 qp.X_bulk_passed_to_magemin_individual,
                                                 cached_melt_fraction,
                                                 false);
                if (qp.active == nullptr)
                  {
                    qp.assemblage = results[result_index];
                    qp.active = &qp.assemblage;
                  }
              }
          }
      }
      template <int dim>
      bool meltMagemin<dim>::try_kdtree_lookup(quad_point_data &qp) const
      {
        if (useKd_tree != 1 || kdtree == nullptr ||
            kdtree_points.points.size() < static_cast<std::size_t>(num_required_points_for_kd_tree))
          return false;

        std::vector<double> query(oxide_names.size() + 2, 0.0);
        query[0] = qp.P_kbar / 10.0 * hash_table_P_scale;
        query[1] = qp.T_C * hash_table_T_scale;
        for (unsigned int i = 0; i < oxide_names.size(); ++i)
          {
            const double scale = have_per_oxide_scales
                                 ? hash_X_scales[i] : hash_table_X_scale;
            query[i+2] = qp.X_bulk_passed_to_magemin_individual[i] * scale;
          }

        KDTree::BoundingBox search_box;
        search_box.resize(query.size());
        search_box[0].low = query[0] - kdtree_pressure_tol_bins;
        search_box[0].high = query[0] + kdtree_pressure_tol_bins;
        search_box[1].low = query[1] - kdtree_temperature_tol_bins;
        search_box[1].high = query[1] + kdtree_temperature_tol_bins;
        for (unsigned int i = 0; i < oxide_names.size(); ++i)
          {
            search_box[i+2].low = query[i+2] - kdtree_oxide_tol_bins;
            search_box[i+2].high = query[i+2] + kdtree_oxide_tol_bins;
          }

        const auto t0 = std::chrono::high_resolution_clock::now();
        std::vector<std::size_t> candidate_indices;
        nanoflann::BoxResultSet<std::size_t> result_set(candidate_indices);
        const std::size_t n_candidates = kdtree->findWithinBox(result_set, search_box);
        (void) n_candidates;
        time_kdtree_us += std::chrono::duration<double, std::micro>(
                            std::chrono::high_resolution_clock::now() - t0).count();

        struct Neighbour
        {
          std::size_t point_index;
          double distance;
          double melt_fraction;
          const stableAssemblage *assemblage;
        };

        std::vector<Neighbour> neighbours;
        neighbours.reserve(candidate_indices.size());
        std::vector<std::size_t> stale_points;
        std::vector<double> composition(oxide_names.size());

        for (const std::size_t point_index : candidate_indices)
          {
            const MAGEMinLookupPoint &point = kdtree_points.points[point_index];
            const double pressure_GPa = point.coordinates[0] / hash_table_P_scale;
            const double temperature_C = point.coordinates[1] / hash_table_T_scale;
            for (unsigned int i = 0; i < oxide_names.size(); ++i)
              {
                const double scale = have_per_oxide_scales
                                     ? hash_X_scales[i] : hash_table_X_scale;
                composition[i] = point.coordinates[i+2] / scale;
              }

            double cached_melt_fraction = 0.0;
            MAGEMin_hash::MAGEMin_Hash_table &cache = point.uses_fine_cache
                                                      ? magemin_cache : magemin_cache_coarse;
            const stableAssemblage *assemblage = cache.lookup(pressure_GPa,
                                                              temperature_C,
                                                              composition,
                                                              cached_melt_fraction,
                                                              false);
            if (assemblage == nullptr)
              {
                stale_points.push_back(point_index);
                continue;
              }

            double distance = 0.0;
            for (unsigned int d = 0; d < query.size(); ++d)
              {
                const double difference = point.coordinates[d] - query[d];
                distance += difference * difference;
              }
            neighbours.push_back({point_index, distance, cached_melt_fraction, assemblage});
          }

        for (const std::size_t point_index : stale_points)
          kdtree->removePoint(point_index);

        if (neighbours.size() < static_cast<std::size_t>(minimum_neighbours_to_trust_kd_tree))
          return false;

        std::sort(neighbours.begin(), neighbours.end(),
                  [](const Neighbour &left, const Neighbour &right)
        {
          return left.distance < right.distance;
        });

        const std::uint64_t reference_mask =
          kdtree_points.points[neighbours[0].point_index].assemblage_mask;
        for (const Neighbour &neighbour : neighbours)
          if (kdtree_points.points[neighbour.point_index].assemblage_mask != reference_mask)
            {
              return false;
            }

        double minimum_melt = std::numeric_limits<double>::max();
        double maximum_melt = -std::numeric_limits<double>::max();
        for (const Neighbour &neighbour : neighbours)
          {
            minimum_melt = std::min(minimum_melt, neighbour.melt_fraction);
            maximum_melt = std::max(maximum_melt, neighbour.melt_fraction);
          }
        if (maximum_melt - minimum_melt > max_melt_range)
          {
            return false;
          }

        const unsigned int n_minerals = neighbours[0].assemblage->mineral_proportions.size();
        for (unsigned int mineral = 0; mineral < n_minerals; ++mineral)
          {
            double minimum = std::numeric_limits<double>::max();
            double maximum = -std::numeric_limits<double>::max();
            for (const Neighbour &neighbour : neighbours)
              {
                const double value = neighbour.assemblage->mineral_proportions[mineral];
                minimum = std::min(minimum, value);
                maximum = std::max(maximum, value);
              }
            if (maximum - minimum > max_mineral_range)
              {
                return false;
              }
          }

        // Mineral proportions can remain smooth while the composition of a
        // phase changes appreciably. Check every oxide of every present phase
        // before reusing a complete assemblage, because these compositions are
        // subsequently used by the reactive-transport coupling.
        for (unsigned int mineral = 0; mineral < n_minerals; ++mineral)
          {
            if ((reference_mask & (std::uint64_t(1) << mineral)) == 0)
              continue;

            for (unsigned int oxide = 0; oxide < oxide_names.size(); ++oxide)
              {
                double minimum = std::numeric_limits<double>::max();
                double maximum = -std::numeric_limits<double>::max();
                for (const Neighbour &neighbour : neighbours)
                  {
                    if (mineral >= neighbour.assemblage->mineral_oxides.size() ||
                        oxide >= neighbour.assemblage->mineral_oxides[mineral].size())
                      {
                        return false;
                      }
                    const double value =
                      neighbour.assemblage->mineral_oxides[mineral][oxide];
                    minimum = std::min(minimum, value);
                    maximum = std::max(maximum, value);
                  }

                if (maximum - minimum > max_oxide_deviation)
                  {
                    return false;
                  }
              }
          }

        qp.fraction_of_melt = neighbours[0].melt_fraction;
        qp.run_type = 12;
        // KD reuse is approximate and read-only. Never insert this result into
        // either hash tier or export it as an exact NN-training state.
        qp.active = neighbours[0].assemblage;
        ++kdtree_hits;
        return true;
      }

      template<int dim> void meltMagemin<dim>::extract_stable_assemblages_and_physical_parameters(quad_point_data &qp)       const
      {
        const stableAssemblage &a = *qp.active;
        qp.fraction_of_melt = std::min(std::max(a.meltFraction, 0.0), 1.0);

        int liq_idx = -1;
        for (std::size_t k = 0; k < a.mineral_names.size(); ++k)
          if (a.mineral_names[k] == "liq")
            {
              liq_idx = static_cast<int>(k);
              break;
            }

        double solid_specific_volume = 0.0;
        double solid_weight_fraction = 0.0;
        for (std::size_t j = 0; j < a.mineral_proportions.size(); ++j)
          if (std::isfinite(a.mineral_proportions[j]) &&
              a.mineral_proportions[j] > 0.0 &&
              j < a.fixed_density.size() &&
              std::isfinite(a.fixed_density[j]) &&
              a.fixed_density[j] > 0.0)
            {
              const double specific_volume =
                a.mineral_proportions[j] / a.fixed_density[j];
              if (static_cast<int>(j) != liq_idx)
                {
                  solid_specific_volume += specific_volume;
                  solid_weight_fraction += a.mineral_proportions[j];
                }
            }

        // Do not convert the liquid phase proportion to a volume fraction using
        // phase densities. Use the phase proportion directly as porosity.
        // const double total_specific_volume =
        //   liquid_specific_volume + solid_specific_volume;
        // qp.fraction_of_melt = total_specific_volume > 0.0
        //                       ? liquid_specific_volume / total_specific_volume
        //                       : 0.0;

        std::fill(qp.X_liq_magemin.begin(),   qp.X_liq_magemin.end(),   0.0);
        std::fill(qp.X_solid_magemin.begin(), qp.X_solid_magemin.end(), 0.0);
        qp.density_solid_new = 0.0;
        qp.density_liquid_new  = 0.0;
        qp.specific_heat_solid_new = 0.0;

        const double inv_solid_weight_fraction = solid_weight_fraction > 0.0
                                                 ? 1.0 / solid_weight_fraction
                                                 : 0.0;

        for (std::size_t j = 0; j < a.mineral_proportions.size(); ++j)
          {
            const double prop = a.mineral_proportions[j];
            if (prop < 1e-7) continue;

            if (static_cast<int>(j) == liq_idx)
              {
                for (int i = 0; i < a.len_oxides; ++i)
                  qp.X_liq_magemin[i] += a.mineral_oxides[j][i];
                qp.density_liquid_new = a.fixed_density[j];
              }
            else
              {
                const double w = prop * inv_solid_weight_fraction;
                for (int i = 0; i < a.len_oxides; ++i)
                  qp.X_solid_magemin[i] += w * a.mineral_oxides[j][i];
                if (j < a.fixed_cp.size())
                  qp.specific_heat_solid_new += w * a.fixed_cp[j];
              }
          }

        if (solid_specific_volume > 0.0)
          qp.density_solid_new = solid_weight_fraction / solid_specific_volume;

        // existing solid-oxide normalization guard
        double s_sum = 0.0;
        for (double x : qp.X_solid_magemin) s_sum += x;
        if (s_sum > 0.01)
          for (double &x : qp.X_solid_magemin) x /= s_sum;

        qp.enthalpy_of_fusion = compute_enthalpy_of_fusion(
                                  a.mineral_names, a.mineral_proportions, a.fixed_entropy,
                                  a.fixed_density, a.fixed_volume_molar, qp.T_K);
      }


      template<int dim>
      void meltMagemin<dim>::resolve_magemin_batch(
        std::vector<quad_point_data> &quad_points,
        const std::vector<unsigned int> &point_indices) const
      {
        total_points_evaluated += point_indices.size();

        for (const unsigned int point_index : point_indices)
          {
            quad_point_data &qp = quad_points[point_index];
            qp.active = nullptr;

            if (snap_queries_to_bins == 1)
              {
                MAGEMin_hash::MAGEMin_Hash_table &cache = pick_tier(qp);
                double pressure_GPa = qp.P_kbar / 10.0;
                cache.snap_query(pressure_GPa,
                                 qp.T_C,
                                 qp.X_bulk_passed_to_magemin_individual);
                qp.P_kbar = pressure_GPa * 10.0;
                qp.T_K = qp.T_C + 273.15;

                std::vector<double> solver_composition =
                  qp.X_bulk_passed_to_magemin_individual;
                double sum = 0.0;
                for (const double value : solver_composition)
                  sum += value;
                if (sum > 0.5)
                  for (double &value : solver_composition)
                    value /= sum;

                // Keep the snapped composition as the cache key. Normalize only
                // the solver copy because MAGEMin expects a unit-sum bulk.
                qp.X_bulk_passed_to_magemin = {solver_composition};
              }
          }

        std::vector<unsigned int> hash_misses;
        hash_misses.reserve(point_indices.size());
        for (const unsigned int point_index : point_indices)
          if (!try_running_hash_table(quad_points[point_index]))
            hash_misses.push_back(point_index);

        // The NN route is self-contained: reuse its own in-memory predictions,
        // then evaluate all new states in one TorchScript batch. It never
        // reaches the KD-tree or MAGEMin path below.
        if (surrogateMode == 1)
          {
            run_neural_network_batch(quad_points, hash_misses);
            for (const unsigned int point_index : point_indices)
              if (quad_points[point_index].run_type == 10)
                extract_stable_assemblages_and_physical_parameters(
                  quad_points[point_index]);
            return;
          }

        std::vector<unsigned int> kdtree_misses;
        kdtree_misses.reserve(hash_misses.size());
        for (const unsigned int point_index : hash_misses)
          if (!try_kdtree_lookup(quad_points[point_index]))
            kdtree_misses.push_back(point_index);

        fallback_to_magemin_batch(quad_points, kdtree_misses);

        for (const unsigned int point_index : point_indices)
          {
            quad_point_data &qp = quad_points[point_index];
            AssertThrow(qp.active != nullptr,
                        ExcMessage("No thermodynamic result was assigned to a quadrature point."));
            extract_stable_assemblages_and_physical_parameters(qp);
          }
      }


      template<int dim>
      void meltMagemin<dim>::write_to_aspect_data_structures(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out, unsigned int q, MeltOutputs<dim> *melt_out, EnthalpyOutputs<dim> *enthalpy_out) const
      {
        const std::shared_ptr<MageDiagnosticOutputs<dim>> diag_out =
                                                         out.template get_additional_output_object<MageDiagnosticOutputs<dim>>();
        if (diag_out != nullptr)
          {
            diag_out->mage_enthalpy_of_fusion[q] = qp.enthalpy_of_fusion;
            diag_out->mage_bulk_specific_enthalpy[q] =
              qp.active->bulk_specific_enthalpy;
            diag_out->mage_log10_fO2[q] = qp.active->fO2;
            diag_out->mage_delta_QFM[q] = qp.active->dQFM;
            diag_out->output_mechanism[q] = qp.run_type;

            const stableAssemblage &assemblage = *qp.active;
            double total_specific_volume = 0.0;
            for (std::size_t phase = 0;
                 phase < assemblage.mineral_proportions.size();
                 ++phase)
              if (phase < assemblage.fixed_density.size() &&
                  assemblage.fixed_density[phase] > 0.0)
                total_specific_volume += assemblage.mineral_proportions[phase]
                                         / assemblage.fixed_density[phase];

            for (std::size_t mineral = 0;
                 mineral < assemblage.mineral_names.size() &&
                 mineral < assemblage.mineral_proportions.size();
                 ++mineral)
              {
                const double proportion =
                  mineral < assemblage.fixed_density.size() &&
                  assemblage.fixed_density[mineral] > 0.0 &&
                  total_specific_volume > 0.0
                  ? (assemblage.mineral_proportions[mineral]
                     / assemblage.fixed_density[mineral])
                  / total_specific_volume
                  : 0.0;
                if (assemblage.mineral_names[mineral] == "ol")
                  diag_out->proportion_ol[q] = proportion;
                else if (assemblage.mineral_names[mineral] == "cpx")
                  diag_out->proportion_cpx[q] = proportion;
                else if (assemblage.mineral_names[mineral] == "g")
                  diag_out->proportion_g[q] = proportion;
                else if (assemblage.mineral_names[mineral] == "spl")
                  diag_out->proportion_spl[q] = proportion;
                else if (assemblage.mineral_names[mineral] == "opx")
                  diag_out->proportion_opx[q] = proportion;
                else if (assemblage.mineral_names[mineral] == "fsp")
                  diag_out->proportion_fsp[q] = proportion;
              }
          }

        // Apply the equilibrated melt volume fraction to the complete
        // solid-plus-melt system.
        const double equilibrated_volume =
          qp.phi_solid + qp.amount_of_melt_reacting_with_solid;

        // Phase fractions relative to the complete system.
        const double equilibrium_melt_fraction =
          std::isfinite(qp.fraction_of_melt)
          ? std::min(std::max(qp.fraction_of_melt, 0.0), 1.0)
          : 0.0;
        const double new_melt_from_equilibration =
          equilibrium_melt_fraction * equilibrated_volume;

        // Include melt that did not participate in equilibration.
        double total_final_melt =
          std::min(std::max(qp.amount_of_free_melt
                            + new_melt_from_equilibration,
                            0.0),
                   1.0);
        double total_final_solid = 1.0 - total_final_melt;

        // // If the advected melt composition was invalid and equilibration produces
        // // no melt, return that free melt mass to the solid. Otherwise the MAGEMin
        // // phase compositions define the updated state below.
        // if (!qp.liquid_valid && new_melt_from_equilibration < 1e-5)
        //   {
        //     total_final_melt =0;
        //     total_final_solid = 1.0;
        //   }


        // Changing some physical parameters
        //  1. density
        if (melt_out!=nullptr && in.requests_property(MaterialProperties::additional_outputs) && use_MAGEMin_density==1)
          {
            const double fluid_compressibility = melt_compressibility /(1.0 + in.pressure[q] * melt_bulk_modulus_derivative * melt_compressibility);
            melt_out->fluid_densities[q] = qp.density_liquid_new;
            melt_out->fluid_density_gradients[q] = qp.density_liquid_new * qp.density_liquid_new * fluid_compressibility * this->get_gravity_model().gravity_vector(in.position[q]);
            out.densities[q] = qp.density_solid_new;
          }

        if (use_MAGEMin_specific_heat == 1 &&
            in.requests_property(MaterialProperties::specific_heat) &&
            std::isfinite(qp.specific_heat_solid_new) &&
            qp.specific_heat_solid_new > 0.0)
          out.specific_heat[q] = qp.specific_heat_solid_new;

        if (enthalpy_out != nullptr && std::isfinite(qp.enthalpy_of_fusion))
          enthalpy_out->enthalpies_of_fusion[q] = qp.enthalpy_of_fusion;


        // Now changing the reaction terms after melting
        // 1. liquid
        if (total_final_melt < 1e-5)
          {
            // This scenario means that all melt is used up and converted to solid
            for (int i=0; i<qp.oxideLength; i++)
              {
                out.reaction_terms[q][qp.idx_X_melt[i]] = -in.composition[q][qp.idx_X_melt[i]];
              }
            total_final_melt=0;
            total_final_solid=1.0;
          }
        else
          {
            for (int i=0; i<qp.oxideLength; i++)
              {
                double final_val = 0.0;
                // if (qp.liquid_valid)
                //   {
                const double old_liquid_mass = qp.amount_of_free_melt;
                const double new_liquid_mass = new_melt_from_equilibration;
                final_val = std::min(std::max(
                                       (old_liquid_mass * qp.X_melt[i]
                                        + new_liquid_mass * qp.X_liq_magemin[i])
                                       / (old_liquid_mass + new_liquid_mass),
                                       0.0), 1.0);
                // }
                // else
                //   {
                //     final_val = qp.X_liq_magemin[i];
                //   }

                double change = final_val - in.composition[q][qp.idx_X_melt[i]];
                change = std::max(change, -in.composition[q][qp.idx_X_melt[i]]);
                out.reaction_terms[q][qp.idx_X_melt[i]] = change;
              }
          }
        // 2. Solid
        for (int i = 0; i < qp.oxideLength; i++)
          {
            double final_val = std::max(qp.X_solid_magemin[i], 0.0);
            double change    = final_val - in.composition[q][qp.idx_X_solid[i]];
            change = std::max(change, -in.composition[q][qp.idx_X_solid[i]]);
            out.reaction_terms[q][qp.idx_X_solid[i]] = change;
          }
        // AssertThrow(std::abs(total_final_solid + total_final_melt - 1.0) < 1e-8,
        //             ExcMessage("MAGEMin reaction did not conserve solid-plus-melt volume fraction."));
        qp.fraction_of_melt = total_final_melt;
      }

      template<int dim> void meltMagemin<dim>::remove_melt_below_cutoff_depth(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out, unsigned q, MeltOutputs<dim> *melt_out, const double reference_T) const
      {

        if (melt_out != nullptr && in.requests_property(MaterialProperties::additional_outputs) && use_MAGEMin_density==1)
          {
            double temperature_dependence = 1.0;
            if (this->include_adiabatic_heating())
              {
                temperature_dependence -= (in.temperature[q] -
                                           this->get_adiabatic_conditions().temperature(in.position[q])) *
                                          out.thermal_expansion_coefficients[q];

              }
            else
              {
                temperature_dependence -= (in.temperature[q] - reference_T) *
                                          out.thermal_expansion_coefficients[q];
              }

            const double fluid_compressibility =melt_compressibility /(1.0 + in.pressure[q] * melt_bulk_modulus_derivative * melt_compressibility);
            melt_out->fluid_densities[q] = reference_rho_fluid *std::exp(fluid_compressibility * (in.pressure[q] - this->get_surface_pressure())) * temperature_dependence;
            melt_out->fluid_density_gradients[q] = melt_out->fluid_densities[q] * melt_out->fluid_densities[q] * fluid_compressibility * this->get_gravity_model().gravity_vector(in.position[q]);
          }

        // Clear finite-element liquid fields, but retain the normalized dormant
        // composition carried by liquid particles.
        for (int i = 0; i < qp.oxideLength; i++)
          {
            out.reaction_terms[q][qp.idx_X_melt[i]] = -in.composition[q][qp.idx_X_melt[i]];
          }
      }




      template <int dim>
      MAGEMin_hash::MAGEMin_Hash_table &
      meltMagemin<dim>::pick_tier(const quad_point_data &qp) const
      {
        if (surrogateMode == 1 || enable_two_tier == 0)
          return magemin_cache;
        return qp.near_solidus ? magemin_cache : magemin_cache_coarse;
      }

      template <int dim>
      void meltMagemin<dim>::setup_hash_tiers()
      {
        if (surrogateMode == 1)
          cache_context_signature = "solver=neural-network;model=" + surrogate_directory;
        else
          {
            // Bump the solver signature whenever changes can alter exact labels.
            cache_context_signature =
              "solver=MAGEMin-1.9.3-aspect-exact-1;database=" + database_name +
              ";included=";
            for (const std::string &phase : included_magemin_phases)
              cache_context_signature += phase + ",";
            cache_context_signature += ";excluded=";
            for (const std::string &phase : excluded_magemin_phases)
              cache_context_signature += phase + ",";
          }

        // Fine tier
        magemin_cache.set_oxide_names(oxide_names);
        magemin_cache.set_context_signature(cache_context_signature);
        magemin_cache.set_scales(hash_table_P_scale, hash_table_T_scale, hash_table_X_scale);
        if (have_per_oxide_scales)
          magemin_cache.set_oxide_scales(hash_X_scales);
        magemin_cache.resize_hash_table(size_of_hash_table);

        // Coarse tier (12): same machinery, fatter bins
        if (surrogateMode == 0 && enable_two_tier == 1)
          {
            magemin_cache_coarse.set_oxide_names(oxide_names);
            magemin_cache_coarse.set_context_signature(cache_context_signature);
            magemin_cache_coarse.set_scales(hash_table_P_scale / coarse_tier_PT_factor,
                                            hash_table_T_scale / coarse_tier_PT_factor,
                                            hash_table_X_scale / coarse_tier_X_factor);
            if (have_per_oxide_scales)
              {
                std::vector<double> coarse_x(hash_X_scales.size());
                for (unsigned int i = 0; i < hash_X_scales.size(); ++i)
                  coarse_x[i] = hash_X_scales[i] / coarse_tier_X_factor;
                magemin_cache_coarse.set_oxide_scales(coarse_x);
              }
            magemin_cache_coarse.resize_hash_table(size_of_coarse_hash_table);

          }
      }
      template<int dim>
      void meltMagemin<dim>::compute_solidus_state(quad_point_data &qp) const
      {
        // Hydrous Katz solidus, shifted by depletion.
        const double T_solidus_dry =
          A1 + 273.15 + A2 * qp.P_Pa + A3 * qp.P_Pa * qp.P_Pa;

        const double X_H2O_saturated =
          Xie1 * std::pow(std::max(0.0, qp.P_Pa), Lambda)
          + Xie2 * std::max(0.0, qp.P_Pa);
        const double deltaT_saturated =
          K * std::pow(std::max(0.0, X_H2O_saturated), gamma);
        const double bulk_water = water_oxide_index < qp.X_solid.size()
                                  ? qp.X_solid[water_oxide_index]*100 : 0.0;
        const double deltaT =
          std::min(K * std::pow(std::max(0.0, bulk_water), gamma), deltaT_saturated);

        const double drained_depletion =
          std::max(qp.depletion - qp.phi_melt, 0.0);
        const double dT_depl = depletion_solidus_change
                               * drained_depletion;

        qp.T_solidus_K = T_solidus_dry - deltaT + dT_depl;

        const double dT_sol = qp.T_K - qp.T_solidus_K;
        qp.near_solidus = (dT_sol >= -fine_band_below_solidus) &&
                          (dT_sol <=  fine_band_above_solidus);
      }



      /**
      * The function which prepares the magemin point for evaluation. It loads the variables
      * from the present time step, normalizes the drifted compositions, and checks if the
      * point is below the cutoff depth. If it is, it removes the melt below the cutoff depth
      * and returns false. It also computes the solidus state and checks if the point is above
      * solidus and should be evaluated. If not, it returns false. Finally, it calculates the
      * bulk composition and returns true if the point is ready for evaluation.
      */
      template<int dim>
      bool meltMagemin<dim>::prepare_magemin_point(
        const typename Interface<dim>::MaterialModelInputs &in,
        typename Interface<dim>::MaterialModelOutputs &out,
        const unsigned int q,
        const double reference_T,
        quad_point_data &qp_dat,
        double &equilibrium_melt_fraction) const
      {
        const std::shared_ptr<MeltOutputs<dim>> melt_out = out.template get_additional_output_object<MeltOutputs<dim>>();
        const std::shared_ptr<EnthalpyOutputs<dim>> enthalpy_out =
                                                   out.template get_additional_output_object<EnthalpyOutputs<dim>>();
        equilibrium_melt_fraction = 0.0;

        // Points that do not require an equilibrium calculation have no
        // melting reaction here. Initializing the additional output also
        // prevents signaling NaNs in latent_heat_melt.cc for those points.
        if (enthalpy_out != nullptr)
          enthalpy_out->enthalpies_of_fusion[q] = 0.0;

        qp_dat.database = const_cast<char *>(database_name.c_str());

        // loading the variables from the present time step
        qp_dat.density_solid_current =
          std::isfinite(out.densities[q]) && out.densities[q] > 0.0
          ? out.densities[q] : 3300.0;
        qp_dat.density_liquid_current = reference_rho_fluid;
        load_variables_from_present_time_step(qp_dat, in, q);
        // normalize_drifted_compositions(qp_dat);
        // normalize_solid_compositions(qp_dat);
        // normalize_liquid_compositions(qp_dat);
        normalize_compositions(qp_dat,in, out, q);



        // if (qp_dat.depth_m > cutOff_depth)
        //   {
        //     remove_melt_below_cutoff_depth(qp_dat, in, out,q, melt_out.get(), reference_T);
        //     return false;
        //   }

        // computing solidus values
        compute_solidus_state(qp_dat);

        if (qp_dat.phi_melt > 1e-5 && qp_dat.T_K < fractional_crystallization_start)
          {
            const double normal_porosity_value =
              freeze_melt_below_certain_tenperature(qp_dat, in, out, q);
            equilibrium_melt_fraction = normal_porosity_value;
            return false;
          }

        // spurious melt fraction at low temperatures, which is not physical. This is a temporary fix to avoid this issue.
        if (qp_dat.phi_melt > 1e-5 && (qp_dat.T_K < qp_dat.T_solidus_K - 400))
          {
            for (int i = 0; i < qp_dat.oxideLength; ++i)
              {
                // const double X_new =(qp_dat.X_solid[i] * (1 - qp_dat.phi_melt) + qp_dat.X_melt[i] * qp_dat.phi_melt);
                out.reaction_terms[q][qp_dat.idx_X_solid[i]] = 0.0; //X_new - std::max(in.composition[q][qp_dat.idx_X_solid[i]], 0.0);
                out.reaction_terms[q][qp_dat.idx_X_melt[i]]  = -in.composition[q][qp_dat.idx_X_melt[i]];
              }
            equilibrium_melt_fraction = 0.0;
            return false;
          }


        // if (qp_dat.phi_solid<1e-4)
        //   {
        //     handle_cases_when_there_is_no_solid_left(qp_dat, out, q);
        //     equilibrium_melt_fraction = qp_dat.phi_melt;
        //     return false;
        //   }

        // Check if this point is above solidus and should we evaluate it?
        if (skip_magemin_evaluation(qp_dat, in, out, q))
          {
            return false;
          }

        calculate_bulk_composition(qp_dat);
        return true;
      }



      template<int dim>
      double meltMagemin<dim>::finalize_magemin_point(
        const typename Interface<dim>::MaterialModelInputs &in,
        typename Interface<dim>::MaterialModelOutputs &out,
        const unsigned int q,
        quad_point_data &qp_dat) const
      {
        const std::shared_ptr<MeltOutputs<dim>> melt_out =
                                               out.template get_additional_output_object<MeltOutputs<dim>>();
        const std::shared_ptr<EnthalpyOutputs<dim>> enthalpy_out =
                                                   out.template get_additional_output_object<EnthalpyOutputs<dim>>();

        if (qp_dat.run_type == 20)
          {
            for (int i = 0; i < qp_dat.oxideLength; ++i)
              {
                out.reaction_terms[q][qp_dat.idx_X_solid[i]] = 0.0;
                out.reaction_terms[q][qp_dat.idx_X_melt[i]] = 0.0;
              }
            return qp_dat.phi_melt;
          }

        write_to_aspect_data_structures(qp_dat, in,out, q, melt_out.get(), enthalpy_out.get());
        return qp_dat.fraction_of_melt;
      }



      template <int dim>
      void meltMagemin<dim>::calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                                             typename Interface<dim>::MaterialModelOutputs &out, const double reference_T) const
      {

        const std::shared_ptr<ReactionRateOutputs<dim>> reaction_rate_out = out.template get_additional_output_object<ReactionRateOutputs<dim>>();

        // if melt transport is not included, then we don't need to compute reaction rates
        if (!this->include_melt_transport())
          {
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                if (reaction_rate_out != nullptr)
                  reaction_rate_out->reaction_rates[i][c] = 0.0;
            return;
          }

        const unsigned int n_points = in.n_evaluation_points();

        // creatue a vector of quad_point_data to store the data for all the points
        std::vector<quad_point_data> quad_points;
        quad_points.reserve(n_points);
        for (unsigned int i = 0; i < n_points; ++i)
          quad_points.emplace_back(oxide_names.size());
        std::vector<double> equilibrium_melt_fractions(n_points, 0.0);
        std::vector<unsigned int> points_requiring_thermodynamics;
        points_requiring_thermodynamics.reserve(n_points);

        for (unsigned int i = 0; i < n_points; ++i)
          if (prepare_magemin_point(in,
                                    out,
                                    i,
                                    reference_T,
                                    quad_points[i],
                                    equilibrium_melt_fractions[i]))
            points_requiring_thermodynamics.push_back(i);

        resolve_magemin_batch(quad_points, points_requiring_thermodynamics);

        for (const unsigned int point_index : points_requiring_thermodynamics)
          equilibrium_melt_fractions[point_index] =
            finalize_magemin_point(in, out, point_index, quad_points[point_index]);

        const unsigned int porosity_idx =
          this->introspection().compositional_index_for_name("porosity");
        for (unsigned int i = 0; i < n_points; ++i)
          {
            quad_point_data &qp_dat = quad_points[i];

            // The porosity from aspect is in volume fraction
            const double old_volume_fraction_porosity = std::clamp(in.composition[i][porosity_idx], 0.0, 1.0);
            const double old_mass_fraction_porosity = compute_mass_fraction(old_volume_fraction_porosity, qp_dat.density_liquid_current, qp_dat.bulk_density_current);


            const double new_mass_fraction_porosity = std::clamp(equilibrium_melt_fractions[i], 0.0, 1.0);
            const double liquid_density_for_conversion =
              use_MAGEMin_density==1 && std::isfinite(qp_dat.density_liquid_new) &&
              qp_dat.density_liquid_new > 0.0
              ? qp_dat.density_liquid_new
              : qp_dat.density_liquid_current;
            const double new_volume_fraction_porosity = std::clamp(new_mass_fraction_porosity * qp_dat.bulk_density_current / liquid_density_for_conversion, 0.0, 1.0);

            const double volume_fraction_porosity_change =std::clamp(new_volume_fraction_porosity - old_volume_fraction_porosity, -old_volume_fraction_porosity, 1.0 - old_volume_fraction_porosity);
            const double mass_fraction_porosity_change = std::clamp(new_mass_fraction_porosity - old_mass_fraction_porosity, -old_mass_fraction_porosity, 1.0 - old_mass_fraction_porosity);

            // Match ASPECT's operator-split reaction convention: timestep zero
            // establishes the initial state but does not apply a reaction.
            if (this->get_timestep_number() == 0)
              {
                for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                  {
                    out.reaction_terms[i][c] = 0.0;
                    if (reaction_rate_out != nullptr &&
                        in.requests_property(MaterialProperties::reaction_rates))
                      reaction_rate_out->reaction_rates[i][c] = 0.0;
                  }
                continue;
              }

            double depletion_change = 0.0;
            double crust_change = 0.0;

            if (qp_dat.idx_depletion !=
                std::numeric_limits<unsigned int>::max())
              {
                const double generated_melt = std::max(new_mass_fraction_porosity - old_mass_fraction_porosity, 0.0);
                const double old_remaining_solid = 1.0 - old_mass_fraction_porosity;
                if (old_remaining_solid > 1e-12)
                  depletion_change = generated_melt * (1.0 - qp_dat.depletion) / old_remaining_solid;

                depletion_change = std::clamp(depletion_change, 0.0, 1.0 - qp_dat.depletion);
                out.reaction_terms[i][qp_dat.idx_depletion] = 0.0;
              }

            if (qp_dat.idx_phi_crust != std::numeric_limits<unsigned int>::max() &&
                std::isfinite(qp_dat.phi_crust))
              {
                const double transported_crust = in.composition[i][qp_dat.idx_phi_crust];
                crust_change = std::clamp(qp_dat.phi_crust - transported_crust, -transported_crust, 1.0 - transported_crust);
              }


            const bool use_operator_split_rates =
              reaction_rate_out != nullptr &&
              in.requests_property(MaterialProperties::reaction_rates);

            for (unsigned int c = 0; c < in.composition[i].size(); ++c)
              if (use_operator_split_rates)
                {
                  if (c == porosity_idx && this->get_timestep_number() > 0)
                    reaction_rate_out->reaction_rates[i][c] =
                      volume_fraction_porosity_change / melting_time_scale;
                  else if (c == qp_dat.idx_depletion &&
                           this->get_timestep_number() > 0)
                    reaction_rate_out->reaction_rates[i][c] =
                      depletion_change / melting_time_scale;
                  else if (c == qp_dat.idx_phi_crust &&
                           this->get_timestep_number() > 0)
                    reaction_rate_out->reaction_rates[i][c] =
                      crust_change / melting_time_scale;
                  else
                    reaction_rate_out->reaction_rates[i][c] = 0.0;
                }
          }
      }

      template <int dim>
      void meltMagemin<dim>::calculate_fluid_outputs( const typename Interface<dim>::MaterialModelInputs &in,
                                                      typename Interface<dim>::MaterialModelOutputs &out,
                                                      const double reference_T) const
      {
        const std::shared_ptr<MeltOutputs<dim>> melt_out = out.template get_additional_output_object<MeltOutputs<dim>>();
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        if (melt_out != nullptr && in.requests_property(MaterialProperties::additional_outputs))
          {
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                double porosity = std::max(in.composition[i][porosity_idx], 0.0);

                melt_out->fluid_viscosities[i] = viscosity_fluid;
                melt_out->permeabilities[i] = reference_permeability *
                                              Utilities::fixed_power<3>(porosity) *
                                              Utilities::fixed_power<2>(1.0 - porosity);

                // Calculate dependence of fluid density on temperature if adiabatic
                // heating in switched on first, calculate temperature dependence of
                // density
                double temperature_dependence = 1.0;
                if (this->include_adiabatic_heating())
                  {
                    // temperature dependence is 1 - alpha * (T - T(adiabatic))
                    temperature_dependence -=
                      (in.temperature[i] -
                       this->get_adiabatic_conditions().temperature(in.position[i])) *
                      out.thermal_expansion_coefficients[i];
                  }
                else
                  temperature_dependence -= (in.temperature[i] - reference_T) *
                                            out.thermal_expansion_coefficients[i];

                // the fluid compressibility includes two parts, a constant
                // compressibility, and a pressure-dependent one this is a simplified
                // formulation, experimental data are often fit to the Birch-Murnaghan
                // equation of state

                // Provide a valid density for points where thermodynamics is
                // skipped. An NN/hash result overrides it later for evaluated
                // points when predicted density is enabled.
                const double fluid_compressibility = melt_compressibility /
                                                     (1.0 + in.pressure[i] * melt_bulk_modulus_derivative * melt_compressibility);
                if (use_MAGEMin_density==0 || !std::isfinite(melt_out->fluid_densities[i]) || melt_out->fluid_densities[i] <= 0.0)
                  {
                    melt_out->fluid_densities[i] = reference_rho_fluid *
                                                   std::exp(fluid_compressibility *
                                                            (in.pressure[i] - this->get_surface_pressure())) *
                                                   temperature_dependence;
                  }
                melt_out->fluid_density_gradients[i] =
                  melt_out->fluid_densities[i] * melt_out->fluid_densities[i] *
                  fluid_compressibility *
                  this->get_gravity_model().gravity_vector(in.position[i]);
                const double phi_0 = 0.05;
                porosity = std::max(std::min(porosity, 0.995), 1e-4);
                melt_out->compaction_viscosities[i] = xi_0 * phi_0 / porosity;

                double visc_temperature_dependence = 1.0;
                if (this->include_adiabatic_heating())
                  {
                    const double delta_temp = in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]);
                    visc_temperature_dependence = std::max( std::min(std::exp(-thermal_bulk_viscosity_exponent * delta_temp / this->get_adiabatic_conditions().temperature( in.position[i])), 1e4), 1e-4);
                  }
                else
                  {
                    const double delta_temp = in.temperature[i] - reference_T;
                    const double T_dependence =
                      (thermal_bulk_viscosity_exponent == 0.0
                       ? 0.0
                       : thermal_bulk_viscosity_exponent * delta_temp / reference_T);
                    visc_temperature_dependence = std::max(std::min(std::exp(-T_dependence), 1e4), 1e-4);
                  }
                melt_out->compaction_viscosities[i] *= visc_temperature_dependence;
              }
          }

        if (this->include_melt_transport() && in.requests_property(MaterialProperties::viscosity))
          {
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx], 0.0));
                out.viscosities[i] *= std::exp(-alpha_phi * porosity);
              }
          }
      }

      template <int dim>
      void meltMagemin<dim>::declare_parameters(ParameterHandler &prm)
      {

        prm.declare_entry("Use surrogate instead of MAGEMin", "0", Patterns::Integer(0, 1),
                          "Whether to use the independent neural-network route. It starts with "
                          "an empty in-memory hash, evaluates all new states with one batched "
                          "neural-network call, and caches accepted predictions for reuse. This "
                          "route does not initialize or call MAGEMin and does not use the KD-tree.");
        prm.declare_entry("Surrogate Dir", "$ASPECT_SOURCE_DIR/data/material-model/xgBoost/",
                          Patterns::DirectoryName(), "Directory containing the neural-network models.");
        prm.declare_entry("MAGEMin database", "ig", Patterns::Anything(),
                          "MAGEMin thermodynamic database name. The native database metadata "
                          "determines the required oxide order and available phases.");
        prm.declare_entry("Solid oxide field names", "", Patterns::Anything(),
                          "Optional comma-separated ASPECT compositional field names mapped to "
                          "the native MAGEMin oxide order. Empty uses the native oxide names.");
        prm.declare_entry("Liquid oxide field names", "", Patterns::Anything(),
                          "Optional comma-separated ASPECT compositional field names mapped to "
                          "the native MAGEMin oxide order. Empty appends 'l' to each native name.");

        prm.declare_entry("Reference melt density", "2500.", Patterns::Double(0.),
                          "Reference density of the melt/fluid$\\rho_{f,0}$. "
                          "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
        prm.declare_entry(
          "Reference bulk viscosity", "1e22", Patterns::Double(0.),
          "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
          "This viscosity may be modified by both temperature and porosity "
          "dependencies. Units: \\si{\\pascal\\second}.");
        prm.declare_entry("Reference melt viscosity", "10.", Patterns::Double(0.),
                          "The value of the constant melt viscosity $\\eta_f$. "
                          "Units: \\si{\\pascal\\second}.");
        prm.declare_entry(
          "Exponential melt weakening factor", "27.", Patterns::Double(0.),
          "The porosity dependence of the viscosity. Units: dimensionless.");
        prm.declare_entry(
          "Thermal bulk viscosity exponent", "0.0", Patterns::Double(0.),
          "The temperature dependence of the bulk viscosity. Dimensionless "
          "exponent. "
          "See the general documentation "
          "of this model for a formula that states the dependence of the "
          "viscosity on this factor, which is called $\\beta$ there.");
        prm.declare_entry("Melt compressibility", "0.0", Patterns::Double(0.),
                          "The value of the compressibility of the melt. "
                          "Units: \\si{\\per\\pascal}.");
        prm.declare_entry("Melt bulk modulus derivative", "0.0", Patterns::Double(0.),
                          "The value of the pressure derivative of the melt bulk "
                          "modulus. "
                          "Units: None.");
        prm.declare_entry("Melting time scale for operator splitting", "1e3",
                          Patterns::Double(0.),
                          "Because the operator splitting scheme is used, the "
                          "porosity field can not "
                          "be set to a new equilibrium melt fraction instantly, but "
                          "the model has to "
                          "provide a melting time scale instead. This time scale "
                          "defines how fast melting "
                          "happens, or more specifically, the parameter defines the "
                          "time after which "
                          "the deviation of the porosity from the equilibrium melt "
                          "fraction will be "
                          "reduced to a fraction of $1/e$. So if the melting time "
                          "scale is small compared "
                          "to the time step size, the reaction will be so fast that "
                          "the porosity is very "
                          "close to the equilibrium melt fraction after reactions "
                          "are computed. Conversely, "
                          "if the melting time scale is large compared to the time "
                          "step size, almost no "
                          "melting and freezing will occur."
                          "\n\n"
                          "Also note that the melting time scale has to be larger "
                          "than or equal to the reaction "
                          "time step used in the operator splitting scheme, "
                          "otherwise reactions can not be "
                          "computed. "
                          "Units: yr or s, depending on the ``Use years in output "
                          "instead of seconds'' parameter.");
        prm.declare_entry(
          "Depletion solidus change", "200.0", Patterns::Double(0.),
          "The solidus temperature change for a depletion of 100\\%. For positive "
          "values, the solidus gets increased for a positive depletion field "
          "(depletion) and lowered for a negative depletion field (enrichment). "
          "Scaling with depletion is linear. "
          "Units: \\si{\\kelvin}.");
        prm.declare_entry("Reference permeability", "1e-8", Patterns::Double(0.),
                          "Reference permeability of the solid host rock."
                          "Units: \\si{\\meter\\squared}.");
        prm.declare_entry(
          "Cutoff depth for magemin", "200e3", Patterns::Double(0.),
          "Cutoff depth: Melting will only be considered above this depth"
          "Units: \\si{\\meter}.");

        prm.declare_entry("Size of hash table", "65536", Patterns::Integer(1),
                          "Capacity of the fine cache. Must be a power of two.");

        prm.declare_entry("Save hash table as binary files", "0", Patterns::Integer(0, 1),
                          "Whether to save and reload per-rank MAGEMin cache files at checkpoints. "
                          "Ignored by the independent neural-network route.");
        prm.declare_entry("Store timestep in MAGEMin binary files", "1",
                          Patterns::Integer(0, 1),
                          "Record the ASPECT solve timestep for each exact cache entry. "
                          "Synthetic preseed states and disabled recording use distinct "
                          "sentinels in the NN extractor.");
        prm.declare_entry("Where in the middle earth shall I save the binary files", "./",
                          Patterns::DirectoryName(), "Directory for per-rank binary cache files.");
        prm.declare_entry("Hash table P scale", "50", Patterns::Double(0.),
                          "Reciprocal pressure-bin width in 1/GPa.");
        prm.declare_entry("Hash table T scale", "1", Patterns::Double(0.),
                          "Reciprocal temperature-bin width in 1/C.");
        prm.declare_entry("Hash table X scale", "200", Patterns::Double(0.),
                          "Reciprocal composition-bin width for all oxides unless overridden.");



        prm.declare_entry(
          "Number of points before kd tree starts", "500", Patterns::Integer(1),
          "Minimum number of cached equilibria before kd-tree queries begin.");
        prm.declare_entry(
          "Minimum number of neighbours to trust kd tree", "3", Patterns::Integer(1),
          "Minimum number of nearby cached equilibria required to validate a "
          "nearest-neighbour result.");
        prm.declare_entry(
          "Maximum melt range in nearest neighbours", "0.001", Patterns::Double(0.),
          "Maximum accepted melt-fraction range across neighbouring equilibria.");
        prm.declare_entry(
          "Use KD Tree", "1", Patterns::Integer(0, 1),
          "Whether to reuse nearby equilibria through the nanoflann kd-tree. "
          "A value of 1 enables the kd-tree and 0 disables it.");

        prm.declare_entry("Use MAGEMin density", "0",
                          Patterns::Integer(0, 1),
                          "Whether to use the solid and liquid densities returned by MAGEMin.");
        prm.declare_entry("Use MAGEMin specific heat", "0",
                          Patterns::Integer(0, 1),
                          "1: use the mass-weighted solid specific heat returned by "
                          "MAGEMin in J/(kg K); 0: retain the host material model value.");
        prm.declare_entry("Excluded MAGEMin phases", "",
                          Patterns::Anything(),
                          "Comma-separated MAGEMin phase names to remove from every "
                          "equilibrium calculation. Empty keeps the full database. "
                          "Use only deliberate physical exclusions; automatic removal "
                          "from a previous assemblage can suppress phase nucleation.");
        prm.declare_entry("Included MAGEMin phases", "",
                          Patterns::Anything(),
                          "Comma-separated MAGEMin phase names to retain while removing "
                          "all other phases. Empty keeps the full database. This option "
                          "is mutually exclusive with Excluded MAGEMin phases and should "
                          "only be used for a physically justified restricted system.");
        prm.declare_entry("Max mineral range in nearest neighbours", "0.10", Patterns::Double(0.),
                          "Max range of any mineral proportion among neighbors (0.10 = 10%)");
        prm.declare_entry("Max oxide deviation from neighbours", "0.01", Patterns::Double(0.),
                          "Maximum range of a mineral's oxide mass fraction among "
                          "neighbours before kd-tree reuse is rejected.");
        prm.declare_entry("Preseed KD tree", "0",
                          Patterns::Integer(0, 1),
                          "1: Pre-seed the KD-tree at initialization by running MAGEMin "
                          "on a P-T grid covering the local domain. Maps phase boundaries "
                          "before the simulation starts. 0: disabled (default).");
        prm.declare_entry("Preseed P spacing", "0.5",
                          Patterns::Double(0.),
                          "Pressure spacing in GPa for the pre-seeding grid.");
        prm.declare_entry("Preseed T spacing", "5.0",
                          Patterns::Double(0.),
                          "Temperature spacing in Celsius for the pre-seeding grid.");
        prm.declare_entry("Preseed T min", "800.0",
                          Patterns::Double(0.),
                          "Minimum temperature in Celsius for the pre-seeding grid. "
                          "Set this near the solidus to avoid wasting calls on cold rock.");
        prm.declare_entry("Preseed T max", "1600.0",
                          Patterns::Double(0.),
                          "Maximum temperature in Celsius for the pre-seeding grid.");
        prm.declare_entry("Preseed bulk composition",
                          "0.4486, 0.03512, 0.0307, 0.39524, 0.08202, 0.000183, 0.003, 0.00155, 0.000298, 0.00321, 0.0",
                          Patterns::List(Patterns::Double(0.)),
                          "Initial bulk composition in the native oxide order of the selected "
                          "MAGEMin database for pre-seeding. "
                          "The composition is normalized before it is passed to MAGEMin.");
        prm.declare_entry("Times hash bin size for kd-tree oxides", "1", Patterns::Integer(1),
                          "KD-tree composition search half-width in hash bins.");
        prm.declare_entry("Times hash bin size for kd-tree pressure", "1", Patterns::Integer(1),
                          "KD-tree pressure search half-width in hash bins.");
        prm.declare_entry("Times hash bin size for kd-tree temperature", "1", Patterns::Integer(1),
                          "KD-tree temperature search half-width in hash bins.");

        prm.declare_entry("Snap queries to hash bin centers", "1", Patterns::Integer(0, 1),
                          "1: quantize (P,T,X) to bin centers before lookup AND before "
                          "calling MAGEMin, so cached results are solved exactly at their "
                          "bin center (max error = half bin, symmetric). 0: legacy behavior "
                          "(cached result corresponds to an arbitrary point in the bin; "
                          "max error = full bin width).");

        prm.declare_entry("Hash table X scales per oxide", "",
                          Patterns::Anything(),
                          "Optional comma-separated per-oxide scales in the native order of "
                          "the selected database (scale = 1/bin width in weight fraction). "
                          "Empty: use the global 'Hash table X scale' for all oxides.");

        prm.declare_entry("Enable two tier hash table", "1", Patterns::Integer(0, 1),
                          "1: fine bins near the Katz solidus, coarse bins elsewhere.");
        prm.declare_entry("Size of coarse hash table", "32768", Patterns::Integer(1),
                          "Capacity of the coarse tier (power of 2).");
        prm.declare_entry("Fine tier band below solidus", "30.0", Patterns::Double(0.),
                          "Kelvin below the Katz solidus that still uses fine bins.");
        prm.declare_entry("Fine tier band above solidus", "100.0", Patterns::Double(0.),
                          "Kelvin above the Katz solidus that still uses fine bins. "
                          "dF/dT is steepest just above the solidus (beta power law), "
                          "so keep this generous.");
        prm.declare_entry("Coarse tier PT factor", "2.0", Patterns::Double(1.),
                          "Coarse tier P,T scales = fine scales / this factor.");
        prm.declare_entry("Coarse tier X factor", "2.0", Patterns::Double(1.),
                          "Coarse tier composition scales = fine scales / this factor.");


        prm.declare_entry("Melt retention threshold", "0.01", Patterns::Double(0., 1.),
                          "Retained porosity entering the MAGEMin bulk; melt above this advects "
                          "and never re-equilibrates. Cite McKenzie 1984; von Bargen & Waff 1986.");
        prm.declare_entry("Basalt gate temperature", "600.0", Patterns::Double(0.),
                          "Coarse T gate (C) below which melt freezes to crust.");
      }

      template <int dim>
      void meltMagemin<dim>::parse_parameters(ParameterHandler &prm)
      {

        surrogate_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get("Surrogate Dir"));
        database_name = prm.get("MAGEMin database");
        configured_solid_oxide_field_names =
          Utilities::split_string_list(prm.get("Solid oxide field names"));
        configured_liquid_oxide_field_names =
          Utilities::split_string_list(prm.get("Liquid oxide field names"));
        AssertThrow(!database_name.empty(),
                    ExcMessage("'MAGEMin database' must not be empty."));

        reference_rho_fluid = prm.get_double("Reference melt density");
        xi_0 = prm.get_double("Reference bulk viscosity");
        viscosity_fluid = prm.get_double("Reference melt viscosity");
        thermal_bulk_viscosity_exponent = prm.get_double("Thermal bulk viscosity exponent");
        alpha_phi = prm.get_double("Exponential melt weakening factor");
        melt_compressibility = prm.get_double("Melt compressibility");
        melting_time_scale =prm.get_double("Melting time scale for operator splitting");
        melt_bulk_modulus_derivative = prm.get_double("Melt bulk modulus derivative");
        depletion_solidus_change = prm.get_double("Depletion solidus change");
        reference_permeability = prm.get_double("Reference permeability");
        cutOff_depth = prm.get_double("Cutoff depth for magemin");
        surrogateMode = prm.get_integer("Use surrogate instead of MAGEMin");

        // MAGEMin physical parameters.
        use_MAGEMin_density = prm.get_integer("Use MAGEMin density");
        use_MAGEMin_specific_heat = prm.get_integer("Use MAGEMin specific heat");
        included_magemin_phases =
          Utilities::split_string_list(prm.get("Included MAGEMin phases"));
        excluded_magemin_phases =
          Utilities::split_string_list(prm.get("Excluded MAGEMin phases"));
        AssertThrow(included_magemin_phases.empty() || excluded_magemin_phases.empty(),
                    ExcMessage("'Included MAGEMin phases' and 'Excluded MAGEMin phases' "
                               "cannot both be set."));
        if (surrogateMode == 0)
          {
            wrap.setIncludedPhases(included_magemin_phases);
            wrap.setExcludedPhases(excluded_magemin_phases);
          }
        size_of_hash_table = prm.get_integer("Size of hash table");
        save_hash_table = prm.get_integer("Save hash table as binary files");
        store_cache_timestep =
          (prm.get_integer("Store timestep in MAGEMin binary files") == 1);
        hash_storage_location = Utilities::expand_ASPECT_SOURCE_DIR(prm.get("Where in the middle earth shall I save the binary files"));
        hash_table_P_scale=prm.get_double("Hash table P scale");
        hash_table_T_scale=prm.get_double("Hash table T scale");
        hash_table_X_scale=prm.get_double("Hash table X scale");
        AssertThrow(hash_table_P_scale > 0.0 &&
                    hash_table_T_scale > 0.0 &&
                    hash_table_X_scale > 0.0,
                    ExcMessage("Hash table P, T, and X scales must be positive."));
        num_required_points_for_kd_tree = prm.get_integer("Number of points before kd tree starts");
        minimum_neighbours_to_trust_kd_tree = prm.get_integer("Minimum number of neighbours to trust kd tree");
        useKd_tree = prm.get_integer("Use KD Tree");
        max_melt_range  = prm.get_double("Maximum melt range in nearest neighbours");
        max_mineral_range = prm.get_double("Max mineral range in nearest neighbours");
        max_oxide_deviation = prm.get_double("Max oxide deviation from neighbours");


        preseed_kdtree_enabled = prm.get_integer("Preseed KD tree");
        preseed_P_spacing = prm.get_double("Preseed P spacing");
        preseed_T_spacing = prm.get_double("Preseed T spacing");
        preseed_T_min = prm.get_double("Preseed T min");
        preseed_T_max = prm.get_double("Preseed T max");
        AssertThrow(preseed_P_spacing > 0.0 && preseed_T_spacing > 0.0,
                    ExcMessage("Preseed P and T spacing must be positive."));
        kdtree_oxide_tol_bins = prm.get_integer("Times hash bin size for kd-tree oxides");
        kdtree_pressure_tol_bins = prm.get_integer("Times hash bin size for kd-tree pressure");
        kdtree_temperature_tol_bins = prm.get_integer("Times hash bin size for kd-tree temperature");



        snap_queries_to_bins      = prm.get_integer("Snap queries to hash bin centers");
        enable_two_tier           = prm.get_integer("Enable two tier hash table");
        size_of_coarse_hash_table = prm.get_integer("Size of coarse hash table");
        fine_band_below_solidus   = prm.get_double("Fine tier band below solidus");
        fine_band_above_solidus   = prm.get_double("Fine tier band above solidus");
        coarse_tier_PT_factor     = prm.get_double("Coarse tier PT factor");
        coarse_tier_X_factor      = prm.get_double("Coarse tier X factor");

        melt_retention_threshold   = prm.get_double("Melt retention threshold");
        fractional_crystallization_start            = prm.get_double("Basalt gate temperature") + 273.15;

        have_per_oxide_scales = false;
        const std::string xs = prm.get("Hash table X scales per oxide");
        if (!xs.empty())
          {
            const std::vector<std::string> xv = Utilities::split_string_list(xs);
            hash_X_scales.resize(xv.size());
            for (unsigned int i = 0; i < xv.size(); ++i)
              {
                hash_X_scales[i] = Utilities::string_to_double(xv[i]);
                AssertThrow(hash_X_scales[i] > 0.0,
                            ExcMessage("Every per-oxide hash scale must be positive."));
              }
            have_per_oxide_scales = true;
          }

        const std::vector<std::string> comp_strs = Utilities::split_string_list(prm.get("Preseed bulk composition"));
        preseed_bulk_composition.resize(comp_strs.size());
        for (size_t i = 0; i < comp_strs.size(); ++i)
          preseed_bulk_composition[i] = Utilities::string_to_double(comp_strs[i]);

        AssertThrow(preseed_T_min <= preseed_T_max,
                    ExcMessage("'Preseed T min' must not exceed 'Preseed T max'."));

        if (this->convert_output_to_years() == true)
          melting_time_scale *= year_in_seconds;

        AssertThrow(melting_time_scale > 0,
                    ExcMessage("The Melting time scale for operator splitting must "
                               "be larger than 0!"));

        if (this->get_parameters().reaction_solver_type ==
            Parameters<dim>::ReactionSolverType::fixed_step)
          {
            AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                        ExcMessage("The reaction time step " +
                                   Utilities::to_string(
                                     this->get_parameters().reaction_time_step) +
                                   " in the operator splitting scheme is too large to "
                                   "compute melting rates! "
                                   "You have to choose it in such a way that it is "
                                   "smaller than the 'Melting time scale for "
                                   "operator splitting' chosen in the material model, "
                                   "which is currently " +
                                   Utilities::to_string(melting_time_scale) + "."));

          }
      }

    } // namespace ReactionModel
  } // namespace MaterialModel
} // namespace aspect

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim)                                                       \
  namespace ReactionModel {                                                    \
    template class meltMagemin<dim>;                                           \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  } // namespace MaterialModel
} // namespace aspect












