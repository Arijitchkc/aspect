/**
 * ASPECT reaction model coupled to MAGEMin.
 */

#ifndef _aspect_material_model_reaction_model_meltMagemin_h
#define _aspect_material_model_reaction_model_meltMagemin_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/reaction_model/magemin_hash.h>
#include <aspect/material_model/reaction_model/magemin_point_cloud.h>
#include <aspect/material_model/reaction_model/padawan.h>
#include <aspect/melt.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/simulator_access.h>

#include <nanoflann.hpp>

#include <deal.II/base/config.h>

#include <chrono>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim> class meltMagemin : public ::aspect::SimulatorAccess<dim>
      {
        public:
          meltMagemin();
          mutable stableAssemblage sAssemblage;
          mutable MAGEMin_wrapper wrap;

          struct quad_point_data
          {
            char *database = nullptr;
            int oxideLength = 0;
            double P_Pa = 0.0, P_kbar = 0.0, T_K = 0.0, T_C = 0.0;
            double depth_m = 0.0;
            double phi_solid = 0.0;
            double phi_melt = 0.0;
            double depletion = 0.0;
            std::vector<double> X_solid, X_melt;
            double sum_solid = 0.0, sum_liquid = 0.0;
            unsigned int idx_phi_porosity = std::numeric_limits<unsigned int>::max();
            unsigned int idx_depletion = std::numeric_limits<unsigned int>::max();
            std::vector<unsigned int> idx_X_melt, idx_X_solid;

            double amount_of_melt_reacting_with_solid = 0.0;
            double amount_of_free_melt = 0.0;
            double fraction_of_melt = 0.0;
            std::vector<double> X_liq_magemin, X_solid_magemin;
            std::vector<double> X_bulk_passed_to_magemin_individual;
            std::vector<std::vector<double>> X_bulk_passed_to_magemin;
            double density_solid_current = 0.0, density_liquid_current = 0.0;
            double density_liquid_new = 0.0, density_solid_new = 0.0;
            double bulk_density_current=0.0;
            double bulk_density_new=0.0;
            double specific_heat_solid_new = 0.0;
            double enthalpy_of_fusion = 0.0;
            bool liquid_valid = false;
            int run_type = 0; // 5=nn; 10=cache; 12=kd-tree; 15=MAGEMin; 20=NN reject
            mutable stableAssemblage assemblage;
            const stableAssemblage *active=nullptr;
            double T_solidus_K = 0.0;
            bool near_solidus = true;
            double phi_crust = 0.0;
            unsigned int idx_phi_crust = std::numeric_limits<unsigned int>::max();

            explicit quad_point_data(const std::size_t n_oxides)
              : X_solid(n_oxides,0.0),
                X_melt(n_oxides,0.0),
                idx_X_melt(n_oxides, 0u),
                idx_X_solid(n_oxides, 0u),
                X_liq_magemin(n_oxides, 0.0),
                X_solid_magemin(n_oxides, 0.0),
                X_bulk_passed_to_magemin_individual(n_oxides, 0.0) {}
          };
          /**
          * Declare the parameters this function takes through input files.
          */
          static void declare_parameters(ParameterHandler &prm);

          /**
          * Read the parameters from the parameter file.
          */
          void parse_parameters(ParameterHandler &prm);

          double compute_enthalpy_of_fusion(const std::vector<std::string> &names,
                                            const std::vector<double> &props,
                                            const std::vector<double> &entropies,
                                            const std::vector<double> &densities,
                                            const std::vector<double> &volumes,
                                            const double temperature) const;

          bool prepare_magemin_point(const typename Interface<dim>::MaterialModelInputs &in,
                                     typename Interface<dim>::MaterialModelOutputs &out,
                                     unsigned int q,
                                     const double reference_T,
                                     quad_point_data &qp_dat,
                                     double &equilibrium_melt_fraction) const;

          double finalize_magemin_point(const typename Interface<dim>::MaterialModelInputs &in,
                                        typename Interface<dim>::MaterialModelOutputs &out,
                                        unsigned int q,
                                        quad_point_data &qp_dat) const;
          double guess_MeltFraction(double pressure, double temperature, double bulk_water) const;

          void calculate_reaction_rate_outputs(
            const typename Interface<dim>::MaterialModelInputs &in,
            typename Interface<dim>::MaterialModelOutputs &out, const double reference_T) const;

          void calculate_fluid_outputs(
            const typename Interface<dim>::MaterialModelInputs &in,
            typename Interface<dim>::MaterialModelOutputs &out,
            const double reference_T) const;

          void initializeNewandModern();
          void connect_to_signals();
          void preseed_kdtree();
          void create_additional_named_outputs(typename aspect::MaterialModel::Interface<dim>::MaterialModelOutputs &out) const;
          int size_of_hash_table;
          int num_required_points_for_kd_tree;
          int minimum_neighbours_to_trust_kd_tree;
          double max_melt_range;
          // Pre-seeding parameters
          int preseed_kdtree_enabled;
          double preseed_P_spacing;       // GPa
          double preseed_T_spacing;       // Celsius
          double preseed_T_min;           // Celsius
          double preseed_T_max;           // Celsius
          std::vector<double> preseed_bulk_composition;
          std::vector<std::string> included_magemin_phases;
          std::vector<std::string> excluded_magemin_phases;

          // Flag for pre-seeding.
          bool preseed_done;

        private:
          mutable padawan nn;
          mutable MAGEMin_hash::MAGEMin_Hash_table magemin_cache;
          mutable size_t total_points_evaluated;
          mutable size_t total_actual_magemin_calls;


          // ---- hashing variables
          mutable MAGEMin_hash::MAGEMin_Hash_table magemin_cache_coarse;
          int    enable_two_tier;
          int    size_of_coarse_hash_table;
          double fine_band_below_solidus;   // K below Katz solidus that still uses fine bins
          double fine_band_above_solidus;   // K above
          double coarse_tier_PT_factor;     // coarse P/T scales = fine / factor
          double coarse_tier_X_factor;      // coarse X scales   = fine / factor

          int  snap_queries_to_bins;
          bool have_per_oxide_scales;
          std::vector<double> hash_X_scales;

          MAGEMin_hash::MAGEMin_Hash_table &pick_tier(const quad_point_data &qp) const;
          void setup_hash_tiers();

          // mass and volume fraction conversion functions
          double compute_bulk_density (const double porosity,
                                       const double solid_density,
                                       const double fluid_density) const;
          double compute_mass_fraction (const double volume_frac,
                                        const double material_density,
                                        const double bulk_density) const;


          // Helpers for the individual stages of the batched MAGEMin pipeline.
          void load_variables_from_present_time_step(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in, unsigned int q) const;
          // void normalize_drifted_compositions(quad_point_data &qp) const;
          // void normalize_solid_compositions(quad_point_data &qp) const;
          // void normalize_liquid_compositions(quad_point_data &qp) const;
          void normalize_compositions(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in,
                                      typename Interface<dim>::MaterialModelOutputs &out,
                                      const unsigned int q) const;
          void remove_melt_below_cutoff_depth(quad_point_data &qp,const typename Interface<dim>::MaterialModelInputs &in,  typename Interface<dim>::MaterialModelOutputs &out, unsigned q,MeltOutputs<dim> *melt_out, const double reference_T) const;
          double freeze_melt_below_certain_tenperature(quad_point_data &qp,const typename Interface<dim>::MaterialModelInputs &in ,typename Interface<dim>::MaterialModelOutputs &out, unsigned q) const;
          void handle_cases_when_there_is_no_solid_left(quad_point_data &qp, typename Interface<dim>::MaterialModelOutputs &out, unsigned q) const;
          bool skip_magemin_evaluation(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out, unsigned int q) const;
          void calculate_bulk_composition(quad_point_data &qp) const;
          void resolve_magemin_batch(std::vector<quad_point_data> &quad_points,
                                     const std::vector<unsigned int> &point_indices) const;
          void run_neural_network_batch(std::vector<quad_point_data> &quad_points,
                                        const std::vector<unsigned int> &point_indices) const;
          void fallback_to_magemin_batch(std::vector<quad_point_data> &quad_points,
                                         const std::vector<unsigned int> &point_indices) const;

          void fill_stableAssemblage_from_nn(const padawan::Prediction &pred,
                                             std::size_t prediction_index,
                                             quad_point_data &qp) const;

          bool try_running_hash_table(quad_point_data &qp) const;
          bool try_kdtree_lookup(quad_point_data &qp) const;

          void add_kdtree_point(double pressure_GPa,
                                double temperature_C,
                                const std::vector<double> &composition,
                                std::uint64_t assemblage_mask,
                                bool uses_fine_cache) const;
          void clear_kdtree() const;

          void extract_stable_assemblages_and_physical_parameters(quad_point_data &qp) const;

          void write_to_aspect_data_structures(quad_point_data &qp, const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out, unsigned int q, MeltOutputs<dim> *melt_out, EnthalpyOutputs<dim> *enthalpy_out) const;
          mutable double time_hash_us = 0;
          mutable double time_kdtree_us = 0;
          mutable double time_magemin_us = 0;


          // Neural-network timing and outcome counters.
          mutable double time_nn_us        = 0;
          mutable size_t count_nn_calls    = 0;  // total NN.predict() invocations
          mutable size_t count_nn_ood      = 0;  // NN flagged OOD

          void clear_cache_on_timestep(const SimulatorAccess<dim> &simulator_access);


          // Private helper methods
          void save_caches_on_checkpoint(
            typename parallel::distributed::Triangulation<dim> &);
          void load_caches_on_resume(
            typename parallel::distributed::Triangulation<dim> &);
          bool cache_loaded_before_resume = false;

          // Incremental KD-tree over cached thermodynamic states.
          using KDTreeDistance = nanoflann::L2_Simple_Adaptor<double, MAGEMinPointCloud>;
          using KDTree = nanoflann::KDTreeSingleIndexIncrementalAdaptor<
                         KDTreeDistance, MAGEMinPointCloud, -1, std::size_t>;

          mutable MAGEMinPointCloud kdtree_points;
          mutable std::unique_ptr<KDTree> kdtree;
          mutable size_t kdtree_hits;

          int surrogateMode;
          std::string database_name;
          std::vector<std::string> oxide_names;
          std::vector<std::string> solid_oxide_field_names;
          std::vector<std::string> liquid_oxide_field_names;
          std::vector<std::string> configured_solid_oxide_field_names;
          std::vector<std::string> configured_liquid_oxide_field_names;
          unsigned int water_oxide_index = std::numeric_limits<unsigned int>::max();

          double reference_rho_fluid;
          double xi_0;
          double viscosity_fluid;
          double thermal_bulk_viscosity_exponent;
          double alpha_phi;
          double melt_compressibility;
          double melt_bulk_modulus_derivative;
          double depletion_solidus_change;
          double melting_time_scale;
          double reference_permeability;
          double cutOff_depth;
          int useKd_tree;
          int use_MAGEMin_density;
          int use_MAGEMin_specific_heat;
          double max_mineral_range;
          double max_oxide_deviation;

          // Hash bin scales (inverse bin widths).
          double hash_table_P_scale;
          double hash_table_T_scale;
          double hash_table_X_scale;


          int kdtree_oxide_tol_bins;
          int kdtree_temperature_tol_bins;
          int kdtree_pressure_tol_bins;

          std::string surrogate_directory;
          std::string hash_storage_location;
          std::string cache_context_signature;
          int save_hash_table;
          bool store_cache_timestep = true;

          double A1 = 1085.7;
          double A2 = 1.329e-7;
          double A3 = -5.1e-18;
          double B1 = 1475;
          double B2 = 8.0e-8;
          double B3 = -3.2e-18;
          double C1 = 1780.0;
          double C2 = 4.50e-8;
          double C3 = -2.0e-18;
          double r1 = 0.5;
          double r2 = 8e-11;
          double beta = 1.5;
          double M_cpx = 0.15;

          // Hydrous Katz melting parameters.
          double K = 43.0;
          double gamma = 0.75;
          double D_H2O = 0.01;

          double Xie1 = 4.777286e-5;
          double Xie2 = 1.0e-9;
          double Lambda = 0.6;



          double melt_retention_threshold;
          double fractional_crystallization_start;
          void compute_solidus_state(quad_point_data &qp) const;


      };
    } // namespace ReactionModel
  } // namespace MaterialModel
} // namespace aspect

#endif














