/**
 * Initial implementation of magemin using a lookup table
 */

#ifndef _aspect_material_model_reaction_model_magemin_lookup_h
#define _aspect_material_model_reaction_model_magemin_lookup_h


#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {   
        struct magelookupFile
        {
            std::vector<double> Pressure;
            std::vector<double> Temperature;
            std::vector<std::vector<double>> materialProperties;
        };
        
        namespace internal
        {
            class readLookUpTable
            {
                public:
                    readLookUpTable(const std::string &filename);

                    std::vector<std::vector<double>> getData() const;
                    std::vector<std::vector<double>> data_magemin;
            // private:
            };
        }

        namespace internal2
        {
            class readLookUpTable
            {
                public:
                    // readLookUpTable(const std::string &filename);

                    readLookUpTable(std::string data_directory_magemin, std::vector<std::string>& compositionalLookupFileNames, std::vector<std::unique_ptr<magelookupFile>>& magelookupFiles);
                    magelookupFile getData() const;
                    // std::vector<std::vector<double>> data_magemin;
                private:
                    magelookupFile lF;
            // private:
            };
        }


        template <int dim>
        class mageminLookup : public ::aspect::SimulatorAccess<dim>
        {
            public:
                // constructor
                mageminLookup();

                /**
                 * Declare the parameters this function takes through input files.
                 */
                static
                void
                declare_parameters (ParameterHandler &prm);

                /**
                 * Read the parameters from the parameter file.
                 */
                void
                parse_parameters (ParameterHandler &prm);


                /**
                 * Function which for now finds the closest PT datapoint from the lookup table and fetches required variables
                 */
                double
                melt_fraction (const MaterialModel::MaterialModelInputs<dim> &in, unsigned int q) const;

                void calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                                typename Interface<dim>::MaterialModelOutputs &out) const;

                void calculate_fluid_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                       typename Interface<dim>::MaterialModelOutputs &out,
                                       const double reference_T) const;

                void initialize();
                void initializeNewandModern();

                // create a pointer to the structure holding magemin lookuptable properties
                std::vector<std::unique_ptr<magelookupFile>> magelookupFiles;
                std::vector<std::string> fileNames;
                std::vector<std::string> compositionalLookupNames;
                std::vector<std::string> compositionalLookupFileNames;

            private:
                double reference_rho_fluid;
                double xi_0;
                double viscosity_fluid;
                double thermal_bulk_viscosity_exponent;
                double alpha_phi;
                double extraction_depth;
                double melt_compressibility;
                double melt_bulk_modulus_derivative;
                double depletion_solidus_change;
                bool fractional_melting;
                double freezing_rate;
                double melting_time_scale;
                double reference_permeability;

                /**
                 * Magemin lookup table Parameters
                 */
                double crust_min_P;
                double crust_max_P; 
                double crust_min_T; 
                double crust_max_T;

                double mantle_min_P; 
                double mantle_max_P;
                double mantle_min_T; 
                double mantle_max_T;

                /**
                 * Magemin filepaths
                 */
                std::string data_directory_magemin;
                std::string mantle_data_file_name;
                std::string crust_data_file_name;


                std::unique_ptr<internal::readLookUpTable> mantleLookupData;
                std::unique_ptr<internal::readLookUpTable> crustLookupData;

                std::vector<std::vector<double>> mantle_data;
                std::vector<std::vector<double>> crust_data;

                std::vector<std::vector<double>> compositionalLookupNames;
                std::vector<std::vector<double>> compositionalLookupFileNames;


        };
    }
  }
}

#endif