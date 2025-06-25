/**
 * Initial implementation of magemin using a lookup table
 */

#ifndef _aspect_material_model_reaction_model_magemin_lookup_h
#define _aspect_material_model_reaction_model_magemin_lookup_h

#include <aspect/material_model/interface.h>
#include <aspect/melt.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/simulator_access.h>

// #ifdef ASPECT_WITH_MAGEMin
#include <MAGEMin_cpp.h>
// #endif

namespace aspect {
namespace MaterialModel {
namespace ReactionModel {
struct magelookupFile {
  std::string compositionName;
  std::vector<double> Pressure;
  std::vector<double> Temperature;
  std::vector<std::vector<double>> materialProperties;
};

namespace internal {
class readLookUpTable {
public:
  readLookUpTable(const std::string &filename, const MPI_Comm comm);

  std::vector<std::vector<double>> getData() const;
  std::vector<std::vector<double>> data_magemin;
  // private:
};
} // namespace internal

namespace internal2 {
class readLookUpTableFunky {
public:
  // readLookUpTable(const std::string &filename);

  readLookUpTableFunky(
      std::string data_directory_magemin,
      std::vector<std::string> &compositionalLookupFileNames,
      std::vector<std::string> &compositionalLookupNames,
      std::vector<std::unique_ptr<magelookupFile>> &magelookupFiles,
      const MPI_Comm comm);
  magelookupFile getData() const;
  // std::vector<std::vector<double>> data_magemin;
private:
  magelookupFile magelF;
  // private:
};
} // namespace internal2

template <int dim> class mageminLookup : public ::aspect::SimulatorAccess<dim> {
public:
  // constructor
  mageminLookup();
  // Make object for MAGEMin class
  mutable MAGEMin_wrapper wrap;

  // create a pointer to the structure holding magemin lookuptable properties

  std::vector<std::unique_ptr<magelookupFile>> magelookupFiles;
  std::unique_ptr<internal2::readLookUpTableFunky> allDataPtr;

  std::vector<std::string> fileNames;

  std::vector<std::string> compositionalLookupNames;
  std::vector<std::string> compositionalLookupFileNames;

  /**
   * Declare the parameters this function takes through input files.
   */
  static void declare_parameters(ParameterHandler &prm);

  /**
   * Read the parameters from the parameter file.
   */
  void parse_parameters(ParameterHandler &prm);

  /**
   * Function which for now finds the closest PT datapoint from the lookup table
   * and fetches required variables
   */
  double melt_fraction(const typename Interface<dim>::MaterialModelInputs &in,
                       unsigned int q) const;

  /**
   * Function which directly calls MAGEMin
   */
  double
  melt_fractionMAGEMin(const typename Interface<dim>::MaterialModelInputs &in,
                       typename Interface<dim>::MaterialModelOutputs &out,
                       unsigned int q) const;

  /**
   * Function which uses Katz parameterization of solidus lines to predict melt;
   * just to reduce number of points where we call MAGEMin_cpp
   */
  bool guess_MeltFraction(double pressure, double temperature) const;

  void calculate_reaction_rate_outputs(
      const typename Interface<dim>::MaterialModelInputs &in,
      typename Interface<dim>::MaterialModelOutputs &out) const;

  void calculate_fluid_outputs(
      const typename Interface<dim>::MaterialModelInputs &in,
      typename Interface<dim>::MaterialModelOutputs &out,
      const double reference_T) const;

  void initialize();

  void initializeNewandModern();

  int get_closest_index(
      float T, float P,
      const std::vector<std::unique_ptr<magelookupFile>> &magelookupFiles,
      int largestCompIndex) const;

  // struct magelookupFile
  // {
  //     std::vector<double> Pressure;
  //     std::vector<double> Temperature;
  //     std::vector<std::vector<double>> materialProperties;
  // };

private:
  // void readLookUpTableFunky(std::string data_directory_magemin,
  // std::vector<std::string>& compositionalLookupFileNames,
  // std::vector<std::unique_ptr<magelookupFile>>& magelookupFiles);
  // magelookupFile lF;

  // void readLookUpTable(const std::string &filename);
  // std::vector<std::vector<double>> getData() const;
  // std::vector<std::vector<double>> data_magemin;

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
  double cutOff_depth;

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

  // Katz parameters for dummy calculations;
  double A1 = 1085.7;
  double A2 = 1.329e-7;
  double A3 = -5.1e-18;

  // std::vector<std::vector<double>> compositionalLookupNames;
  // std::vector<std::vector<double>> compositionalLookupFileNames;
};
} // namespace ReactionModel
} // namespace MaterialModel
} // namespace aspect

#endif
