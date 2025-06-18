
/**
 * Initial implementation of magemin using a lookup table
 */

#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/reaction_model/magemin_lookup.h>
#include <aspect/utilities.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect {
namespace MaterialModel {
namespace ReactionModel {
namespace internal {
readLookUpTable::readLookUpTable(const std::string &filename,
                                 const MPI_Comm comm) {

  // std::ifstream file(filename);
  std::string line;

  // Read data from disk and distribute among processes
  std::istringstream file(
      Utilities::read_and_distribute_file_content(filename, comm));

  if (!file) {
    std::cerr << "Error: Unable to open file " << filename << "\n";
  }

  while (std::getline(file, line)) {

    // Skip header lines
    if (line[0] == '#')
      continue;

    std::stringstream ss(line);
    double value;
    std::vector<double> row;

    while (ss >> value) {
      row.push_back(value);
    }

    // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<" "<<row[10]<<std::endl;
    // Ensure there are enough columns before accessing indices
    if (row.size() >= 17) {
      data_magemin.push_back(
          {row[0], row[1], row[2],
           row[10]}); // 5th and 6th columns (zero-based index 4 and 5)
    }
  }
  // file.close();
}

std::vector<std::vector<double>> readLookUpTable::getData() const {
  return data_magemin;
}
} // namespace internal

namespace internal2 {
readLookUpTableFunky::readLookUpTableFunky(
    std::string data_directory_magemin,
    std::vector<std::string> &compositionalLookupFileNames,
    std::vector<std::string> &compositionalLookupNames,
    std::vector<std::unique_ptr<magelookupFile>> &magelookupFiles,
    const MPI_Comm comm) {

  for (size_t i = 0; i < compositionalLookupFileNames.size(); i++) {
    // std::cout<<"\nReading"<<compositionalLookupFileNames[i]<<"\n";
    // std::ifstream file(compositionalLookupFileNames[i]);
    std::string line;

    // Read data from disk and distribute among processes
    std::istringstream file(Utilities::read_and_distribute_file_content(
        data_directory_magemin + compositionalLookupFileNames[i], comm));

    if (!file) {
      std::cerr << "Error: Unable to open file "
                << compositionalLookupFileNames[i] << "\n";
    }

    while (std::getline(file, line)) {

      // Skip header lines
      if (line[0] == '#')
        continue;

      std::stringstream ss(line);
      double value;
      std::vector<double> row;

      while (ss >> value) {
        row.push_back(value);
      }

      // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<" "<<row[10]<<std::endl;
      // Ensure there are enough columns before accessing indices
      if (row.size() >= 17) {
        magelookupFiles[i]->compositionName = compositionalLookupNames[i];
        magelookupFiles[i]->Pressure.push_back(row[0]);
        magelookupFiles[i]->Temperature.push_back(row[1]);
        magelookupFiles[i]->materialProperties.push_back({row[2], row[10]});
        // data_magemin.push_back({row[0],row[1],row[2],row[10]}); // 5th and
        // 6th columns (zero-based index 4 and 5)
      }
    }
    // file.close();
  }
}
magelookupFile readLookUpTableFunky::getData() const {
  // dummy=1;
  return magelF;
}
} // namespace internal2

template <int dim> mageminLookup<dim>::mageminLookup() = default;

template <int dim> void mageminLookup<dim>::initializeNewandModern() {

  // Create unique_ptr for each magelookupFile structure instance for each
  // material
  for (size_t i = 0; i < compositionalLookupFileNames.size(); i++) {
    magelookupFiles.push_back(std::make_unique<magelookupFile>());
  }
  allDataPtr = std::make_unique<internal2::readLookUpTableFunky>(
      data_directory_magemin, compositionalLookupFileNames,
      compositionalLookupNames, magelookupFiles, this->get_mpi_communicator());

  // std::cout<<"Initialize new and modern"<<"\n";
  // // Verifying Data
  // for (size_t i = 0; i < magelookupFiles.size(); i++)
  // {
  //     std::cout<<"Printing "<<compositionalLookupFileNames[i]<<"\n";
  //     std::cout<<"Printing "<<magelookupFiles[i]->Pressure.size()<<"\n";
  //     for (int j = 0; j < 2; j++)
  //     {
  //         std::cout << magelookupFiles[i]->Pressure[j]<<"
  //         "<<magelookupFiles[i]->Temperature[j]<<"
  //         "<<magelookupFiles[i]->materialProperties[j][0] << "\n";
  //     }
  // }
}

template <int dim> void mageminLookup<dim>::initialize() {

  crustLookupData = std::make_unique<internal::readLookUpTable>(
      data_directory_magemin + crust_data_file_name,
      this->get_mpi_communicator());
  mantleLookupData = std::make_unique<internal::readLookUpTable>(
      data_directory_magemin + mantle_data_file_name,
      this->get_mpi_communicator());

  mantle_data = mantleLookupData->getData();
  crust_data = crustLookupData->getData();

  // // Iterating using indexes
  // std::cout<<"initialize"<<"\n";
  // for (int i = 0; i < 2; i++)
  // {
  //     for (size_t j = 0; j < mantle_data[i].size(); j++)
  //     {
  //         std::cout << mantle_data[i][j] << " ";
  //     }
  //     std::cout << "\n";
  // }
}

template <int dim>
int mageminLookup<dim>::get_closest_index(
    float T, float P,
    const std::vector<std::unique_ptr<magelookupFile>> &magelookupFiles,
    int largestCompIndex) const {
  int closestIndex = -1;
  double minDiff = 999999999999;
  std::cout << "Checking if this works or not"
            << magelookupFiles[largestCompIndex]->compositionName;
  for (size_t i = 0; i < magelookupFiles[largestCompIndex]->Pressure.size();
       ++i) {
    double diff =
        std::abs(magelookupFiles[largestCompIndex]->Pressure[i] - P) +
        std::abs(magelookupFiles[largestCompIndex]->Temperature[i] - T);
    if (diff < minDiff) {
      minDiff = diff;
      closestIndex = i;
    }
  }

  if (minDiff == 999999999999) {
    return -999;
  }

  return closestIndex;
}

template <int dim>
double mageminLookup<dim>::melt_fraction(
    const MaterialModel::MaterialModelInputs<dim> &in, unsigned int q) const {
  double meltFraction = 0.0;

  // const unsigned int bound_fluid_idx =
  // this->introspection().compositional_index_for_name("bound_fluid"); // Needs
  // to go

  float currP = this->get_adiabatic_conditions().pressure(in.position[q]);
  ; // this->get_adiabatic_conditions().pressure(in.position[q]);
  const float currT = in.temperature[q]; // in.temperature[q];
  const float currDepth = this->get_geometry_model().depth(
      in.position[q]); // this->get_geometry_model().depth(in.position[q]);

  // Put in a check for composition; This is where I call the function to
  // calculate closest P and T conditions;
  float compVal = 0.0;
  int largestCompIndex = 0;
  int closestMaterialPropertiesIndex = 0;

  if (currDepth < cutOff_depth) {
    // Loop through all the composition to find the maximum composition
    for (size_t comp = 0; comp < compositionalLookupNames.size(); comp++) {
      const unsigned int comp_idx =
          this->introspection().compositional_index_for_name(
              compositionalLookupNames[comp]);
      if (in.composition[q][comp_idx] > compVal) {
        std::cout << "\n"
                  << in.composition[q][comp_idx] << "\n"
                  << compositionalLookupNames[comp];
        largestCompIndex = comp;
        compVal = in.composition[q][comp_idx];
      }
    }

    closestMaterialPropertiesIndex =
        get_closest_index(currT, currP, magelookupFiles, largestCompIndex);

    if (closestMaterialPropertiesIndex == -999) {
      meltFraction = 0.0;
    } else {
      meltFraction =
          magelookupFiles[largestCompIndex]
              ->materialProperties[closestMaterialPropertiesIndex][0];
    }
    std::cout << "\n"
              << compositionalLookupNames[largestCompIndex] << " " << currP
              << "  " << currT << " " << meltFraction;
  } else {
    meltFraction = 0.0;
  }

  return meltFraction;
}

template <int dim>
void mageminLookup<dim>::calculate_fluid_outputs(
    const typename Interface<dim>::MaterialModelInputs &in,
    typename Interface<dim>::MaterialModelOutputs &out,
    const double reference_T) const {
  MeltOutputs<dim> *melt_out =
      out.template get_additional_output<MeltOutputs<dim>>();
  if (melt_out != nullptr) {
    for (unsigned int i = 0; i < in.n_evaluation_points(); ++i) {
      const unsigned int porosity_idx =
          this->introspection().compositional_index_for_name("porosity");
      double porosity = std::max(in.composition[i][porosity_idx], 0.0);

      melt_out->fluid_viscosities[i] = viscosity_fluid;
      melt_out->permeabilities[i] = reference_permeability *
                                    Utilities::fixed_power<3>(porosity) *
                                    Utilities::fixed_power<2>(1.0 - porosity);

      /// need to find nearest element for density

      // Calculate dependence of fluid density on temperature if adiabatic
      // heating in switched on first, calculate temperature dependence of
      // density
      double temperature_dependence = 1.0;
      if (this->include_adiabatic_heating()) {
        // temperature dependence is 1 - alpha * (T - T(adiabatic))
        temperature_dependence -=
            (in.temperature[i] -
             this->get_adiabatic_conditions().temperature(in.position[i])) *
            out.thermal_expansion_coefficients[i];
      } else
        temperature_dependence -= (in.temperature[i] - reference_T) *
                                  out.thermal_expansion_coefficients[i];

      // the fluid compressibility includes two parts, a constant
      // compressibility, and a pressure-dependent one
      // this is a simplified formulation, experimental data are often fit to
      // the Birch-Murnaghan equation of state
      const double fluid_compressibility =
          melt_compressibility /
          (1.0 + in.pressure[i] * melt_bulk_modulus_derivative *
                     melt_compressibility);

      melt_out->fluid_densities[i] =
          reference_rho_fluid *
          std::exp(fluid_compressibility *
                   (in.pressure[i] - this->get_surface_pressure())) *
          temperature_dependence;

      melt_out->fluid_density_gradients[i] =
          melt_out->fluid_densities[i] * melt_out->fluid_densities[i] *
          fluid_compressibility *
          this->get_gravity_model().gravity_vector(in.position[i]);

      const double phi_0 = 0.05;
      porosity = std::max(std::min(porosity, 0.995), 1e-4);
      melt_out->compaction_viscosities[i] = xi_0 * phi_0 / porosity;

      double visc_temperature_dependence = 1.0;
      if (this->include_adiabatic_heating()) {
        const double delta_temp =
            in.temperature[i] -
            this->get_adiabatic_conditions().temperature(in.position[i]);
        visc_temperature_dependence = std::max(
            std::min(std::exp(-thermal_bulk_viscosity_exponent * delta_temp /
                              this->get_adiabatic_conditions().temperature(
                                  in.position[i])),
                     1e4),
            1e-4);
      } else {
        const double delta_temp = in.temperature[i] - reference_T;
        const double T_dependence =
            (thermal_bulk_viscosity_exponent == 0.0
                 ? 0.0
                 : thermal_bulk_viscosity_exponent * delta_temp / reference_T);
        visc_temperature_dependence =
            std::max(std::min(std::exp(-T_dependence), 1e4), 1e-4);
      }
      melt_out->compaction_viscosities[i] *= visc_temperature_dependence;
    }
  }

  if (this->include_melt_transport() &&
      in.requests_property(MaterialProperties::viscosity)) {
    for (unsigned int i = 0; i < in.n_evaluation_points(); ++i) {
      const double porosity =
          std::min(1.0, std::max(in.composition[i][porosity_idx], 0.0));
      out.viscosities[i] *= std::exp(-alpha_phi * porosity);
    }
  }
}

/**
 * Function to modify material properties based on minimization using MAGEMin
 */
template <int dim>
void mageminLookup<dim>::calculate_reaction_rate_outputs(
    const typename Interface<dim>::MaterialModelInputs &in,
    typename Interface<dim>::MaterialModelOutputs &out) const {
  ReactionRateOutputs<dim> *reaction_rate_out =
      out.template get_additional_output<ReactionRateOutputs<dim>>();

  // Fill reaction rate outputs if the model uses operator splitting.
  // Specifically, change the porosity (representing the amount of free fluid)
  // based on the water solubility and the fluid content.
  if (this->get_parameters().use_operator_splitting &&
      reaction_rate_out != nullptr) {
    std::vector<double> eq_free_fluid_fractions(out.n_evaluation_points());
    melt_fractions(in, eq_free_fluid_fractions);

    for (unsigned int q = 0; q < out.n_evaluation_points(); ++q)
      for (unsigned int c = 0; c < in.composition[q].size(); ++c) {
        double porosity_change =
            eq_free_fluid_fractions[q] - in.composition[q][porosity_idx];
        // do not allow negative porosity
        if (in.composition[q][porosity_idx] + porosity_change < 0)
          porosity_change = -in.composition[q][porosity_idx];

        const unsigned int bound_fluid_idx =
            this->introspection().compositional_index_for_name("bound_fluid");
        if (c == bound_fluid_idx && this->get_timestep_number() > 0)
          reaction_rate_out->reaction_rates[q][c] =
              -porosity_change / fluid_reaction_time_scale;
        else if (c == porosity_idx && this->get_timestep_number() > 0)
          reaction_rate_out->reaction_rates[q][c] =
              porosity_change / fluid_reaction_time_scale;
        else
          reaction_rate_out->reaction_rates[q][c] = 0.0;
      }
  }
}

template <int dim>
void mageminLookup<dim>::declare_parameters(ParameterHandler &prm) {

  prm.declare_entry("Crust Minimum Pressure", "1e8", Patterns::Double(0.),
                    "Crust Minimum Pressure. Units: \\si{\\Pa}.");
  prm.declare_entry("Crust Maximum Pressure", "0.2e9", Patterns::Double(0.),
                    "Crust Maximum Pressure. Units: \\si{\\Pa}.");
  prm.declare_entry("Crust Minimum Temperature", "973.", Patterns::Double(0.),
                    "Crust Minimum Temperature. Units: \\si{\\kelvin}.");
  prm.declare_entry("Crust Maximum Temperature", "1273.", Patterns::Double(0.),
                    "Crust Maximum Tempetrature. Units: \\si{\\kelvin}.");

  prm.declare_entry("Mantle Minimum Pressure", "1e8", Patterns::Double(0.),
                    "Mantle Minimum Pressure. Units: \\si{\\Pa}.");
  prm.declare_entry("Mantle Maximum Pressure", "2.0e9", Patterns::Double(0.),
                    "Mantle Maximum Pressure. Units: \\si{\\Pa}.");
  prm.declare_entry("Mantle Minimum Temperature", "1000", Patterns::Double(0.),
                    "Mantle Minimum Temperature. Units: \\si{\\kelvin}.");
  prm.declare_entry("Mantle Maximum Temperature", "1673.", Patterns::Double(0.),
                    "Mantle Maximum Tempetrature. Units: \\si{\\kelvin}.");

  prm.declare_entry("Data directory Lookup Table Magemin",
                    "$ASPECT_SOURCE_DIR/data/material-model/",
                    Patterns::DirectoryName(), "magemin loopup table");

  prm.declare_entry("Composition names", "upper_crust, lower_crust, mantle",
                    Patterns::Anything(),
                    "List of material model names, ideally same to all "
                    "composition field you used. ");

  prm.declare_entry("File names for compositions",
                    "upper_crust.dat, lower_crust.dat, mantle.dat",
                    Patterns::Anything(),
                    "List of file names of corresponding composition models. ");

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
  prm.declare_entry(
      "Melt extraction depth", "1000.0", Patterns::Double(0.),
      "Depth above that melt will be extracted from the model, "
      "which is done by a negative reaction term proportional to the "
      "porosity field. "
      "Units: \\si{\\meter}.");
  prm.declare_entry("Melt compressibility", "0.0", Patterns::Double(0.),
                    "The value of the compressibility of the melt. "
                    "Units: \\si{\\per\\pascal}.");
  prm.declare_entry("Melt bulk modulus derivative", "0.0", Patterns::Double(0.),
                    "The value of the pressure derivative of the melt bulk "
                    "modulus. "
                    "Units: None.");
  prm.declare_entry(
      "Use fractional melting", "false", Patterns::Bool(),
      "If fractional melting should be used (if true), including a solidus "
      "change based on depletion (in this case, the amount of melt that has "
      "migrated away from its origin), and freezing of melt when it has moved "
      "to a region with temperatures lower than the solidus; or if batch "
      "melting should be used (if false), assuming that the melt fraction only "
      "depends on temperature and pressure, and how much melt has already been "
      "generated at a given point, but not considering movement of melt in "
      "the melting parameterization."
      "\n\n"
      "Note that melt does not freeze unless the 'Freezing rate' parameter is "
      "set "
      "to a value larger than 0.");
  prm.declare_entry(
      "Freezing rate", "0.0", Patterns::Double(0.),
      "Freezing rate of melt when in subsolidus regions. "
      "If this parameter is set to a number larger than 0.0, it specifies the "
      "fraction of melt that will freeze per year (or per second, depending on "
      "the "
      "``Use years in output instead of seconds'' parameter), as soon as the "
      "porosity "
      "exceeds the equilibrium melt fraction, and the equilibrium melt "
      "fraction "
      "falls below the depletion. In this case, melt will freeze according to "
      "the "
      "given rate until one of those conditions is not fulfilled anymore. The "
      "reasoning behind this is that there should not be more melt present "
      "than "
      "the equilibrium melt fraction, as melt production decreases with "
      "increasing "
      "depletion, but the freezing process of melt also reduces the depletion "
      "by "
      "the same amount, and as soon as the depletion falls below the "
      "equilibrium "
      "melt fraction, we expect that material should melt again (no matter how "
      "much melt is present). This is quite a simplification and not a "
      "realistic "
      "freezing parameterization, but without tracking the melt composition, "
      "there "
      "is no way to compute freezing rates accurately. "
      "If this parameter is set to zero, no freezing will occur. "
      "Note that freezing can never be faster than determined by the "
      "``Melting time scale for operator splitting''. The product of the "
      "``Freezing rate'' and the ``Melting time scale for operator splitting'' "
      "defines how fast freezing occurs with respect to melting (if the "
      "product is 0.5, melting will occur twice as fast as freezing). "
      "Units: 1/yr or 1/s, depending on the ``Use years "
      "in output instead of seconds'' parameter.");
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
      "values, the solidus gets increased for a positive peridotite field "
      "(depletion) and lowered for a negative peridotite field (enrichment). "
      "Scaling with depletion is linear. Only active when fractional melting "
      "is used. "
      "Units: \\si{\\kelvin}.");
  prm.declare_entry("Reference permeability", "1e-8", Patterns::Double(),
                    "Reference permeability of the solid host rock."
                    "Units: \\si{\\meter\\squared}.");
  prm.declare_entry(
      "Cutoff depth for magemin", "300e3", Patterns::Double(),
      "Cutoff depth: Melting will only be considered above this depth"
      "Units: \\si{\\meter}.");
}

template <int dim>
void mageminLookup<dim>::parse_parameters(ParameterHandler &prm) {
  // New Parameters
  crust_min_P = prm.get_double("Crust Minimum Pressure");
  crust_max_P = prm.get_double("Crust Maximum Pressure");
  crust_min_T = prm.get_double("Crust Minimum Temperature");
  crust_max_T = prm.get_double("Crust Maximum Temperature");

  mantle_min_P = prm.get_double("Mantle Minimum Pressure");
  mantle_max_P = prm.get_double("Mantle Maximum Pressure");
  mantle_min_T = prm.get_double("Mantle Minimum Temperature");
  mantle_max_T = prm.get_double("Mantle Maximum Temperature");

  mantle_data_file_name = "arranged_output_mantle.dat";
  crust_data_file_name = "arranged_output_crust.dat";

  // Read my models Names
  data_directory_magemin = Utilities::expand_ASPECT_SOURCE_DIR(
      prm.get("Data directory Lookup Table Magemin"));

  compositionalLookupNames = {"c1", "c2"};
  // Utilities::split_string_list(prm.get("Composition names"), ',');
  compositionalLookupFileNames = {"c1.dat", "c2.dat"};
  // Utilities::split_string_list(prm.get("File names for compositions"), ',');

  reference_rho_fluid = prm.get_double("Reference melt density");
  xi_0 = prm.get_double("Reference bulk viscosity");
  viscosity_fluid = prm.get_double("Reference melt viscosity");
  thermal_bulk_viscosity_exponent =
      prm.get_double("Thermal bulk viscosity exponent");
  alpha_phi = prm.get_double("Exponential melt weakening factor");
  extraction_depth = prm.get_double("Melt extraction depth");
  melt_compressibility = prm.get_double("Melt compressibility");
  fractional_melting = prm.get_bool("Use fractional melting");
  freezing_rate = prm.get_double("Freezing rate");
  melting_time_scale =
      prm.get_double("Melting time scale for operator splitting");
  melt_bulk_modulus_derivative = prm.get_double("Melt bulk modulus derivative");
  depletion_solidus_change = prm.get_double("Depletion solidus change");
  reference_permeability = prm.get_double("Reference permeability");
  cutOff_depth = prm.get_double("Cutoff depth for magemin");

  if (this->convert_output_to_years() == true) {
    melting_time_scale *= year_in_seconds;
    freezing_rate /= year_in_seconds;
  }

  AssertThrow(melting_time_scale > 0,
              ExcMessage("The Melting time scale for operator splitting must "
                         "be larger than 0!"));

  if (this->get_parameters().reaction_solver_type ==
      Parameters<dim>::ReactionSolverType::fixed_step) {
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

    AssertThrow(
        freezing_rate * this->get_parameters().reaction_time_step <= 1.0,
        ExcMessage(
            "The reaction time step " +
            Utilities::to_string(this->get_parameters().reaction_time_step) +
            " in the operator splitting scheme is too large to compute "
            "freezing rates! "
            "You have to choose it in such a way that it is smaller than the "
            "inverse of the "
            "'Freezing rate' chosen in the material model, which is "
            "currently " +
            Utilities::to_string(1.0 / freezing_rate) + "."));
  }
}

} // namespace ReactionModel
} // namespace MaterialModel
} // namespace aspect

// explicit instantiations
namespace aspect {
namespace MaterialModel {
#define INSTANTIATE(dim)                                                       \
  namespace ReactionModel {                                                    \
  template class mageminLookup<dim>;                                           \
  }

ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
} // namespace MaterialModel
} // namespace aspect
