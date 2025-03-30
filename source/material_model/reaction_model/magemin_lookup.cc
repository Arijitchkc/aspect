/**
 * Initial implementation of magemin using a lookup table
 */

#include <aspect/material_model/reaction_model/katz2003_mantle_melting.h>
#include <aspect/utilities.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
    namespace MaterialModel
    {
        namespace ReactionModel
        {
            namespace internal
            {
                readLookUpTable::readLookUpTable(const std::string &filename) 
                {

                    std::ifstream file(filename);
                    std::string line;

                    if (!file) 
                    {
                        std::cerr << "Error: Unable to open file " << filename << std::endl;
                    }

                    while (std::getline(file, line)) 
                    {
                        
                        // Skip header lines
                        if (line[0] == '#') continue;
                        
                        std::stringstream ss(line);
                        double value;
                        std::vector<double> row;
                        
                        while (ss >> value) 
                        {
                            row.push_back(value);
                        }

                        // std::cout<<row[0]<<" "<<row[1]<<" "<<row[2]<<" "<<row[10]<<std::endl;
                        // Ensure there are enough columns before accessing indices
                        if (row.size() >= 17) 
                        {
                            data_magemin.push_back({row[0],row[1],row[2],row[10]}); // 5th and 6th columns (zero-based index 4 and 5)
                        }
                    } 
                    file.close();
                }

                std::vector<std::vector<double>> readLookUpTable::getData() const
                {
                    return data_magemin;
                }
            }



            template <int dim>
            mageminLookup<dim>::mageminLookup()
                = default;
            
            template <int dim>
            mageminLookup<dim>::initialize()
            {
                crustLookupData = std::make_unique<internal::readLookUpTable>(data_directory_magemin+crust_data_file_name);
                mantleLookupData = std::make_unique<internal::readLookUpTable>(data_directory_magemin+mantle_data_file_name);
 
                mantle_data = mantleLookupData->getData(); 
                crust_data = mantleLookupData->getData();

                // Iterating using indexes
                std::cout<<"initialize"<<std::endl;
                for (int i = 0; i < 2; i++) 
                {
                    for (int j = 0; j < mantle_data[i].size(); j++) 
                    {
                        std::cout << mantle_data[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
            }

            template <int dim>
            double
            mageminLookup<dim>::
            melt_fraction (const MaterialModel::MaterialModelInputs<dim> &in, unsigned int q) const
            {
                double meltFraction;


                const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid"); // Needs to go


                float currP=this->get_adiabatic_conditions().pressure(in.position[q]);; //this->get_adiabatic_conditions().pressure(in.position[q]);
                const float currT=in.temperature[q]; //in.temperature[q];
                const float currDepth=this->get_geometry_model().depth(in.position[q]); //this->get_geometry_model().depth(in.position[q]);

                int reqd_index=0;
                float melt_fraction=0.0;
                
                // // Check which File to use -> crust or mantle
                // if((currDepth>=0) && (currDepth<=5e3))
                // {
                    
                //     float check_min_T=abs(crust_data[0][1]-currT), check_min_P=abs(crust_data[0][0]*1e+8-currP);  
                //     if((currT<=crust_max_T) && (currT>=crust_min_T-200) && (currP>=crust_min_P) && (currP<=crust_max_P))
                //     {
                //         for(int i=1;i<crust_data.size();i++)
                //         {
                //             if((check_min_T>=abs(crust_data[i][1]-currT)) && (check_min_P>=abs(crust_data[i][0]*1e+8-currP)))
                //             {
                //                 check_min_T=abs(crust_data[i][1]-currT);
                //                 check_min_P=abs(crust_data[i][0]*1e+8-currP);
                //                 reqd_index=i;
                //             }
                //         }
                //         melt_fraction=crust_data[reqd_index][2] + crust_data[reqd_index][3];
                //     }
                //     else
                //     {
                //         melt_fraction=0.0;
                //     }
                //     melt_fractions[q]=melt_fraction;
                // }

                std::cout<<"inside melt_fraction"<<std::endl;
                for (int i = 0; i < 2; i++) 
                {
                    for (int j = 0; j < mantle_data[i].size(); j++) 
                    {
                        std::cout << mantle_data[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                reqd_index=0;
                // double calcPressure=3000*9.8*currDepth;
                // currP=calcPressure;
                if ((currDepth>=5e3) && (currDepth<80e3))
                {
                    if(currP<0)
                    {
                    currP=0;
                    }

                    // Magemin file in kbar and Kelvin
                    float check_min_T=abs(mantle_data[0][1]-currT);
                    float check_min_P=abs(mantle_data[0][0]*1e8-currP);

                    if((currT<=mantle_max_T+20) && (currT>=mantle_min_T-20) && (currP>=mantle_min_P) && (currP<=mantle_max_P))
                    {    
                        int checkFound=999;
                        for(int i=0;i<mantle_data.size();i++)
                        {
                            if((check_min_T>=abs(mantle_data[i][1]-currT)) && (check_min_P>=abs(mantle_data[i][0]*1e+8-currP)))
                            {
                                check_min_T=abs(mantle_data[i][1]-currT);
                                check_min_P=abs(mantle_data[i][0]*1e+8-currP);
                                reqd_index=i;
                                checkFound=0;
                            }
                        }
                        
                        if(checkFound==0)
                        {
                            meltFraction=mantle_data[reqd_index][2]; //+ mantle_data[reqd_index][3];
                        }
                        else
                        {
                            meltFraction=0.0;
                        }
                        
                    }
                    else
                    {
                        meltFraction=0.0;
                    }

                    // melt_fractions[q]=melt_fraction;
                }
                else
                {
                    meltFraction=0.0;
                }

                return meltFraction;
            }



            /**
             * Function to modify material properties based on melting
             */
            template <int dim>
            void
            Katz2003MantleMelting<dim>::
            calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                            typename Interface<dim>::MaterialModelOutputs &out) const
            {
                
            }


            template <int dim>
            void
            mageminLookup<dim>::declare_parameters (ParameterHandler &prm)
            {


                prm.declare_entry ("Crust Minimum Pressure", "1e8",
                            Patterns::Double (0.),
                            "Crust Minimum Pressure. Units: \\si{\\Pa}.");
                prm.declare_entry ("Crust Maximum Pressure", "0.2e9",
                                    Patterns::Double (0.),
                                    "Crust Maximum Pressure. Units: \\si{\\Pa}.");
                prm.declare_entry ("Crust Minimum Temperature", "973.",
                                    Patterns::Double (0.),
                                    "Crust Minimum Temperature. Units: \\si{\\kelvin}.");
                prm.declare_entry ("Crust Maximum Temperature", "1273.",
                                    Patterns::Double (0.),
                                    "Crust Maximum Tempetrature. Units: \\si{\\kelvin}.");
                        

                prm.declare_entry ("Mantle Minimum Pressure", "1e8",
                                    Patterns::Double (0.),
                                    "Mantle Minimum Pressure. Units: \\si{\\Pa}.");
                prm.declare_entry ("Mantle Maximum Pressure", "2.0e9",
                                    Patterns::Double (0.),
                                    "Mantle Maximum Pressure. Units: \\si{\\Pa}.");
                prm.declare_entry ("Mantle Minimum Temperature", "1000",
                                    Patterns::Double (0.),
                                    "Mantle Minimum Temperature. Units: \\si{\\kelvin}.");
                prm.declare_entry ("Mantle Maximum Temperature", "1673.",
                                    Patterns::Double (0.),
                                    "Mantle Maximum Tempetrature. Units: \\si{\\kelvin}.");


                prm.declare_entry ("Data directory Lookup Table Magemin", "$ASPECT_SOURCE_DIR/data/material-model/",
                            Patterns::DirectoryName (),
                            "magemin loopup table");
                prm.declare_entry ("Mantle data file name", "arranged_output_mantle.dat",
                                    Patterns::Anything (),
                                    "The file name of Mantle data. ");
                prm.declare_entry ("Crust data file name", "arranged_output_crust.dat",
                                    Patterns::Anything (),
                                    "The file name of the Crustal data. ");


                prm.declare_entry ("Reference melt density", "2500.",
                        Patterns::Double (0.),
                        "Reference density of the melt/fluid$\\rho_{f,0}$. "
                        "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
                prm.declare_entry ("Reference bulk viscosity", "1e22",
                                Patterns::Double (0.),
                                "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                                "This viscosity may be modified by both temperature and porosity "
                                "dependencies. Units: \\si{\\pascal\\second}.");
                prm.declare_entry ("Reference melt viscosity", "10.",
                                Patterns::Double (0.),
                                "The value of the constant melt viscosity $\\eta_f$. Units: \\si{\\pascal\\second}.");
                prm.declare_entry ("Exponential melt weakening factor", "27.",
                                Patterns::Double (0.),
                                "The porosity dependence of the viscosity. Units: dimensionless.");
                prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                                Patterns::Double (0.),
                                "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                                "See the general documentation "
                                "of this model for a formula that states the dependence of the "
                                "viscosity on this factor, which is called $\\beta$ there.");
                prm.declare_entry ("Melt extraction depth", "1000.0",
                                Patterns::Double (0.),
                                "Depth above that melt will be extracted from the model, "
                                "which is done by a negative reaction term proportional to the "
                                "porosity field. "
                                "Units: \\si{\\meter}.");
                prm.declare_entry ("Melt compressibility", "0.0",
                                Patterns::Double (0.),
                                "The value of the compressibility of the melt. "
                                "Units: \\si{\\per\\pascal}.");
                prm.declare_entry ("Melt bulk modulus derivative", "0.0",
                                Patterns::Double (0.),
                                "The value of the pressure derivative of the melt bulk "
                                "modulus. "
                                "Units: None.");
                prm.declare_entry ("Use fractional melting", "false",
                                Patterns::Bool (),
                                "If fractional melting should be used (if true), including a solidus "
                                "change based on depletion (in this case, the amount of melt that has "
                                "migrated away from its origin), and freezing of melt when it has moved "
                                "to a region with temperatures lower than the solidus; or if batch "
                                "melting should be used (if false), assuming that the melt fraction only "
                                "depends on temperature and pressure, and how much melt has already been "
                                "generated at a given point, but not considering movement of melt in "
                                "the melting parameterization."
                                "\n\n"
                                "Note that melt does not freeze unless the 'Freezing rate' parameter is set "
                                "to a value larger than 0.");
                prm.declare_entry ("Freezing rate", "0.0",
                                Patterns::Double (0.),
                                "Freezing rate of melt when in subsolidus regions. "
                                "If this parameter is set to a number larger than 0.0, it specifies the "
                                "fraction of melt that will freeze per year (or per second, depending on the "
                                "``Use years in output instead of seconds'' parameter), as soon as the porosity "
                                "exceeds the equilibrium melt fraction, and the equilibrium melt fraction "
                                "falls below the depletion. In this case, melt will freeze according to the "
                                "given rate until one of those conditions is not fulfilled anymore. The "
                                "reasoning behind this is that there should not be more melt present than "
                                "the equilibrium melt fraction, as melt production decreases with increasing "
                                "depletion, but the freezing process of melt also reduces the depletion by "
                                "the same amount, and as soon as the depletion falls below the equilibrium "
                                "melt fraction, we expect that material should melt again (no matter how "
                                "much melt is present). This is quite a simplification and not a realistic "
                                "freezing parameterization, but without tracking the melt composition, there "
                                "is no way to compute freezing rates accurately. "
                                "If this parameter is set to zero, no freezing will occur. "
                                "Note that freezing can never be faster than determined by the "
                                "``Melting time scale for operator splitting''. The product of the "
                                "``Freezing rate'' and the ``Melting time scale for operator splitting'' "
                                "defines how fast freezing occurs with respect to melting (if the "
                                "product is 0.5, melting will occur twice as fast as freezing). "
                                "Units: 1/yr or 1/s, depending on the ``Use years "
                                "in output instead of seconds'' parameter.");
                prm.declare_entry ("Melting time scale for operator splitting", "1e3",
                                Patterns::Double (0.),
                                "Because the operator splitting scheme is used, the porosity field can not "
                                "be set to a new equilibrium melt fraction instantly, but the model has to "
                                "provide a melting time scale instead. This time scale defines how fast melting "
                                "happens, or more specifically, the parameter defines the time after which "
                                "the deviation of the porosity from the equilibrium melt fraction will be "
                                "reduced to a fraction of $1/e$. So if the melting time scale is small compared "
                                "to the time step size, the reaction will be so fast that the porosity is very "
                                "close to the equilibrium melt fraction after reactions are computed. Conversely, "
                                "if the melting time scale is large compared to the time step size, almost no "
                                "melting and freezing will occur."
                                "\n\n"
                                "Also note that the melting time scale has to be larger than or equal to the reaction "
                                "time step used in the operator splitting scheme, otherwise reactions can not be "
                                "computed. "
                                "Units: yr or s, depending on the ``Use years in output instead of seconds'' parameter.");
                prm.declare_entry ("Depletion solidus change", "200.0",
                                Patterns::Double (0.),
                                "The solidus temperature change for a depletion of 100\\%. For positive "
                                "values, the solidus gets increased for a positive peridotite field "
                                "(depletion) and lowered for a negative peridotite field (enrichment). "
                                "Scaling with depletion is linear. Only active when fractional melting "
                                "is used. "
                                "Units: \\si{\\kelvin}.");
                prm.declare_entry ("Reference permeability", "1e-8",
                                Patterns::Double(),
                                "Reference permeability of the solid host rock."
                                "Units: \\si{\\meter\\squared}.");
            }

            template <int dim>
            void
            mageminLookup<dim>::parse_parameters (ParameterHandler &prm)
            {
                // New Parameters 
                crust_min_P=prm.get_double ("Crust Minimum Pressure");
                crust_max_P=prm.get_double ("Crust Maximum Pressure");
                crust_min_T=prm.get_double ("Crust Minimum Temperature");
                crust_max_T=prm.get_double ("Crust Maximum Temperature");

                mantle_min_P=prm.get_double ("Mantle Minimum Pressure");
                mantle_max_P=prm.get_double ("Mantle Maximum Pressure");
                mantle_min_T=prm.get_double ("Mantle Minimum Temperature");
                mantle_max_T=prm.get_double ("Mantle Maximum Temperature");


                // Read my models Names
                data_directory_magemin                  = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory Lookup Table Magemin"));
                mantle_data_file_name      = prm.get ("Mantle data file name");
                crust_data_file_name     = prm.get ("Crust data file name");



                reference_rho_fluid        = prm.get_double ("Reference melt density");
                xi_0                       = prm.get_double ("Reference bulk viscosity");
                viscosity_fluid            = prm.get_double ("Reference melt viscosity");
                thermal_bulk_viscosity_exponent = prm.get_double ("Thermal bulk viscosity exponent");
                alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
                extraction_depth           = prm.get_double ("Melt extraction depth");
                melt_compressibility       = prm.get_double ("Melt compressibility");
                fractional_melting         = prm.get_bool ("Use fractional melting");
                freezing_rate              = prm.get_double ("Freezing rate");
                melting_time_scale         = prm.get_double ("Melting time scale for operator splitting");
                melt_bulk_modulus_derivative = prm.get_double ("Melt bulk modulus derivative");
                depletion_solidus_change   = prm.get_double ("Depletion solidus change");
                reference_permeability     = prm.get_double ("Reference permeability");


                if (this->convert_output_to_years() == true)
                {
                    melting_time_scale *= year_in_seconds;
                    freezing_rate /= year_in_seconds;
                }

                AssertThrow(melting_time_scale > 0,
                            ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));

                if (this->get_parameters().reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step)
                {
                    AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                                ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                        + " in the operator splitting scheme is too large to compute melting rates! "
                                        "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                        "operator splitting' chosen in the material model, which is currently "
                                        + Utilities::to_string(melting_time_scale) + "."));

                    AssertThrow(freezing_rate * this->get_parameters().reaction_time_step <= 1.0,
                                ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                        + " in the operator splitting scheme is too large to compute freezing rates! "
                                        "You have to choose it in such a way that it is smaller than the inverse of the "
                                        "'Freezing rate' chosen in the material model, which is currently "
                                        + Utilities::to_string(1.0/freezing_rate) + "."));
                }
            }
        

        }
    }
}
