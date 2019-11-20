/***********************************************************************
/
/ Set default parameter values
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

int grackle_verbose = FALSE;

chemistry_data *grackle_data = NULL;
chemistry_data_storage grackle_rates;

chemistry_data _set_default_chemistry_parameters(void);

int set_default_chemistry_parameters(chemistry_data *my_grackle)
{
  *my_grackle = _set_default_chemistry_parameters();
  grackle_data = my_grackle;
  return SUCCESS;
}

chemistry_data _set_default_chemistry_parameters(void)
{

  chemistry_data my_chemistry;

  my_chemistry.Gamma                          = 5./3.;
  my_chemistry.use_grackle                    = FALSE;
  my_chemistry.with_radiative_cooling         = TRUE;
  my_chemistry.primordial_chemistry           = 0;
  my_chemistry.dust_chemistry                 = 0;
  my_chemistry.metal_cooling                  = FALSE;
  my_chemistry.h2_on_dust                     = FALSE;
  my_chemistry.use_dust_density_field         = FALSE;

  my_chemistry.cmb_temperature_floor          = TRUE;
  my_chemistry.grackle_data_file              = "";

  my_chemistry.three_body_rate                = 0;
  my_chemistry.cie_cooling                    = 0;
  my_chemistry.h2_optical_depth_approximation = 0;

  my_chemistry.photoelectric_heating          = 0;
  // epsilon=0.05, G_0=1.7 (in erg s^-1 cm^-3)
  my_chemistry.photoelectric_heating_rate     = 8.5e-26;
  my_chemistry.interstellar_radiation_field   = 1.7;

  my_chemistry.use_volumetric_heating_rate    = 0;
  my_chemistry.use_specific_heating_rate      = 0;

  my_chemistry.UVbackground                   = 0;

  my_chemistry.UVbackground_redshift_on      = FLOAT_UNDEFINED;
  my_chemistry.UVbackground_redshift_off     = FLOAT_UNDEFINED;
  my_chemistry.UVbackground_redshift_fullon  = FLOAT_UNDEFINED;
  my_chemistry.UVbackground_redshift_drop    = FLOAT_UNDEFINED;

  my_chemistry.Compton_xray_heating   = 0;

  my_chemistry.LWbackground_intensity = 0.0;   // [in units of 10^21 erg/s/cm^2/Hz/sr]
  my_chemistry.LWbackground_sawtooth_suppression = 0;

  my_chemistry.HydrogenFractionByMass       = 0.76;
  /* The DToHRatio is by mass in the code, so multiply by 2. */
  my_chemistry.DeuteriumToHydrogenRatio     = 2.0*3.4e-5; // Burles & Tytler 1998

  /*
     Up until version 2.1, the solar metal mass fraction was 0.02041.  
     This is close to 0.0194 of Anders & Grevesse (1989), but significantly 
     higher than the more recent value of 0.0122 from Asplund et al. (2005).
     As of version 2.1, the solar metal mass fraction has been set to 0.01295, 
     which is consistent with the abundances used in Cloudy when generating the 
     cooling tables.
  */
  my_chemistry.SolarMetalFractionByMass     = 0.01295; // Cloudy v13 abundances

  /*
    The dust to gas ratio in local molecular clouds.
     Table 2 from Pollack et al. (1994).
  */
  my_chemistry.local_dust_to_gas_ratio      = 0.009387;

  my_chemistry.NumberOfTemperatureBins      = 600;
  my_chemistry.ih2co                        = 1;
  my_chemistry.ipiht                        = 1;
  my_chemistry.TemperatureStart             = 1.0;
  my_chemistry.TemperatureEnd               = 1.0e9;
  my_chemistry.CaseBRecombination           = 0; // default to case A rates
  my_chemistry.NumberOfDustTemperatureBins  = 250;
  my_chemistry.DustTemperatureStart         = 1.0;
  my_chemistry.DustTemperatureEnd           = 1500.0;

  my_chemistry.cloudy_electron_fraction_factor = 9.153959e-3; // Cloudy 07.02 abundances

  /* radiative transfer parameters */
  my_chemistry.use_radiative_transfer                 = 0;
  my_chemistry.radiative_transfer_coupled_rate_solver = 0;
  my_chemistry.radiative_transfer_intermediate_step   = 0;
  my_chemistry.radiative_transfer_hydrogen_only       = 0;

  /* approximate self-shielding */
  my_chemistry.self_shielding_method                  = 0;
  my_chemistry.H2_self_shielding                      = 0;

  my_chemistry.use_dust_evol       = 0;
  my_chemistry.SolarAbundances[0]  = 0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction . Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
  my_chemistry.SolarAbundances[1]  = 3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3)
  my_chemistry.SolarAbundances[2]  = 1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3)
  my_chemistry.SolarAbundances[3]  = 8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3)
  my_chemistry.SolarAbundances[4]  = 2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3)
  my_chemistry.SolarAbundances[5]  = 9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4)
  my_chemistry.SolarAbundances[6]  = 1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4)
  my_chemistry.SolarAbundances[7]  = 6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4)
  my_chemistry.SolarAbundances[8]  = 1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4)
  my_chemistry.SolarAbundances[9] = 1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3)

//number of OpenMP threads
# ifdef _OPENMP
  my_chemistry.omp_nthreads = omp_get_max_threads(); // maximum allowed number
# endif

  return my_chemistry;
}
