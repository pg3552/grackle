#include <stdio.h>
#include <math.h>

#include "grackle_macros.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#include "dust_evol.h"
/* Routines for dust growth and destruction */
/* Qi Li, 06/17/2019 */



int initialize_dust_data(chemistry_data *my_chemistry){
	my_chemistry->SolarAbundances[0]=0.0134; // Asplund (2009)
	my_chemistry->SolarAbundances[1]=0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
	my_chemistry->SolarAbundances[2]=3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3)
	my_chemistry->SolarAbundances[3]=1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3)
	my_chemistry->SolarAbundances[4]=8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3)
	my_chemistry->SolarAbundances[5]=2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3)
	my_chemistry->SolarAbundances[6]=9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4)
	my_chemistry->SolarAbundances[7]=1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4)
	my_chemistry->SolarAbundances[8]=6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4)
	my_chemistry->SolarAbundances[9]=1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4)
	my_chemistry->SolarAbundances[10]=1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3)
}
