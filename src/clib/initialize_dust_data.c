#include <stdio.h>
#include <math.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "phys_constants.h"
#include "dust_evol.h"
/* Routines for dust growth and destruction */
/* Qi Li, 06/17/2019 */



int initialize_dust_data(chemistry_data *my_chemistry);
	my_chemistry->SolarAbundances[0]=0.0134; // Asplund (2009)
	SolarAbundances[1]=0.28;    // He  (10.93 in units where log[H]=12, so photospheric mass fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse Y=0.2485, X=0.7314)
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


/* Density must be in [CGS], if not, covert the unit before calling routines here*/
/* Mass must be in code units */
/* Timescales are in code units */

/* calculate characteristic timescales used for subcycling*/
void dust_tau_for_timestep(gr_float *Density, gr_float *Tgas,
		gr_float *Mass, gr_float *dust_Mass,
		gr_float (*Metallicity)[11], gr_float *SNe_ThisTimeStep,
		double *tau_accr, double *tau_dest, double *tau_sput, int i, int j, int k,
		double *SolarAbundances,
		double time_units, double dtime){ //CHECK: can dtime be passed easily?

	int n;
	gr_float Mgas;
	gr_float tau_accr0;
	Mgas = Mass[i][j][k] - dust_Mass[i][j][k];

	/* growth */
	tau_accr0 = tau_ref * (DUST_GROWTH_DENSREF/Density[i][j][k]) * pow(DUST_GROWTH_TREF/Tgas[i][j][k],0.5);
	tau_accr[i][j][k] = INFINITY;
	for (n=1;n<NUM_METAL_SPECIES;n++)
		tau_accr[i][j][k] = fmin(tau_accr0 * (SolarAbundances[n] / Metallicity[n]),tau_accr[i][j][k]);

	/* destruction by SN shocks*/
	if (SNe_ThisTimeStep<=0.)
		tau_dest[i][j][k] = INFINITY;
	else 
		tau_dest[i][j][k] = Mgas/(Ms100*SNe_ThisTimeStep[i][j][k]*GALSF_DUST_DESTRUCTION_EFF) * dtime;

	/* destruction by thermal sputtering */
	tau_sput[i][j][k] = 1.7e8 * SEC_PER_YEAR / time_units
					  * (GALSF_DUST_GRAINSIZE/0.1) * (1e-27/Density[i][j][k])
					  * (pow(2e6/Tgas[i][j][k],2.5)+1.0); // Tsai & Mathews (1995)
}

/* main routine for dust growth and destruction */
void dust_growth_and_destruction(gr_float *Density, gr_float *Tgas,
		gr_float (*Metallicity)[11], gr_float (*dust_Metallicity)[11],
		gr_float *Mass, gr_float *dust_Mass, gr_float *SNe_ThisTimeStep,
		double *tau_dest, double *tau_sput, int i, int j, int k,
		double Ms100, double tau_ref, double *SolarAbundances, double dtime){
		
	int n;
	gr_float Mgas;
	gr_float dM[NUM_METAL_SPECIES]={0.},dMs, //dMs: dust mass loss via shock and sputtering
	gr_float Mmet[NUM_METAL_SPECIES]={0.},Mdust[NUM_METAL_SPECIES];
	gr_float tau_accr0,tau_accr;
	
	Mgas = Mass[i][j][k] - dust_Mass[i][j][k];

	for (n=0;n<NUM_METAL_SPECIES;n++){ 
		Mdust[n] = dust_Metallicity[n] * dust_Mass[i][j][k];
		dM[n] = 0.;
	}

	/* metals accreted from gas to dust */
	tau_accr0 = tau_ref * (DUST_GROWTH_DENSREF/Density[i][j][k]) * pow(DUST_GROWTH_TREF/Tgas[i][j][k],0.5);
	for (n=1;n<NUM_METAL_SPECIES;n++){
		Mmet[n] = Metallicity[n] * Mgas + Mdust[n];
		tau_accr = tau_accr0 * (SolarAbundances[n] / Metallicity[n]);
		dM[n] += fmin((1 - Mdust[n] / Mmet[n]) * (Mdust[n]/tau_accr) * dtime_sub, Mmet[n] - Mdust[n]); //growth capped by gas in a single cell
		dM[0] += dM[n];
	}

	/* SNe shock and thermal sputtering */
	if (SNe_ThisTimeStep[i][j][k] <= 0.0)
		dMs = 0.0;
	else
		dMs = fmin(dust_Mass[i][j][k]/tau_dest[i][j][k]*dtime,dust_Mass[i][j][k]);
	if (dMs >= dust_Mass[i][j][k]){
		printf("WARNING: dMs >= Mdust; SNe shock destruction:%f %g\n",SNe_ThisTimeStep[i][j][k],tau_dest[i][j][k]);	
	}
	else{
		dMs += dust_Mass[i][j][k] / tau_sput[i][j][k] * 3.0 * dtime;
		dMs = fmin(dMs,dust_Mass[i][j][k]);
	}
	for (n = 0; n < NUM_METAL_SPECIES; n++)
		dM[n] -= dust_Metallicity[n] * dMs;

	/* recalculate metallicity */
	dust_Mass[i][j][k] += dM[0];
	for (n=0;n<NUM_METAL_SPECIES;n++){
		dust_Metallicity[n] = dust_Mass[i][j][k] > 0 ?  (dM[n] + Mdust[n]) / dust_Mass[i][j][k] : 0.;  // SNe shocks on neighbors are ignored
		Metallicity[n] = (Metallicity[n] * Mgas - dM[n])/(Mgas - dM[0]);
	}
} // end of dust growth and destruction module
