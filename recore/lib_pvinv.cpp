#include <math.h>
#include <cmath>
#include <limits>
#include "lib_pvinv.h"



partload_inverter_t::partload_inverter_t( )
{
	Paco = Pdco = std::numeric_limits<double>::quiet_NaN();
}

bool partload_inverter_t::acpower(
	/* inputs */
	double Pdc,     /* Input power to inverter (Wdc) */

	/* outputs */
	double *Pac,    /* AC output power (Wac) */
	double *Plr,    /* Part load ratio (Pdc_in/Pdc_rated, 0..1) */
	double *Eff	    /* Conversion efficiency (0..1) */
	)
{

	// linear interpolation based on Pdc/Pdco and *Partload and *Efficiency arrays


	// clipping loss
	if ( *Pac > Paco ) *Pac = Paco;

	*Plr = Pdc / Pdco;
	*Eff = *Pac / Pdc;
	if ( *Eff < 0.0 ) *Eff = 0.0;

	return true;
}



