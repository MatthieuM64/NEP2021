/*LIBRARY FOR RANDOM NUMBER GENERATOR*/

#ifndef DEF_RANDOM_MANGEAT_CPP
#define DEF_RANDOM_MANGEAT_CPP

#include <gsl/gsl_randist.h>

gsl_rng *GSL_r;

void init_gsl_ran()
{
	GSL_r=gsl_rng_alloc(gsl_rng_mt19937);
}

double ran()
{
	double r=gsl_rng_uniform(GSL_r);
	if (r!=0)
	{
		return r;
	}
	else
	{
		return ran();
	}
}

double gaussian(const double &sig)
{
	double phi=ran()*2*M_PI;
	double psi=ran();
	return sqrt(-2*sig*log(psi))*cos(phi);
}

#endif


