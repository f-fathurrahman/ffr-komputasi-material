/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
*/


// Modified by ffr

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include "system.h"

double bessel(int *n, double *x)
{
  return gsl_sf_bessel_Jn(*n, *x);
}

double bessel_(int *n, double *x)
{
  return gsl_sf_bessel_Jn(*n, *x);
}

double besseli(int *n, double *x)
{
  return gsl_sf_bessel_In(*n, *x);
}

double besseli_(int *n, double *x)
{
  return gsl_sf_bessel_In(*n, *x);
}
