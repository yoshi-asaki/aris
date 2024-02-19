#include <stdio.h>
#include <math.h>
#include "astrotools.h"


double LST(int  *TimUTC,  double dT, double  UT1_UTC, double Lambda)
{
  double theta;

  theta = GST(TimUTC, dT, UT1_UTC) + Lambda;

  return (theta);
}
