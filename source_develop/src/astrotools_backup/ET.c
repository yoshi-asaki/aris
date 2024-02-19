#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"


double ET(int *TimUT, double UT1_UTC)
{
  double t;

  t = (MJD(TimUT[0], TimUT[1], TimUT[2], TimUT[3], TimUT[4], TimUT[5], UT1_UTC)
     - MJD(1975, 1, 0, 0, 0, 0, 0.0)) / 365.25;
  t += (0.0317*t + 1.43)*1.0e-6;

  return t;
}
