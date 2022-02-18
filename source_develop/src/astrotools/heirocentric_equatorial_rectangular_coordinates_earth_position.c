#include <stdio.h>
#include <math.h>
#include "astrotools.h"

void  heirocentric_equatorial_rectangular_coordinates_earth_position
         (int *TimUT, double UT1_UTC, double *earth_vec)
{
  double RA, DC;
  double R = 1.49597870e+11;

/****************************/
  R = 1.0;
/****************************/

  sun_position(TimUT, UT1_UTC, &RA, &DC);
  earth_vec[0] = cos(DC) * cos(RA);
  earth_vec[1] = cos(DC) * sin(RA);
  earth_vec[2] = sin(DC);

  earth_vec[0] *= -R;
  earth_vec[1] *= -R;
  earth_vec[2] *= -R;

  return;
}
