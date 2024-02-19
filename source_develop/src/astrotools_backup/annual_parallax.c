#include <stdio.h>
#include <math.h>
#include "astrotools.h"

void  annual_parallax(double *ap, double RA, double DC,
                      int *TimUT, double UT1_UTC)
{
  double vector[3];
  double e[3], s[3], q[4], A[3][3];
  double dpi = 3.141592653589793238462643;

  heirocentric_equatorial_rectangular_coordinates_earth_position
         (TimUT, UT1_UTC, vector);
  e[0] = 0.0;
  e[1] = 0.0;
  e[2] = 1.0;
  coordinate_rotation(e, 3.0/2.0*dpi-RA, q, A[0], vector, s);
  e[0] = 1.0;
  e[1] = 0.0;
  e[2] = 0.0;
  coordinate_rotation(e, DC-dpi/2.0,     q, A[0], s, vector);

/****
  ap[0] = -1.0 * vector[0] / D;
  ap[1] = -1.0 * vector[1] / D;
****/

  ap[0] = -vector[0];
  ap[1] = -vector[1];
}
