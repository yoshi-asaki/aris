#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <aris.h>


void    spherical_geometry(double H,    double  el, double *z,
                           double *rho, double *EL, double *Z,
                           double *dz)
{
  *z   = 0.5 * dpi - el;
  *Z   = asin(earth_radius * sin(*z) / (earth_radius + H));
  *dz  = *z - *Z;
  *EL  = 0.5 * dpi - *Z;
  *rho = (earth_radius + H) * *dz;

  return;
}
