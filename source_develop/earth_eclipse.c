#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

int   earth_eclipse(double *s, double *P, double *R, double *L,
                    double sep_ang_lim)
{
  *L = s[0]*P[0] + s[1]*P[1] + s[2]*P[2];
  *R = sqrt(P[0]*P[0] + P[1]*P[1] + P[2]*P[2] - *L* *L);

  if (*R - fabs(*L) * tan(sep_ang_lim) <= earth_radius && *L < 0.0) {
    return NIGHT;
  } else {
    return DAY;
  }
}
