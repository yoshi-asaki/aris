#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"


void  earth_rotation(int    *TimUTC,  double UT1_UTC,
                     double *s0,      double *s1)
{
  double s_tmp[3];

  s_tmp[0] = s0[0];
  s_tmp[1] = s0[1];
  s_tmp[2] = s0[2];

  drotate(s_tmp, GST(TimUTC, 0.0, UT1_UTC), "z");
  s1[0] = s_tmp[0];
  s1[1] = s_tmp[1];
  s1[2] = s_tmp[2];

  return;
}
