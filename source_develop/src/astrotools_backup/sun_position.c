#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"

/****
#define __DEBUG__
****/

void sun_position(int  *TimUT, double UT1_UTC, double *RA, double *DC)
{
  double e[3], Q[4], A[3][3], s[3];
  double lambda, lambda_dash, rs, q, t;
  double ti, alpha, sun_vec[3];
  double dpi = 3.141592653589793238462643;
  double DPI;

/*
--------
*/

  DPI = dpi / 180.0;

  t = ET(TimUT, UT1_UTC);
  lambda_dash =279.0358 + 360.00769*t
        + (1.9159 - 0.00005*t) * sin(DPI*(356.531 + 359.991*t))
        + 0.0200               * sin(DPI*(353.06  + 719.981*t))
        - 0.0048               * sin(DPI*(248.64  -  19.341*t))
        + 0.0020               * sin(DPI*(285.0   +  329.64*t))
        + 0.0018               * sin(DPI*(334.2   - 4452.67*t))

        + 0.0018               * sin(DPI*(293.7   -    0.20*t))
        + 0.0015               * sin(DPI*(242.4   +  450.37*t))
        + 0.0013               * sin(DPI*(211.1   +  225.18*t))
        + 0.0008               * sin(DPI*(208.0   +  659.29*t))
        + 0.0007               * sin(DPI*( 53.5   +   90.38*t))

        + 0.0007               * sin(DPI*( 12.1   -   30.35*t))
        + 0.0006               * sin(DPI*(239.1   +  337.18*t))
        + 0.0005               * sin(DPI*( 10.1   -    1.50*t))
        + 0.0005               * sin(DPI*( 99.1   -   22.81*t))
        + 0.0004               * sin(DPI*(264.8   +  315.56*t))

        + 0.0004               * sin(DPI*(233.8   +  299.30*t))
        - 0.0004               * sin(DPI*(198.1   +  720.02*t))
        + 0.0003               * sin(DPI*(349.6   + 1079.97*t))
        + 0.0003               * sin(DPI*(241.2   -   44.43*t));
  lambda_dash *= DPI;

  q = (-0.007261 + 0.0000002*t) * cos(DPI*(356.53 + 359.991*t))+0.000030
      - 0.000091                * cos(DPI*(353.1  + 719.98 *t))
      + 0.000013                * cos(DPI*(205.8  +4452.67 *t))
      + 0.000007                * cos(DPI*( 62.0  + 450.40 *t))
      + 0.000007                * cos(DPI*(105.0  + 329.6  *t));

  rs = pow((double)10.0, q);
  lambda = lambda_dash + DPI*0.0057;

  s[0] = cos(lambda_dash);
  s[1] = sin(lambda_dash);
  s[2] = 0.0;
  e[0] = 1.0;
  e[1] = 0.0;
  e[2] = 0.0;
  coordinate_rotation(e, EPSIRON(t), Q, A[0], s, sun_vec);

  *RA = atan2(sun_vec[1], sun_vec[0]);
  *DC = atan2(sun_vec[2], vlen2(sun_vec));

  return;
}
