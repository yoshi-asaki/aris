#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"

void   B1950toJ2000(double RA_B1950, double DC_B1950,
                    double *RA_J2k,  double *DC_J2k)
{
  int    i;
  double s[3], p[3][3], x[3], eterm[3];
  double ec, cdel;
  double dpi = 3.141592653589793238462643;

/* Following coefficients are brought from JAPANESE EPHEMERIS 1992. */

  p[0][0] =  0.9999256782;
  p[0][1] = -0.011182061;
  p[0][2] = -0.0048579477;

  p[1][0] =  0.0111820609;
  p[1][1] =  0.9999374784;
  p[1][2] = -0.0000271765;

  p[2][0] =  0.0048579479;
  p[2][1] = -0.0000271474;
  p[2][2] =  0.9999881997;

  eterm[0] = -1.62557e-6;
  eterm[1] = -3.1919e-7;
  eterm[2] = -1.3843e-7;

  s[0] = cos(DC_B1950) * cos(RA_B1950);
  s[1] = cos(DC_B1950) * sin(RA_B1950);
  s[2] = sin(DC_B1950);

  ec = 0.0;
  ec += eterm[0] * s[0];
  ec += eterm[1] * s[1];
  ec += eterm[2] * s[2];

  s[0] = (1.0 + ec) * s[0] - eterm[0];
  s[1] = (1.0 + ec) * s[1] - eterm[1];
  s[2] = (1.0 + ec) * s[2] - eterm[2];

  for (i=0; i<3; i++) {
    x[i] = 0.0;
    x[i] += p[i][0] * s[0];
    x[i] += p[i][1] * s[1];
    x[i] += p[i][2] * s[2];
  }

  cdel = vlen2(x);
  *RA_J2k = atan2(x[1], x[0]);
  *DC_J2k = atan2(x[2], cdel);

  return;
}
