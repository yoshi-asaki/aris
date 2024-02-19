#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"

/*
----  MJD = JD - 2400000.5 ----
*/

double MJD(int  Y, int  M, int  D, int  h, int  m, int  s, double d_sec)
{
  int    YP, MP, B;
  double mjd;

  YP = Y;
  MP = M;

  if (M <= 2) {
    M += 12;
    Y--;
  }

  if (YP <  1582 ||              /* until 4 Oct. 1582 (Julian Calendarn) */
     (YP == 1582 && MP < 10) ||
     (YP == 1582 && MP == 10 && D <= 4)) {
    B = -2 + (int)((Y + 4716)/4) - 1179;
  } else if (YP >  1582 ||       /* from 10 Oct. 1582 (Gregorian Calendarn) */
            (YP == 1582 && MP > 10) ||
            (YP == 1582 && MP == 10 && D >= 10)) {
    B = (int)(Y/400) - (int)(Y/100) + (int)(Y/4);
  }

  mjd = 365.0*(double)Y - 679004.0 + (double)B
      + floor(30.6001 * (double)(M+1)) + (double)D;
  mjd += ((double)(3600*h + 60*m + s) + d_sec) / 86400.0;


/****
  mjd = floor(365.25 * (double)Y) + floor((double)Y / 400.0)
      - floor((double)Y / 100.0) + floor(30.59 * (double)(M - 2)) + (double)D
      - 678912.0;
  mjd += ((double)(3600*h + 60*m + s) + d_sec) / 86400.0;
****/

  return mjd;
}
