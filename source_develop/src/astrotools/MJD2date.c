#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"

/*
---- MJD = JD - 2400000.5 ----
*/

int  MJD2date(double mjd,
              int    *Y,    int    *M,    int    *D,
              int    *h,    int    *m,    int    *s)
{
  int    i, MFLG, W, F;
  int    MONTH[2][12] = {
            31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
            31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  *Y = (int)(0.0027379093 * mjd + 1858.877);
  *D = (int)(mjd - MJD(*Y, 1, 0, 0, 0, 0, 0.0));

  if ((*Y % 4 == 0 && *Y % 400 == 0) || (*Y % 4 == 0 && *Y % 100 != 0)) {
    MFLG = 0;
  } else {
    MFLG = 1;
  }

  for (i=0; i<12; i++) {
    *D -= MONTH[MFLG][i];
    if (*D <= 0) {
      *D += MONTH[MFLG][i];
      *M = ++i;
      break;
    }
  }
  W = ((int)mjd-5) % 7;


  mjd -= MJD(*Y, *M, *D, 0, 0, 0, 0.0);
  mjd *= 86400.0;
  F = (int)mjd;
  if (fabs(mjd-floor(mjd)) >= 0.5) {
    F++;
  }
  *h = F / 3600;
  *m = (F % 3600) / 60;
  *s = (F % 3600) % 60;

  return W;

/****
  double q;
  int    a, b, c, d, e, f, F, W;

  a = (int)mjd + 2400001;
  q = mjd - (int)mjd;

  if (a < 2299161) {
    b = 0;
    c = a + 1524;
  } else {
    b = (int)(((double)a - (double)1867216.25) / (double)36524.25); 
    c = a + b - (int)(b/4) + 1525;
  }

  d = (int)(((double)c - (double)121.1) / (double)365.25);
  e = (int)((double)365.25 * (double)d);
  f = (int)(((double)c - (double)e) / (double)30.6001);

  *D = c - e - (int)(30.6001 * (double)f);
  *M = f - 1 - 12*(int)(f/14);
  *Y = d - 4715 - (int)((7 + *M) / 10);

  q *= 86400.0;
  F = (int)rint(q);
  *h = F / 3600;
  *m = (F % 3600) / 60;
  *s = (F % 3600) % 60;

  W = ((int)mjd-5) % 7;

  return W;
****/
}
