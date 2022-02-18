#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include "./pgtools.h"

void pg_skyplot_trajectory(int ndat, float *az, float *el, int pgclr, int pgmrk)
{
  int    idat;
  float  a;
  float  pgx, pgy;
  double DPI, dpi=3.141592653589793238462643;

/*
----
*/

  DPI = dpi / 180.0;
  cpgsci(pgclr);

  for (idat=0; idat<ndat; idat++) {
    a = 1.0 - 2.0 * (float)el[idat] / (float)dpi;
    pgx = a * cos(az[idat]);
    pgy = a * sin(az[idat]);
    if (el[idat] >= 10.0 * DPI) {
      cpgpt(1, &pgx, &pgy, pgmrk);
    }
  }
}
