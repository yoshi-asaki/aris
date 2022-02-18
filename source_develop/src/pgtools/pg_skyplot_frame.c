#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cpgplot.h>
#include <math.h>
#include "./pgtools.h"

void pg_skyplot_frame()
{
  int    i, j;
  float  pgx[361], pgy[361];
  float  a, b;
  double DPI, dpi=3.141592653589793238462643;

/*
--------
*/

  DPI = dpi / 180.0;

  cpgpap(5.0, 1.0);
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);
  cpgsch(1.2);
  cpgscf(2);
  cpgsvp(0.10, 0.90, 0.10, 0.90);
  cpgswin(-1.0, 1.0, -1.0, 1.0);
  cpgsci(1);

/*
--------
*/

  cpgsls(1);
  for (i=0; i<3; i++) {
    b = (float)(i + 1) * 30.0 / 90.0;
    for (j=0; j<361; j++) {
      a = (float)j * DPI;
      pgx[j] = b * cos(a);
      pgy[j] = b * sin(a);
    }
    cpgline(361, pgx, pgy);
  }

  cpgsls(2);
  b = 10.0 / 90.0;
  for (i=0; i<361; i++) {
    a = (float)i * DPI;
    pgx[i] = b * cos(a);
    pgy[i] = b * sin(a);
  }
  cpgline(361, pgx, pgy);
  b = 80.0 / 90.0;
  for (i=0; i<361; i++) {
    a = (float)i * DPI;
    pgx[i] = b * cos(a);
    pgy[i] = b * sin(a);
  }
  cpgline(361, pgx, pgy);
  cpgsls(1);
  pgx[0] = -1.0;
  pgx[1] =  1.0;
  pgy[0] =  0.0;
  pgy[1] =  0.0;
  cpgline(2, pgx, pgy);
  pgx[0] =  0.0;
  pgx[1] =  0.0;
  pgy[0] = -1.0;
  pgy[1] =  1.0;
  cpgline(2, pgx, pgy);
  cpgtext(0.0, 1.0, "N");
  cpgtext(1.0, 0.0, "E");
  cpgtext(0.0, -1.00,  "0");
  cpgtext(0.0, -0.67, "30");
  cpgtext(0.0, -0.33, "60");
  cpgtext(0.0,  0.00, "Zenith");

  return;
}
