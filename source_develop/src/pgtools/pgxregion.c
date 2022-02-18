#include <stdio.h>
#include <stdlib.h>
#include <cpgplot.h>
#include "./pgtools.h"


int   pgxregion(float pgxmin, float pgxmax, float *pgx, float *pgy)
{
  char   c;
  float  ftmp;

  cpgband(4, 1, 0.5*(pgxmin+pgxmax), 0.0, pgx,   pgy,   &c);
  cpgband(4, 1, pgx[0],              0.0, pgx+1, pgy+1, &c);

  if (pgx[0] > pgx[1]) {
    ftmp = pgx[0];
    pgx[0] = pgx[1];
    pgx[1] = ftmp;
  }

  return 1;
}
