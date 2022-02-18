#include <stdio.h>
#include <stdlib.h>
#include <cpgplot.h>
#include "./pgtools.h"


int   pgarea_set(float *pgx, float *pgy, float *cx, float *cy)
{
  char   string[5];
  float  x[5], y[5];

  cpgband(2, 1, pgx[0], pgy[0], &pgx[0], &pgy[0], string);
  pgx[1] = pgx[0];
  pgy[1] = pgy[0];
  cpgband(2, 1, pgx[0], pgy[0], &pgx[1], &pgy[1], string);
  if (pgx[0] < pgx[1]) {
    cx[0] = pgx[0];
    cx[1] = pgx[1];
  } else {
    cx[0] = pgx[1];
    cx[1] = pgx[0];
  }
  if (pgy[0] < pgy[1]) {
    cy[0] = pgy[0];
    cy[1] = pgy[1];
  } else {
    cy[0] = pgy[1];
    cy[1] = pgy[0];
  }

  x[0] = cx[0];
  x[1] = cx[0];
  x[2] = cx[1];
  x[3] = cx[1];
  x[4] = cx[0];

  y[0] = cy[0];
  y[1] = cy[1];
  y[2] = cy[1];
  y[3] = cy[0];
  y[4] = cy[0];
  cpgline(5, x, y);

  return 1;
}
