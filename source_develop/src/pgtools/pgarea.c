#include <stdio.h>
#include <stdlib.h>
#include <cpgplot.h>
#include "./pgtools.h"


int   pgarea(float pgcenter_x, float pgcenter_y, float *cx, float *cy)
{
  int    i;
  char   string[10];
  float  x[5], y[5];

  for (i=0; i<2; i++) {
    cx[i] = pgcenter_x;
    cy[i] = pgcenter_y;
  }
  cpgband(2, 1, pgcenter_x, pgcenter_y, &cx[0], &cy[0], string);
  cpgband(2, 1, cx[0],      cy[0],      &cx[1], &cy[1], string);
  if (cx[0] > cx[1]) {
    x[0]  = cx[0];
    cx[0] = cx[1];
    cx[1] = x[0];
  }
  if (cy[0] > cy[1]) {
    y[0]  = cy[0];
    cy[0] = cy[1];
    cy[1] = y[0];
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
