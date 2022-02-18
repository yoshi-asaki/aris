#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


float  ordinate_axis(int ndat, float *y, float *pgmin, float *pgmax,
                     char *string)
{
  int    i;
  float  f;

  if (*pgmax >= 0.1) {
    sprintf(string, "[s");
    f = 1.0;
  } else if (*pgmax >= 1.0e-4 && *pgmax < 0.1) {
    sprintf(string, "[milli s");
    f = 1.0e3;
  } else if (*pgmax >= 1.0e-7 && *pgmax < 1.0e-4) {
    sprintf(string, "[micro s");
    f = 1.0e6;
  } else if (*pgmax >= 1.0e-10 && *pgmax < 1.0e-7) {
    sprintf(string, "[ns");
    f = 1.0e9;
  } else if (*pgmax >= 1.0e-13 && *pgmax < 1.0e-10) {
    sprintf(string, "[ps");
    f = 1.0e12;
  } else if (*pgmax >= 1.0e-16 && *pgmax < 1.0e-13) {
    sprintf(string, "[fs");
    f = 1.0e15;
  } else {
    sprintf(string, "[s");
    f = 1.0;
  }

  *pgmin *= f;
  *pgmax *= f;
  for (i=0; i<ndat; i++) {
    y[i] *= f;
  }

  return (f);
}
