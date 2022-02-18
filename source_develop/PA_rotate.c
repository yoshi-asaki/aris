#include <stdio.h>
#include <math.h>
#include <aris.h>


void PA_rotate(int n, float *x, float *y, double theta)
{
  int    i;
  double xr, yr, ar, br;

  xr = cos (theta);
  yr = sin (theta);
  for (i=0; i<n; i++) {
    ar = x[i];
    br = y[i];
    x[i] = ar * xr - br * yr;
    y[i] = ar * yr + br * xr;
  }
}
