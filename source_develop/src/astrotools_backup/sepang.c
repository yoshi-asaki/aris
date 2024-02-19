#include <stdio.h>
#include <math.h>
#include "astrotools.h"

double sepang(double *a, double *b)
{
  double sa;

  sa = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]) / vlen3(a) / vlen3(b);
  if (sa > 1.0) {
    return (0.0);
  }
  return (acos(sa));
}
