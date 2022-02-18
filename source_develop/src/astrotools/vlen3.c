#include <stdio.h>
#include <math.h>
#include "astrotools.h"

double vlen3(double *v)
{
  return (sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
}
