#include <stdio.h>
#include <math.h>
#include "astrotools.h"

double vlen2(double *v)
{
  return (sqrt(v[0]*v[0] + v[1]*v[1]));
}
