#include <math.h>
#include "mathtools.h"

void  dvector_calc(double M[][3], double *a)
{
  double b[3];

  b[0] = M[0][0]*a[0] + M[0][1]*a[1] + M[0][2]*a[2];
  b[1] = M[1][0]*a[0] + M[1][1]*a[1] + M[1][2]*a[2];
  b[2] = M[2][0]*a[0] + M[2][1]*a[1] + M[2][2]*a[2];

  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];

  return;
}
