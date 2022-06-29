#include "mathtools.h"

void cross_product3(double *a1, double *a2, double *A)
{
  double  t;

  A[0] = a1[1]*a2[2] - a1[2]*a2[1];
  A[1] = a1[2]*a2[0] - a1[0]*a2[2];
  A[2] = a1[0]*a2[1] - a1[1]*a2[0];

  return;
}
