#include "mathtools.h"


double    inner_product3(double *a, double *b)
{
  double ipro;

  ipro = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

  return ipro;
}
