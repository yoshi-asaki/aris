#include "mathtools.h"

double dvec2scalar(double *t, int    dim)
{
  int    i;
  double y = 0.0;

  for (i=0; i<dim; i++) {
    y += pow(t[i], 2.0);
  }

  return (sqrt(y));
}
