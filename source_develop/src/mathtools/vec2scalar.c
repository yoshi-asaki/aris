#include "mathtools.h"

float  vec2scalar(float  *t, int    dim)
{
  int    i;
  float  y = 0.0;

  for (i=0; i<dim; i++) {
    y += pow(t[i], 2.0);
  }

  return (sqrt(y));
}
