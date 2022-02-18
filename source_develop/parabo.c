#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <aris.h>


void   parabo(float *xmax,    float *ymax,
              float x_center, float delta_x, float *y,
              float *a, float *b, float *c)
{
  float  u, v, w; 

/*
--------
*/

  u = y[2] -y[0];
  v = y[0] +y[2] -(float)2.0 *y[1];
  w = u / v / (float)2.0;
  *xmax = x_center - w * delta_x;
  *ymax = y[1] - u*u / 8.0 / v;

  *a = v / 2.0 / delta_x / delta_x;
  *b = u / 2.0 / delta_x - 2.0 * *a * x_center;
  *c = y[1] - *a * *a * x_center * x_center - *b * x_center;

  return;
}
