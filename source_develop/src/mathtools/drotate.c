#include "mathtools.h"

int    drotate(double *x, double e, char *axis)
{
  double u, v, w;

  u = x[0];
  v = x[1];
  w = x[2];

  if (*axis == 'x' || *axis == 'X') {
    x[1] = v*cos(e) - w*sin(e);
    x[2] = v*sin(e) + w*cos(e);
    return 1;
  }

  else if (*axis == 'y' || *axis == 'Y') {
    x[2] = w*cos(e) - u*sin(e);
    x[0] = w*sin(e) + u*cos(e);
    return 1;
  }

  else if (*axis == 'z' || *axis == 'Z') {
    x[0] = u*cos(e) - v*sin(e);
    x[1] = u*sin(e) + v*cos(e);
    return 1;
  }

  else {
    printf("drotate: bad flag to rotate %s\n", axis);
    return 0;
  }
}
