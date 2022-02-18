#include "mathtools.h"

void   dmaxmin(int  ndat, double *y, double *ymin, double *ymax,
                                   int   *nmin, int   *nmax)
{
  int    i;
  double tmp;

  *ymin = y[0];
  *ymax = y[0];
  *nmin = 0;
  *nmax = 0;

  for(i=1; i<ndat; i++) {
    if (y[i] > *ymax) {
      *ymax = y[i];
      *nmax = i;
    }
    if (y[i] < *ymin) {
      *ymin = y[i];
      *nmin = i;
    }
  }
  return;
}
