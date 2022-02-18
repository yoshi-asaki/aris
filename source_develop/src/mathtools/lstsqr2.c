/* LeaST SQuaRe fitting                          */
/*                                               */
/*                                               */

#include "mathtools.h"

void   lstsqr2(int n,
               float *x, float *y, float *xe, float *ye,
               float *A, float *a, float *variance)
{
  int    i;
  float  xm, ym, sx, sy, ey;
  float  weight;

  xm = 0.0; ym = 0.0;
  sx = 0.0; sy = 0.0;
  for (i=0; i<n; i++) {
    weight = 1.0 / pow(ye[i], 2.0);
    xm += x[i] *weight;
    ym += y[i] *weight;
    sx += x[i] * x[i] *weight;
    sy += x[i] * y[i] *weight;
  }
  xm /= (float)n;
  ym /= (float)n;
  sx /= (float)n;
  sy /= (float)n;
  A[0] = (sy -xm*ym) / (sx -xm*xm);
  A[1] = (ym*sx -xm*sy) / (sx -xm*xm);


  xm = 0.0; ym = 0.0;
  sx = 0.0; sy = 0.0;
  ey = 0.0;
  for (i=0; i<n; i++) {
    weight = 1.0 / pow(ye[i]+fabs(A[0])*xe[i], 2.0);
    xm += x[i] /weight;
    ym += y[i] /weight;
    sx += x[i] * x[i] /weight;
    sy += x[i] * y[i] /weight;
    ey += weight;
  }
  xm /= (float)n;
  ym /= (float)n;
  sx /= (float)n;
  sy /= (float)n;
  A[0] = (sy -xm*ym) / (sx -xm*xm);
  A[1] = (ym*sx -xm*sy) / (sx -xm*xm);
  a[0] = weight/ (sx*ey-xm*xm);
  a[1] = sx    / (sx*ey-xm*xm);

  for (i=0; i<n; i++) {
    variance[0] += pow((y[i] -A[0]*x[i] -A[1]), 2.0);
    variance[1] += pow((y[i] -A[0]*x[i] -A[1]), 2.0) /
                   pow(ye[i]+fabs(A[0])*xe[i], 2.0);
  }
  variance[0] /= (float)(n -1);
  variance[1] /= (float)(n -1);

}
