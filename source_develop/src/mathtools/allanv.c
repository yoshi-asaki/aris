#include "mathtools.h"


void   allanv(int    nsample,            int    ndata,
              float  *sampling_interval, float  *variance,
              float  *y,                 float  parttm)
{
  int    i, j, N;
  float  ftmp, tau;

  for(i=1; i<=nsample; i++) {
    tau = parttm * (float)i;
    sampling_interval[i-1] = tau;

    ftmp = (float)0.0;
    N = ndata - 2*i;
    for (j=0; j<N; j++) {
      ftmp += pow((y[j+2*i] - 2.0*y[j+i] + y[j]), 2.0);
    }
    variance[i-1] = ftmp / pow(tau, 2.0) / (float)N / 2.0;
  }
}
