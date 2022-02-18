#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <aris.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])

void   Gaussian_Noise1D(int N, double sigma, double *gn)
{
  int    i, NDIV, M;
  double fai, amp;
  double *wx;

  M = 1;
  while (M < N) {
    M *= 2;
  }
  if ((wx = (double *)calloc(2*M, sizeof(double))) == NULL) {
    printf("fail in Gaussian_Noise1D for memory allocation of wx.\n");;
    exit (-1);
  }

  for (i=0; i<M; i++) {
    fai = dpi * random_val1();
    wx[2*i  ] = (double)cos(fai);
    wx[2*i+1] = (double)sin(fai);
  }

  gsl_fft_complex_radix2_forward(wx, 1, M);
  amp = 1.0 / sqrt((double)(M/2)) * sigma;
  for (i=0; i<N; i++) {
    gn[i] = amp * wx[i];
  }

  free (wx);
}
