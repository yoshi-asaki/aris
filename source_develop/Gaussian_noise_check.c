#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])

void  Gaussian_noise_check()
{
  int    i, j, k, I, J;
  int    N1 = 200, N2 = 50;
  int    NMTRX;
  double *gnoise;
  double sigma = 0.5;
  double amp;
  float  *pgx, *pgy;
  float  dx;

  NMTRX = NOISE_MTRX;

  cpgsubp(2, 2);

/*
------------------------------------------------
*/

  gnoise = (double *)calloc(NMTRX, sizeof(double));
  pgx    = (float  *)calloc(NMTRX, sizeof(float));
  pgy    = (float  *)calloc(NMTRX, sizeof(float));

  for (I=0; I<N2; I++) {
    pgx[I] = (float)I;

    pgy[I] = 0.0;
    for (J=0; J<N1; J++) {
      Gaussian_Noise1D(NMTRX, sigma, gnoise);
      for (i=I; i<NMTRX; i++) {
        pgy[I] += gnoise[i] * gnoise[i-I];
      }
    }
    pgy[I] /= (float)(N1 * (NMTRX - I));
    pgy[I] = sqrt(fabsf(pgy[I]));
    printf("data:  %d   %f \n", I, pgy[I]);
  }

  cpgenv(0.0, (float)N2, 0.0, 1.0, 0, 0);
  cpgline(N2, pgx, pgy);

/*
------------------------------------------------
*/

  Gaussian_Noise1D(NMTRX, sigma, gnoise);
  for (i=NMTRX/2-1; i>=0; i--) {
    gnoise[2*i  ] = gnoise[i];
    gnoise[2*i+1] = 0.0;
  }
  gsl_fft_complex_radix2_inverse(gnoise, 1, NMTRX/2);
  for (i=0; i<NMTRX/2; i++) {
    pgy[i] = sqrt(gnoise[2*i]*gnoise[2*i] + gnoise[2*i+1]*gnoise[2*i+1]);
    pgx[i] = (float)i;
  }

  cpgsci(1);
  cpgenv(0.0, (float)(NMTRX/2), 0.0, 1.0e-1, 0, 0);
  cpgsci(2);
  cpgline(NMTRX/2, pgx, pgy);

/*
------------------------------------------------
*/

  Gaussian_Noise1D(NMTRX, sigma, gnoise);
  for (i=0; i<NMTRX; i++) {
    pgx[i] = (float)i;
    pgy[i] = (float)gnoise[i];
  }

  cpgsci(1);
  cpgenv(0.0, (float)NMTRX, -4.0, 4.0, 0, 0);
  cpgline(NMTRX, pgx, pgy);

/*
------------------------------------------------
*/

  dx = 0.10;
  for (i=0; i<N2; i++) {
    pgy[i] = 0.0;
  }
  for (J=0; J<N1; J++) {
    Gaussian_Noise1D(NMTRX, sigma, gnoise);
    for (i=0; i<N2; i++) {
      for (j=0; j<NMTRX; j++) {
        if (gnoise[j] >= dx * ((float)(i - N2/2 - 1) + 0.5) &&
            gnoise[j] <  dx * ((float)(i - N2/2    ) + 0.5)) {
          pgy[i] += 1.0;
        }
      }
    }
  }
  for (i=0; i<N2; i++) {
    pgx[i] = dx * (float)(i - N2/2);
    pgy[i] /= (float)(NMTRX * N1);
    pgy[i] /= dx;
  }

  amp = 1.0 / sqrt(2.0 * dpi) / sigma;
  cpgsci(1);
  cpgenv(pgx[0], pgx[N2-1], 0.0, 1.2*amp, 0, 0);
  cpgline(N2, pgx, pgy);

  for (i=0; i<N2; i++) {
    pgy[i] = amp * exp(-pow(pgx[i], 2.0) / 2.0 / pow(sigma, 2.0));
  }
  cpgsls(4);
  cpgsci(2);
  cpgline(N2, pgx, pgy);

/*
------------------------------------------------
*/

  free (gnoise);
  free (pgx);
  free (pgy);
  cpgsubp(1, 1);
}
