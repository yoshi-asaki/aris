#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <aris.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])

void  wfft2d(float *mapr, float *mapi, float *weight, int  fmax, int nflag)
{
  int    i, j;
  int    I, J, pnt_ad;
  int    IREF, JREF;
  int    LEN1, LEN2;
  float  wt;
  double *wx;

/*
--------
*/

  wx = (double *)calloc(2*fmax, sizeof(double));

  IREF = fmax / 2;
  JREF = fmax / 2;

  LEN1 = sizeof(float) * fmax;
  LEN2 = sizeof(float) * JREF;
  for (i=0; i<fmax; i++) {
    I = i * fmax;
    for (j=0; j<fmax; j++) {
      wt = *(weight + i*fmax + j);
      wx[2*j  ] = wt * mapr[I+j];
      wx[2*j+1] = wt * mapi[I+j];
    }
    if (nflag == 1) {
      gsl_fft_complex_radix2_backward(wx, 1, fmax);
    } else {
      gsl_fft_complex_radix2_inverse(wx, 1, fmax);
    }
    for (j=0; j<JREF; j++) {
      mapr[I+JREF+j] = wx[2*j  ];
      mapi[I+JREF+j] = wx[2*j+1];
      mapr[I     +j] = wx[2*(j+JREF)  ];
      mapi[I     +j] = wx[2*(j+JREF)+1];
    }
  }

  for (j=0; j<fmax; j++) {
    for (i=0; i<fmax; i++) {
      pnt_ad = i*fmax + j;
      wx[2*i  ] = *(mapr + pnt_ad);
      wx[2*i+1] = *(mapi + pnt_ad);
    }
    gsl_fft_complex_radix2_inverse(wx, 1, fmax);
    for (i=0; i<IREF; i++) {
      pnt_ad = (i+IREF)*fmax + j;
      *(mapr + pnt_ad) = wx[2*i  ];
      *(mapi + pnt_ad) = wx[2*i+1];
    }
    for (i=IREF; i<fmax; i++) {
      pnt_ad = (i-IREF)*fmax + j;
      *(mapr + pnt_ad) = wx[2*i  ];
      *(mapi + pnt_ad) = wx[2*i+1];
    }
  }

  free (wx);
}
