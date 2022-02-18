#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <aris.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])


void  fft2d(float *mapr, float *mapi, int  fmax, int nflag)
{
  int    i, j;
  int    IREF, JREF, LEN1, LEN2;
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
    for (j=0; j<fmax; j++) {
      wx[2*j  ] = mapr[i*fmax+j];
      wx[2*j+1] = mapi[i*fmax+j];
    }
    gsl_fft_complex_radix2_forward(wx, 1, fmax);
    for (j=0; j<JREF; j++) {
      mapr[i*fmax+j+JREF] = wx[2*j  ];
      mapi[i*fmax+j+JREF] = wx[2*j+1];
      mapr[i*fmax+j     ] = wx[2*(j+JREF)  ];
      mapi[i*fmax+j     ] = wx[2*(j+JREF)+1];
    }
  }

/*
----------------------------------
*/

  for (j=0; j<fmax; j++) {
    for (i=0; i<fmax; i++) {
      wx[2*i  ] = *(mapr + i*fmax + j);
      wx[2*i+1] = *(mapi + i*fmax + j);
    }
    gsl_fft_complex_radix2_forward(wx, 1, fmax);
    for (i=0; i<IREF; i++) {
      *(mapr + (i+IREF)*fmax + j) = wx[2*i  ];
      *(mapi + (i+IREF)*fmax + j) = wx[2*i+1];
    }
    for (i=IREF; i<fmax; i++) {
      *(mapr + (i-IREF)*fmax + j) = wx[2*i];
      *(mapi + (i-IREF)*fmax + j) = wx[2*i+1];
    }
  }

  free (wx);
}
