#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <aris.h>

#include <gsl/gsl_multifit_nlin.h>

/****
#include <Crecipes.h>
#include <nrutil.h>
****/

/****
#define __Sgr_A_star__
****/
#define __Sgr_A_star__
#ifndef __Sgr_A_star__
  #define __M87__
#endif

void gaussian_floar(float , float [], float *, float [], int );

#define  NBH        242
#define  MDL_RSLT   0.250e-3

#define  ALAMDA_LIM 1.0e-10

int   bhs_model(char  *ch_bhs_file, float *dist, int dmax, int fmax,
                double *pix_uvl, double *pix_mas, double wave_length,
                double uv_max, double *uv_factor)
{
  FILE   *bhsm_fp;
  int    i, j, k, l;
  int    N_OFF, M, M2;
  float  fdum1, fdum2, fdum3;
  char   string[100];
  int    hdmax;
  int    scale_factor;
  float  scale_factor_sq;
  float  sf;
  float  *bhm;

  float  x[NBH], y[NBH], sig[NBH];
  int    ma = 4;
  int    ia[10];
  float  chisq[10], dyda[10];
  float  A[10], B[10];
  float  A_tmp[10], B_tmp[10];
  float  alamda, alamda_tmp;
  float  **alpha, **covar;
  float  r;
  int    ndat;
  double data_pix_scale, t_scale;

  int    k1, k2, l1, l2;
  float  *ddd;

  hdmax = dmax / 2;

#ifdef __Sgr_A_star__
  data_pix_scale = 4.8 * (double)MDL_RSLT;   /* micro arcsec */
#elif defined __M87__
  data_pix_scale = 2.0 * (double)MDL_RSLT;   /* micro arcsec */
#endif

/*
-----------------------------------------
*/

  for (i=0; i<dmax*dmax; i++) {
    dist[i] =  0.0;
  }

  if ((bhm = (float *)calloc(NBH*NBH, sizeof(float))) == NULL) {
    printf("ERROR: data array for Blackhole model cannot be allocated.\n");
    return (__NG__);
  }

  if ((bhsm_fp=fopen(ch_bhs_file, "r")) == NULL) {
    printf("ERROR: data file cannot be found.\n");
    free (bhm);
    return (__NG__);
  }

  for (i=0; i<NBH; i++) {
    for (j=0; j<NBH; j++) {
      if (fgets(string, sizeof(string), bhsm_fp) == NULL) {
        printf("ERROR: BHS_MODEL: Invalid input.\n");
        free (bhm);
        return (__NG__);
      }
      sscanf(string, "%f %f %f\n", &fdum1, &fdum2, &fdum3);
      bhm[NBH*i + j] = fdum3;
    }
    if (fgets(string, sizeof(string), bhsm_fp) == NULL) {
      printf("ERROR: BHS_MODEL: Invalid input.\n");
      free (bhm);
      return (__NG__);
    }
  }
  fclose (bhsm_fp);

/*
----------------------------------------------
*/

/********
  sf = *pix_mas / data_pix_scale;
  aaa = 0.5 * log(sf) / log(2.0);
********/

/****
  if (sf > 1.0 && sf <= 4.0) {
****/
    *pix_mas = data_pix_scale;
    pixel_calc(pix_uvl, pix_mas, wave_length, uv_max, uv_factor, fmax, -1);
    scale_factor = 1;
/****
  }
****/



  if (scale_factor == 0) {
    printf("ERROR: Model resolution is not enough. \n");
    exit (0);
  }

  scale_factor_sq = (float)(scale_factor * scale_factor);

  scale_factor_sq *= 600.0;

/*
----------------------------------------------
*/

  M     =  NBH / 2 / scale_factor;
  N_OFF = (NBH / 2) % M;
  M2    = 2 * M;

  for (i=0; i<M2; i++) {
    for (j=0; j<M2; j++) {
      fdum1 = 0.0;
      for (k=0; k<scale_factor; k++) {
        for (l=0; l<scale_factor; l++) {
          fdum1 += bhm[NBH * (N_OFF + i*scale_factor + k)
                           + (N_OFF + j*scale_factor + l)];
        }
      }
      dist[dmax*(hdmax-M+i) + hdmax-M+j] = fdum1 / scale_factor_sq;
    }
  }

/*
------------------------------
*/

/**
  alpha = matrix(1, ma, 1, ma);
  covar = matrix(1, ma, 1, ma);
**/
  N_OFF = 50;

  fdum1 = 0.0;
  for (i=0; i<NBH; i++) {
    for (j=0; j<NBH; j++) {
      if (bhm[NBH*i + j] > fdum1) {
        fdum1 = bhm[NBH*i + j];
      }
    }
  }
  fdum1 *= 0.3;

  ndat = 1;
  for (i=0; i<NBH/2-N_OFF; i++) {
    x[ndat]   = (float)(i - NBH/2);
    y[ndat]   = bhm[NBH*i + NBH/2];
    sig[ndat] = 0.2 * y[1];
    ndat++;
  }
  for (i=NBH/2+N_OFF; i<NBH; i++) {
    x[ndat] = (float)(i - NBH/2);
    y[ndat] = bhm[NBH*i + NBH/2];
    sig[ndat] = 0.2 * y[1];
    ndat++;
  }
  ndat -= 1;

  A[1] = fdum1;
  A[2] = -log(y[1] / A[1]) / x[1] / x[1];
  A[3] = 0.0;
  A[4] = 0.20 * y[1];
  alamda = -1.0;
  alamda_tmp = 0.1;

/****
  for (i=1; i<=ma; i++) {
    A_tmp[i] = A[i];
  }
  while (1) {
    mrqmin(x, y, sig, ndat, A_tmp, ia, ma, covar, alpha, chisq,
           gaussian_floar, &alamda);
    if (alamda > alamda_tmp || alamda <= ALAMDA_LIM) {
      break;
    } else {
      alamda_tmp = alamda;
      for (i=1; i<=ma; i++) {
        A[i] = A_tmp[i];
      }
    }
  }
****/

  ndat = 1;
  for (i=0; i<NBH/2-N_OFF; i++) {
    x[ndat]   = (float)(i - NBH/2);
    y[ndat]   = bhm[NBH*NBH/2 + i];
    sig[ndat] = 0.2 * y[1];
    ndat++;
  }
  for (i=NBH/2+N_OFF; i<NBH; i++) {
    x[ndat] = (float)(i - NBH/2);
    y[ndat]   = bhm[NBH*NBH/2 + i];
    sig[ndat] = 0.2 * y[1];
    ndat++;
  }
  ndat -= 1;

  B[1] = fdum1;
  B[2] = -log(y[1] / B[1]) / x[1] / x[1];
  B[3] = 0.0;
  B[4] = 0.20 * y[1];
  B[4] = 0.0;
  alamda = -1.0;
  alamda_tmp = 0.1;

/****
  for (i=1; i<=ma; i++) {
    B_tmp[i] = B[i];
  }
  while (1) {
    mrqmin(x, y, sig, ndat, A_tmp, ia, ma, covar, alpha, chisq,
           gaussian_floar, &alamda);
    if (alamda > alamda_tmp || alamda <= ALAMDA_LIM) {
      break;
    } else {
      alamda_tmp = alamda;
      for (i=1; i<=ma; i++) {
        B[i] = B_tmp[i];
      }
    }
  }
****/

/**
  free_matrix(alpha, 1, ma, 1, ma);
  free_matrix(covar, 1, ma, 1, ma);
**/

  A[1] = sqrt(A[1] * B[1]);
  A[4] = 0.5 * (A[4] + B[4]);

/****
  for (i=0; i<dmax; i++) {
    if (i < hdmax - M || i >= hdmax + M) {
      for (j=0; j<dmax; j++) {
          r = pow((float)((i-hdmax)*scale_factor)-A[3], 2.0) * A[2]
            + pow((float)((j-hdmax)*scale_factor)-B[3], 2.0) * B[2];
          dist[dmax*i + j] = A[1] * exp(-r) + A[4];
      }
    }
  }

  for (j=0; j<dmax; j++) {
    if (j < hdmax - M || j >= hdmax + M) {
      for (i=0; i<dmax; i++) {
          r = pow((float)((i-hdmax)*scale_factor)-A[3], 2.0) * A[2]
            + pow((float)((j-hdmax)*scale_factor)-B[3], 2.0) * B[2];
          dist[dmax*i + j] = A[1] * exp(-r) + A[4];
      }
    }
  }
****/

/*
----------------
*/

  free (bhm);

/*
----------------
*/

  fdum1 = 0.0;
  for (i=0; i<dmax; i++) {
    for (j=0; j<dmax; j++) {
      fdum1 += dist[dmax*i + j];
    }
  }

  for (i=0; i<dmax; i++) {
    for (j=0; j<dmax; j++) {
      dist[dmax*i + j] *= (3.0 / fdum1);
    }
  }

/*
----------------
*/

/****
  ddd = calloc(dmax*dmax, sizeof(float));
  memcpy(ddd, dist, sizeof(float) * dmax * dmax);

-- 100GHZ --
  sf = -log(2.0) * sqrt(*pix_mas / 0.050);
-- 350GHZ --
  sf = -log(2.0) * sqrt(*pix_mas / 0.005);
  sf = -log(2.0) * sqrt(*pix_mas / 0.005);

  for (i=0; i<dmax; i++) {
    for (j=0; j<dmax; j++) {

      dist[dmax*i+j] = 0.0;
      k1 = i - 50;
      if (k1 < 0) k1 = 0;
      l1 = j - 50;
      if (l1 < 0) l1 = 0;
      k2 = i + 50;
      if (k2 > dmax) k2 = dmax;
      l2 = j + 50;
      if (l2 > dmax) l2 = dmax;

      for (k=k1; k<k2; k++) {
        for (l=l1; l<l2; l++) {
          dist[dmax*i+j] += ddd[dmax*k+l]
                *exp(sf*(float)((k-i)*(k-i)+(l-j)*(l-j)));
        }
      }
    }
    printf("%d  %d\n", i, j);
  }
  sf = 0.0;
  for (i=0; i<dmax; i++) {
    for (j=0; j<dmax; j++) {
      sf += dist[dmax*i+j];
    }
  }
  sf = 3.0 / sf;
  for (i=0; i<dmax; i++) {
    for (j=0; j<dmax; j++) {
      dist[dmax*i+j] *= sf;
    }
  }
  free (ddd);
****/

/*
----------------
*/

  return (__GO__);
}


void gaussian_floar(float x, float a[], float *y, float dyda[], int na)
{
  float E;

  E = exp(-a[2]*(x-a[3])*(x-a[3]));

  *y = a[1] * E + a[4];

  dyda[1] = E;
  dyda[2] =  -(x-a[3]) * (x-a[3]) * *y;
  dyda[3] = 2.0 * a[2] * (x-a[3]) * *y;
  dyda[4] = 0.0;

  return;
}
