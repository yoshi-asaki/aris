#include <stdio.h>
#include <math.h>
#include <cpgplot.h>
#include <mathtools.h>
#include <aris.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])


/**
#define __LANDSCAPE__
**/
#define __LANDSCAPE__
#ifndef __LANDSCAPE__
  #define __PORTRATE__
#endif

int   spectl_disp(int  NX, int  NY, int  nobs,
                  struct antenna_parameter *ant_prm, double *nu,
                  int  *TimUTC, struct st_observable *int_obs[])
{
  int    i, j, n, iobs, ndat;
  int    ix, iy;
  int    nmin, nmax;
  float  *pgx, *pgy;
  float  pgxmin, pgxmax;
  float  pgymin, pgymax;
  float  src_flag;
  float  f;
  char   p_title[50], y_lab[50];
  float  F2;
  double *wx;

/****
  float  PGY[86400];
  float  a1, b1, c1;
  float  RM;
****/

/*
-----------------------------
*/

#ifdef __LANDSCAPE__
  cpgpap(1.5*pgpap_prm, 0.4);
  cpgsch(1.4);
#elif defined __PORTRATE__
  cpgpap(0.75*pgpap_prm, 1.5);
  cpgsch(1.4);
#endif

/*
-----------------------------
*/

  if ((pgx = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for pgx.\n");
    return (__NG__);
  }
  if ((pgy = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for pgy.\n");
    free (pgx);
    return (__NG__);
  }

  F2 = 1.0 / (float)(nu[0] * nu[0]);

/*
-----------------------------
*/

  ndat = 0;
  ix = NX * nobs;
  iy = NY * nobs;
  src_flag = 100.0;
  for (iobs=0; iobs<nobs; iobs++) {
    if (weight_flag_chk(src_flag, &(int_obs[0][ix].wt),
                                      GRT_ELEVATION_LIMIT) == false &&
        weight_flag_chk(src_flag, &(int_obs[0][iy].wt),
                                      GRT_ELEVATION_LIMIT) == false) {
      pgx[ndat] = (float)(TimUTC[3]*3600 + TimUTC[4]*60 + TimUTC[5] + iobs);
      pgy[ndat] = diff(int_obs[0][ix].grp_delay, int_obs[0][iy].grp_delay)
           + F2 * diff(int_obs[0][ix].ion_delay, int_obs[0][iy].ion_delay);
      ndat++;
    }
    ix++;
    iy++;
  }
  fmaxmin(ndat, pgy, &pgymin, &pgymax, &nmin, &nmax);
  if (fabs(pgymin) > fabs(pgymax)) {
    pgymax = fabs(pgymin);
  } else {
    pgymax = fabs(pgymax);
  }
  pgymin = -pgymax;
  if (pgymin == pgymax) {
    pgymin -= 1.0;
    pgymax += 1.0;
  }
  pgxmin = (float)(TimUTC[3]*3600 + TimUTC[4]*60 + TimUTC[5]);
  pgxmax = (float)(TimUTC[3]*3600 + TimUTC[4]*60 + TimUTC[5] + nobs);

  sprintf(y_lab,  "Excess group delay ");
  f = ordinate_axis(ndat, pgy, &pgymin, &pgymax,  y_lab+strlen(y_lab));
  sprintf(y_lab+strlen(y_lab),  "]");
  sprintf(p_title, "Baseline: %s-%s", ant_prm[NX].IDC, ant_prm[NY].IDC);

#ifdef __LANDSCAPE__
  cpgsvp(0.05, 0.45, 0.15, 0.90);
#elif defined __PORTRATE__
/**
  cpgsvp(0.15, 0.90, 0.60, 0.98);
**/
  cpgsvp(0.15, 0.90, 0.70, 0.98);
#endif
  cpgswin(pgxmin, pgxmax, 1.2*pgymin, 1.2*pgymax);
  cpgsci(1);
  cpgtbox("BCSTNZHP", 0.0, 0, "BCNTS", 0.0, 10);
  cpglab("Time (UTC)", y_lab, p_title);
  cpgsci(1);
  cpgline(ndat, pgx, pgy);

/*
-----------------------------
*/

/****
  for (i=0; i<ndat-600; i++) {
    RM = 0.0;
    for (j=i; j<i+600; j++) {
      RM += pgy[j];
    }
    RM /= (float)600;
    PGY[i] = pgy[i+300] - RM;
  }
  ndat -= 600;

  lstsqr(ndat, pgx, PGY, &a1, &b1, &c1, 1);
  printf("#### Strandard Deviation ####  %f\n", sqrt(c1)/1.0e16);
****/

/*
-----------------------------
*/

/****
  cpgsci(2);
  cpgline(ndat, pgx, PGY);
  cpgsci(1);
****/

/*
-------------------------------
*/

  for (i=0; i<ndat; i++) {
    pgx[i] = 0.0;
    pgy[i] /= f;
/**
    pgy[i] /= 1.0e16;
    pgy[i] /= f;
    pgy[i] = PGY[i] / 1.0e16;
**/
  }

  n = 1;
  while (1) {
    n *= 2;
    if (n > ndat) {
      break;
    }
  }
  n *= 2;
  if ((wx = (double *)calloc(n, sizeof(double))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for wx.\n");
    return (__NG__);
  }
  for (i=0; i<ndat; i++) {
    wx[2*i  ] = pgy[i] / f;
    wx[2*i+1] = 0.0;
  }
  gsl_fft_complex_radix2_forward(wx, 1, n/2);
  for (i=1; i<n/4; i++) {
    pgy[i] =  log10(sqrt(wx[2*i]*wx[2*i]));
    pgx[i] =  log10((float)i / (float)(n/2));
  }
  free (wx);

#ifdef __LANDSCAPE__
  cpgsvp(0.55, 0.95, 0.15, 0.90);
#elif defined __PORTRATE__
/**
  cpgsvp(0.15, 0.90, 0.10, 0.48);
**/
  cpgsvp(0.15, 0.90, 0.30, 0.58);
#endif

  fmaxmin(n/4-1, pgy+1, &pgymin, &pgymax, &nmin, &nmax);
  cpgswin(log10(1.0/(float)(n/2)), log10(0.5), pgymin, pgymax);
  cpgbox("BCLNTS", 0, 0, "BCLNTS", 0, 0);

  cpglab("Frequency [Hz]", "Power [sec/Hz\\u1/2\\d]", "");
  cpgline(n/4-1, pgx+1, pgy+1);

/*
-------------------------------
*/

  free (pgx);
  free (pgy);

  return (__GO__);
}
