#include <stdio.h>
#include <math.h>
#include <cpgplot.h>
#include <mathtools.h>
#include <aris.h>

/**
#define __LANDSCAPE__
**/
#define __LANDSCAPE__
#ifndef __LANDSCAPE__
  #define __PORTRATE__
#endif

/****
#define __FILE_OUT__
#define __DEBUG__
****/
#define __DEBUG__


int   allanv_disp(int  NX, int  NY, int  nobs,
                  struct antenna_parameter *ant_prm, double *nu,
                  int  *TimUTC, struct st_observable *int_obs[])
{
  int    i, j, n, iobs, ndat;
  int    ix, iy;
  int    nmin, nmax;
  float  *pgx, *pgy, *pgz;
  float  pgxmin, pgxmax;
  float  pgymin, pgymax;
  float  src_flag;
  float  f;
  float  F2;
  char   p_title[50], y_lab[50];

#ifdef __FILE_OUT__
  FILE   *log_fp;
#endif /* __FILE_OUT__ */

/*
-----------------------------
*/

  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);

#ifdef __LANDSCAPE__
  cpgpap(1.5*pgpap_prm, 0.4);
  cpgsch(1.0);
#elif defined __PORTRATE__
/****
  cpgpap(0.50*pgpap_prm, 1.0);
****/
  cpgpap(1.00*pgpap_prm, 1.0);
  cpgsch(1.0);
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
  if ((pgz = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for pgz.\n");
    free (pgx);
    free (pgy);
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
      pgy[ndat] = int_obs[0][ix].grp_delay;
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

  sprintf(y_lab,  "Excess delay ");
  f = ordinate_axis(ndat, pgy, &pgymin, &pgymax, y_lab+strlen(y_lab));
  sprintf(y_lab+strlen(y_lab),  "]");
  sprintf(p_title, "Baseline: %s-%s", ant_prm[NX].IDC, ant_prm[NY].IDC);



#ifdef __LANDSCAPE__
  cpgsvp(0.05, 0.45, 0.15, 0.90);
#elif defined __PORTRATE__
  cpgsvp(0.20, 0.80, 0.60, 0.90);
#endif
  cpgswin(pgxmin, pgxmax, 1.2*pgymin, 1.2*pgymax);
  cpgsci(1);
  cpgtbox("BCSTNZHP", 0.0, 0, "BCNTS", 0.0, 0);
/********
  cpglab("Time (UTC)", y_lab, p_title);
********/
  cpglab("Time (UTC)", y_lab, "");
  cpgsci(1);
  cpgline(ndat, pgx, pgy);

/*
-----------------------------
*/

#ifdef __FILE_OUT__
  if ((log_fp = fopen("aris_log/delay_err.log", "w")) == NULL) {
    printf("Warning; ORBIT_DISP: ./aris_log/delay_err.log cannot be made.\n");
  } else {
    fprintf(log_fp, "time [s]  %s\n", y_lab);
    for (i=0; i<ndat; i++) {
      fprintf(log_fp, "%6d,   %12.8f\n", iobs, pgy[i]);
    }
    fclose (log_fp);
  }
#endif /* __FILE_OUT__ */

/*
-----------------------------
*/

  for (i=0; i<ndat; i++) {
    pgy[i] /= f;
  }
  allanv(ndat/2, ndat, pgx, pgz, pgy, 1.0);
  n = ndat / 2 - 1;
  n--;

  for (i=0; i<n; i++) {
    pgx[i] = log10(pgx[i]);
    pgy[i] = log10(sqrt(pgz[i]));
  }
  pgxmin = log10(1.0);
  pgxmax = log10(1.0e5);
  if (pgx[n-1] > log10(1.0e5)) {
    pgxmax = pgy[n-1];
  }

#ifdef __LANDSCAPE__
  cpgsvp(0.55, 0.95, 0.15, 0.90);
#elif defined __PORTRATE__
  cpgsvp(0.20, 0.80, 0.15, 0.45);
#endif
  fmaxmin(n, pgy, &pgymin, &pgymax, &nmin, &nmax);
  if (pgymax < pgymin+0.75*log10(4.0e3)) {
    pgymax = pgymin + 0.75*log10(4.0e3);
  }
  pgymin -= 0.3;
  pgymax += 0.3;

  cpgswin(pgxmin, pgxmax, pgymin, pgymax);
  cpgbox("BCLNTS", 0, 0, "BCLNTS", 0, 0);
  cpglab("Sampling interval [s]", "Allan standard deviation [s/s]", "");
  cpgsci(1);
  cpgline(n, pgx, pgy);

#ifdef __FILE_OUT__
  if ((log_fp = fopen("aris_log/allanv.log", "w")) == NULL) {
    printf("Warning; ORBIT_DISP: ./aris_log/allanv.log cannot be made.\n");
  } else {
    fprintf(log_fp, "Sampling interval [s]    Stability [s/s]\n");
    fprintf(log_fp, "(log scale)              (log scale)    \n");
    for (i=0; i<n; i++) {
      fprintf(log_fp, "%12f,   %12.8f\n", pgx[i], pgy[i]);
    }
    fclose (log_fp);
  }
#endif /* __FILE_OUT__ */

/*
-----------------------------
*/

#ifdef __DEBUG__
  for (i=0; i<n; i++) {
    pgx[i] = pow(1.2, (double)i);
    pgy[i] = (float)asd_model((double)pgx[i], CSO_10);

    pgx[i] = log10(pgx[i]);
    pgy[i] = log10(pgy[i]);
  }
  cpgsci(2);
  cpgline(n, pgx, pgy);
  cpgsci(1);

  for (i=0; i<n; i++) {
    pgx[i] = pow(1.2, (double)i);
    pgy[i] = (float)asd_model((double)pgx[i], CSO_100);

    pgx[i] = log10(pgx[i]);
    pgy[i] = log10(pgy[i]);
  }
  cpgsci(3);
  cpgline(n, pgx, pgy);
  cpgsci(1);

  for (i=0; i<n; i++) {
    pgx[i] = pow(1.2, (double)i);
    pgy[i] = (float)asd_model((double)pgx[i], H_M);

    pgx[i] = log10(pgx[i]);
    pgy[i] = log10(pgy[i]);
  }
  cpgsci(4);
  cpgline(n, pgx, pgy);
  cpgsci(1);

  for (i=0; i<n; i++) {
    pgx[i] = pow(1.2, (double)i);
    pgy[i] = (float)asd_model((double)pgx[i], TRP3);

    pgx[i] = log10(pgx[i]);
    pgy[i] = log10(pgy[i]);
  }
  cpgsci(5);
  cpgline(n, pgx, pgy);
  cpgsci(1);
#endif

/*
-----------------------------
*/

  free (pgx);
  free (pgy);
  free (pgz);

  return (__GO__);
}
