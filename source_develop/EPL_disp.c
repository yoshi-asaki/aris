#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <cpgplot.h>
#include <aris.h>

/****
#define __FILE_OUT__
****/

void EPL_disp(int NX, int NY,  int  nobs,
              struct antenna_parameter     *ant_prm, double *nu,
              int *TimUTC, struct st_observable *int_obs[])
{
  int    i, IX, IY, ix, iy, iobs, ns;
  int    I, J, k, l;
  int    timUTC[6];
  char   string[100];
  int    ndat, nmax, nmin;
  float  *tim;
  float  *epl[2];
  float  pgxmin,  pgxmax;
  float  pgymin1, pgymax1, pgymin2, pgymax2, pgymin, pgymax;
  double ELX, ELY, rho, z, Z, dz;
  float  src_flag;
  float  F2[SRC_NUM];

#ifdef __FILE_OUT__
  FILE   *log_fp;
  char   filename[20];
#endif /*__FILE_OUT__*/

  for (i=0; i<6; i++) {
    timUTC[i] = TimUTC[i];
  }

/*
----------------------------------------------------
*/

  IX = NX * nobs;
  IY = NY * nobs;

  if ((tim = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("EPL_disp: fail in calloc for tim.\n");
    exit (-1);
  }
  for (ns=0; ns<SRC_NUM; ns++) {
    if ((epl[ns] = (float *)calloc(nobs, sizeof(float))) == NULL) {
      printf("EPL_disp: fail in calloc for epl[%d].\n", ns);
      exit (-1);
    }
  }

  for (ns=0; ns<SRC_NUM; ns++) {
    F2[ns] = 1.0 / (float)(nu[ns] * nu[ns]);
  }

/*
----------------------------------------------------
*/

  ndat = 0;
  ix = NX * nobs;
  for (iobs=0; iobs<nobs; iobs++) {
    tim[ndat] = (float)(timUTC[3]*3600 + timUTC[4]*60 + timUTC[5] + iobs);
    for (ns=0; ns<SRC_NUM; ns++) {
      epl[ns][ndat] = int_obs[ns][ix].grp_delay
           + F2[ns] * int_obs[ns][ix].ion_delay;
      epl[ns][ndat] *= 1.0e+12;
    }
    ndat++;
    ix++;
  }

/*
----
*/

  pgxmin = (float)(timUTC[3]*3600 + timUTC[4]*60 + timUTC[5]);
  pgxmax = (float)(timUTC[3]*3600 + timUTC[4]*60 + timUTC[5] + nobs - 1);

  cpgpap(pgpap_prm, 1.0);
  cpgsch(0.8);

  fmaxmin(ndat, epl[0], &pgymin1, &pgymax1, &nmin, &nmax);
  fmaxmin(ndat, epl[1], &pgymin2, &pgymax2, &nmin, &nmax);

  if (fabs(pgymin1) <= fabs(pgymax1)) {
    pgymax1 = fabs(pgymax1);
  } else {
    pgymax1 = fabs(pgymin1);
  }
  if (fabs(pgymin2) <= fabs(pgymax2)) {
    pgymax2 = fabs(pgymax2);
  } else {
    pgymax2 = fabs(pgymin2);
  }
  if (pgymax1 >= pgymax2) {
    pgymax = pgymax1;
  } else {
    pgymax = pgymax2;
  }
  pgymin = -pgymax;
  if (pgymin == pgymax) {
    pgymin -= 1.0;
    pgymax += 1.0;
  }

  cpgsvp(0.2, 0.8, 0.74, 0.94);
  cpgswin(pgxmin, pgxmax, 1.2*pgymin, 1.2*pgymax);
  cpgsci(1);
  cpgtbox("BCSTNZHI", 0.0, 0, "BCNTS", 0.0, 10);
  sprintf(string, "Baseline: %s", ant_prm[NX].IDC);
  cpglab("", "", string);

  for (ns=0; ns<SRC_NUM; ns++) {
    src_flag = 100.0 * (float)(ns + 1);
    cpgsci(ns+1);
    ix = NX * nobs;
    iobs = 0;
    while (1) {
      l = 0;
      for (k=iobs; k<nobs; k++) {
        I = ix + k;
        if (! (int_obs[ns][I].wt >= src_flag)) {
          break;
        } else {
          l++;
        }
      }
      if (l != 0) {
        cpgline(l, tim+iobs, epl[ns]+iobs);
        iobs += l;
      } else {
        iobs++;
      }
      if (k >= nobs - 1) {
        break;
      }
    }
  }

#ifdef __FILE_OUT__
  if (NX < 9) {
    sprintf(filename, "aris_log/epl0%1d.dat", NX+1);
  } else {
    sprintf(filename, "aris_log/epl%2d.dat",  NX+1);
  }
  log_fp = fopen(filename, "w");
  fprintf(log_fp, "%d\n", ndat);
  for (i=0; i<ndat; i++) {
    fprintf(log_fp, "%f  %f  %f  %f\n", tim[i], epl[0][i], epl[1][i], epl[0][i]-epl[1][i]);
  }
  fclose (log_fp);
#endif /*__FILE_OUT__*/

/*
----------------------------------------------------------
*/

  ndat = 0;
  iy = NY * nobs;
  for (iobs=0; iobs<nobs; iobs++) {
    tim[ndat] = (float)(timUTC[3]*3600 + timUTC[4]*60 + timUTC[5] + iobs);
    for (ns=0; ns<SRC_NUM; ns++) {
      epl[ns][ndat] = int_obs[ns][iy].grp_delay
           + F2[ns] * int_obs[ns][iy].ion_delay;
      epl[ns][ndat] *= 1.0e+12;
    }
    ndat++;
    iy++;
  }

/*
----
*/

  fmaxmin(ndat, epl[0], &pgymin1, &pgymax1, &nmin, &nmax);
  fmaxmin(ndat, epl[1], &pgymin2, &pgymax2, &nmin, &nmax);

  if (fabs(pgymin1) <= fabs(pgymax1)) {
    pgymax1 = fabs(pgymax1);
  } else {
    pgymax1 = fabs(pgymin1);
  }
  if (fabs(pgymin2) <= fabs(pgymax2)) {
    pgymax2 = fabs(pgymax2);
  } else {
    pgymax2 = fabs(pgymin2);
  }
  if (pgymax1 >= pgymax2) {
    pgymax = pgymax1;
  } else {
    pgymax = pgymax2;
  }
  pgymin = -pgymax;
  if (pgymin == pgymax) {
    pgymin -= 1.0;
    pgymax += 1.0;
  }

  cpgsci(1);
  cpgsvp(0.2, 0.8, 0.42, 0.62);
  cpgswin(pgxmin, pgxmax, 1.2*pgymin, 1.2*pgymax);
  cpgtbox("BCSTNZHI", 0.0, 0, "BCNTS", 0.0, 0);
  sprintf(string, "Baseline: %s", ant_prm[NY].IDC);
  cpglab("", "Excess delay [ps]", string);

  for (ns=0; ns<SRC_NUM; ns++) {
    src_flag = 100.0 * (float)(ns + 1);
    cpgsci(ns+1);
    iy = NY * nobs;
    iobs = 0;
    while (1) {
      l = 0;
      for (k=iobs; k<nobs; k++) {
        J = iy + k;
        if (! (int_obs[ns][J].wt >= src_flag)) {
          break;
        } else {
          l++;
        }
      }
      if (l != 0) {
        cpgline(l, tim+iobs, epl[ns]+iobs);
        iobs += l;
      } else {
        iobs++;
      }
      if (k >= nobs - 1) {
        break;
      }
    }
  }

#ifdef __FILE_OUT__
  if (NY < 9) {
    sprintf(filename, "aris_log/epl0%1d.dat", NY+1);
  } else {
    sprintf(filename, "aris_log/epl%2d.dat",  NY+1);
  }
  log_fp = fopen(filename, "w");
  fprintf(log_fp, "%d\n", ndat);
  for (i=0; i<ndat; i++) {
    fprintf(log_fp, "%f  %f  %f  %f\n", tim[i], epl[0][i], epl[1][i], epl[0][i]-epl[1][i]);
  }
  fclose (log_fp);
#endif /*__FILE_OUT__*/

/*
----------------------------------------------------------
*/

  ndat = 0;
  ix = NX * nobs;
  iy = NY * nobs;
  for (iobs=0; iobs<nobs; iobs++) {
    tim[ndat] = (float)(timUTC[3]*3600 + timUTC[4]*60 + timUTC[5] + iobs);
    for (ns=0; ns<SRC_NUM; ns++) {
      epl[ns][ndat] =
                  diff(int_obs[ns][ix].grp_delay, int_obs[ns][iy].grp_delay)
       + F2[ns] * diff(int_obs[ns][ix].ion_delay, int_obs[ns][iy].ion_delay);
      epl[ns][ndat] *= 1.0e+12;
    }
    ndat++;
    ix++;
    iy++;
  }

  fmaxmin(ndat, epl[0], &pgymin1, &pgymax1, &nmin, &nmax);
  fmaxmin(ndat, epl[1], &pgymin2, &pgymax2, &nmin, &nmax);
  if (fabs(pgymin1) <= fabs(pgymax1)) {
    pgymax1 = fabs(pgymax1);
  } else {
    pgymax1 = fabs(pgymin1);
  }
  if (fabs(pgymin2) <= fabs(pgymax2)) {
    pgymax2 = fabs(pgymax2);
  } else {
    pgymax2 = fabs(pgymin2);
  }
  if (pgymax1 >= pgymax2) {
    pgymax = pgymax1;
  } else {
    pgymax = pgymax2;
  }
  pgymin = -pgymax;
  if (pgymin == pgymax) {
    pgymin -= 1.0;
    pgymax += 1.0;
  }

  cpgsci(1);
  cpgsvp(0.2, 0.8, 0.10, 0.30);
  cpgswin(pgxmin, pgxmax, 1.2*pgymin, 1.2*pgymax);
  cpgtbox("BCSTNZHI", 0.0, 0, "BCNTS", 0.0, 0);
  sprintf(string, "Baseline %s - %s", ant_prm[NX].IDC, ant_prm[NY].IDC);
  cpglab("Time (UT)", "", string);

  for (ns=0; ns<SRC_NUM; ns++) {
    src_flag = 100.0 * (float)(ns + 1);
    cpgsci(ns+1);
    ix = NX * nobs;
    iy = NY * nobs;
    iobs = 0;
    while (1) {
      l = 0;
      for (k=iobs; k<nobs; k++) {
        I = ix + k;
        J = iy + k;
        if (int_obs[ns][I].wt < src_flag || int_obs[ns][J].wt < src_flag) {
          break;
        } else {
          l++;
        }
      }
      if (l != 0) {
        cpgline(l, tim+iobs, epl[ns]+iobs);
        iobs += l;
      } else {
        iobs++;
      }
      if (k >= nobs - 1) {
        break;
      }
    }
  }

/*
----------------------------------------------------------
*/

  free (tim);
  for (ns=0; ns<SRC_NUM; ns++) {
    free (epl[ns]);
  }

  return;
}
