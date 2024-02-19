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
#define __RMS_PHASE__
****/
#ifndef __RMS_PHASE__
  #define __STD__DVTN__
#endif


int   SSF_disp(int    ANT_NUM,
               int    GRT_NUM,
               int    BGN_ANT_I, int END_ANT_I,
               int    BGN_ANT_J, int END_ANT_J,
               int    nobs,      int nfrq,
               int    *TimUT,    double UT1_UTC,
               double *nu,
               struct antenna_parameter *ant_prm,
               struct fringe            *frng[],
               float  *fringe_weight[],
               struct st_observable *int_obs[],
               char   *ascii_out,
               int    PROC_MODE)
{
  int    i, j, k, iobs, ifrq, ns, NS;
  int    ibase, nbase[2];
  int    I, IF;
  int    ndat;
  int    iant, jant;
  int    idat, jdat;
  int    nmin, nmax;
  int    ix, iy;
  float  *pgx[2], *pgy[2];
  float  ftmp1,  ftmp2;
  float  pgxmin, pgxmax;
  float  pgymin, pgymax;
  float  dphs;
  float  a, b, c;
  double ar, ai;
  float  *phs, *tim;
  char   p_title[50], y_lab[50];
  FILE   *log_fp=NULL;
  float  F2[SRC_NUM];

/*
-----------------------------
*/

  cpgbbuf();
  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);
  cpgpap(1.0*pgpap_prm, 1.0);
  cpgsch(1.0);

/*
-----------------------------
*/

  ibase = 0;
  for (iant=BGN_ANT_I; iant<END_ANT_I; iant++) {
    for (jant=BGN_ANT_J; jant<END_ANT_J; jant++) {
      ibase++;
    }
  }

/*
-----------------------------
*/

  if ((pgx[0] = (float *)calloc(ibase, sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for pgx.\n");
    return (__NG__);
  }
  if ((pgx[1] = (float *)calloc(ibase, sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for pgx.\n");
    free (pgx[0]);
    return (__NG__);
  }
  if ((pgy[0] = (float *)calloc(ibase, sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for pgy.\n");
    free (pgx[0]);
    free (pgx[1]);
    return (__NG__);
  }
  if ((pgy[1] = (float *)calloc(ibase, sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for pgy.\n");
    free (pgx[0]);
    free (pgx[1]);
    free (pgy[0]);
    return (__NG__);
  }

  if ((tim = (float *)calloc(nobs,  sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for tim.\n");
    free (pgx[0]);
    free (pgx[1]);
    free (pgy[0]);
    free (pgy[1]);
    return (__NG__);
  }
  if ((phs = (float *)calloc(nobs,  sizeof(float))) == NULL) {
    printf("ERROR: allanv_disp: memory allocation for rms.\n");
    free (pgx[0]);
    free (pgx[1]);
    free (pgy[0]);
    free (pgy[1]);
    free (tim);
    return (__NG__);
  }

/*
-----------------------------
*/

  for (ns=0; ns<SRC_NUM; ns++) {
    F2[ns] = 1.0 / (float)(nu[ns] * nu[ns]);
  }

  if (ascii_out[0] == '!') {
    if ((log_fp=fopen(ascii_out+1, "w")) == NULL) {
      printf("ERROR: SSF_DISP\n");
      return (__NG__);
    }
  }

  for (IF=0; IF<2; IF++) {
    if (IF == 0) {
      NS = 0;
    } else {
      NS = 2;
    }

    nbase[IF] = 0;
    for (iant=0; iant<ANT_NUM; iant++) {
      for (jant=iant+1; jant<ANT_NUM; jant++) {
        if (iant >= BGN_ANT_I && iant < END_ANT_I &&
            jant >= BGN_ANT_J && jant < END_ANT_J) {

          if (PROC_MODE == SSF_PHASE) {
            i = 0;
            for (iobs=0; iobs<nobs; iobs++) {
              j = nbase[IF] * nobs + iobs;
              if (fringe_weight[NS][j] > 0.0) {
                tim[i] = (float)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5] + iobs);
                ar = 0.0;
                ai = 0.0;
                for (ifrq=0; ifrq<nfrq; ifrq++) {
                  k = j * nfrq + ifrq;
                  ar += frng[NS][k].rl;
                  ai += frng[NS][k].im;
                }
                phs[i] = atan2(ai, ar);
                i++;
              }
            }

            if (i != 0) {
              ndat = i;
              for (i=0; i<ndat-1; i++) {
                dphs = phs[i+1] - phs[i];
                if (fabs(dphs) > (float)dpi) {
                  phs[i+1] -= 2.0 * dpi * lrint(dphs / 2.0 / dpi);
                }
              }
#ifdef __RMS_PHASE__
              lstsqr(ndat, tim, phs, &a, &b, &c, 2);
#elif defined __STD__DVTN__
              lstsqr(ndat, tim, phs, &a, &b, &c, 1);
#endif
/**
              pgy[IF][nbase[IF]] = sqrt(c) / 2.0 / dpi / nu[NS] * speed_of_light;
**/
              pgy[IF][nbase[IF]] = sqrt(c);
              pgx[IF][nbase[IF]] = sqrt(
                   pow(diff(ant_prm[iant].XYZ[0], ant_prm[jant].XYZ[0]), 2.0) 
                 + pow(diff(ant_prm[iant].XYZ[1], ant_prm[jant].XYZ[1]), 2.0) 
                 + pow(diff(ant_prm[iant].XYZ[2], ant_prm[jant].XYZ[2]), 2.0));

              if (ascii_out[0] == '!') {
                fprintf(log_fp, "%s  %s: %f, %f\n",
                     ant_prm[iant].IDC, ant_prm[jant].IDC,
                     pgx[IF][nbase[IF]], pgy[IF][nbase[IF]]);
              }
              nbase[IF]++;
            }

/*
--------
*/

          } else if (PROC_MODE == SSF_DELAY) {

            i = 0;
            for (iobs=0; iobs<nobs; iobs++) {
              j = nbase[IF] * nobs + iobs;
              if (fringe_weight[NS][j] > 0.0) {
                ix = iant * nobs + iobs;
                iy = jant * nobs + iobs;
                tim[i] = (float)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5] + iobs);
                phs[i] = diff(int_obs[0][ix].grp_delay, int_obs[0][iy].grp_delay)
               + F2[0] * diff(int_obs[0][ix].ion_delay, int_obs[0][iy].ion_delay);
                phs[i] = diff(int_obs[0][ix].grp_delay, int_obs[0][iy].grp_delay);
                if (IF == 1) {
                  phs[i] -= diff(int_obs[1][ix].grp_delay, int_obs[1][iy].grp_delay)
                  - F2[1] * diff(int_obs[1][ix].ion_delay, int_obs[1][iy].ion_delay);
                  phs[i] -= diff(int_obs[1][ix].grp_delay, int_obs[1][iy].grp_delay);
                }
                phs[i] *= 2.0 * dpi * nu[0];
                i++;
              }
            }

            if (i != 0) {
              ndat = i;
#ifdef __RMS_PHASE__
              lstsqr(ndat, tim, phs, &a, &b, &c, 2);
#elif defined __STD__DVTN__
              lstsqr(ndat, tim, phs, &a, &b, &c, 1);
#endif
              pgy[IF][nbase[IF]] = sqrt(c);
              pgx[IF][nbase[IF]] = sqrt(
                   pow(diff(ant_prm[iant].XYZ[0], ant_prm[jant].XYZ[0]), 2.0) 
                 + pow(diff(ant_prm[iant].XYZ[1], ant_prm[jant].XYZ[1]), 2.0) 
                 + pow(diff(ant_prm[iant].XYZ[2], ant_prm[jant].XYZ[2]), 2.0));

              if (ascii_out[0] == '!') {
                fprintf(log_fp, "%s  %s: %f, %f\n",
                     ant_prm[iant].IDC, ant_prm[jant].IDC,
                     pgx[IF][nbase[IF]], pgy[IF][nbase[IF]]);
              }
              nbase[IF]++;
            }
          }

/*
--------
*/

        }
      }
    }
  }
  if (ascii_out[0] == '!') {
    fclose (log_fp);
  }

/*
-----------------------------
*/

  fmaxmin(nbase[0], pgx[0], &ftmp1, &ftmp2, &nmin, &nmax);
  pgxmin = ftmp1;
  pgxmax = ftmp2;
  fmaxmin(nbase[0], pgy[0], &ftmp1, &ftmp2, &nmin, &nmax);
  pgymin = ftmp1;
  pgymax = ftmp2;

  fmaxmin(nbase[1], pgx[1], &ftmp1, &ftmp2, &nmin, &nmax);
  if (pgxmin > ftmp1) {
    pgxmin = ftmp1;
  }
  if (pgxmax < ftmp2) {
    pgxmax = ftmp2;
  }
  fmaxmin(nbase[1], pgy[1], &ftmp1, &ftmp2, &nmin, &nmax);
  if (pgymin > ftmp1) {
    pgymin = ftmp1;
  }
  if (pgymax < ftmp2) {
    pgymax = ftmp2;
  }

  pgxmin = log10(pgxmin);
  pgxmax = log10(pgxmax);
  pgymin = log10(pgymin);
  pgymax = log10(pgymax);

  ftmp1 = pgxmax - pgxmin;
  ftmp2 = pgymax - pgymin;
  if (ftmp1 > ftmp2) {
    ftmp1 = 0.5 * (ftmp1 - ftmp2);
    pgymin -= ftmp1;
    pgymax += ftmp1;
  } else {
    ftmp1 = 0.5 * (ftmp2 - ftmp1);
    pgxmin -= ftmp1;
    pgxmax += ftmp1;
  }
  ftmp1 = 0.05 * (pgxmax - pgxmin);
  pgxmin -= ftmp1;
  pgxmax += ftmp1;
  pgymin -= ftmp1;
  pgymax += ftmp1;

  for (IF=0; IF<2; IF++) {
    for (ibase=0; ibase<nbase[IF]; ibase++) {
/********
      printf("Baseline length [m] : %7d  \n", (int)lrint(pgx[IF][ibase]));
********/
      pgx[IF][ibase] = log10(pgx[IF][ibase]);
      pgy[IF][ibase] = log10(pgy[IF][ibase]);
    }
  }

  cpgscf(2);
  cpgsvp(0.15, 0.90, 0.15, 0.90);
  cpgswin(pgxmin, pgxmax, pgymin, pgymax);
  cpgsci(1);
  cpgtbox("BCNLTS", 0.0, 0, "BCNLTS", 0.0, 0);
  if (PROC_MODE == SSF_PHASE) {
/**
    cpglab("\\frBaseline [m]", "\\frRMS Phase (EPL) [m]", "");
**/
    cpglab("\\frBaseline [m]", "\\frRMS Phase (rad)", "");
  } else if (PROC_MODE == SSF_DELAY) {
    cpglab("\\frBaseline [m]", "\\frRMS Phase (rad)", "");
  }

  cpgsch(3.0);
  cpgsci(2);
  cpgpt(nbase[0], pgx[0], pgy[0], 22);
/**
  cpgsci(4);
  cpgpt(nbase[1], pgx[1], pgy[1],  2);
**/

  cpgebuf();

/*
-----------------------------
*/

  for (IF=0; IF<2; IF++) {
    i = 0;
    for (ibase=0; ibase<nbase[IF]; ibase++) {
      if (pgx[IF][ibase] >= 3.0) {
        pgx[IF][i] = pgx[IF][ibase];
        pgy[IF][i] = pgy[IF][ibase];
        i++;
      }
    }
    lstsqr(i, pgx[IF], pgy[IF], &a, &b, &c, 0);
    printf("RMS values [micron] (1km, 10km)     %f    %f\n",
            pow(10.0, a*3.0+b) * 1.0e6, pow(10.0, a*4.0+b) * 1.0e6);
  }



/*
-----------------------------
*/

  free (pgx[0]);
  free (pgx[1]);
  free (pgy[0]);
  free (pgy[1]);
  free (tim);
  free (phs);

  return (__GO__);
}
