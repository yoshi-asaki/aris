#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


/****
#define __COLOR__
****/
#define __COLOR__
#ifndef __COLOR__
  #define __BLACK_WHITE__
#endif

/****
#define __DEBUG_FRINGE_DISP__
****/

/****
#define __AMPLITUDE__
#define __REF_SOURCE__
#define __LOG_FILE__
****/
#define __REF_SOURCE__
#define __AMPLITUDE__



int  fringe_disp(int FRNG_NUM, int ANT_NUM, int NX, int NY,
                 struct antenna_parameter *ant_prm,
                 int nobs, int nfrq,
                 int    *TimUT,  double UT1_UTC,
                 struct fringe *frng[], float *fringe_weight[],
                 float  *rms_phase, float *coherence)
{
  int    i, j, k, IBASE, IF, iobs, ifrq;
  int    ndat;
  char   p_title[100];
  float  *pgx, *pgy;
  float  ar, ai;
  float  tr, ti;
  float  pgxmin, pgxmax, pgymax;
  float  a, b, c;
  int    mrk[2];
#ifdef __LOG_FILE__
  FILE   *log_fp;
#endif

/*
--------
*/

  mrk[0] = 1;
  mrk[1] = 1;

  IBASE = baseline_number(ANT_NUM, NX,        NY);

  if ((pgx = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("FRINGE_DISP: ERROR; memory allocation.\n");
    return (-1);
  }
  if ((pgy = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("FRINGE_DISP: ERROR; memory allocation.\n");
    free (pgx);
    return (-1);
  }

  pgxmin = (float)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5]);
  pgxmax = (float)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5] + nobs);

  cpgbbuf();
  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);
/****
  cpgpap(0.6*pgpap_prm, 1.0);
****/
  cpgpap(0.75 * pgpap_prm, 1.0);
  cpgsch(1.0);

  cpgsci(1);
  cpgsvp(0.20, 0.80, 0.60, 0.90);
  cpgswin(pgxmin, pgxmax, -179.9, 180.0);
  cpgtbox("BCSTNZHI", 0.0, 0, "BCNTS", 90.0, 3);

#ifdef __DEBUG_FRINGE_DISP__
  printf("%2d/%2d", IBASE+1, (ANT_NUM * (ANT_NUM - 1)) / 2);
#endif /* __DEBUG_FRINGE_DISP__ */

  sprintf(p_title, "%s - %s", ant_prm[NX].IDC, ant_prm[NY].IDC);
  cpglab("Time (UT)", "Fringe phase [deg]", p_title);

#ifndef __REF_SOURCE__
  FRNG_NUM = 1;
#endif /* __REF_SOURCE__ */

  for (IF=0; IF<FRNG_NUM; IF++) {
    i = 0;
    for (iobs=0; iobs<nobs; iobs++) {
      j = IBASE * nobs + iobs;
      if (fringe_weight[IF][j] > 0.0) {
        pgx[i] = (float)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5] + iobs);
        ar = 0.0;
        ai = 0.0;
        for (ifrq=0; ifrq<nfrq; ifrq++) {
          k = j * nfrq + ifrq;
          ar += frng[IF][k].rl;
          ai += frng[IF][k].im;
        }
        pgy[i] = 180.0 / dpi * atan2(ai, ar);
        i++;
      }
    }

#ifdef __COLOR__
    cpgsci(IF+1);
#endif /* __COLOR__ */
    cpgsch(3.0);
    cpgpt(i, pgx, pgy, mrk[IF]);
    cpgsch(1.0);
    cpgsci(1);

    lstsqr(i, pgx, pgy, &a, &b, &c, 1);
    rms_phase[IF] = sqrt(c);
  }

/*
---------------------------------------------------
*/

#ifdef __AMPLITUDE__

  pgymax = 0.0;
  ndat = 0;
  tr   = 0.0;
  ti   = 0.0;
  for (IF=0; IF<FRNG_NUM; IF++) {
    i = 0;
    for (iobs=0; iobs<nobs; iobs++) {
      j = IBASE * nobs + iobs;
      if (fringe_weight[IF][j] > 0.0) {
        ar = 0.0;
        ai = 0.0;
        for (ifrq=0; ifrq<nfrq; ifrq++) {
          k = j * nfrq + ifrq;
          ar += frng[IF][k].rl;
          ai += frng[IF][k].im;
        }
        pgy[i] = sqrt(ar*ar + ai*ai) / (float)nfrq;
        tr   += ar;
        ti   += ai;
        ndat += nfrq;
        if (pgy[i] > pgymax) {
          pgymax = pgy[i];
        }
        i++;
      }
    }
  }
  if (pgymax == 0.0) {
    pgymax = 1.0;
  }
  *coherence = sqrt(tr * tr + ti * ti) / (float)ndat;

  cpgsci(1);
  cpgsvp(0.20, 0.80, 0.15, 0.45);
  cpgswin(pgxmin, pgxmax, 0.0, 1.2*pgymax);
  cpgtbox("BCSTNZHI", 0.0, 0, "BCNTS", 0.0, 0);
  cpglab("Time (UT)", "Fringe amplitude [Jy]", "");

  for (IF=0; IF<FRNG_NUM; IF++) {
    i = 0;
    for (iobs=0; iobs<nobs; iobs++) {
      j = IBASE * nobs + iobs;
      if (fringe_weight[IF][j] > 0.0) {
        pgx[i] = (float)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5] + iobs);
        ar = 0.0;
        ai = 0.0;
        for (ifrq=0; ifrq<nfrq; ifrq++) {
          k = j * nfrq + ifrq;
          ar += frng[IF][k].rl;
          ai += frng[IF][k].im;
        }
        pgy[i] = sqrt(ar*ar + ai*ai) / (float)nfrq;
        i++;
      }
    }
#ifdef __COLOR__
    cpgsci(IF+1);
#endif /* __COLOR__ */
    cpgpt(i, pgx, pgy, mrk[IF]);
    cpgsci(1);
  }

#endif /* __AMPLITUDE__ */

/*
---------------------------------------------------
*/

  cpgebuf();

/*
---------------------------------------------------
*/

#ifdef __LOG_FILE__
  if ((log_fp=fopen("aris_log/fringe.dat", "w")) == NULL) {
    printf("$$ WARNING: fringe_disp: log file cannot be opened.\n");
  } else {
    for (IF=0; IF<FRNG_NUM; IF++) {
      i = 0;
      for (iobs=0; iobs<nobs; iobs++) {
        j = IBASE * nobs + iobs;
        if (fringe_weight[IF][j] > 0.0) {
          ar = 0.0;
          ai = 0.0;
          for (ifrq=0; ifrq<nfrq; ifrq++) {
            k = j * nfrq + ifrq;
            ar += frng[IF][k].rl;
            ai += frng[IF][k].im;
          }
          fprintf(log_fp, "%6d  %10.1lf   %lf  %lf\n",
              i, (double)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5] + iobs),
              ar, ai);
          i++;
        }
      }
    }
    fclose(log_fp);
  }
#endif

/*
---------------------------------------------------
*/

  free (pgx);
  free (pgy);

  return (1);
}
