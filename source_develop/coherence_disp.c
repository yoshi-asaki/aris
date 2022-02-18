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
#define __AMPLITUDE__
****/
#define __AMPLITUDE__

/****
#define __DEBUG_FRINGE_DISP__
****/
#define __DEBUG_FRINGE_DISP__

/****
#define __REF_SOURCE__
****/


int  coherence_disp(int FRNG_NUM, int ANT_NUM, int NX, int NY,
                    struct antenna_parameter *ant_prm,
                    int nobs, int nfrq,
                    int    *TimUT,  double UT1_UTC,
                    struct fringe *frng[], float *fringe_weight[],
                    float  *rms_phase)
{
  int    i, j, k, IBASE, IF, iobs, ifrq;
  int    N, M, NT, I, J, NOBS, NDIV=10;
  int    nt;
  char   p_title[100];
  float  *pgx, *pgy, *pgz;
  float  *ar, *ai;
  float  pgxmin, pgxmax, pgymax;
  float  a, b, c;
  float  amp;
  float  xr, xi, phs;
  int    mrk[2];

/*
--------
*/

  mrk[0] = 17;
  mrk[1] = 16;

  IBASE = baseline_number(ANT_NUM, NX,        NY);

  if ((pgx = (float *)calloc(NDIV, sizeof(float))) == NULL) {
    printf("FRINGE_DISP: ERROR; memory allocation.\n");
    return (-1);
  }
  if ((pgy = (float *)calloc(NDIV, sizeof(float))) == NULL) {
    printf("FRINGE_DISP: ERROR; memory allocation.\n");
    free (pgx);
    return (-1);
  }
  if ((pgz = (float *)calloc(NDIV, sizeof(float))) == NULL) {
    printf("FRINGE_DISP: ERROR; memory allocation.\n");
    free (pgx);
    free (pgy);
    return (-1);
  }
  if ((ar  = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("FRINGE_DISP: ERROR; memory allocation.\n");
    free (pgx);
    free (pgy);
    free (pgz);
    return (-1);
  }
  if ((ai  = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("FRINGE_DISP: ERROR; memory allocation.\n");
    free (pgx);
    free (pgy);
    free (pgz);
    free (ar);
    return (-1);
  }

  pgxmin = log10((float)(nobs/NDIV));
  pgxmax = log10((float)nobs);

  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);
/****
  cpgpap(0.6*pgpap_prm, 1.0);
****/
  cpgpap(0.75 * pgpap_prm, 1.0);
  cpgsch(1.0);

  cpgsci(1);
  cpgsvp(0.20, 0.80, 0.60, 0.90);
  cpgswin(pgxmin, pgxmax, 0.0, 1.5);
  cpgbox("BCNLTS", 0, 0, "BCNTS", 0.0, 0);

#ifdef __DEBUG_FRINGE_DISP__
  printf("%2d/%2d", IBASE+1, (ANT_NUM * (ANT_NUM - 1)) / 2);
#endif /* __DEBUG_FRINGE_DISP__ */

  sprintf(p_title, "%s - %s", ant_prm[NX].IDC, ant_prm[NY].IDC);
  cpglab("Integration Time [sec]", "Coherence", p_title);

#ifndef __REF_SOURCE__
  FRNG_NUM = 1;
#endif /* __REF_SOURCE__ */

  for (IF=0; IF<FRNG_NUM; IF++) {
    for (iobs=0; iobs<nobs; iobs++) {
      j = IBASE * nobs + iobs;
      ar[iobs] = 0.0;
      ai[iobs] = 0.0;
      for (ifrq=0; ifrq<nfrq; ifrq++) {
        k = j * nfrq + ifrq;
        ar[iobs] += frng[IF][k].rl;
        ai[iobs] += frng[IF][k].im;
      }
      phs = atan2(ai[iobs], ar[iobs]);
      ar[iobs] = cos(phs);
      ai[iobs] = sin(phs);
    }

#ifdef __COLOR__
    cpgsci(IF+1);
#endif /* __COLOR__ */

    for (N=0; N<NDIV; N++) {
      NOBS = nobs / (N + 1);
      I = 0;
      for (M=0; M<=N; M++) {
        nt = 0;
        xr = 0.0;
        xi = 0.0;
        for (iobs=NOBS*M; iobs<NOBS*(M+1); iobs++) {
          j = IBASE * nobs + iobs;
          if (fringe_weight[IF][j] > 0.0) {
            xr += ar[iobs];
            xi += ai[iobs];
            nt++;
          } 
        }
        if (nt > 0) {
          pgy[I] = sqrt(xr*xr + xi*xi) / (float)nt;
          I++;
        }
      }
      if (I == 0) {
        pgz[N] = pgy[0];
      } else {
        lstsqr(I, pgx, pgy, &a, &pgz[N], &c, 1);
      }
    }

    for (N=0; N<NDIV; N++) {
      pgx[N] = log10((float)(nobs / (N + 1)));
    }

    cpgsch(3.0);
    cpgpt(NDIV, pgx, pgz, mrk[IF]);
    cpgline(NDIV, pgx, pgz);
    cpgsch(1.0);
    cpgsci(1);
  }

/*
---------------------------------------------------
*/

  free (pgx);
  free (pgy);
  free (pgz);
  free (ar );
  free (ai );

  return (1);
}
