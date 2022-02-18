#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


#define NCOM      24
#define MENU_NUM  20

#define REAL(z, i) ((z)[2*(i)])
#define IMAG(z, i) ((z)[2*(i)+1])

#define ANTNUM    12
#define F1      2048

#define NFDAT     32
#define NRDAT     16
#define NDAT      64


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

/****
#define __REF_SOURCE__
****/
#define __REF_SOURCE__



int  fringe_fitting(int FRNG_NUM, int ANT_NUM, int NX, int NY,
                    struct antenna_parameter *ant_prm,
                    int  nobs, int  nfrq,
                    int  *TimUT,   double UT1_UTC,
                    struct fringe *frng[],
                    float  *fringe_weight[],
                    float  *rms_phase)
{
  int    i, j, k, IBASE, IF, iobs, ifrq;
  char   p_title[100];
  float  *pgx, *pgy;
  float  a, b, c;
  int    mrk[2];

/*
--------
*/

  int    I, J, itmp, jtmp, icount;
  int    ICNT, IF_NUM;
  int    IFREF, JFREF;
  int    refant;
  int    ndat, idat, iant, i_if, ncnt;
  FILE   *fp;
  char   string[1000];
  long   firstelem=1;
  double ar, ai, v_amp;
  double t_tmp0, timtmp;
  double DPI;
  double D_TIM;
  _Bool  PROC_SWT;
  float  pgxmin, pgxmax, pgymin, pgymax;

  float  *dist;
  int    dmax;

  double tgt_frq, *ref_frq;
  double tgt_dfq, ref_dfq;
  int    tgt_fcn, ref_fcn;
  double *diff_frq;

  float  pmax, p_tmp;
  float  rl, im, uu, vv, ww, wt;
  double *wx;
  double *vis_rl, *vis_im;
  float  *stim, *frn_phs;
  double *tim;

  double pix_mas = 1.0;
  double uv_max, pix_uvl, uv_factor, uvl;
  double DELA[ANTNUM][3], AMPL[ANTNUM][10];
  double delay, phs, dphs;

  int    all_hdunum, hdunum, ihd;
  int    status=0;
  int    hdutype;
  char   comment[80], value[80];
  double obs_freq,    if_freq;

  float  fr[NRDAT][NFDAT], fi[NRDAT][NFDAT];
  float  aar, aai;
  float  *dela, *rate;
  float  tmin, tmax, noise, err_x, err_y, delta_x, delta_y, resx, resy;
  float  peak_val, peak_tmp, peak_x, peak_y;
  double *srchwr, *srchwi;
  double solint, solsub;
  int    solstep;

/*
--------
*/

  solint = 1024.0;
  solsub = 4.0;

/*
--------
*/

  pgxmin = (float)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5]);
  pgxmax = (float)(3600*TimUT[3] + 60*TimUT[4] + TimUT[5] + nobs);

  mrk[0] = 4;
  mrk[1] = 2;

  mrk[0] = 1;
  mrk[1] = 1;

  IBASE = baseline_number(ANT_NUM, NX,        NY);

  if ((pgx = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("FRINGE_FITTING: ERROR: memory allocation.\n");
    return (-1);
  }
  if ((pgy = (float *)calloc(nobs, sizeof(float))) == NULL) {
    printf("FRINGE_FITTING: ERROR: memory allocation.\n");
    free (pgx);
    return (-1);
  }
  if ((wx = (double *)calloc(2*nfrq, sizeof(double))) == NULL) {
    printf("FRINGE_FITTING: ERROR: memory allocation.\n");
    free (pgx);
    free (pgy);
    return (-1);
  }
  if ((srchwr = (double *)calloc(4*nobs*nfrq, sizeof(double))) == NULL) {
    printf("FRINGE_FITTING: ERROR: memory allocation.\n");
    free (pgx);
    free (pgy);
    free (wx);
    return (-1);
  }
  if ((srchwi = (double *)calloc(4*nobs*nfrq, sizeof(double))) == NULL) {
    printf("FRINGE_FITTING: ERROR: memory allocation.\n");
    free (pgx);
    free (pgy);
    free (wx);
    free (srchwr);
    return (-1);
  }

  for (iobs=0; iobs<nobs; iobs++) {
    j = IBASE * nobs + iobs;
    for (ifrq=0; ifrq<4*nfrq; ifrq++) {
      wx[ifrq] = 0.0;
    }
    for (ifrq=0; ifrq<nfrq; ifrq++) {
      k = j * nfrq + ifrq;
      wx[2*ifrq  ] = frng[0][k].rl;
      wx[2*ifrq+1] = frng[0][k].im;
    }

    gsl_fft_complex_radix2_backward(wx, 1, 2*nfrq);
    for (ifrq=0; ifrq<nfrq; ifrq++) {
      k = j * 4*nfrq + nfrq + ifrq;
      srchwr[k] = wx[2*ifrq  ] / (float)nfrq;
      srchwi[k] = wx[2*ifrq+1] / (float)nfrq;
    }
    for (ifrq=nfrq; ifrq<2*nfrq; ifrq++) {
      k = j * 4*nfrq + ifrq;
      srchwr[k] = wx[2*ifrq  ] / (float)nfrq;
      srchwi[k] = wx[2*ifrq+1] / (float)nfrq;
    }
  }

  free (wx);

/*
--------
*/

  solstep = (int)ceil(solint / solsub);
  if ((dela = (float *)calloc(ceil(nobs/solstep), sizeof(float))) == NULL) {
    printf("FRINGE_FITTING: ERROR: memory allocation for dela.\n");
    free (pgx);
    free (pgy);
    free (srchwr);
    free (srchwi);
    return (-1);
  }
  if ((rate = (float *)calloc(ceil(nobs/solstep), sizeof(float))) == NULL) {
    printf("FRINGE_FITTING: ERROR: memory allocation for dela.\n");
    free (pgx);
    free (pgy);
    free (srchwr);
    free (srchwi);
    free (dela);
    return (-1);
  }

/*
--------
*/

  if ((wx = (double *)calloc(2*F1, sizeof(double))) == NULL) {
    printf("FRINGE_FITTING: ERROR: memory allocation.\n");
    free (pgx);
    free (pgy);
    free (srchwr);
    free (srchwi);
    free (dela);
    free (rate);
    return (-1);
  }

  icount = 0;
  ICNT   = 0;
  PROC_SWT = true;

  itmp = 0;
  jtmp = 0;
  ifrq = 0;
  for (iobs=0; iobs<nobs; iobs+=solstep) {
    for (i=0; i<2*nfrq; i++) {
      for (j=0; j<F1; j++) {
        wx[2*j  ] = srchwr[(icount+j)*4*nfrq+ifrq];
        wx[2*j+1] = srchwi[(icount+j)*4*nfrq+ifrq];
      }
      gsl_fft_complex_radix2_forward(wx, 1, F1);
      peak_tmp = 0.0;
      peak_val = 0.0;
      peak_x = 0;
      peak_y = 0;
      for (j=0; j<F1; j++) {
        peak_tmp = wx[2*j  ]*wx[2*j  ] + wx[2*j+1]*wx[2*j+1];
        if (peak_tmp > peak_val) {
          peak_val = peak_tmp;
          peak_x = i;
          peak_y = j;
        }
      }
    }
    if (peak_y > F1/2) {
      peak_y -= F1;
    }

    pgx[iobs]  = pgxmin + (float)(icount / 2);
    dela[iobs] = (float)itmp / (float)NFDAT;
    rate[iobs] = (float)jtmp / (float)NRDAT;
    icount += F1;

    if (iobs + icount >= nobs) {
      break;
    }
  }

/*
-------------------------------------------
*/

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

  free (pgx);
  free (pgy);

  return (1);
}
