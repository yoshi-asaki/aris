#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __DEBUG__
****/

#define FLDMAXA  1024
#define LFDMAXA   128


int qlook_imager(int ANT_NUM,
                 int BGN_ANT_I, int END_ANT_I, int BGN_ANT_J, int END_ANT_J,
                 int   nswt,   int nobs,      int nfrq,
                 struct baseline_uvw *bluvw[],
                 struct fringe *frng[],
                 float *fringe_weight[], double *wave_length)
{
  int    NSWT;
  int    N1, N2, N3;
  int    i, j, I, J, iobs, ibase, iant, jant, ifrq;
  int    TSAMPLE;
  char   *title[2];
  float  *dist[2];
  float  pmin[2], pmax[2], noise[2];
  float  err_x[2],   err_y[2];
  float  delta_x[2], delta_y[2];
  float  s_x_cnt[2], s_y_cnt[2], s_x_w[2], s_y_w[2];
  double ar, ai, wt;
  double pix_mas[2];
  struct baseline_uvw  *im_uvw;
  struct fringe        *imfrng;

  int    IDAT, NDAT, iave;
  _Bool  AVE_SWT;
  double ar_tmp, ai_tmp, wt_tmp;
  double uu, vv, ww;

#ifdef __DEBUG__
  char   string[100];
#endif

/*
--------------------------
*/

  for (N1=0; N1<2; N1++) {
    if ((dist[N1]=(float *)calloc(LFDMAXA*LFDMAXA, sizeof(float))) == NULL) {
      printf("ERROR: qlook_iamger: calloc error.\n");
      exit (1);
    }
  }

/*
--------------------------
*/

  if (nswt == 0) {
    NSWT = 10;
  } else {
    NSWT = nswt;
  }

  ibase = 0;
  for (iant=0; iant<ANT_NUM; iant++) {
    for (jant=iant+1; jant<ANT_NUM; jant++) {
      if (iant >= BGN_ANT_I && iant < END_ANT_I &&
          jant >= BGN_ANT_J && jant < END_ANT_J) {
        ibase++;
      }
    }
  }
  NDAT = ibase * (nobs / NSWT + 1);

  if ((im_uvw = (struct baseline_uvw *)
            calloc(NDAT, sizeof(struct baseline_uvw))) == NULL) {
    printf("ERROR: qlook_imager: calloc error\n");
    free (dist[0]);
    free (dist[1]);
    exit (1);
  }
  if ((imfrng = (struct fringe *)calloc(NDAT, sizeof(struct fringe)))
       == NULL) {
    printf("ERROR: qlook_imager: calloc error\n");
    free (dist[0]);
    free (dist[1]);
    free (im_uvw);
    exit (1);
  }

/*
----------------
*/

  for (N1=0; N1<2; N1++) {
    if        (N1 == 0) {
      N2 = 0;
      N3 = 0;
    } else if (N1 == 1) {
      N2 = 2;
      N3 = 0;
    }

    IDAT = 0;
    for (iant=0; iant<ANT_NUM; iant++) {
      for (jant=iant+1; jant<ANT_NUM; jant++) {
        if (iant >= BGN_ANT_I && iant < END_ANT_I &&
            jant >= BGN_ANT_J && jant < END_ANT_J) {
          ibase = baseline_number(ANT_NUM, iant, jant);

          iobs = 0;
          while (iobs < nobs) {
            AVE_SWT = true;
            ar = (double)0.0;
            ai = (double)0.0;
            wt = (double)0.0;
            uu = (double)0.0;
            vv = (double)0.0;
            ww = (double)0.0;
            iave = 0;

            while (1) {
              ar_tmp = 0.0;
              ai_tmp = 0.0;
              wt_tmp = 0.0;
              J = ibase * nobs + iobs;
              if (fringe_weight[N2][J] > 0.0) {
                for (ifrq=0; ifrq<nfrq; ifrq++) {
                  I = (ibase * nobs + iobs) * nfrq + ifrq;
                  ar_tmp += frng[N2][I].rl;
                  ai_tmp += frng[N2][I].im;
                  wt_tmp += frng[N2][I].wt;
                }
                uu += bluvw[N3][J].u;
                vv += bluvw[N3][J].v;
                ww += bluvw[N3][J].w;
              } else {
                break;
              }

              ar += ar_tmp / (double)nfrq;
              ai += ai_tmp / (double)nfrq;
              wt += wt_tmp / (double)nfrq;
              iave++;
              iobs++;
              if (iobs >= nobs - 1) {
                AVE_SWT = false;
                break;
              }
              if (iave >= NSWT) {
                break;
              }
            }

            if (AVE_SWT == false) {
              break;
            } else {
              if (iave != 0) {
                if (uu <= 0.0) {
                  imfrng[IDAT].rl =  ar / (double)iave;
                  imfrng[IDAT].im =  ai / (double)iave;
                  imfrng[IDAT].wt =  wt / (double)iave;
                  im_uvw[IDAT].u  =  uu / (double)iave / wave_length[N2];
                  im_uvw[IDAT].v  =  vv / (double)iave / wave_length[N2];
                  im_uvw[IDAT].w  =  ww / (double)iave / wave_length[N2];
                } else {
                  imfrng[IDAT].rl =  ar / (double)iave;
                  imfrng[IDAT].im = -ai / (double)iave;
                  imfrng[IDAT].wt =  wt / (double)iave;
                  im_uvw[IDAT].u  = -uu / (double)iave / wave_length[N2];
                  im_uvw[IDAT].v  = -vv / (double)iave / wave_length[N2];
                  im_uvw[IDAT].w  = -ww / (double)iave / wave_length[N2];
                }

/*
    The following part should be verified. 
*/
                imfrng[IDAT].im *= -1.0;
                IDAT++;
              }

              while (1) {
                J = ibase * nobs + iobs;
                if (fringe_weight[N2][J] <= 0.0) {
                  iobs++;
                  if (iobs >= nobs - 1) {
                    break;
                  }
                } else {
                  break;
                }
              }
            }
          }
        }
      }
    }
    imager(IDAT, im_uvw, &pix_mas[N1], FLDMAXA, imfrng, dist[N1], LFDMAXA);
  }

  free (im_uvw);
  free (imfrng);

/*
----------------
*/

  for (N1=0; N1<2; N1++) {
    if ((title[N1]=(char *)calloc(20, sizeof(char))) == NULL) {
      printf("ERROR: qlook_iamger: calloc error.\n");
      exit (1);
    }
  }
  sprintf(title[0], "Target (before P-R)");
  sprintf(title[1], "Target (after P-R)");
  brightness_disp(2, LFDMAXA, LFDMAXA/2, LFDMAXA/2,
                  pix_mas[0], pix_mas[0],
                  100.0*pix_mas[0], 100.0*pix_mas[0], 0.0, 0.0,
                  false, s_x_cnt, s_y_cnt, s_x_w, s_y_w, dist,
                  "[mas]",  title,
                  true, true, true,  true, true, 64, "clr",
                  pmin, pmax, noise, err_x, err_y, delta_x, delta_y);

#ifdef __DEBUG__
  for (N1=0; N1<2; N1++) {
    sprintf(string, "Maximum Peak: %7.2E\n", pmax[0]);
/**
    if (TV_SWT == true) {
      cpgsvp(0.0, 1.0, 0.0, 1.0);
      cpgswin(0.0, 1.0, 0.0, 1.0);
      cpgsch(0.65);
      comment_disp(cmnt, comment, string, ON);
    } else {
**/
      printf("%s", string);
/**
    }
**/
    sprintf(string,  "(X, Y)      : (%f, %f)\n", delta_x[0], delta_y[0]);
/**
    if (TV_SWT == true) {
      comment_disp(cmnt, comment, string, ON);
    } else {
**/
      printf("%s", string);
/**
    }
**/
  }
#endif

  for (N1=0; N1<2; N1++) {
    free (dist[N1]);
    free (title[N1]);
  }

  return 1;
}
