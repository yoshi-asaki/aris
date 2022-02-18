#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <aris.h>

/**
#define __DRAW__
**/
#define __DRAW__
#ifndef __DRAW__
  #define __CONT__
#endif


int   modified_distribution(
             int       , float    *, int       , int       ,
             float     , float     ,
             int       , float    *, int      *, int      *,
             float    *, float    *, int      *);
int   modified_distribution_XY(
             int       , int       , float    *,
             int       , int       , int       , float    *);


int  brightness_disp(int  NP,             int N,
                     int   center_nx,     int   center_ny,
                     float resx,          float resy,
                     float widthx,        float widthy,
                     float px,            float py,
                     _Bool srch_box_set_swt,
                     float *srch_center_x, float  *srch_center_y,
                     float *srch_width_x,  float  *srch_width_y,
                     float *dist[],        char   *label, char  *title[],
                     _Bool TV_SWT,         _Bool  BAR_SWT,
                     _Bool COMPARA_SWT,    _Bool  NEGATIVE_COMP_SWT,
                     _Bool X_REVERSE,      int    MAP_DISP_SIZE,
                     char  *clr_mode,
                     float *pmin,         float *pmax,
                     float *noise,
                     float *err_x,        float *err_y,
                     float *delta_x,      float *delta_y)
{
  int    i, n, ntmp;
  int    nx, ny, NW;
  _Bool  SMOOTH_SWT;
  int    center_mx, center_my;
  int    x0, x1, y0, y1;
  float  pmin_tmp, pmax_tmp;
  float  *m_dist=NULL;
  float  posi_pmax, snr;
  float  mrx, mry;
  float  bias[4],    width[4];
  float  barxmin[4], barxmax[4], barymin[4], barymax[4];
  float  mapxmin[4], mapxmax[4], mapymin[4], mapymax[4];
  char   *labx[4], *laby[4], *null_lab="\0";

/*
---------------------------------
*/

  for (i=0; i<NP; i++) {
    peak_normalize(N, pmin+i, pmax+i, noise+i, err_x+i, err_y+i,
                      delta_x+i, delta_y+i, resx, resy,
                      0, N, 0, N, dist[i]);
  }

/*
--------
*/

  if (NP >= 2 && COMPARA_SWT == true) {
    pmin_tmp = pmin[0];
    pmax_tmp = pmax[0];

    for (i=1; i<NP; i++) {
      if (pmin[i] < pmin_tmp) {
        pmin_tmp = pmin[i];
      }
      if (pmax[i] > pmax_tmp) {
        pmax_tmp = pmax[i];
      }
    }

    for (i=0; i<NP; i++) {
      pmin[i] = pmin_tmp;
      pmax[i] = pmax_tmp;
    }
  }

  for (i=0; i<NP; i++) {
    if (pmin[i] != pmax[i]) {
      bias[i]  = pmin[i];
      width[i] = pmax[i] - pmin[i];
    } else {
      bias[i]  = 0.0;
      width[i] = pmax[i];
    }
  }

  if (TV_SWT == false) {
    return (__GO__);
  }

/*
--------
*/

  SMOOTH_SWT = false;
  if (MAP_DISP_SIZE > 0) {
    nx = (int)(widthx / resx);
    ny = (int)(widthy / resy);
    if (nx >= ny) {
      NW = nx;
    } else {
      NW = ny;
    }
    if (N > MAP_DISP_SIZE && NW > MAP_DISP_SIZE) {
      n = NW / MAP_DISP_SIZE;
      if (NW % MAP_DISP_SIZE != 0) {
        n++;
      }
      if ((m_dist =
          (float *)calloc(MAP_DISP_SIZE*MAP_DISP_SIZE, sizeof(float))) == NULL) {
        printf("ERROR: brightness_disp: memory allocation for m_dist.\n");
        return (__NG__);
      }
      SMOOTH_SWT = true;
    }
  }

/*
--------
*/

  if (COMPARA_SWT == false) {
    if (NP <= 2) {
      for (i=0; i<NP; i++) {
        labx[i] = label;
        laby[i] = label;
      }
    } else {
      for (i=0; i<NP; i++) {
        labx[i] = null_lab;
        laby[i] = label;
      }
    }
  } else if (COMPARA_SWT == true) {
    for (i=0; i<NP; i++) {
      labx[i] = null_lab;
      laby[i] = label;
    }
  }

  if (BAR_SWT == true) {
    if (COMPARA_SWT == false) {
      if (NP == 1) {
        barxmin[0] = 0.30;
        barxmax[0] = 0.70;
        barymin[0] = 0.20;
        barymax[0] = 0.25;
/*aaaaaaa*/
        barymin[0] = 0.30;
        barymax[0] = 0.35;
/*aaaaaaa*/
      } else if (NP == 2) {
        barxmin[0] = 0.08;
        barxmax[0] = 0.48;
        barymin[0] = 0.30;
        barymax[0] = 0.35;
        barxmin[1] = barxmin[0] + 0.50;
        barxmax[1] = barxmax[0] + 0.50;
        barymin[1] = barymin[0];
        barymax[1] = barymax[0];
      } else if (NP == 3 || NP == 4) {
        barxmin[0] = 0.08;
        barxmax[0] = 0.48;
        barymin[0] = 0.57;
        barymax[0] = 0.59;
        barxmin[1] = barxmin[0] + 0.50;
        barxmax[1] = barxmax[0] + 0.50;
        barymin[1] = barymin[0];
        barymax[1] = barymax[0];
        barxmin[2] = barxmin[0];
        barxmax[2] = barxmax[0];
        barymin[2] = barymin[0] - 0.49;
        barymax[2] = barymax[0] - 0.49;
        barxmin[3] = barxmin[1];
        barxmax[3] = barxmax[1];
        barymin[3] = barymin[2];
        barymax[3] = barymax[2];
      }
    } else if (COMPARA_SWT == true) {
      if (NP <= 2) {
        barxmin[0] = 0.20;
        barxmax[0] = 0.80;
        barymin[0] = 0.35;
        barymax[0] = 0.40;
/****
        barymin[0] = 0.15;
        barymax[0] = 0.20;
****/
      } else {
        barxmin[0] = 0.30;
        barxmax[0] = 0.70;
        barymin[0] = 0.05;
        barymax[0] = 0.07;
      }
    }

    if (NP == 1) {
      mapxmin[0] = 0.25;
      mapxmax[0] = 0.75;
      mapymin[0] = 0.40;
      mapymax[0] = 0.90;
    } else if (NP == 2) {
      mapxmin[0] = 0.08;
      mapxmax[0] = 0.48;
      mapymin[0] = 0.50;
      mapymax[0] = 0.90;
      mapxmin[1] = mapxmin[0] + 0.50;
      mapxmax[1] = mapxmax[0] + 0.50;
      mapymin[1] = mapymin[0];
      mapymax[1] = mapymax[0];
    } else if (NP == 3 || NP == 4) {
      mapxmin[0] = 0.14;
      mapxmax[0] = 0.44;
      mapymin[0] = 0.63;
      mapymax[0] = 0.93;
      mapxmin[1] = mapxmin[0] + 0.50;
      mapxmax[1] = mapxmax[0] + 0.50;
      mapymin[1] = mapymin[0];
      mapymax[1] = mapymax[0];
      mapxmin[2] = mapxmin[0];
      mapxmax[2] = mapxmax[0];
      mapymin[2] = mapymin[0] - 0.49;
      mapymax[2] = mapymax[0] - 0.49;
      mapxmin[3] = mapxmin[1];
      mapxmax[3] = mapxmax[1];
      mapymin[3] = mapymin[2];
      mapymax[3] = mapymax[2];
    }
  } else if (BAR_SWT == false) {
    if (NP == 1) {
      mapxmin[0] = 0.15;
      mapxmax[0] = 0.90;
      mapymin[0] = 0.15;
      mapymax[0] = 0.90;
    } else if (NP == 2) {
      mapxmin[0] = 0.08;
      mapxmax[0] = 0.48;
      mapymin[0] = 0.50;
      mapymax[0] = 0.90;
      mapxmin[1] = mapxmin[0] + 0.50;
      mapxmax[1] = mapxmax[0] + 0.50;
      mapymin[1] = mapymin[0];
      mapymax[1] = mapymax[0];
    } else if (NP == 3 || NP == 4) {
      mapxmin[0] = 0.14;
      mapxmax[0] = 0.44;
      mapymin[0] = 0.63;
      mapymax[0] = 0.93;
      mapxmin[1] = mapxmin[0] + 0.50;
      mapxmax[1] = mapxmax[0] + 0.50;
      mapymin[1] = mapymin[0];
      mapymax[1] = mapymax[0];
      mapxmin[2] = mapxmin[0];
      mapxmax[2] = mapxmax[0];
      mapymin[2] = mapymin[0] - 0.49;
      mapymax[2] = mapymax[0] - 0.49;
      mapxmin[3] = mapxmin[1];
      mapxmax[3] = mapxmax[1];
      mapymin[3] = mapymin[2];
      mapymax[3] = mapymax[2];
    }
  }

/*
--------
*/

#ifdef __DRAW__
  if (BAR_SWT == true) {
    if (COMPARA_SWT == false) {
      for (i=0; i<NP; i++) {
/***
        pg_color_bar(bias[i], width[i],
                     barxmin[i], barxmax[i], barymin[i], barymax[i],
                     NEGATIVE_COMP_SWT, clr_mode);
***/
      }
    } else if (COMPARA_SWT == true) {
/****
      pg_color_bar(bias[0], width[0],
                   barxmin[0], barxmax[0], barymin[0], barymax[0],
                   NEGATIVE_COMP_SWT, clr_mode);
****/
    }
  }
#endif /* __DRAW__ */

/*
--
*/

  for (i=0; i<NP; i++) {
    if (SMOOTH_SWT == false) {
#ifdef __DRAW__
      pg_color_map(N, center_nx, center_ny, resx, resy, bias[i], width[i],
                   widthx, widthy, px, py, dist[i],
                   mapxmin[i], mapxmax[i], mapymin[i], mapymax[i],
                   labx[i], laby[i], title[i], NEGATIVE_COMP_SWT, X_REVERSE,
                   clr_mode);
#else
      pgcont_map  (N, center_nx, center_ny, resx, resy, bias[i], width[i],
                   widthx, widthy, px, py, 0.05, dist[i],
                   mapxmin[i], mapxmax[i], mapymin[i], mapymax[i],
                   labx[i], laby[i], title[i], NEGATIVE_COMP_SWT, X_REVERSE,
                   clr_mode);
#endif
    } else if (SMOOTH_SWT == true) {
      modified_distribution(N, dist[i], center_nx, center_ny, resx, resy,
                    MAP_DISP_SIZE, m_dist, &center_mx, &center_my,
                    &mrx, &mry, &n);
#ifdef __DRAW__
      pg_color_map(MAP_DISP_SIZE, center_mx, center_my, mrx, mry,
                   bias[i], width[i],
                   widthx, widthy, px, py, m_dist,
                   mapxmin[i], mapxmax[i], mapymin[i], mapymax[i],
                   labx[i], laby[i], title[i], NEGATIVE_COMP_SWT, X_REVERSE,
                   clr_mode);
#else
      pgcont_map  (MAP_DISP_SIZE, center_mx, center_my, mrx, mry,
                   bias[i], width[i],
                   widthx, widthy, px, py, 0.05, m_dist,
                   mapxmin[i], mapxmax[i], mapymin[i], mapymax[i],
                   labx[i], laby[i], title[i], NEGATIVE_COMP_SWT, X_REVERSE,
                   clr_mode);
#endif
    }
  }

/*
----
*/

  if (SMOOTH_SWT == true) {
    free (m_dist);
  }

/*
--------
*/

  for (i=0; i<NP; i++) {
    if (srch_box_set_swt == true) {
      x0 = N/2+1+(int)lrint((srch_center_x[i] - 0.5 * srch_width_x[i]) / resx);
      x1 = N/2+1+(int)lrint((srch_center_x[i] + 0.5 * srch_width_x[i]) / resx);
      if (x0 < 0) {
        x0 = 0;
      } else if (x0 > N) {
        x0 = N;
      }
      if (x1 < 0) {
        x1 = 0;
      } else if (x1 > N) {
        x1 = N;
      }
      if (x0 > x1) {
        ntmp = x0;
        x0 = x1;
        x1 = ntmp;
      }
      y0 = N/2+1+(int)lrint((srch_center_y[i] - 0.5 * srch_width_y[i]) / resy);
      y1 = N/2+1+(int)lrint((srch_center_y[i] + 0.5 * srch_width_y[i]) / resy);
      if (y0 < 0) {
        y0 = 0;
      } else if (y0 > N) {
        y0 = N;
      }
      if (y1 < 0) {
        y1 = 0;
      } else if (y1 > N) {
        y1 = N;
      }
      if (y0 > y1) {
        ntmp = y0;
        y0 = y1;
        y1 = ntmp;
      }
    } else {
      x0 = 0;
      x1 = N;
      y0 = 0;
      y1 = N;
    }
    peak_normalize(N, pmin+i, pmax+i, noise+i, err_x+i, err_y+i,
                      delta_x+i, delta_y+i, resx, resy,
                      x0, x1, y0, y1,
                      dist[i]);
/********
    printf(  "Maximum Peak: %7.2E\n",    pmax[i]);
    if (noise[i] == 0.0) {
      printf("SNR         : Infinity\n");
    } else {
      printf("SNR         : %lf\n", pmax[i] / noise[i]);
    }
    printf(  "(X, Y)      : (%f, %f)\n", delta_x[i], delta_y[i]);
********/
  }

/*
--------
*/

  return (__GO__);
}



int   modified_distribution(
             int      N, float   *x, int    ocx, int    ocy,
             float  orx, float  ory,
             int      M, float   *y, int   *mcx, int   *mcy,
             float *mrx, float *mry, int     *n)
{

  int    i, j, k, l, m, n2;
  int    imax;
  float  *a, fmax;

  *n = N / M;
  n2 = *n * *n;
  if ((a = (float  *)calloc(n2, sizeof(float))) == NULL) {
    printf("ERROR: modified_distribution: memory allocation error.\n");
    return (__NG__);
  }

  for (i=0; i<M; i++) {
    for (j=0; j<M; j++) {
      m = 0;
      for (k=0; k<*n; k++) {
        for (l=0; l<*n; l++) {
          a[m] = *(x + N * (*n * i + k) + *n * j + l);
          m++;
        }
      }
      imax = 0;
      fmax = fabsf(a[0]);
      for (k=1; k<n2; k++) {
        if (fabsf(a[k]) > fmax) {
          fmax = fabsf(a[k]);
          imax = k;
        }
      }
      *(y + M * i + j) = a[imax];
    }
  }

  free (a);

  *mrx = orx * *n;
  *mry = ory * *n;
  *mcx = ocx / *n;
  *mcy = ocy / *n;

  return (__GO__);
}




int   modified_distribution_XY(
             int     NX, int     NY, float   *A,
             int     MX, int     MY, int    n_bandle,
             float   *B)
{

  int    i, j, k, l, nxy;
  int    imax;
  float  fmax, ftmp;

  nxy = n_bandle * n_bandle;

  for (i=0; i<MX; i++) {
    for (j=0; j<MY; j++) {
      fmax = 0.0;
      for (k=0; k<n_bandle; k++) {
        for (l=0; l<n_bandle; l++) {
          imax = NY * (n_bandle * i + k) + n_bandle * j + l;
          if (imax < NX*NY) {
            ftmp = fabsf(*(A + imax));
            if (ftmp > fmax) {
              fmax = ftmp;
            }
          }
        }
      }
      *(B + MY * i + j) = fmax;
    }
  }

  return (__GO__);
}
