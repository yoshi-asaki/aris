#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


double weight_calc(double , int   );
int    visibility_pnt_ad_calc(int  ,  int  ,  int   );


int   imager(int     NDAT,   struct baseline_uvw *bluvw, double  *pix_mas,
             int     fmax,   struct fringe *frng,
             float   *dist,  int    dmax)
{
  int    i, j, n, CNT, I, N;
  int    pnt_ad;
  int    nuvmax;
  int    IFREF, JFREF;
  int    IDREF, JDREF;
  int    LEN;
  int    *nx, *ny;
  _Bool  *USED;
  float  *mapr, *mapi, *wght;
  float  taper;
  double uvmax, uvtmp, pix_uvl;
  double mr, mi, wt, MAPR, MAPI, WGHT, w;

/*
-------------------
*/

  IFREF = fmax / 2;
  JFREF = fmax / 2;
  IDREF = dmax / 2;
  JDREF = dmax / 2;

/*
--------
*/

  if ((USED  = (_Bool  *)calloc(NDAT,   sizeof(_Bool)  )) == NULL) {
    printf("ERROR : memory allocation of USED in IMAGER.\n");
    return (__NG__);
  }
  if ((nx    = (int    *)calloc(NDAT,   sizeof(int)  )) == NULL) {
    printf("ERROR : memory allocation of nx   in IMAGER.\n");
    free (USED);
    return (__NG__);
  }
  if ((ny    = (int    *)calloc(NDAT,   sizeof(int)  )) == NULL) {
    printf("ERROR : memory allocation of ny   in IMAGER.\n");
    free (USED);
    free (nx);
    return (__NG__);
  }

/*
---------------------------------------
*/

  if ((mapr = (float  *)calloc(fmax*fmax, sizeof(float))) == NULL) {
    printf("ERROR : memory allocation of mapr in IMAGER.\n");
    free (USED);
    free (nx);
    free (ny);
    return (__NG__);
  }

  if ((mapi = (float  *)calloc(fmax*fmax, sizeof(float))) == NULL) {
    printf("ERROR : memory allocation of mapi in IMAGER.\n");
    free (USED);
    free (nx);
    free (ny);
    free (mapr);
    return (__NG__);
  }

  if ((wght = (float  *)calloc(fmax*fmax, sizeof(float))) == NULL) {
    printf("ERROR : memory allocation of mapi in IMAGER.\n");
    free (USED);
    free (nx);
    free (ny);
    free (mapr);
    free (mapi);
    return (__NG__);
  }

/*
---------------------------------------
*/

  uvmax = 0.0;
  for (i=0; i<NDAT; i++) {
    uvtmp = pow(bluvw[i].u, 2.0) + pow(bluvw[i].v, 2.0);
    if (uvtmp > uvmax) {
      uvmax = uvtmp;
    }
  }
  uvmax = ceil(sqrt(uvmax));
  pix_uvl = uvmax / (double)(fmax / 2);
  *pix_mas = 180.0 / dpi * 3.6e6 / (pix_uvl * (double)fmax);

  for (i=0; i<NDAT; i++) {
    nx[i] = (int)lrint(bluvw[i].u / pix_uvl);
    ny[i] = (int)lrint(bluvw[i].v / pix_uvl);
  }

/*
---------------------------------------
*/

  n = fmax * fmax;
  for (i=0; i<n; i++) {
    *(mapr + i) = 0.0;
    *(mapi + i) = 0.0;
    *(wght + i) = 0.0;
  }

  for (i=0; i<NDAT; i++) {
    USED[i] = false;
  }

  N = 0;
  for (i=0; i<NDAT; i++) {
    if (USED[i] == false) {
      MAPR = frng[i].rl;
      MAPI = frng[i].im;
      WGHT = frng[i].wt;
      CNT = 1;
      for (j=i+1; j<NDAT; j++) {
        if (USED[j] == false) {
          if        (nx[j] ==  nx[i] && ny[j] ==  ny[i]) {
            MAPR += frng[j].rl;
            MAPI += frng[j].im;
            WGHT += frng[j].wt;
            USED[j] = true;
            CNT++;
          } else if (nx[j] == -nx[i] && ny[j] == -ny[i]) {
            MAPR += frng[j].rl;
            MAPI -= frng[j].im;
            WGHT += frng[j].wt;
            USED[j] = true;
            CNT++;
          }
        }
      }

      mr = MAPR / (float)CNT;
      mi = MAPI / (float)CNT;
      wt = weight_calc(WGHT, CNT);

      pnt_ad = visibility_pnt_ad_calc(nx[i], ny[i],  fmax);
      *(mapr + pnt_ad) =  mr;
      *(mapi + pnt_ad) =  mi;
      *(wght + pnt_ad) =  wt;

/**** 2012.02.13
      if (! ((I == 0 && J == 0) ||
             nx[i] + IFREF == 0 ||
             ny[i] + JFREF == 0)) {
        pnt_ad = visibility_pnt_ad_calc(-nx[i], -ny[i],  fmax);
        *(mapr + pnt_ad) =  mr;
        *(mapi + pnt_ad) = -mi;
        *(wght + pnt_ad) =  wt;
      }
****/
      pnt_ad = visibility_pnt_ad_calc(-nx[i], -ny[i],  fmax);
      *(mapr + pnt_ad) =  mr;
      *(mapi + pnt_ad) = -mi;
      *(wght + pnt_ad) =  wt;

      N++;
    }
  }
  free (nx);
  free (ny);

/*
---------------------------------------
*/

/****
  for (i=0; i<fmax; i++) {
    for (j=0; j<fmax; j++) {
      taper = (float)((i-IFREF)*(i-IFREF) + (j-JFREF)*(j-JFREF))
                    / 2.0 / 0.42 / (float)nuvmax;
      taper = exp(taper);
      pnt_ad = i * fmax + j;
      *(mapr + pnt_ad) *= taper;
      *(mapi + pnt_ad) *= taper;
    }
  }
****/

/*
---------------------------------------
*/

  for (i=0; i<fmax; i++) {
    w = 0.0;
    for (j=0; j<fmax; j++) {
      w += *(wght + i * fmax + j);
    }
    for (j=0; j<fmax; j++) {
      if (w != 0.0) {
        *(wght + i * fmax + j) = 1.0 / w;
      }  else {
        *(wght + i * fmax + j) = 0.0;
      }
    }
  }
  wfft2d(mapr, mapi, wght, fmax, 1);

/*
---------------------------------------
*/

  LEN = sizeof(float) * dmax;
  for (i=0; i<dmax; i++) {
    I = IFREF - IDREF + i;
    memcpy(dist + i*dmax, mapr + I*fmax + JFREF - JDREF, LEN);
  }

  free (USED);
  free (mapr);
  free (mapi);

  return (__GO__);
}


double  weight_calc(double  wght,   int  n)
{
/****
  return (wght / (double)n);
****/
  return 1.0;
}



int    visibility_pnt_ad_calc(int  nx,  int  ny,  int fmax)
{
  int     pnt_ad=0, I, J;

  if        (nx >= 0 && ny >= 0) {
    I = nx;
    J = ny;
    pnt_ad = I*fmax + J;
  } else if (nx <  0 && ny >= 0) {
    I = fmax + nx;
    J = ny;
    pnt_ad = I*fmax + J;
  } else if (nx >= 0 && ny <  0) {
    I = nx;
    J = fmax + ny;
    pnt_ad = I*fmax + J;
  } else if (nx <  0 && ny <  0) {
    I = fmax + nx;
    J = fmax + ny;
    pnt_ad = I*fmax + J;
  }

  return pnt_ad;
}
