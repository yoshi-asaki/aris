#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


void vis_smoothing(double *mr, double *mi, float *mapr, float *mapi,
                   double  su, double  sv,
                   int NFMAX, int IFREF, int JFREF, double pix_uvl)
{
  int    nxtmp, nytmp, ntmp;
  float  xtmp, ytmp, pxtmp, pytmp, ptmp;
  float  c0, cp1=0.0, cm1=0.0, cx, cy;
  float  *ad_tmp;


  xtmp = (float)su / (float)pix_uvl + (float)IFREF;
  ytmp = (float)sv / (float)pix_uvl + (float)JFREF;
  nxtmp = (int)xtmp;
  nytmp = (int)ytmp;
  pxtmp = xtmp - (float)nxtmp;
  pytmp = ytmp - (float)nytmp;

/*
--------------------------
*/

  ad_tmp = mapr;
  c0  = *(ad_tmp + nxtmp*NFMAX + nytmp);

  ntmp   = nxtmp;
  ptmp   = pxtmp;
  if (ntmp != 0) {
    cm1 = *(ad_tmp + (nxtmp-1)*NFMAX + nytmp);
  }
  if (ntmp != NFMAX - 1) {
    cp1 = *(ad_tmp + (nxtmp+1)*NFMAX + nytmp);
  }
  if (ntmp != 0 && ntmp != NFMAX - 1) {
    cx = c0 + 0.5 * ptmp * (cp1 - cm1);
  } else if (ntmp == 0) {
    cx = c0 + ptmp * (cp1 - c0);
  } else if (ntmp == NFMAX - 1) {
    cx = c0 + ptmp * (c0 - cm1);
  }

  ntmp   = nytmp;
  ptmp   = pytmp;
  if (ntmp != 0) {
    cm1 = *(ad_tmp + nxtmp*NFMAX + nytmp-1);
  }
  if (ntmp != NFMAX - 1) {
    cp1 = *(ad_tmp + nxtmp*NFMAX + nytmp+1);
  }
  if (ntmp != 0 && ntmp != NFMAX - 1) {
    cy = c0 + 0.5 * ptmp * (cp1 - cm1);
  } else if (ntmp == 0) {
    cy = c0 + ptmp * (cp1 - c0);
  } else if (ntmp == NFMAX - 1) {
    cy = c0 + ptmp * (c0 - cm1);
  }

  *mr = 0.5 * (cx + cy);

/*
--------------------------
*/

  ad_tmp = mapi;
  c0  = *(ad_tmp + nxtmp*NFMAX + nytmp);

  ntmp   = nxtmp;
  ptmp   = pxtmp;
  if (ntmp != 0) {
    cm1 = *(ad_tmp + (nxtmp-1)*NFMAX + nytmp);
  }
  if (ntmp != NFMAX - 1) {
    cp1 = *(ad_tmp + (nxtmp+1)*NFMAX + nytmp);
  }
  if (ntmp != 0 && ntmp != NFMAX - 1) {
    cx = c0 + 0.5 * ptmp * (cp1 - cm1);
  } else if (ntmp == 0) {
    cx = c0 + ptmp * (cp1 - c0);
  } else if (ntmp == NFMAX - 1) {
    cx = c0 + ptmp * (c0 - cm1);
  }

  ntmp   = nytmp;
  ptmp   = pytmp;
  if (ntmp != 0) {
    cm1 = *(ad_tmp + nxtmp*NFMAX + nytmp-1);
  }
  if (ntmp != NFMAX - 1) {
    cp1 = *(ad_tmp + nxtmp*NFMAX + nytmp+1);
  }
  if (ntmp != 0 && ntmp != NFMAX - 1) {
    cy = c0 + 0.5 * ptmp * (cp1 - cm1);
  } else if (ntmp == 0) {
    cy = c0 + ptmp * (cp1 - c0);
  } else if (ntmp == NFMAX - 1) {
    cy = c0 + ptmp * (c0 - cm1);
  }

  *mi = 0.5 * (cx + cy);

/****
The following part should be verified. 
****/
  *mi *= -1.0;

  return;
}
