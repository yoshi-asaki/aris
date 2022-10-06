#include <stdio.h>
#include <cpgplot.h>
#include "./pgtools.h"

void  pg_panel_tile(int NPLANE, int IPLANE, float *fs,
                    int   *NX,    int   *NY,    int   *nx,    int   *ny,
                    float *sxmin, float *sxmax, float *symin, float *symax,
                    float *bxmin, float *bxmax, float *bymin, float *bymax)
{
  int    NTMP;
  float  d1, d2;

/*
----
*/

  if (*NX == 0 && *NY == 0) {
    if (NPLANE == 1) {
      *NX = 1;
    } else if (NPLANE == 2 || NPLANE == 4) {
      *NX = 2;
    } else if (NPLANE == 3 || NPLANE > 4 && NPLANE <= 9) {
      *NX = 3;
    } else if (NPLANE > 9 && NPLANE <= 16) {
      *NX = 4;
    } else if (NPLANE > 16 && NPLANE <= 25) {
      *NX = 5;
    } else if (NPLANE > 25 && NPLANE <= 36) {
      *NX = 6;
    } else if (NPLANE > 36 && NPLANE <= 49) {
      *NX = 7;
    } else if (NPLANE > 49 && NPLANE <= 64) {
      *NX = 8;
    } else if (NPLANE > 64 && NPLANE <= 90) {
      *NX = 9;
    } else if (NPLANE > 91) {
      *NX = 10;
    }
    *NY = NPLANE / *NX;
    if (NPLANE % *NX != 0) {
      *NY++;
    }
  } else {
    if (*NX == 0) {
      *NX = NPLANE / *NY;
      if (NPLANE % *NX != 0) {
        *NY++;
      }
    } else if (*NY == 0) {
      *NY = NPLANE / *NX;
      if (NPLANE % *NX != 0) {
        *NY++;
      }
    }
  }

  if (*NX >= *NY) {
    NTMP = *NX;
  } else {
    NTMP = *NY;
  }

  if (NTMP == 1) {
    d1 = 0.845;
    d2 = 0.840;
    *bxmin = 0.119;
    d1 = 0.795;
    d2 = 0.790;
    *bxmin = 0.120;
    *bymax = 0.925;
    *fs = 0.7;
  } else if (NTMP == 2) {
    d1 = 0.395;
    d2 = 0.390;
    *bxmin = 0.119;
    *bymax = 0.945;
    *fs = 0.7;
  } else if (NTMP == 3) {
    d1 = 0.295;
    d2 = 0.290;
    *bxmin = 0.119;
    *bymax = 0.945;
    *fs = 0.7;
  } else if (NTMP == 4) {
    d1 = 0.205;
    d2 = 0.200;
    *bxmin = 0.119;
    *bymax = 0.945;
    *fs = 0.6;
  } else if (NTMP == 5) {
    d1 = 0.185;
    d2 = 0.180;
    *bxmin = 0.119;
    *bymax = 0.945;
    *fs = 0.6;
  } else if (NTMP == 6) {
    d1 = 0.145;
    d2 = 0.140;
    *bxmin = 0.119;
    *bymax = 0.945;
    *fs = 0.5;
  } else if (NTMP == 7) {
    d1 = 0.120;
    d2 = 0.115;
    *bxmin = 0.119;
    *bymax = 0.945;
    *fs = 0.4;
  } else if (NTMP == 8) {
    d1 = 0.105;
    d2 = 0.100;
    *bxmin = 0.119;
    *bymax = 0.945;
    *fs = 0.4;
  } else if (NTMP == 9) {
    d1 = 0.090;
    d2 = 0.088;
    *bxmin = 0.102;
    *bymax = 0.945;
    *fs = 0.3;
  } else if (NTMP == 10) {
    d1 = 0.085;
    d2 = 0.083;
    *bxmin = 0.102;
    *bymax = 0.945;
    *fs = 0.3;
  }

  cpgsch(*fs);

  *nx = IPLANE % *NX;
  *ny = IPLANE / *NX;

  *bymin = *bymax - d1 * (float)(*NY);

  *bxmax = *bxmin + d1 * (float)(*NX-1) + d2;
  *sxmin = *bxmin + d1 * (float)(*nx);
  *sxmax = *sxmin + d2;
  *symax = *bymax - d1 * (float)(*ny);
  *symin = *symax - d2;

  return;
}
