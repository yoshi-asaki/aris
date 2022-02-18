#include <mathtools.h>
#include <cpgplot.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <aris.h>


void  get_vis(float *mapr, float *mapi, int  fmax,
              float *dist, int dmax)
{
  int   i, j;
  int   I, J;
  int   IFREF, JFREF, IDREF, JDREF, LEN;

/*
------------------------------------
*/

  for (i=0; i<fmax*fmax; i++) {
    *(mapr + i) = 0.0;
    *(mapi + i) = 0.0;
  }

/*
------------------------------------
*/

  IFREF = fmax / 2;
  JFREF = fmax / 2;
  IDREF = dmax / 2;
  JDREF = dmax / 2;

  LEN = sizeof(float) * JDREF;
  for (i=0; i<IDREF; i++) {
    I = fmax - IDREF + i;
    memcpy(mapr+I*fmax+fmax-JDREF, dist+i*dmax,       LEN);
    memcpy(mapr+I*fmax,            dist+i*dmax+JDREF, LEN);
  }
  for (i=IDREF; i<dmax; i++) {
    I = i - IDREF;
    memcpy(mapr+I*fmax+fmax-JDREF, dist+i*dmax,       LEN);
    memcpy(mapr+I*fmax,            dist+i*dmax+JDREF, LEN);
  }

/*
------------------------------------
*/

/********
  IFREF = fmax / 2 + 1;
  JFREF = fmax / 2 + 1;
  IDREF = dmax / 2 + 1;
  JDREF = dmax / 2 + 1;

  LEN = sizeof(float) * JDREF;
  for (i=0; i<IDREF; i++) {
    I = fmax - IDREF + i;
    memcpy(mapr+I*fmax+fmax-JDREF, dist+i*dmax,       LEN);
    memcpy(mapr+I*fmax,            dist+i*dmax+JDREF, LEN);
  }
  for (i=IDREF; i<dmax; i++) {
    I = i - IDREF;
    memcpy(mapr+I*fmax+fmax-JDREF, dist+i*dmax,       LEN);
    memcpy(mapr+I*fmax,            dist+i*dmax+JDREF, LEN);
  }
********/

/*
------------------------------------
*/

/****
  IFREF = fmax / 2 + 1;
  JFREF = fmax / 2 + 1;
  IDREF = dmax / 2 + 1;
  JDREF = dmax / 2 + 1;

  for (i=0; i<IDREF; i++) {
    for (j=0; j<JDREF; j++) {
      *(mapr + (IDREF-1-i)*fmax + (JDREF-1-j))
                          = *(dist + i*dmax + j);
    }
  }

  for (i=0; i<IDREF; i++) {
    for (j=JDREF; j<dmax; j++) {
      *(mapr + (IDREF-1-i)*fmax + (fmax-1-(j-JDREF)))
                          = *(dist + i*dmax + j);
    }
  }

  for (i=IDREF; i<dmax; i++) {
    for (j=0; j<JDREF; j++) {
      *(mapr + (fmax-1-(i-IDREF))*fmax + (JDREF-1-j))
                          = *(dist + i*dmax + j);
    }
  }

  for (i=IDREF; i<dmax; i++) {
    for (j=JDREF; j<dmax; j++) {
      *(mapr + (fmax-1-(i-IDREF))*fmax + (fmax-1-(j-JDREF)))
                          = *(dist + i*dmax + j);
    }
  }
****/

/*
----------------------------------------------------------
*/

  fft2d(mapr, mapi, fmax, -1);
  return;
}
