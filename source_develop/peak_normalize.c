#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <aris.h>


void   peak_normalize(int  N,                    
                      float *tmin,      float *tmax,
                      float *noise,
                      float *err_x,     float *err_y,
                      float *delta_x,   float *delta_y,
                      float resx,       float resy,
                      int   x0,         int   x1,
                      int   y0,         int   y1,
                      float *dist)
{
  int   i, j, k;
  int   NX, NY;
  float XMAX, YMAX, y[3];
  float ax, bx, cx;
  float ay, by, cy;

  NX = 0;
  NY = 0;
  *tmin = *dist;
  *tmax = *dist;
  for (i=x0; i<x1; i++) {
    for (j=y0; j<y1; j++) {
      k = i*N + j;
      if (*(dist + k) > *tmax) {
        *tmax = *(dist + k);
        NX = i;
        NY = j;
      }
      if (*(dist + k) < *tmin) {
        *tmin = *(dist + k);
      }
    }
  }

/*
--------------------------------------------------
*/

  if (NX != 0 && NX != N - 1) {
    y[0] = *(dist +N*(NX - 1) + NY);
    y[1] = *(dist +N* NX      + NY);
    y[2] = *(dist +N*(NX + 1) + NY);
    parabo(delta_x, &XMAX, (float)(NX - N/2)*resx, resx, y, &ax, &bx, &cx);
  } else {
    XMAX = *(dist +N*NX + NY);
    *delta_x = (float)(NX - N/2) * resx;
  }

  if (NY != 0 && NY != N - 1) {
    y[0] = *(dist +N* NX      + NY - 1);
    y[1] = *(dist +N* NX      + NY    );
    y[2] = *(dist +N* NX      + NY + 1);
    parabo(delta_y, &YMAX, (float)(NY - N/2)*resy, resy, y, &ay, &by, &cy);
  } else {
    YMAX = *(dist +N*NX + NY);
    *delta_y = (float)(NY - N/2) * resy;
  }

  *tmax  = 0.5 * (XMAX + YMAX);
  *err_x = sqrt( fabs((XMAX - *tmax) / ax) );
  *err_y = sqrt( fabs((YMAX - *tmax) / ay) );

  *noise = 0.0;
  for (i=0; i<N/4; i++) {
    for (j=0; j<N/4; j++) {
      *noise += pow(*(dist +i*N         + j),       2.0);
      *noise += pow(*(dist +(i+3/4*N)*N + j),       2.0);
      *noise += pow(*(dist +i*N         + j+3/4*N), 2.0);
      *noise += pow(*(dist +(i+3/4*N)*N + j+3/4*N), 2.0);
    }
  }
  *noise /= (float)(N*N/4);
  *noise = sqrt(*noise);

  return;
}
