#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <aris.h>


void source_model3(float *dist, int dmax, float res,
                   int   n_source, float *tamp,
                   float *radius_maj, float *radius_min, float *pa,
                   float *xoff, float *yoff, char  *s_type)
{
  int    i, j, dref, I, N;
  float  *DIST;
  float  AMP, x, y, ar, ai;
  float  alpha, r, phs, u, v;
  float  RA2, RB2, ln2;
  int    NN, NCNT;
  int    ncenx, nceny;
  float  factor, sigma_factor=3.0, sigma_factor2=9.0;

/*
--------------------
*/

  dref = dmax / 2;
  DIST = (float *)calloc(dmax*dmax, sizeof(float));
  ln2 = log(2.0);

  for (N=0; N<n_source; N++) {

    for (i=0; i<dmax*dmax; i++) {
      DIST[i] = 0.0;
    }

    ncenx = (int)lrint(xoff[N] / res);
    nceny = (int)lrint(yoff[N] / res);


    if (s_type[N] == 'P') {
      *(DIST + dmax * (dref + ncenx) + dref + nceny) = tamp[N];


    } else if (s_type[N] == 'S') {
      NN = (int)lrint(radius_maj[N] / res);
      phs = -(0.5*dpi - pa[N]);
      ar = cos(phs);
      ai = sin(phs);
      RA2 = radius_maj[N] * radius_maj[N];
      RB2 = radius_min[N] * radius_min[N];

      NCNT = 0;
      for (i=ncenx+dref-NN; i<ncenx+dref+NN; i++) {
        for (j=nceny+dref-NN; j<nceny+dref+NN; j++) {
          u = res * (float)(i - dref) - xoff[N];
          v = res * (float)(j - dref) - yoff[N];
          x = ar * u - ai * v;
          y = ai * u + ar * v;
          if ((x*x/RA2 + y*y/RB2) <= 1.0) {
            *(DIST + dmax * i + j) = tamp[N];
            NCNT++;
          }
        }
      }
      factor = 1.0 / (float)NCNT;
      for (i=ncenx+dref-NN; i<ncenx+dref+NN; i++) {
        for (j=nceny+dref-NN; j<nceny+dref+NN; j++) {
          *(DIST + dmax * i + j) *= factor;
        }
      }


    } else if (s_type[N] == 'G') {
      NN = (int)lrint(sigma_factor * radius_maj[N] / res);
      phs = -(0.5*dpi - pa[N]);
      ar = cos(phs);
      ai = sin(phs);
      RA2 = radius_maj[N] * radius_maj[N];
      RB2 = radius_min[N] * radius_min[N];

      AMP = 0.0;
      for (i=ncenx+dref-NN; i<ncenx+dref+NN; i++) {
        for (j=nceny+dref-NN; j<nceny+dref+NN; j++) {
          u = res * (float)(i - dref) - xoff[N];
          v = res * (float)(j - dref) - yoff[N];
          x = ar * u - ai * v;
          y = ai * u + ar * v;
          if ((x*x/RA2 + y*y/RB2) <= sigma_factor2) {
            *(DIST + dmax * i + j) = exp(-ln2 * (x*x/RA2 + y*y/RB2));
            AMP += *(DIST + dmax * i + j);
          }
        }
      }
      factor = tamp[N] / AMP;
      for (i=ncenx+dref-NN; i<ncenx+dref+NN; i++) {
        for (j=nceny+dref-NN; j<nceny+dref+NN; j++) {
          *(DIST + dmax * i + j) *= factor;
        }
      }
    }

    for (i=0; i<dmax*dmax; i++) {
      *(dist + i) += *(DIST + i);
    }
  }

  free (DIST);

  return;
}
