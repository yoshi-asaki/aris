#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <aris.h>


void source_model(float *dist, int dmax, float res,
                  int   n_source,
                  double *amp, float *radius,
                  float *xoff, float *yoff, char  *s_type)
{
  int    i, j, dref, N;
  int    ix0, ix1, iy0, iy1;
  float  x, y, R, R2;
  float  r_factor=1.0, r;

  dref = dmax / 2;

  for (N=0; N<n_source; N++) {
    if (radius[N] == 0.0) {
      *(dist + dmax * (dref + (int)lrint(xoff[N]/res))
                     + dref + (int)lrint(yoff[N]/res)) += (float)amp[N];
    } else {
      if (n_source > 5) {
        if (s_type[N] == 'G') {
          r_factor = 8.0;
        } else if (s_type[N] == 'S') {
          r_factor = 2.0;
        }

        r = r_factor * radius[N];
        ix0 = dref + (int)lrint((xoff[N] - r) / res);
        if (ix0 < 0) {
          ix0 = 0;
        }
        ix1 = dref + (int)lrint((xoff[N] + r) / res);
        if (ix1 > dmax) {
          ix1 = dmax;
        }
        iy0 = dref + (int)lrint((yoff[N] - r) / res);
        if (iy0 < 0) {
          iy0 = 0;
        }
        iy1 = dref + (int)lrint((yoff[N] + r) / res);
        if (iy1 > dmax) {
          iy1 = dmax;
        }
      } else {
        ix0 = 0;
        ix1 = dmax;
        iy0 = 0;
        iy1 = dmax;
      }

      for (i=ix0; i<ix1; i++) {
        for (j=iy0; j<iy1; j++) {
          x = res * (float)(i - dref) - xoff[N];
          y = res * (float)(j - dref) - yoff[N];
          R2 = x*x + y*y;
          R  = sqrt(R2);

          if (s_type[N] == 'G') {
            *(dist + i*dmax + j) += (float)amp[N]
                                  * exp(-log(2.0)/radius[N]/radius[N]*R2);
          } else if (s_type[N] == 'S') {
            if (R <= radius[N]) {
              *(dist + i*dmax + j) += (float)amp[N];
            }
          }
        }
      }
    }
  }

  return;
}
