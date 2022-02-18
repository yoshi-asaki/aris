#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <cpgplot.h>
#include <aris.h>


void pgcont_map(int   N,
                  int   cenx,        int   ceny,
                  float resx,        float resy,
                  float bias,        float width,
                  float xwidth,      float ywidth,
                  float px,          float py,
                  float min_lev,
                  float *dist,
                  float x1, float x2, float y1, float y2,
                  char *labelx, char *labely, char *title,
                  _Bool NEGATIVE_COMP_SWT, _Bool X_REVERSE)
{
  int    i, j;
  int    nxs, nxe, nys, nye;
  float  xmin, xmax, ymin, ymax;
  float  e, f, f_max, tx[2], ty[2];
  float  CR, CB, CG;
  float  lev, tr[6];

/*
---------------------------
*/

  nxs = cenx - (int)(0.5 * xwidth / resx);
  nxe = cenx + (int)(0.5 * xwidth / resx);
  nys = ceny - (int)(0.5 * ywidth / resy);
  nye = ceny + (int)(0.5 * ywidth / resy);

  xmin = (float)(nxs - N/2) * resx + px;
  ymin = (float)(nys - N/2) * resy + py;
  xmax = (float)(nxe - N/2) * resx + px;
  ymax = (float)(nye - N/2) * resy + py;

  cpgsvp(x1, x2, y1, y2);
  if (X_REVERSE == false) {
    cpgswin(xmin, xmax, ymin, ymax);
  } else if (X_REVERSE == true) {
    cpgswin(xmax, xmin, ymin, ymax);
  }

  if (nxs < 0) {
    nxs = 0;
  }
  if (nxe > N) {
    nxe = N;
  }
  if (nys < 0) {
    nys = 0;
  }
  if (nye > N) {
    nye = N;
  }

/*
-------------------------------------------
*/

  for (i=nxs; i<nxe; i++) {
    for (j=nys; j<nye; j++) {
      f = *(dist + N*i + j);
      if (NEGATIVE_COMP_SWT == true) {
        if (bias >= 0.0) {
          e = (f - bias) / width;
        } else if (bias < 0.0) {
          if (bias + width < 0.0) {
            e = (f - bias) / width;
          } else {
            if (f < 0.0) {
              e = f / fabs(bias);
            } else {
              e = f / (width - fabs(bias));
            }
          }
        }
      } else if (NEGATIVE_COMP_SWT == false) {
        e = (f - bias) / width;
      }
    }
  }

/*
-------------------------------------------
*/

  cpgsci(1);
  cpgbox("BCNT", 0.0, 0, "BCNT", 0.0, 0);
/****
  cpgbox("BCT", 0.0, 0, "BCT", 0.0, 0);
****/
  cpglab(labelx, labely, title);

  f = 0.0;
  f_max = 0.0;
  for (i=0; i<N*N; i++) {
    f += pow(*(dist + i), 2.0);
    if (*(dist + i) > f_max) {
      f_max = *(dist + i);
    }
  }
  f = sqrt(f) / (float)(N * N);

  tr[3] =  -resx * (float)(N / 2 + 1) + px;
  tr[4] =   resx;
  tr[5] =   0.0;

  tr[0] =  -resy * (float)(N / 2 + 1) + py;
  tr[1] =   0.0;
  tr[2] =   resy;

  if (min_lev < 0.0) {
    cpgscr(20, 0.6, 0.6, 0.6);
    cpgsci(20);
    lev = f_max * min_lev;
    cpgcont(dist, N, N, 1, N, 1, N, &lev, -1, tr);
    min_lev = fabs(min_lev);
    cpgsci(1);
  }

  i = 0;
  while (1) {
    lev = f_max * min_lev * (float)pow(2.0, (double)i);
    if (lev >= f_max) {
      break;
    }
    cpgcont(dist, N, N, 1, N, 1, N, &lev, -1, tr);
    i++;
  }

  return;
}
