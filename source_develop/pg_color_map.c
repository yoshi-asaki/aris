#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <aris.h>


void pg_color_map(int   N,
                  int   cenx,        int   ceny,
                  float resx,        float resy,
                  float bias,        float width,
                  float xwidth,      float ywidth,
                  float px,          float py,
                  float *dist,
                  float x1, float x2, float y1, float y2,
                  char *labelx, char *labely, char *title,
                  _Bool  NEGATIVE_COMP_SWT, _Bool X_REVERSE,
                  char *mode)
{
  int    i, j;
  int    nxs, nxe, nys, nye;
  float  xmin, xmax, ymin, ymax;
  float  e, f, tx[2], ty[2];
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

  e = 0.0;
  if (NEGATIVE_COMP_SWT == true) {
    if (bias >= 0.0) {
      e = - bias / width;
    } else if (bias < 0.0) {
      if (bias + width < 0.0) {
        e = - bias / width;
      } else {
        e = 0.0;
      }
    }
  } else if (NEGATIVE_COMP_SWT == false) {
    e =  bias / width;
  }

  set_color(e, &CR, &CG, &CB, NEGATIVE_COMP_SWT, mode);
  cpgscr(10, CR, CG, CB);
  cpgsci(10);
  cpgrect(xmin, xmax, ymin, ymax);

/*
-------------------------------------------
*/

  tr[0] =  -resy * (float)(N / 2 + 1) + py;
  tr[1] =   0.0;
  tr[2] =   resy;

  tr[3] =  -resx * (float)(N / 2 + 1) + px;
  tr[4] =   resx;
  tr[5] =   0.0;

  e = *dist;
  f = *dist;
  for (i=1; i<N*N; i++) {
    if (e > *(dist+i)) {
      e = *(dist + i);
    }
    if (f < *(dist+i)) {
      f = *(dist + i);
    }
  }
  if (strncmp(mode, "clr", 3) == 0) {
    palett(2, 1.0, 0.5);
  } else if (strncmp(mode, "b2w", 3) == 0) {
    palett(1, 1.0, 0.5);
  } else if (strncmp(mode, "w2b", 3) == 0) {
    palett(1, 1.0, 0.5);
  }

  cpgimag(dist, N, N, 1, N, 1, N, e, f, tr);
  cpgsci(1);
/****
  cpgwedg("BI", 4.0, 5.0, e, f, "");
  cpgbox("BCNT", 0.0, 0, "BCNT", 0.0, 0);
  cpglab(labelx, labely, title);
****/
  cpgbox("BC", 0.0, 0, "BC", 0.0, 0);

  return;
}
