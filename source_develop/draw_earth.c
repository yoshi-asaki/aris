#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __ORBD_ANIMATION__
****/
#ifdef __ORBD_ANIMATION__
  #define __ANIMA_CLOCK__
  #define __TRACKING_STATUS__
#endif

void  draw_earth(double Re, 
                 double epsiron,
                 double GST_theta,
                 struct source_parameter src,
                 struct source_parameter sun,
                 struct earth_shape      *earth_shape,
                 int *clrplt,
                 _Bool SHADOW_SWT,
                 _Bool COLOR_SWT,
                 _Bool AXIS_SWT)
{
  int    i, j, k;
  int    ns, ne;
  int    ndiv;
  float  theta, fai, r, ang1, ang2, ang3, ang4;
  float  fpos[3], pgx[721], pgy[721], pgz[721];
  double DPI;
  double epos[3], spos[3], A[9], q[4];

/*
--------
*/

  DPI = dpi / 180.0;
  r = Re * 1.0e-3;

  ang1 = -src.RA - dpi/2.0;
  ang2 =  src.DC - dpi/2.0;

  earth_shape->n_day = 0;
  earth_shape->n_ngt = 0;
  earth_shape->n_shd = 0;

#ifdef __ORBD_ANIMATION__
  ndiv = 4;
#else
  ndiv = 1;
#endif /* __ORBD_ANIMATION__ */

/*
---- SHADOW ------------------------
*/

  epos[0] =  0.0;
  epos[1] = -sin(epsiron);
  epos[2] =  cos(epsiron);

  for (i=0; i<=360; i++) {
    coordinate_rotation(sun.s, (double)i*DPI, q, A, epos, spos);
    drotate(spos, (double)ang1, "z");
    drotate(spos, (double)ang2, "x");
    pgx[i] = (float)(r * spos[0]);
    pgy[i] = (float)(r * spos[1]);
    pgz[i] = (float)(r * spos[2]);
  }

  ns = 0;
  ne = 0;
  for (i=0; i<360; i++) {
    if (pgz[i] <  0.0 && pgz[i+1] >= 0.0 ||
        pgz[i] <= 0.0 && pgz[i+1] >  0.0) {
      ns = i;
    } else if (pgz[i] >  0.0 && pgz[i+1] <= 0.0 ||
               pgz[i] >= 0.0 && pgz[i+1] <  0.0) {
      ne = i;
    }
  }

  if (ns <= ne) {
    j = 0;
    for (i=ns; i<=ne; i++) {
      pgx[j] = pgx[i];
      pgy[j] = pgy[i];
      j++;
    }
  } else {
    j = 0;
    for (i=ns; i<360; i++) {
      pgz[j] = pgx[i];
      j++;
    }
    for (i=0; i<ne; i++) {
      pgz[j] = pgx[i];
      j++;
    }
    for (i=0; i<j; i++) {
      pgx[i] = pgz[i];
    }

    j = 0;
    for (i=ns; i<360; i++) {
      pgz[j] = pgy[i];
      j++;
    }
    for (i=0; i<ne; i++) {
      pgz[j] = pgy[i];
      j++;
    }
    for (i=0; i<j; i++) {
      pgy[i] = pgz[i];
    }
  }

  ang4 = atan2(pgy[0],   pgx[0]);
  k = (int)lrint(atan2(pgy[j-1], pgx[j-1]) / DPI);
  while (1) {
    fai = (float)k * DPI;
    fpos[0] = cos(fai);
    fpos[1] = sin(fai);
    fpos[2] = 0.0;
    frotate(fpos, -ang2, "x");
    frotate(fpos, -ang1, "z");
    ang3 = fpos[0]*sun.s[0] + fpos[1]*sun.s[1] + fpos[2]*sun.s[2];
    if (ang3 < 0.0) {
      pgx[j] = r * cos(fai);
      pgy[j] = r * sin(fai);
      j++;
    }
    if (fabs(fai - ang4) <= DPI) {
      break;
    } else {
      k--;
      if (k <= -180) {
        k += 360;
      }
    }
  }

  earth_shape->n_shd = j;
  for (i=0; i<earth_shape->n_shd; i++) {
    earth_shape->s_shd[0][i] = pgx[i];
    earth_shape->s_shd[1][i] = pgy[i];
  }

  if (SHADOW_SWT == true) {
    cpgscr(20, 0.5, 0.5, 0.5);
    cpgsci(20);
    cpgpoly(earth_shape->n_shd, earth_shape->s_shd[0], earth_shape->s_shd[1]);
    cpgsci(1);
  }

/*
---- LONGITUDE ---------------------
*/

  for (i=0; i<12; i++) {
    fai = 30.0 * (float)i * DPI + GST_theta;
    for (j=0; j<360/ndiv; j++) {
      theta = (float)(ndiv*j) * DPI;
      fpos[0] = cos(theta) * cos(fai);
      fpos[1] = cos(theta) * sin(fai);
      fpos[2] = sin(theta);
      ang3 = fpos[0]*sun.s[0] + fpos[1]*sun.s[1] + fpos[2]*sun.s[2];
      frotate(fpos, ang1, "z");
      frotate(fpos, ang2, "x");
      if (fpos[2] >= 0.0) {
        fpos[0] *= r;
        fpos[1] *= r;
        if (ang3 >= 0.0) {
          cpgsci(clrplt[8]);
          cpgsci(clrplt[4]);
          earth_shape->s_day[0][earth_shape->n_day] = fpos[0];
          earth_shape->s_day[1][earth_shape->n_day] = fpos[1];
          (earth_shape->n_day)++;
        } else {
          if (COLOR_SWT == true) {
            cpgsci(clrplt[4]);
          } else {
            cpgsci(clrplt[1]);
          }
          earth_shape->s_ngt[0][earth_shape->n_ngt] = fpos[0];
          earth_shape->s_ngt[1][earth_shape->n_ngt] = fpos[1];
          (earth_shape->n_ngt)++;
        }
        cpgpt(1, fpos, fpos+1, 1);
      }
    }
  }

/*
---- LATITUDE ----------------------
*/

  for (i=-2; i<3; i++) {
    theta = 30.0 * (float)i * DPI;
    for (j=0; j<360/ndiv; j++) {
      fai = (float)(ndiv*j) * DPI;
      fpos[0] = cos(theta) * cos(fai);
      fpos[1] = cos(theta) * sin(fai);
      fpos[2] = sin(theta);
      ang3 = fpos[0]*sun.s[0] + fpos[1]*sun.s[1] + fpos[2]*sun.s[2];
      frotate(fpos, ang1, "z");
      frotate(fpos, ang2, "x");
      if (fpos[2] >= 0.0) {
        fpos[0] *= r;
        fpos[1] *= r;
        if (ang3 >= 0.0) {
          cpgsci(clrplt[8]);
          cpgsci(clrplt[4]);
          earth_shape->s_day[0][earth_shape->n_day] = fpos[0];
          earth_shape->s_day[1][earth_shape->n_day] = fpos[1];
          (earth_shape->n_day)++;
        } else {
          if (COLOR_SWT == true) {
            cpgsci(clrplt[4]);
          } else {
            cpgsci(clrplt[1]);
          }
          earth_shape->s_ngt[0][earth_shape->n_ngt] = fpos[0];
          earth_shape->s_ngt[1][earth_shape->n_ngt] = fpos[1];
          (earth_shape->n_ngt)++;
        }
        cpgpt(1, fpos, fpos+1, 1);
      }
    }
  }

/*
---- LIM ---------------------------
*/

  for (i=0; i<=360; i++) {
    fai = (float)i * DPI;
    fpos[0] = cos(fai);
    fpos[1] = sin(fai);
    fpos[2] = 0.0;
    frotate(fpos, -ang2, "x");
    frotate(fpos, -ang1, "z");
    ang3 = fpos[0]*sun.s[0] + fpos[1]*sun.s[1] + fpos[2]*sun.s[2];
    if (ang3 >= 0.0) {
      cpgsci(clrplt[8]);
      cpgsci(clrplt[4]);
      earth_shape->s_day[0][earth_shape->n_day] = r * cos(fai);
      earth_shape->s_day[1][earth_shape->n_day] = r * sin(fai);
      (earth_shape->n_day)++;
    } else {
      if (COLOR_SWT == true) {
        cpgsci(clrplt[4]);
      } else {
        cpgsci(1);
      }
      earth_shape->s_ngt[0][earth_shape->n_ngt] = r * cos(fai);
      earth_shape->s_ngt[1][earth_shape->n_ngt] = r * sin(fai);
      (earth_shape->n_ngt)++;
    }
    pgx[0] = r * cos(fai);
    pgy[0] = r * sin(fai);
    cpgpt(1, pgx, pgy, 1);
  }

/*
---- AXIS --------------------------
*/

  if (AXIS_SWT == true) {
    cpgsci(1);
    for (i=0; i<3; i++) {
      fpos[0] = 0.0;
      fpos[1] = 0.0;
      fpos[2] = 0.0;
      if (i == 0) {
        fpos[0] = 3.0 * Re;
      } else if (i == 1) {
        fpos[1] = 3.0 * Re;
      } else if (i == 2) {
        if (src.DC >= 0.0) {
          fpos[2] = 3.0 * Re;
        } else {
          fpos[2] = -3.0 * Re;
        }
      }

      frotate(fpos, ang1, "z");
      frotate(fpos, ang2, "x");
      if (fpos[2] >= 0.0 ||
          fpos[2] < 0.0 && sqrt(fpos[0]*fpos[0]+fpos[1]*fpos[1]) > Re) {
        pgx[1] = fpos[0] * 1.0e-3;
        pgy[1] = fpos[1] * 1.0e-3;

        fpos[0] = 0.0;
        fpos[1] = 0.0;
        fpos[2] = 0.0;
        if (i == 0) {
          fpos[0] = Re;
        } else if (i == 1) {
          fpos[1] = Re;
        } else if (i == 2) {
          fpos[2] = Re;
        }
        frotate(fpos, ang1, "z");
        frotate(fpos, ang2, "x");
        if (fpos[2] >= 0.0) {
          pgx[0] = fpos[0] * 1.0e-3;
          pgy[0] = fpos[1] * 1.0e-3;
        } else {
          fai = atan2(fpos[1], fpos[0]);
          pgx[0] = r * cos(fai);
          pgy[0] = r * sin(fai);
        }
        cpgsls(4);
        cpgline(2, pgx, pgy);
        cpgsls(1);
        if (i == 0) {
          cpgptxt(pgx[1], pgy[1], 0.0, 0.5, "Vernal Equinox\0");
        } else if (i == 1) {
          cpgtext(pgx[1], pgy[1], "6h\0");
        } else if (i == 2) {
          if (src.DC >= 0.0) {
            cpgtext(pgx[1], pgy[1], "North Pole\0");
          } else {
            cpgtext(pgx[1], pgy[1], "South Pole\0");
          }
        }
      }
    }
  }

  cpgsci(1);
  cpgsch(1.0);

  return;
}
