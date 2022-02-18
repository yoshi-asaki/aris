#include <stdio.h>
#include <math.h>
#include <cpgplot.h>
#include <aris.h>

#define NA     80000
#define NDISK    300

void ADAF_disk(float *dist, int ND, double THETA,
                double pix_mas, double scale_f)
{
  int    i, j, k, l;
  int    I, NDAT;
  int    nc[10];
  float  X[NA], Y[NA];
  float  R[NA];
  float  pix_t;
  double T[NA];
  char   s_type[NA];

  float  ti, ang1, ang2, theta;
  float  A_mas, B_mas;
  float  tmp, pitch, psi;
  double gnoise[NDISK];
  float  ra, rb;
  float  xoff, yoff;
  float  f;

  pix_t = scale_f * 1.0e9;
  f = 2.0 * 1.38e-23 * pow(pix_mas / 3600.0e3 / 180.0 * dpi, 2.0)
                     / pow(7.0e-3, 2.0) * 1.0e26;
  pix_t *= f;

/*
=======================================================
*/

  A_mas   = 0.020;         /* mas */
  B_mas   = 0.060;         /* mas */
  ti      = 0.05;

  for (i=0; i<NDISK; i++) {
    gnoise[i] = 0.5 * (1.0 + gauss_dev());
  }

  I = 0;
  for (i=0; i<NDISK; i++) {
    ang1 = 2.0 * dpi * random_val1();
    ang2 = 2.0 * dpi * random_val1();
    if (ang1 > ang2) {
      tmp = ang1;
      ang1 = ang2;
      ang2 = tmp;
    }

    xoff = 0.4 * A_mas * random_val1();
    yoff = 0.2 * B_mas * random_val1();

    ra = A_mas * gnoise[i];
    rb = B_mas * gnoise[i];

    j = 0;
    while (1) {
      theta = ang1 + ti * (float)j;
      if (theta >= ang2) {
        break;
      }
      X[I]      = ra * cos(theta) + xoff;
      Y[I]      = rb * sin(theta) + yoff;
      T[I]      = pix_t;
      R[I]      = 0.010;         /* mas */
      s_type[I] = 'G';
      I++;
      j++;
    }
  }
  NDAT = I;

  PA_rotate(NDAT, X, Y, THETA);

  source_model(dist, ND, pix_mas, NDAT, T, R, X, Y, s_type);
  for (i=0; i<NDAT; i++) {
    s_type[i] = 'S';
    R[i] *= 0.13;
  }
  source_model(dist, ND, pix_mas, NDAT, T, R, X, Y, s_type);

/*
=======================================================
*/

}
