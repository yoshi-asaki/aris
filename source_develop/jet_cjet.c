#include <stdio.h>
#include <math.h>
#include <cpgplot.h>
#include <aris.h>

#define NA    30000

void jet_cjet(float *dist, int dmax, double THETA,
              double pix_mas, double scale_f)
{
  int    i, j, I, NDAT;
  float  theta, psi, pitch;
  float  x, y;
  float  X[NA], Y[NA], R[NA];
  double T[NA];
  char   s_type[NA];
  float  f;
  int    IAJET, IDJET;
  float  pix_t;

/*
--------
*/

  IAJET = 100;
  IDJET =  10;
  pix_t = scale_f * 1.0e7;
  f = 2.0 * 1.38e-23 * pow(pix_mas / 3600.0e3 / 180.0 * dpi, 2.0)
                     / pow(7.0e-3, 2.0) * 1.0e26;
  pix_t *= f;

/*
===============================================================
*/

  I = 0;

/*
-------- JET --------
*/

  for (j=0; j<IAJET; j++) {
    psi   = random_val1();
    pitch = random_val1();
    for (i=0; i<230; i++) {
      x = ((float)i + 50.0) * 1.0e-3;        /* mas */
      y = -3.0 * x;
      theta = 2.0 * dpi * pitch * (float)i + psi;
      X[I] = x * sin (theta) + y;
      Y[I] = x * cos (theta)    ;
      T[I]  = 2.0 * pix_t;
      R[I]  = 0.010;         /* mas */
      s_type[I] = 'G';
      I++;
    }
  }

/*
-------- COUNTER JET --------
*/

  for (j=0; j<IDJET; j++) {
    psi   = random_val1();
    pitch = random_val1();
    for (i=0; i<140; i++) {
      x = ((float)i + 100.0) * 1.0e-3;        /* mas */
      y =  3.0 * x;
      theta = 2.0 * dpi * pitch * (float)i + psi;
      X[I] = x * sin (theta) + y;
      Y[I] = x * cos (theta)    ;
      T[I]  = 2.0 * pix_t;
      R[I]  = 0.010;         /* mas */
      s_type[I] = 'G';
      I++;
    }
  }

  NDAT = I;

  for (i=0; i<NDAT; i++) {
    X[i] *= 5.0;
    Y[i] *= 5.0;
    R[i] *= 5.0;
  }
  PA_rotate(NDAT, X, Y, THETA);

#ifdef __DEBUG__
  for (I=0; I<NDAT; I++) {
    printf("  X[%5d]=%9.6f; Y[%5d]=%9.6f; T[%5d]=%fe-3;\n",
           I, X[I], I, Y[I], I, 1000.0*T[I]);
    printf("  R[%5d]=%f; s_type[%5d]='G';\n", I, R[I], I);
  }
#endif

/*
===============================================================
*/

  source_model(dist, dmax, pix_mas, NDAT, T, R, X, Y, s_type);
  for (i=0; i<NDAT; i++) {
    s_type[i] = 'S';
    R[i] *= 0.5;
  }
  source_model(dist, dmax, pix_mas, NDAT, T, R, X, Y, s_type);

}
