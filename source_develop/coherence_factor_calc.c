#include <stdio.h>
#include <math.h>
#include <aris.h>

#define N     10000
#define NORDER    8

double coherence_factor_calc(double nu, double T, int FRQSTD_1, int FRQSTD_2)
{
  int    i, j, M;
  double a, cf;
  double dt, tau, w;

/*
--------
*/

  M = (int)ceil(T) * N;
  dt = T / (double)M;
  w  = pow(nu*dpi, 2.0);

  cf = 0.0;
  for (i=1; i<M; i++) {
    tau = (double)i * dt;
    a = 0.0;
    for (j=0; j<=NORDER; j++) {
      a += (pow(asd_model(pow(2.0, (double)j)*tau, FRQSTD_1), 2.0)
          + pow(asd_model(pow(2.0, (double)j)*tau, FRQSTD_2), 2.0));
    }
    cf += (T - tau) * exp(-w * tau * tau * a) * dt;
  }
  cf *= 2.0;
  cf = sqrt(cf);
  cf /= T;

/*
--------
*/

  return cf;
}
