#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"


void nutation(int  *TimUTC,  double dT, double *s, int swt)
{
  double epsiron, Delta_epsiron, Delta_psi;

  nutation_calc(TimUTC, dT, &epsiron, &Delta_epsiron, &Delta_psi);

/*
---- N^{-1} ( CRS -> TRS ) ----
*/

  if (swt > 0) {
    drotate(s, -epsiron,               "x");
    drotate(s,  Delta_psi,             "z");
    drotate(s,  epsiron+Delta_epsiron, "x");

/*
---- N      ( TRS -> CRS ) ----
*/

  } else {
    drotate(s, -epsiron-Delta_epsiron, "x");
    drotate(s, -Delta_psi,             "z");
    drotate(s,  epsiron,               "x");
  }

  return;
}
