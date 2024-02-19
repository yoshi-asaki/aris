#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"

void precession(int  *TimUTC,  double dT, double *s, int swt)
{
/****
  Base on IAU 2000 Precession-Nutation Model (IERS Conventions, 2003)
****/

  double tueta, theta, z;

  precession_calc(TimUTC, dT, &tueta, &theta, &z);

/*
---- P^{-1} ( CRS -> TRS ) --
*/

  if (swt > 0) {
    drotate(s, -tueta, "z");
    drotate(s, theta,  "y");
    drotate(s, -z,     "z");

/*
---- P      ( TRS -> CRS ) --
*/

  } else {
    drotate(s,  z,     "z");
    drotate(s, -theta, "y");
    drotate(s,  tueta, "z");
  }

  return;
}
