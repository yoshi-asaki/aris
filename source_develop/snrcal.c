#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <aris.h>

#define BOLTZ  1.381e-23

double  snrcal(double flux, double cfact, double d1,  double d2,
               double ae1,  double ae2,   double ts1, double ts2,
               double bw,   double tau)
{
  flux *= 1.0e-26;

  return (
        flux * dpi * cfact * d1 * d2 / (8.0 * BOLTZ) *
            sqrt(ae1 * ae2 / (ts1 * ts2)) * sqrt(2.0 * bw * tau)
         );
}
