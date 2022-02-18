#include <stdio.h>
#include <math.h>
#include <aris.h>

double  SEFD(double Trx, double Tsky, double Dm, double Ae)
{
  return (2.0 * BOLTZ * (Trx + Tsky) / dpi / Ae / pow(0.5*Dm, 2.0));
}
