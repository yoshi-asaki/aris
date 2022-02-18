#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <aris.h>


int    fit_power_two_number
                   (double L, double pix, double *L_m, int *n_power)
{
  int    itmp;
  double lftmp1, lftmp2;

/*
-------------------------------------
*/

  itmp = (int)(log10(L/pix) / log10(2.0));
  lftmp1 = pix * pow(2.0, (double)(itmp  ));
  lftmp2 = 2.0 * lftmp1;

  if (fabs(lftmp2 - L) >= fabs(L - lftmp1)) {
    *n_power = itmp;
    *L_m = lftmp1;
  } else {
    *n_power = itmp + 1;
    *L_m = lftmp2;
  }

  return 1;
}
