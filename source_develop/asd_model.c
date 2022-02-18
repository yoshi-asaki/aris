#include <stdio.h>
#include <math.h>
#include <aris.h>


double asd_model(double tau, int FRQSTD)
{
  double asd_value;

/*
--------
*/

  if (FRQSTD == H_M) {
    asd_value = 1.747e-16 + 8.0123e-14 * pow(tau, -0.7)
              + 8.7699e-21 * tau;
    asd_value *= sqrt(2.0);
    /**** Nand et al., 2011 ****/

  } else if (FRQSTD == CSO_10) {
    asd_value = 5.2e-15 / tau + 3.6e-15 / sqrt(tau) + 4.0e-16
              + 7.0e-20 * tau;
    /**** Nand et al., 2011 ****/

  } else if (FRQSTD == CSO_100) {
    asd_value = 2.1e-15 / sqrt(tau) + 3.0e-16
              + 7.0e-20 * tau;
    /**** Nand et al., 2011 ****/

  } else if (FRQSTD == TRP3) {
    if (tau < 100.0) {
      asd_value = 0.4e-13;
    } else {
      asd_value = 0.4e-12 / sqrt(tau);
    }

  } else {
    asd_value = 0.0;
  }

/*
--------
*/

  return asd_value;
}
