#include <stdio.h>
#include <math.h>
#include <aris.h>


double modified_asd_model(double tau, int FRQSTD)
{
  double asd_value;

/*
--------
*/

  if (FRQSTD == H_M) {
    asd_value = 1.747e-16 + 8.0123e-14 * pow(tau, -0.7)
              + 0.7699e-23 * pow(tau, 1.5);
    asd_value *= sqrt(2.0);
    /**** J. Hartnett 2012, private communication ****/

  } else if (FRQSTD == CSO_10) {
    asd_value = 5.2e-15 / tau + 3.6e-15 / sqrt(tau) + 4.0e-16
              + 0.9e-22 * pow(tau, 1.5);
    asd_value = 5.2e-12 / tau + 3.6e-12 / sqrt(tau) + 4.0e-16
              + 0.9e-22 * pow(tau, 1.5);
    /**** Nand et al., 2011 ****/

  } else if (FRQSTD == CSO_100) {
    asd_value = 2.1e-15 / sqrt(tau) + 3.0e-16
              + 0.9e-22 * pow(tau, 1.5);
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
