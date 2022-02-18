#include <stdio.h>
#include <aris.h>
#include <math.h>

/****
Neill Mapping Function (1996)
A. E. Niell, Journal of Geophysical Research, 101, B2, pp.3227-3246, 1996.
****/

#define T0   28

int  nmf20_coeffi(double lati,   int    DOY,
                  double *h,     double *h_hcr,    double *w)
{
  double nmfh_ave[3][5] = { 1.2769934e-3,      /** a: 15 deg **/
                            1.2683230e-3,      /** a: 30 deg **/
                            1.2465397e-3,      /** a: 45 deg **/
                            1.2196049e-3,      /** a: 60 deg **/
                            1.2045996e-3,      /** a: 75 deg **/

                            2.9153695e-3,      /** b: 15 deg **/
                            2.9152299e-3,      /** b: 30 deg **/
                            2.9288445e-3,      /** b: 45 deg **/
                            2.9022565e-3,      /** b: 60 deg **/
                            2.9024912e-3,      /** b: 75 deg **/

                            62.610505e-3,      /** c: 15 deg **/
                            62.837393e-3,      /** c: 30 deg **/
                            63.721774e-3,      /** c: 45 deg **/
                            63.824265e-3,      /** c: 60 deg **/
                            64.258455e-3       /** c: 75 deg **/
                          };
  double nmfh_amp[3][5] = {          0.0,      /** a: 15 deg **/
                            1.2709626e-5,      /** a: 30 deg **/
                            2.6523662e-5,      /** a: 45 deg **/
                            3.4000452e-5,      /** a: 60 deg **/
                            4.1202191e-5,      /** a: 75 deg **/

                                     0.0,      /** b: 15 deg **/
                            2.1414979e-5,      /** b: 30 deg **/
                            3.0160779e-5,      /** b: 45 deg **/
                            7.2562722e-5,      /** b: 60 deg **/
                            11.723375e-5,      /** b: 75 deg **/

                                     0.0,      /** c: 15 deg **/
                            9.0128400e-5,      /** c: 30 deg **/
                            4.3497037e-5,      /** c: 45 deg **/
                            84.795348e-5,      /** c: 60 deg **/
                            170.37206e-5       /** c: 75 deg **/
                          };
  double nmfh_hcr[3]    = {
                            2.53e-5,    5.49e-3,    1.14e-3
                          };
  double nmfw[3][5]     = { 5.8021897e-4,      /** a: 15 deg **/
                            5.6794847e-4,      /** a: 30 deg **/
                            5.8118019e-4,      /** a: 45 deg **/
                            5.9727542e-4,      /** a: 60 deg **/
                            6.1641693e-4,      /** a: 75 deg **/

                            1.4275268e-3,      /** b: 15 deg **/
                            1.5138625e-3,      /** b: 30 deg **/
                            1.4572752e-3,      /** b: 45 deg **/
                            1.5007428e-3,      /** b: 60 deg **/
                            1.7599082e-3,      /** b: 75 deg **/

                            4.3472961e-2,      /** c: 15 deg **/
                            4.6729510e-2,      /** c: 30 deg **/
                            4.3908931e-2,      /** c: 45 deg **/
                            4.4626982e-2,      /** c: 60 deg **/
                            5.4736038e-2       /** c: 75 deg **/
                          };
  double Lati[5];
  double abs_lati, cos_DOY, d_part;
  double H0, H1;
  int    i, I;

/*
----------------------------
*/

  for (i=0; i<5; i++) {
    Lati[i] = 15.0 * (double)(1 + i) / 180.0 * dpi;
  }
  abs_lati = fabs(lati);

  cos_DOY = cos(2.0 * dpi * (double)(DOY - T0) / 365.25);
  if (lati < 0.0) {
    cos_DOY *= -1.0;
  }

/*
----------------------------
*/

  if (abs_lati <= Lati[0]) {
    for (i=0; i<3; i++) {
      h[i] = nmfh_ave[i][0] + nmfh_amp[i][0] * cos_DOY;
      w[i] = nmfw[i][0];
    }
  } else if (abs_lati >= Lati[4]) {
    for (i=0; i<3; i++) {
      h[i] = nmfh_ave[i][4] + nmfh_amp[i][4] * cos_DOY;
      w[i] = nmfw[i][4];
    }
  } else {
    for (I=0; I<4; I++) {
      if (abs_lati > Lati[I] && abs_lati <= Lati[I+1]) {
        d_part = (abs_lati  - Lati[I]) / (Lati[I+1] - Lati[I]);
        for (i=0; i<3; i++) {
          H0 = nmfh_ave[i][I  ] + nmfh_amp[i][I  ] * cos_DOY;
          H1 = nmfh_ave[i][I+1] + nmfh_amp[i][I+1] * cos_DOY;
          h[i] =         H0 + (          H1 -         H0) * d_part;
          w[i] = nmfw[i][I] + (nmfw[i][I+1] - nmfw[i][I]) * d_part;
        }
        break;
      }
    }
  }

  for (i=0; i<3; i++) {
    h_hcr[i] = nmfh_hcr[i];
  }

/*
----------------------------
*/

  return 1;
}
