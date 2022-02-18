#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrotools.h"

#define L_NUT_NUM   678
#define P_NUT_NUM   687
#define L_NUT_T       5
#define P_NUT_T      14


void   nutation_matrix_set(char *, char *, int  , double *);


void   nutation_calc(int  *TimUTC,  double dT,
             double *epsiron, double *Delta_eps, double *Delta_psi)
/***************************************************************/
/* The IAU 2000 Theory of Precession/Nutation: Nutation        */
/*  IERS Techinical Notes No. 32                               */
/*        (IERS Convensions 2003)                              */
/*                 ed. by McCarthy and Petit (2003)            */
/*                                                             */
/***************************************************************/

{
  double L_NUTATION_A[6];
/*
   L Lm  F  D Om       Period 
                       (days) 
*/

  double L_NUTATION_B[8];
/*
               In Phase                          Out of phase 
       Psi     dPsi/dt   Eps    dEps/dt   Psi   dPsi/dt    Eps  Deps/dt 
       (mas)   (mas/c)    (mas)   (mas/c)  (mas)  (mas/c)  (mas)  (mas/c) 
*/

  double P_NUTATION_A[15];
/*
Planetary nutation terms
  L  L'   F   D  Om  Lm  Lv  Le  LM  Lj  Ls  Lu  Ln  Pa     Period
                                                             (days)
*/

  double P_NUTATION_B[5];
/*
      Longitude          Obliquity        Amplitude 
     In       Out       In       Out      (mas)
    (mas)    (mas)     (mas)    (mas) 
*/

  int    i, j;
  double F[P_NUT_T], ARGUMENT;
  double Delta_TT, Delta_UTC_TAI;
  double t1, t2, t3, t4;
  double sin_a, cos_a;
  double dpi = 3.141592653589793238462643;
  double DPI;

  DPI = dpi / 180.0;

/*
---- Epsiron_A ----
*/

  *epsiron = mean_obliquity_of_the_ecliptic(TimUTC, dT);

/*
---- LUNISOLAR NUTATION ----
*/

  Delta_UTC_TAI = (double)UTC_minus_TAI(TimUTC) + 32.0;
  /* -32 (sec) : UTC - TAI at J2000  */

  Delta_TT   = MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                   TimUTC[3], TimUTC[4], TimUTC[5],  dT)
             - MJD(     2000,         1,         1,
                          12,         0,         0, 0.0)
             + (double)Delta_UTC_TAI / 86400.0;
  /* 32.184 sec to convert from TAI to TT is canceled. */

  t1 = Delta_TT / 36525.0;
  t2 = t1 * t1;
  t3 = t2 * t1;
  t4 = t3 * t1;

  *Delta_psi = 0.0;
  *Delta_eps = 0.0;

/* l     : Mean Anomaly of the Moon                         */
  F[0] = 134.96340251    + (1717915923.2178 * t1
                         +          31.8792 * t2
                         +         0.051635 * t3
                         -       0.00024470 * t4) / 3600.0;

/* l'    : Mean Anomaly of the Sun                          */
  F[1] = 357.52910918    +  (129596581.0481 * t1
                         -           0.5532 * t2
                         +         0.000136 * t3
                         -       0.00001149 * t4) / 3600.0;

/* F     : L - Omega                                        */
  F[2] =  93.27209062    + (1739527262.8478 * t1
                         -          12.7512 * t2
                         -         0.001037 * t3
                         +       0.00000417 * t4) / 3600.0;

/* D     : Mean Elongation of the Moon from the Sun         */
  F[3] = 297.85019547    + (1602961601.2090 * t1
                         -           6.3706 * t2
                         +         0.006593 * t3
                         -       0.00003169 * t4) / 3600.0;

/* Omega : Mean Longitude of the Ascending Node of the Moon */
  F[4] = 125.04455501    - (   6962890.2665 * t1
                         +           7.4722 * t2
                         +         0.007702 * t3
                         -       0.00005939 * t4) / 3600.0;

  for (i=0; i<L_NUT_T; i++) {
    F[i] *= DPI;
  }

  for (i=0; i<L_NUT_NUM; i++) {
    nutation_matrix_set("L", "A", i, L_NUTATION_A);
    nutation_matrix_set("L", "B", i, L_NUTATION_B);
    ARGUMENT = 0.0;
    for (j=0; j<L_NUT_T; j++) {
      ARGUMENT += L_NUTATION_A[j] * F[j];
    }
    cos_a = cos(ARGUMENT);
    sin_a = sin(ARGUMENT);
    *Delta_psi += ((L_NUTATION_B[0] + L_NUTATION_B[1]*t1) * sin_a
                +  (L_NUTATION_B[4] + L_NUTATION_B[5]*t1) * cos_a);
    *Delta_eps += ((L_NUTATION_B[2] + L_NUTATION_B[3]*t1) * cos_a
                +  (L_NUTATION_B[6] + L_NUTATION_B[7]*t1) * sin_a);
  }

/*
---- PLANETARY NUTATION ----
*/

/* lMe [rad]                                                */
  F[ 5] =  4.402608842 + 2608.7903141574 * t1;
/* lVe [rad]                                                */
  F[ 6] =  3.176146697 + 1021.3285546211 * t1;
/* lE  [rad]                                                */
  F[ 7] =  1.753470314 +  628.3075849991 * t1;
/* lMa [rad]                                                */
  F[ 8] =  6.203480913 +  334.0612426700 * t1;
/* lJu [rad]                                                */
  F[ 9] =  0.599546497 +   52.9690962641 * t1;
/* lSa [rad]                                                */
  F[10] =  0.874016757 +   21.3299104960 * t1;
/* lUr [rad]                                                */
  F[11] =  5.481293871 +    7.4781598567 * t1;
/* lNe [rad]                                                */
  F[12] =  5.321159000 +    3.8127774000 * t1;
/* pa  [rad]                                                */
  F[13] =                   0.02438175   * t1
                       +   0.00000538691 * t2;

  for (i=0; i<P_NUT_NUM; i++) {
    nutation_matrix_set("P", "A", i, P_NUTATION_A);
    nutation_matrix_set("P", "B", i, P_NUTATION_B);
    ARGUMENT = 0.0;
    for (j=0; j<P_NUT_T; j++) {
      ARGUMENT += P_NUTATION_A[j] * F[j];
    }
    cos_a = cos(ARGUMENT);
    sin_a = sin(ARGUMENT);
    *Delta_psi += (P_NUTATION_B[0] * sin_a
               +   P_NUTATION_B[1] * cos_a);
    *Delta_eps += (P_NUTATION_B[2] * sin_a
               +   P_NUTATION_B[3] * cos_a);
  }

/*
----------------------------------------------
*/

  *Delta_psi *= (DPI / 3600.0 / 1.0e3);
  *Delta_eps *= (DPI / 3600.0 / 1.0e3);

/*
----------------------------------------------
*/

  return;
}



void   nutation_matrix_set(char *S1, char *S2, int  n, double *nut_data)
{
  if (*S1 == 'L' && *S2 == 'A') {
    if (n ==   0) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =     -6798.3830;
    }
    else if (n ==   1) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       182.6210;
    }
    else if (n ==   2) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.6610;
    }
    else if (n ==   3) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =     -3399.1920;
    }
    else if (n ==   4) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       365.2600;
    }
    else if (n ==   5) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       121.7490;
    }
    else if (n ==   6) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        27.5550;
    }
    else if (n ==   7) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.6330;
    }
    else if (n ==   8) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.1330;
    }
    else if (n ==   9) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       365.2250;
    }
    else if (n ==  10) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       177.8440;
    }
    else if (n ==  11) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        27.0930;
    }
    else if (n ==  12) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        31.8120;
    }
    else if (n ==  13) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        27.6670;
    }
    else if (n ==  14) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -27.4430;
    }
    else if (n ==  15) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.5570;
    }
    else if (n ==  16) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.1210;
    }
    else if (n ==  17) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      1305.4790;
    }
    else if (n ==  18) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.7650;
    }
    else if (n ==  19) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.0960;
    }
    else if (n ==  20) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =   3823589.3370;
    }
    else if (n ==  21) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -205.8920;
    }
    else if (n ==  22) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.8590;
    }
    else if (n ==  23) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        23.9420;
    }
    else if (n ==  24) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        26.9850;
    }
    else if (n ==  25) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        13.7770;
    }
    else if (n ==  26) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        13.6060;
    }
    else if (n ==  27) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       385.9980;
    }
    else if (n ==  28) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        31.9610;
    }
    else if (n ==  29) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        91.3130;
    }
    else if (n ==  30) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -173.3100;
    }
    else if (n ==  31) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -31.6640;
    }
    else if (n ==  32) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -346.6360;
    }
    else if (n ==  33) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.5430;
    }
    else if (n ==  34) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       182.6300;
    }
    else if (n ==  35) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.6430;
    }
    else if (n ==  36) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      1095.1750;
    }
    else if (n ==  37) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.1680;
    }
    else if (n ==  38) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.0880;
    }
    else if (n ==  39) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.1920;
    }
    else if (n ==  40) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.7970;
    }
    else if (n ==  41) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        23.8580;
    }
    else if (n ==  42) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        12.8110;
    }
    else if (n ==  43) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -199.8400;
    }
    else if (n ==  44) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.8520;
    }
    else if (n ==  45) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       346.6040;
    }
    else if (n ==  46) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -14.7330;
    }
    else if (n ==  47) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        34.8470;
    }
    else if (n ==  48) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       212.3230;
    }
    else if (n ==  49) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.6140;
    }
    else if (n ==  50) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       119.6070;
    }
    else if (n ==  51) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        29.8030;
    }
    else if (n ==  52) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      1615.7480;
    }
    else if (n ==  53) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.4920;
    }
    else if (n ==  54) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        15.3870;
    }
    else if (n ==  55) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.3670;
    }
    else if (n ==  56) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        29.5310;
    }
    else if (n ==  57) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.8140;
    }
    else if (n ==  58) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        26.8780;
    }
    else if (n ==  59) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.2360;
    }
    else if (n ==  60) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -13.7490;
    }
    else if (n ==  61) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.9100;
    }
    else if (n ==  62) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.8050;
    }
    else if (n ==  63) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      3232.8620;
    }
    else if (n ==  64) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        25.6220;
    }
    else if (n ==  65) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.1080;
    }
    else if (n ==  66) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -32.6060;
    }
    else if (n ==  67) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        27.7800;
    }
    else if (n ==  68) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -411.7840;
    }
    else if (n ==  69) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.3400;
    }
    else if (n ==  70) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.8020;
    }
    else if (n ==  71) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      6164.1010;
    }
    else if (n ==  72) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      6786.3170;
    }
    else if (n ==  73) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.6380;
    }
    else if (n ==  74) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.6320;
    }
    else if (n ==  75) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -27.3330;
    }
    else if (n ==  76) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        22.4690;
    }
    else if (n ==  77) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.3490;
    }
    else if (n ==  78) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.0570;
    }
    else if (n ==  79) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        12.7870;
    }
    else if (n ==  80) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.6840;
    }
    else if (n ==  81) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.6270;
    }
    else if (n ==  82) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.1850;
    }
    else if (n ==  83) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.7450;
    }
    else if (n ==  84) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        12.6630;
    }
    else if (n ==  85) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.1430;
    }
    else if (n ==  86) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -169.0020;
    }
    else if (n ==  87) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =       187.6620;
    }
    else if (n ==  88) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        10.0850;
    }
    else if (n ==  89) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -943.2270;
    }
    else if (n ==  90) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        15.9060;
    }
    else if (n ==  91) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        35.0260;
    }
    else if (n ==  92) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -388.2670;
    }
    else if (n ==  93) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       409.2340;
    }
    else if (n ==  94) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -13.5790;
    }
    else if (n ==  95) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.1620;
    }
    else if (n ==  96) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        25.4200;
    }
    else if (n ==  97) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.7930;
    }
    else if (n ==  98) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -131.6710;
    }
    else if (n ==  99) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -34.6690;
    }
    else if (n == 100) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        29.2630;
    }
    else if (n == 101) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =     66068.1880;
    }
    else if (n == 102) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        29.9340;
    }
    else if (n == 103) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.7310;
    }
    else if (n == 104) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.3130;
    }
    else if (n == 105) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.4880;
    }
    else if (n == 106) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -329.7910;
    }
    else if (n == 107) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -9.6000;
    }
    else if (n == 108) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.9610;
    }
    else if (n == 109) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.7990;
    }
    else if (n == 110) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -329.8190;
    }
    else if (n == 111) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -38.7420;
    }
    else if (n == 112) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -23.7750;
    }
    else if (n == 113) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.2290;
    }
    else if (n == 114) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.9910;
    }
    else if (n == 115) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.8300;
    }
    else if (n == 116) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.3540;
    }
    else if (n == 117) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        25.2220;
    }
    else if (n == 118) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.1920;
    }
    else if (n == 119) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -117.5390;
    }
    else if (n == 120) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  3.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        73.0510;
    }
    else if (n == 121) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        29.6590;
    }
    else if (n == 122) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.5300;
    }
    else if (n == 123) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.7330;
    }
    else if (n == 124) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        25.7190;
    }
    else if (n == 125) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         8.8980;
    }
    else if (n == 126) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.1270;
    }
    else if (n == 127) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        32.7640;
    }
    else if (n == 128) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        32.1120;
    }
    else if (n == 129) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        27.3220;
    }
    else if (n == 130) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -15.3530;
    }
    else if (n == 131) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -32.4510;
    }
    else if (n == 132) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -29.4030;
    }
    else if (n == 133) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -29.6730;
    }
    else if (n == 134) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.2220;
    }
    else if (n == 135) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.8740;
    }
    else if (n == 136) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.3830;
    }
    else if (n == 137) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.9760;
    }
    else if (n == 138) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.3270;
    }
    else if (n == 139) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -31.5170;
    }
    else if (n == 140) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.7970;
    }
    else if (n == 141) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -26.7720;
    }
    else if (n == 142) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        22.3950;
    }
    else if (n == 143) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.0810;
    }
    else if (n == 144) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       313.0420;
    }
    else if (n == 145) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.6000;
    }
    else if (n == 146) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.5790;
    }
    else if (n == 147) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.3170;
    }
    else if (n == 148) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        12.3770;
    }
    else if (n == 149) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.1070;
    }
    else if (n == 150) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.6760;
    }
    else if (n == 151) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -25.5250;
    }
    else if (n == 152) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.2210;
    }
    else if (n == 153) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.3410;
    }
    else if (n == 154) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         6.8460;
    }
    else if (n == 155) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.2540;
    }
    else if (n == 156) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        10.1000;
    }
    else if (n == 157) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.0450;
    }
    else if (n == 158) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.6800;
    }
    else if (n == 159) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -35.2270;
    }
    else if (n == 160) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        38.5220;
    }
    else if (n == 161) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        13.2760;
    }
    else if (n == 162) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.8170;
    }
    else if (n == 163) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =     -2266.1280;
    }
    else if (n == 164) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  3.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       121.7530;
    }
    else if (n == 165) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -16.1020;
    }
    else if (n == 166) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        15.4220;
    }
    else if (n == 167) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.3910;
    }
    else if (n == 168) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.8950;
    }
    else if (n == 169) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.0830;
    }
    else if (n == 170) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      1616.4310;
    }
    else if (n == 171) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -507.1570;
    }
    else if (n == 172) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -194.1330;
    }
    else if (n == 173) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        29.2630;
    }
    else if (n == 174) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        12.6390;
    }
    else if (n == 175) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         8.7340;
    }
    else if (n == 176) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -129.1690;
    }
    else if (n == 177) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       438.3350;
    }
    else if (n == 178) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        16.0640;
    }
    else if (n == 179) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        15.9430;
    }
    else if (n == 180) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -9.1720;
    }
    else if (n == 181) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.5570;
    }
    else if (n == 182) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.7890;
    }
    else if (n == 183) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.0030;
    }
    else if (n == 184) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -35.8030;
    }
    else if (n == 185) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -15.8690;
    }
    else if (n == 186) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -14.7010;
    }
    else if (n == 187) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        95.4210;
    }
    else if (n == 188) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        29.3900;
    }
    else if (n == 189) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        25.3250;
    }
    else if (n == 190) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.3830;
    }
    else if (n == 191) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.1350;
    }
    else if (n == 192) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.6380;
    }
    else if (n == 193) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       219.1670;
    }
    else if (n == 194) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       187.6710;
    }
    else if (n == 195) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -10.0700;
    }
    else if (n == 196) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        90.1030;
    }
    else if (n == 197) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        37.6250;
    }
    else if (n == 198) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        25.1290;
    }
    else if (n == 199) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        10.3710;
    }
    else if (n == 200) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        10.0850;
    }
    else if (n == 201) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.5000;
    }
    else if (n == 202) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.7260;
    }
    else if (n == 203) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -471.8910;
    }
    else if (n == 204) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       285.4060;
    }
    else if (n == 205) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       134.2720;
    }
    else if (n == 206) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -28.1490;
    }
    else if (n == 207) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -14.9340;
    }
    else if (n == 208) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        29.5170;
    }
    else if (n == 209) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        27.3210;
    }
    else if (n == 210) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        23.4300;
    }
    else if (n == 211) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.7650;
    }
    else if (n == 212) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.1300;
    }
    else if (n == 213) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.8330;
    }
    else if (n == 214) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =        13.6880;
    }
    else if (n == 215) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.4920;
    }
    else if (n == 216) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -9.0960;
    }
    else if (n == 217) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.3000;
    }
    else if (n == 218) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.1970;
    }
    else if (n == 219) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.2200;
    }
    else if (n == 220) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.9830;
    }
    else if (n == 221) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.9530;
    }
    else if (n == 222) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.8560;
    }
    else if (n == 223) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.7440;
    }
    else if (n == 224) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =     -3396.1730;
    }
    else if (n == 225) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      2189.7230;
    }
    else if (n == 226) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -177.8520;
    }
    else if (n == 227) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -38.5220;
    }
    else if (n == 228) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -29.1380;
    }
    else if (n == 229) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        29.1380;
    }
    else if (n == 230) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        25.6220;
    }
    else if (n == 231) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        16.6300;
    }
    else if (n == 232) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        15.3140;
    }
    else if (n == 233) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        15.2420;
    }
    else if (n == 234) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        12.7100;
    }
    else if (n == 235) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.3670;
    }
    else if (n == 236) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.9350;
    }
    else if (n == 237) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.9680;
    }
    else if (n == 238) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         6.8890;
    }
    else if (n == 239) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.7260;
    }
    else if (n == 240) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.5760;
    }
    else if (n == 241) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      2120.6520;
    }
    else if (n == 242) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       124.2030;
    }
    else if (n == 243) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -43.3390;
    }
    else if (n == 244) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -29.9340;
    }
    else if (n == 245) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -27.2120;
    }
    else if (n == 246) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -14.1620;
    }
    else if (n == 247) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -13.7220;
    }
    else if (n == 248) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -12.7630;
    }
    else if (n == 249) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        38.9640;
    }
    else if (n == 250) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        32.4510;
    }
    else if (n == 251) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        27.4320;
    }
    else if (n == 252) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        27.0930;
    }
    else if (n == 253) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        25.5250;
    }
    else if (n == 254) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        21.1670;
    }
    else if (n == 255) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.2760;
    }
    else if (n == 256) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.1960;
    }
    else if (n == 257) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        12.3540;
    }
    else if (n == 258) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -7.1200;
    }
    else if (n == 259) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.6140;
    }
    else if (n == 260) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.0950;
    }
    else if (n == 261) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         8.6650;
    }
    else if (n == 262) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.0500;
    }
    else if (n == 263) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.8450;
    }
    else if (n == 264) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.8230;
    }
    else if (n == 265) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.6330;
    }
    else if (n == 266) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.5660;
    }
    else if (n == 267) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.4110;
    }
    else if (n == 268) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.5760;
    }
    else if (n == 269) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -548.0410;
    }
    else if (n == 270) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -314.5330;
    }
    else if (n == 271) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -299.2620;
    }
    else if (n == 272) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -115.5420;
    }
    else if (n == 273) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       507.0900;
    }
    else if (n == 274) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       328.1530;
    }
    else if (n == 275) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       126.5140;
    }
    else if (n == 276) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =       123.9690;
    }
    else if (n == 277) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -35.4100;
    }
    else if (n == 278) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -23.6920;
    }
    else if (n == 279) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -16.0640;
    }
    else if (n == 280) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -14.2870;
    }
    else if (n == 281) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -14.2240;
    }
    else if (n == 282) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        96.7800;
    }
    else if (n == 283) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  4.0;
      nut_data[ 5] =        91.3110;
    }
    else if (n == 284) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        88.9220;
    }
    else if (n == 285) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        38.7420;
    }
    else if (n == 286) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        34.4750;
    }
    else if (n == 287) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        32.9220;
    }
    else if (n == 288) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        31.1980;
    }
    else if (n == 289) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        31.0550;
    }
    else if (n == 290) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.5690;
    }
    else if (n == 291) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.3480;
    }
    else if (n == 292) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.0680;
    }
    else if (n == 293) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        13.7190;
    }
    else if (n == 294) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        12.2380;
    }
    else if (n == 295) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -9.8590;
    }
    else if (n == 296) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -7.3750;
    }
    else if (n == 297) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.8880;
    }
    else if (n == 298) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.3800;
    }
    else if (n == 299) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.0820;
    }
    else if (n == 300) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.5410;
    }
    else if (n == 301) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.5350;
    }
    else if (n == 302) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.2690;
    }
    else if (n == 303) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.8100;
    }
    else if (n == 304) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.5980;
    }
    else if (n == 305) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.8900;
    }
    else if (n == 306) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.8280;
    }
    else if (n == 307) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.8230;
    }
    else if (n == 308) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.7210;
    }
    else if (n == 309) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.7110;
    }
    else if (n == 310) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.6620;
    }
    else if (n == 311) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.6140;
    }
    else if (n == 312) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.5520;
    }
    else if (n == 313) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.9070;
    }
    else if (n == 314) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.6240;
    }
    else if (n == 315) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.1650;
    }
    else if (n == 316) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.0800;
    }
    else if (n == 317) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         3.5560;
    }
    else if (n == 318) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =     -1656.2530;
    }
    else if (n == 319) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      1077.3170;
    }
    else if (n == 320) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -828.3060;
    }
    else if (n == 321) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -367.2900;
    }
    else if (n == 322) {
      nut_data[ 0] = -4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -244.4450;
    }
    else if (n == 323) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -193.5650;
    }
    else if (n == 324) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -164.9020;
    }
    else if (n == 325) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       468.5450;
    }
    else if (n == 326) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =       385.9590;
    }
    else if (n == 327) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       273.9070;
    }
    else if (n == 328) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  4.0;
      nut_data[ 5] =       192.9890;
    }
    else if (n == 329) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -96.7820;
    }
    else if (n == 330) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -94.1010;
    }
    else if (n == 331) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -88.9240;
    }
    else if (n == 332) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -38.3050;
    }
    else if (n == 333) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -37.4180;
    }
    else if (n == 334) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -35.6150;
    }
    else if (n == 335) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -35.0270;
    }
    else if (n == 336) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -34.4930;
    }
    else if (n == 337) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -28.2660;
    }
    else if (n == 338) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -28.0330;
    }
    else if (n == 339) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -25.8260;
    }
    else if (n == 340) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -25.2310;
    }
    else if (n == 341) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -24.2160;
    }
    else if (n == 342) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -22.3220;
    }
    else if (n == 343) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -16.8450;
    }
    else if (n == 344) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -16.5900;
    }
    else if (n == 345) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -16.0260;
    }
    else if (n == 346) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -15.2800;
    }
    else if (n == 347) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -14.9010;
    }
    else if (n == 348) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -13.2510;
    }
    else if (n == 349) {
      nut_data[ 0] = -4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -12.9130;
    }
    else if (n == 350) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -10.3550;
    }
    else if (n == 351) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -10.1630;
    }
    else if (n == 352) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        87.7740;
    }
    else if (n == 353) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  3.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        73.0510;
    }
    else if (n == 354) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        41.9460;
    }
    else if (n == 355) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        37.8350;
    }
    else if (n == 356) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        35.9920;
    }
    else if (n == 357) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        35.2080;
    }
    else if (n == 358) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        32.6060;
    }
    else if (n == 359) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        30.0660;
    }
    else if (n == 360) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        29.7890;
    }
    else if (n == 361) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        29.0130;
    }
    else if (n == 362) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        27.5540;
    }
    else if (n == 363) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        27.2120;
    }
    else if (n == 364) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =        27.2010;
    }
    else if (n == 365) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        25.8160;
    }
    else if (n == 366) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        25.0360;
    }
    else if (n == 367) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        23.9420;
    }
    else if (n == 368) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        23.5930;
    }
    else if (n == 369) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        23.3490;
    }
    else if (n == 370) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        21.4480;
    }
    else if (n == 371) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        21.1010;
    }
    else if (n == 372) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        21.0360;
    }
    else if (n == 373) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        16.6710;
    }
    else if (n == 374) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        15.4570;
    }
    else if (n == 375) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        15.2420;
    }
    else if (n == 376) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        15.2080;
    }
    else if (n == 377) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.9670;
    }
    else if (n == 378) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.2840;
    }
    else if (n == 379) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.2510;
    }
    else if (n == 380) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.7190;
    }
    else if (n == 381) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.4660;
    }
    else if (n == 382) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.3020;
    }
    else if (n == 383) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        13.1170;
    }
    else if (n == 384) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.1170;
    }
    else if (n == 385) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  4.0;
      nut_data[ 5] =        12.7100;
    }
    else if (n == 386) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        12.6860;
    }
    else if (n == 387) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        10.6040;
    }
    else if (n == 388) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        10.3870;
    }
    else if (n == 389) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        10.1480;
    }
    else if (n == 390) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        10.1150;
    }
    else if (n == 391) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        10.0700;
    }
    else if (n == 392) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        -9.5870;
    }
    else if (n == 393) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -9.5170;
    }
    else if (n == 394) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -7.0740;
    }
    else if (n == 395) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -6.8390;
    }
    else if (n == 396) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.8580;
    }
    else if (n == 397) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.8440;
    }
    else if (n == 398) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.8140;
    }
    else if (n == 399) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.7850;
    }
    else if (n == 400) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.7840;
    }
    else if (n == 401) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.6410;
    }
    else if (n == 402) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.5850;
    }
    else if (n == 403) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.4220;
    }
    else if (n == 404) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.3940;
    }
    else if (n == 405) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.3410;
    }
    else if (n == 406) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.3140;
    }
    else if (n == 407) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =         9.1450;
    }
    else if (n == 408) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         8.9600;
    }
    else if (n == 409) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.9600;
    }
    else if (n == 410) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         8.9230;
    }
    else if (n == 411) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         8.8870;
    }
    else if (n == 412) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.8620;
    }
    else if (n == 413) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.6980;
    }
    else if (n == 414) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.6570;
    }
    else if (n == 415) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.5430;
    }
    else if (n == 416) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.4920;
    }
    else if (n == 417) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.3750;
    }
    else if (n == 418) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.2210;
    }
    else if (n == 419) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.2130;
    }
    else if (n == 420) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.2040;
    }
    else if (n == 421) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.1420;
    }
    else if (n == 422) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.1270;
    }
    else if (n == 423) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.0800;
    }
    else if (n == 424) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.7470;
    }
    else if (n == 425) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.6320;
    }
    else if (n == 426) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.9920;
    }
    else if (n == 427) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.9920;
    }
    else if (n == 428) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.9700;
    }
    else if (n == 429) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.7920;
    }
    else if (n == 430) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.6670;
    }
    else if (n == 431) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.5710;
    }
    else if (n == 432) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.4830;
    }
    else if (n == 433) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.4650;
    }
    else if (n == 434) {
      nut_data[ 0] =  5.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.3490;
    }
    else if (n == 435) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.8530;
    }
    else if (n == 436) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.7410;
    }
    else if (n == 437) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.7310;
    }
    else if (n == 438) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.1290;
    }
    else if (n == 439) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.0480;
    }
    else if (n == 440) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.0010;
    }
    else if (n == 441) {
      nut_data[ 0] =  5.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         3.9270;
    }
    else if (n == 442) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         3.6180;
    }
    else if (n == 443) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         3.4950;
    }
    else if (n == 444) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =   7647178.6740;
    }
    else if (n == 445) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =     65502.2780;
    }
    else if (n == 446) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =     -6810.4930;
    }
    else if (n == 447) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =     -1656.6120;
    }
    else if (n == 448) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =     -1305.9250;
    }
    else if (n == 449) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      6159.1360;
    }
    else if (n == 450) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      3231.4960;
    }
    else if (n == 451) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      3230.1310;
    }
    else if (n == 452) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      3082.0500;
    }
    else if (n == 453) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      2190.3500;
    }
    else if (n == 454) {
      nut_data[ 0] = -4.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -666.5670;
    }
    else if (n == 455) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -601.5210;
    }
    else if (n == 456) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -471.9500;
    }
    else if (n == 457) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -286.6440;
    }
    else if (n == 458) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =      -263.2990;
    }
    else if (n == 459) {
      nut_data[ 0] = -4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -253.5620;
    }
    else if (n == 460) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -219.8970;
    }
    else if (n == 461) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -1.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -188.2010;
    }
    else if (n == 462) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -173.3180;
    }
    else if (n == 463) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -173.3180;
    }
    else if (n == 464) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -137.2610;
    }
    else if (n == 465) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      -126.7610;
    }
    else if (n == 466) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -121.9740;
    }
    else if (n == 467) {
      nut_data[ 0] = -4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =      -102.9460;
    }
    else if (n == 468) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       596.0940;
    }
    else if (n == 469) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       547.9620;
    }
    else if (n == 470) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       346.6200;
    }
    else if (n == 471) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       328.1810;
    }
    else if (n == 472) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       297.9130;
    }
    else if (n == 473) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       199.2370;
    }
    else if (n == 474) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       193.5600;
    }
    else if (n == 475) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       192.9990;
    }
    else if (n == 476) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       188.1970;
    }
    else if (n == 477) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       168.5710;
    }
    else if (n == 478) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       156.5210;
    }
    else if (n == 479) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       152.9990;
    }
    else if (n == 480) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -95.4240;
    }
    else if (n == 481) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -94.1010;
    }
    else if (n == 482) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -49.1730;
    }
    else if (n == 483) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -43.0640;
    }
    else if (n == 484) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -39.6930;
    }
    else if (n == 485) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -34.3010;
    }
    else if (n == 486) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -32.2960;
    }
    else if (n == 487) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -14.7650;
    }
    else if (n == 488) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -29.8030;
    }
    else if (n == 489) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -29.5440;
    }
    else if (n == 490) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -29.2760;
    }
    else if (n == 491) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -27.7910;
    }
    else if (n == 492) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -26.9850;
    }
    else if (n == 493) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -26.6670;
    }
    else if (n == 494) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -26.0340;
    }
    else if (n == 495) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -25.4300;
    }
    else if (n == 496) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -25.4300;
    }
    else if (n == 497) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -24.3020;
    }
    else if (n == 498) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -22.7860;
    }
    else if (n == 499) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -6.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -17.4680;
    }
    else if (n == 500) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -16.8030;
    }
    else if (n == 501) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -15.8320;
    }
    else if (n == 502) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -15.5710;
    }
    else if (n == 503) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -15.4220;
    }
    else if (n == 504) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -15.3180;
    }
    else if (n == 505) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        -9.6140;
    }
    else if (n == 506) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -14.5380;
    }
    else if (n == 507) {
      nut_data[ 0] = -4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -13.9240;
    }
    else if (n == 508) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -13.6910;
    }
    else if (n == 509) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -13.5520;
    }
    else if (n == 510) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -13.3310;
    }
    else if (n == 511) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -13.1710;
    }
    else if (n == 512) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -12.7390;
    }
    else if (n == 513) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -4.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -12.6160;
    }
    else if (n == 514) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =       -12.3320;
    }
    else if (n == 515) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -6.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =       -10.6910;
    }
    else if (n == 516) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -10.1480;
    }
    else if (n == 517) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =       -10.0550;
    }
    else if (n == 518) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        76.5080;
    }
    else if (n == 519) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        75.6570;
    }
    else if (n == 520) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  4.0;
      nut_data[ 5] =        73.0490;
    }
    else if (n == 521) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        71.5120;
    }
    else if (n == 522) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        43.6170;
    }
    else if (n == 523) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -3.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        43.0640;
    }
    else if (n == 524) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        39.1890;
    }
    else if (n == 525) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        38.0680;
    }
    else if (n == 526) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        35.5960;
    }
    else if (n == 527) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        34.6510;
    }
    else if (n == 528) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        32.1280;
    }
    else if (n == 529) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        30.9140;
    }
    else if (n == 530) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        30.0670;
    }
    else if (n == 531) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        29.5170;
    }
    else if (n == 532) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        27.4430;
    }
    else if (n == 533) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        27.4320;
    }
    else if (n == 534) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        27.3220;
    }
    else if (n == 535) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        27.2120;
    }
    else if (n == 536) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        27.2010;
    }
    else if (n == 537) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.6610;
    }
    else if (n == 538) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =      1615.7480;
    }
    else if (n == 539) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        25.9250;
    }
    else if (n == 540) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        24.3890;
    }
    else if (n == 541) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =        24.0270;
    }
    else if (n == 542) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        24.0270;
    }
    else if (n == 543) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        22.6270;
    }
    else if (n == 544) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        22.0180;
    }
    else if (n == 545) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        17.4230;
    }
    else if (n == 546) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -3.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        16.8030;
    }
    else if (n == 547) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        16.1400;
    }
    else if (n == 548) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        15.9850;
    }
    else if (n == 549) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        15.9810;
    }
    else if (n == 550) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        15.3490;
    }
    else if (n == 551) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.9010;
    }
    else if (n == 552) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.8330;
    }
    else if (n == 553) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        31.8120;
    }
    else if (n == 554) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.7330;
    }
    else if (n == 555) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.6980;
    }
    else if (n == 556) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.6980;
    }
    else if (n == 557) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.6980;
    }
    else if (n == 558) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        14.1330;
    }
    else if (n == 559) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        14.0390;
    }
    else if (n == 560) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        13.7770;
    }
    else if (n == 561) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.7470;
    }
    else if (n == 562) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.7190;
    }
    else if (n == 563) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  1.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        13.6910;
    }
    else if (n == 564) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        13.6610;
    }
    else if (n == 565) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.2510;
    }
    else if (n == 566) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.0920;
    }
    else if (n == 567) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =        12.6860;
    }
    else if (n == 568) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        12.2160;
    }
    else if (n == 569) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        12.0600;
    }
    else if (n == 570) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        11.9710;
    }
    else if (n == 571) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        11.9290;
    }
    else if (n == 572) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        10.6740;
    }
    else if (n == 573) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -3.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        10.3710;
    }
    else if (n == 574) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        10.0220;
    }
    else if (n == 575) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -9.6850;
    }
    else if (n == 576) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -9.3540;
    }
    else if (n == 577) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        -9.1600;
    }
    else if (n == 578) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] = -2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =        -8.7230;
    }
    else if (n == 579) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -7.5270;
    }
    else if (n == 580) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -7.2610;
    }
    else if (n == 581) {
      nut_data[ 0] = -4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -6.8820;
    }
    else if (n == 582) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -5.8180;
    }
    else if (n == 583) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        -5.6580;
    }
    else if (n == 584) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.8720;
    }
    else if (n == 585) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.8280;
    }
    else if (n == 586) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.6000;
    }
    else if (n == 587) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.5850;
    }
    else if (n == 588) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =         9.5700;
    }
    else if (n == 589) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.0960;
    }
    else if (n == 590) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =        14.6320;
    }
    else if (n == 591) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.2880;
    }
    else if (n == 592) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.2100;
    }
    else if (n == 593) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.1590;
    }
    else if (n == 594) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.1590;
    }
    else if (n == 595) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =        13.6330;
    }
    else if (n == 596) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  3.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =         9.1070;
    }
    else if (n == 597) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  3.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.0950;
    }
    else if (n == 598) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         9.0700;
    }
    else if (n == 599) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         9.0330;
    }
    else if (n == 600) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         8.6870;
    }
    else if (n == 601) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         8.5300;
    }
    else if (n == 602) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         8.4740;
    }
    else if (n == 603) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.8210;
    }
    else if (n == 604) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.6940;
    }
    else if (n == 605) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.6660;
    }
    else if (n == 606) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.6570;
    }
    else if (n == 607) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -3.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.5350;
    }
    else if (n == 608) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.3990;
    }
    else if (n == 609) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.3660;
    }
    else if (n == 610) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.3330;
    }
    else if (n == 611) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.2770;
    }
    else if (n == 612) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.2530;
    }
    else if (n == 613) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.2440;
    }
    else if (n == 614) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         7.2360;
    }
    else if (n == 615) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         7.1110;
    }
    else if (n == 616) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =         7.1030;
    }
    else if (n == 617) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.6430;
    }
    else if (n == 618) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.5570;
    }
    else if (n == 619) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         7.0430;
    }
    else if (n == 620) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.9980;
    }
    else if (n == 621) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         6.9910;
    }
    else if (n == 622) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         6.9760;
    }
    else if (n == 623) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         6.9610;
    }
    else if (n == 624) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         6.9460;
    }
    else if (n == 625) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =         6.8660;
    }
    else if (n == 626) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.4920;
    }
    else if (n == 627) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         9.1330;
    }
    else if (n == 628) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  3.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  3.0;
      nut_data[ 5] =         6.8450;
    }
    else if (n == 629) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.8380;
    }
    else if (n == 630) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.8300;
    }
    else if (n == 631) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         6.7200;
    }
    else if (n == 632) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         6.5920;
    }
    else if (n == 633) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.5200;
    }
    else if (n == 634) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         6.0920;
    }
    else if (n == 635) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         6.0690;
    }
    else if (n == 636) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.9970;
    }
    else if (n == 637) {
      nut_data[ 0] = -3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.9650;
    }
    else if (n == 638) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.9220;
    }
    else if (n == 639) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.9170;
    }
    else if (n == 640) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  5.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.8850;
    }
    else if (n == 641) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.8180;
    }
    else if (n == 642) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.7520;
    }
    else if (n == 643) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.7220;
    }
    else if (n == 644) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.7160;
    }
    else if (n == 645) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.7060;
    }
    else if (n == 646) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.6330;
    }
    else if (n == 647) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.6090;
    }
    else if (n == 648) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.5620;
    }
    else if (n == 649) {
      nut_data[ 0] =  5.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         5.5110;
    }
    else if (n == 650) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.4830;
    }
    else if (n == 651) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.4600;
    }
    else if (n == 652) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         5.4070;
    }
    else if (n == 653) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] = -2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         5.3230;
    }
    else if (n == 654) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.9740;
    }
    else if (n == 655) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         4.9220;
    }
    else if (n == 656) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.9220;
    }
    else if (n == 657) {
      nut_data[ 0] = -2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.9030;
    }
    else if (n == 658) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.8100;
    }
    else if (n == 659) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         4.8070;
    }
    else if (n == 660) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -2.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.8070;
    }
    else if (n == 661) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         4.7860;
    }
    else if (n == 662) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  3.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.7370;
    }
    else if (n == 663) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  0.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         4.6970;
    }
    else if (n == 664) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  0.0;
      nut_data[ 5] =         4.6770;
    }
    else if (n == 665) {
      nut_data[ 0] =  0.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  4.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.6640;
    }
    else if (n == 666) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.6380;
    }
    else if (n == 667) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  1.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.6310;
    }
    else if (n == 668) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.6210;
    }
    else if (n == 669) {
      nut_data[ 0] =  4.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.5230;
    }
    else if (n == 670) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.2130;
    }
    else if (n == 671) {
      nut_data[ 0] = -1.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  6.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.1630;
    }
    else if (n == 672) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         4.1260;
    }
    else if (n == 673) {
      nut_data[ 0] =  1.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         4.0370;
    }
    else if (n == 674) {
      nut_data[ 0] =  3.0;
      nut_data[ 1] =  1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  2.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         3.9600;
    }
    else if (n == 675) {
      nut_data[ 0] =  5.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  0.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         3.9250;
    }
    else if (n == 676) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] = -1.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  2.0;
      nut_data[ 5] =         3.5910;
    }
    else if (n == 677) {
      nut_data[ 0] =  2.0;
      nut_data[ 1] =  0.0;
      nut_data[ 2] =  2.0;
      nut_data[ 3] =  4.0;
      nut_data[ 4] =  1.0;
      nut_data[ 5] =         3.5540;
    }
  }

  else if (*S1 == 'L' && *S2 == 'B') {
    if (n ==   0) {
      nut_data[ 0] =    -17206.4161;
      nut_data[ 1] =       -17.4666;
      nut_data[ 2] =      9205.2331;
      nut_data[ 3] =         0.9086;
      nut_data[ 4] =         3.3386;
      nut_data[ 5] =         0.0029;
      nut_data[ 6] =         1.5377;
      nut_data[ 7] =         0.0002;
    }
    else if (n ==   1) {
      nut_data[ 0] =     -1317.0906;
      nut_data[ 1] =        -0.1675;
      nut_data[ 2] =       573.0336;
      nut_data[ 3] =        -0.3015;
      nut_data[ 4] =        -1.3696;
      nut_data[ 5] =         0.0012;
      nut_data[ 6] =        -0.4587;
      nut_data[ 7] =        -0.0003;
    }
    else if (n ==   2) {
      nut_data[ 0] =      -227.6413;
      nut_data[ 1] =        -0.0234;
      nut_data[ 2] =        97.8459;
      nut_data[ 3] =        -0.0485;
      nut_data[ 4] =         0.2796;
      nut_data[ 5] =         0.0002;
      nut_data[ 6] =         0.1374;
      nut_data[ 7] =        -0.0001;
    }
    else if (n ==   3) {
      nut_data[ 0] =       207.4554;
      nut_data[ 1] =         0.0207;
      nut_data[ 2] =       -89.7492;
      nut_data[ 3] =         0.0470;
      nut_data[ 4] =        -0.0698;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0291;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==   4) {
      nut_data[ 0] =       147.5877;
      nut_data[ 1] =        -0.3633;
      nut_data[ 2] =         7.3871;
      nut_data[ 3] =        -0.0184;
      nut_data[ 4] =         1.1817;
      nut_data[ 5] =        -0.0015;
      nut_data[ 6] =        -0.1924;
      nut_data[ 7] =         0.0005;
    }
    else if (n ==   5) {
      nut_data[ 0] =       -51.6821;
      nut_data[ 1] =         0.1226;
      nut_data[ 2] =        22.4386;
      nut_data[ 3] =        -0.0677;
      nut_data[ 4] =        -0.0524;
      nut_data[ 5] =         0.0002;
      nut_data[ 6] =        -0.0174;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==   6) {
      nut_data[ 0] =        71.1159;
      nut_data[ 1] =         0.0073;
      nut_data[ 2] =        -0.6750;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0872;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0358;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==   7) {
      nut_data[ 0] =       -38.7298;
      nut_data[ 1] =        -0.0367;
      nut_data[ 2] =        20.0728;
      nut_data[ 3] =         0.0018;
      nut_data[ 4] =         0.0380;
      nut_data[ 5] =         0.0001;
      nut_data[ 6] =         0.0318;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==   8) {
      nut_data[ 0] =       -30.1461;
      nut_data[ 1] =        -0.0036;
      nut_data[ 2] =        12.9025;
      nut_data[ 3] =        -0.0063;
      nut_data[ 4] =         0.0816;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0367;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==   9) {
      nut_data[ 0] =        21.5829;
      nut_data[ 1] =        -0.0494;
      nut_data[ 2] =        -9.5929;
      nut_data[ 3] =         0.0299;
      nut_data[ 4] =         0.0111;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0132;
      nut_data[ 7] =        -0.0001;
    }
    else if (n ==  10) {
      nut_data[ 0] =        12.8227;
      nut_data[ 1] =         0.0137;
      nut_data[ 2] =        -6.8982;
      nut_data[ 3] =        -0.0009;
      nut_data[ 4] =         0.0181;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0039;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  11) {
      nut_data[ 0] =        12.3457;
      nut_data[ 1] =         0.0011;
      nut_data[ 2] =        -5.3311;
      nut_data[ 3] =         0.0032;
      nut_data[ 4] =         0.0019;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  12) {
      nut_data[ 0] =        15.6994;
      nut_data[ 1] =         0.0010;
      nut_data[ 2] =        -0.1235;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0168;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0082;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  13) {
      nut_data[ 0] =         6.3110;
      nut_data[ 1] =         0.0063;
      nut_data[ 2] =        -3.3228;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0027;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0009;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  14) {
      nut_data[ 0] =        -5.7976;
      nut_data[ 1] =        -0.0063;
      nut_data[ 2] =         3.1429;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0189;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0075;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  15) {
      nut_data[ 0] =        -5.9641;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =         2.5543;
      nut_data[ 3] =        -0.0011;
      nut_data[ 4] =         0.0149;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0066;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  16) {
      nut_data[ 0] =        -5.1613;
      nut_data[ 1] =        -0.0042;
      nut_data[ 2] =         2.6366;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0129;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0078;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  17) {
      nut_data[ 0] =         4.5893;
      nut_data[ 1] =         0.0050;
      nut_data[ 2] =        -2.4236;
      nut_data[ 3] =        -0.0010;
      nut_data[ 4] =         0.0031;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0020;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  18) {
      nut_data[ 0] =         6.3384;
      nut_data[ 1] =         0.0011;
      nut_data[ 2] =        -0.1220;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0150;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0029;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  19) {
      nut_data[ 0] =        -3.8571;
      nut_data[ 1] =        -0.0001;
      nut_data[ 2] =         1.6452;
      nut_data[ 3] =        -0.0011;
      nut_data[ 4] =         0.0158;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0068;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  20) {
      nut_data[ 0] =         3.2481;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.3870;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  21) {
      nut_data[ 0] =        -4.7722;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0477;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0018;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0025;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  22) {
      nut_data[ 0] =        -3.1046;
      nut_data[ 1] =        -0.0001;
      nut_data[ 2] =         1.3238;
      nut_data[ 3] =        -0.0011;
      nut_data[ 4] =         0.0131;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0059;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  23) {
      nut_data[ 0] =         2.8593;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.2338;
      nut_data[ 3] =         0.0010;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0003;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  24) {
      nut_data[ 0] =         2.0441;
      nut_data[ 1] =         0.0021;
      nut_data[ 2] =        -1.0758;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0010;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0003;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  25) {
      nut_data[ 0] =         2.9243;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0609;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0074;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0013;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  26) {
      nut_data[ 0] =         2.5887;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0550;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0066;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0011;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  27) {
      nut_data[ 0] =        -1.4053;
      nut_data[ 1] =        -0.0025;
      nut_data[ 2] =         0.8551;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0079;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0045;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  28) {
      nut_data[ 0] =         1.5164;
      nut_data[ 1] =         0.0010;
      nut_data[ 2] =        -0.8001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  29) {
      nut_data[ 0] =        -1.5794;
      nut_data[ 1] =         0.0072;
      nut_data[ 2] =         0.6850;
      nut_data[ 3] =        -0.0042;
      nut_data[ 4] =        -0.0016;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0005;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  30) {
      nut_data[ 0] =         2.1783;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0167;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0013;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0013;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  31) {
      nut_data[ 0] =        -1.2873;
      nut_data[ 1] =        -0.0010;
      nut_data[ 2] =         0.6953;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0037;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0014;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  32) {
      nut_data[ 0] =        -1.2654;
      nut_data[ 1] =         0.0011;
      nut_data[ 2] =         0.6415;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0063;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0026;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  33) {
      nut_data[ 0] =        -1.0204;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.5222;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0025;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0015;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  34) {
      nut_data[ 0] =         1.6707;
      nut_data[ 1] =        -0.0085;
      nut_data[ 2] =         0.0168;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =        -0.0010;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0010;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  35) {
      nut_data[ 0] =        -0.7691;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.3268;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0044;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0019;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  36) {
      nut_data[ 0] =        -1.1024;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0104;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0014;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  37) {
      nut_data[ 0] =         0.7566;
      nut_data[ 1] =        -0.0021;
      nut_data[ 2] =        -0.3250;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0011;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0005;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  38) {
      nut_data[ 0] =        -0.6637;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =         0.3353;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0025;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0014;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  39) {
      nut_data[ 0] =        -0.7141;
      nut_data[ 1] =         0.0021;
      nut_data[ 2] =         0.3070;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  40) {
      nut_data[ 0] =        -0.6302;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =         0.3272;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  41) {
      nut_data[ 0] =         0.5800;
      nut_data[ 1] =         0.0010;
      nut_data[ 2] =        -0.3045;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  42) {
      nut_data[ 0] =         0.6443;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.2768;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0007;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  43) {
      nut_data[ 0] =        -0.5774;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =         0.3041;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0015;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0005;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  44) {
      nut_data[ 0] =        -0.5350;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.2695;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0021;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0012;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  45) {
      nut_data[ 0] =        -0.4752;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =         0.2719;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0003;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  46) {
      nut_data[ 0] =        -0.4940;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =         0.2720;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0021;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0009;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  47) {
      nut_data[ 0] =         0.7350;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0051;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0008;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  48) {
      nut_data[ 0] =         0.4065;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.2206;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  49) {
      nut_data[ 0] =         0.6579;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0199;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0024;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  50) {
      nut_data[ 0] =         0.3579;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.1900;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  51) {
      nut_data[ 0] =         0.4725;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0041;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0006;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0003;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  52) {
      nut_data[ 0] =        -0.3075;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.1313;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  53) {
      nut_data[ 0] =        -0.2904;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.1233;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0015;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0007;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  54) {
      nut_data[ 0] =         0.4348;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0081;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0010;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  55) {
      nut_data[ 0] =        -0.2878;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.1232;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  56) {
      nut_data[ 0] =        -0.4230;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0020;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  57) {
      nut_data[ 0] =        -0.2819;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.1207;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0003;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  58) {
      nut_data[ 0] =        -0.4056;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0040;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  59) {
      nut_data[ 0] =        -0.2647;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.1129;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0005;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  60) {
      nut_data[ 0] =        -0.2294;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.1266;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0010;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  61) {
      nut_data[ 0] =         0.2481;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.1062;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0007;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0003;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  62) {
      nut_data[ 0] =         0.2179;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.1129;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  63) {
      nut_data[ 0] =         0.3276;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  64) {
      nut_data[ 0] =        -0.3389;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0035;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  65) {
      nut_data[ 0] =         0.3339;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0107;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0013;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  66) {
      nut_data[ 0] =        -0.1987;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.1073;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0006;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  67) {
      nut_data[ 0] =        -0.1981;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0854;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  68) {
      nut_data[ 0] =         0.4026;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0553;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0353;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0139;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  69) {
      nut_data[ 0] =         0.1660;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0710;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  70) {
      nut_data[ 0] =        -0.1521;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0647;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  71) {
      nut_data[ 0] =         0.1314;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0700;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  72) {
      nut_data[ 0] =        -0.1283;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0672;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  73) {
      nut_data[ 0] =        -0.1331;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0663;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  74) {
      nut_data[ 0] =         0.1383;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0594;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  75) {
      nut_data[ 0] =         0.1405;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0610;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  76) {
      nut_data[ 0] =         0.1290;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0556;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  77) {
      nut_data[ 0] =        -0.1214;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0518;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  78) {
      nut_data[ 0] =         0.1146;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0490;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  79) {
      nut_data[ 0] =         0.1019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0527;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  80) {
      nut_data[ 0] =        -0.1100;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0465;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  81) {
      nut_data[ 0] =        -0.0970;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0496;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  82) {
      nut_data[ 0] =         0.1575;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0050;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0006;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  83) {
      nut_data[ 0] =         0.0934;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0399;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  84) {
      nut_data[ 0] =         0.0922;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0395;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  85) {
      nut_data[ 0] =         0.0815;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0422;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  86) {
      nut_data[ 0] =         0.0834;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0440;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  87) {
      nut_data[ 0] =         0.1248;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0170;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  88) {
      nut_data[ 0] =         0.1338;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0039;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  89) {
      nut_data[ 0] =         0.0716;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0389;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  90) {
      nut_data[ 0] =         0.1282;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0023;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  91) {
      nut_data[ 0] =         0.0742;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0391;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  92) {
      nut_data[ 0] =         0.1020;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0495;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0025;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0010;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  93) {
      nut_data[ 0] =         0.0715;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0326;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0004;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  94) {
      nut_data[ 0] =        -0.0666;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0369;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  95) {
      nut_data[ 0] =        -0.0667;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0346;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  96) {
      nut_data[ 0] =        -0.0704;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0304;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  97) {
      nut_data[ 0] =        -0.0694;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0294;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  98) {
      nut_data[ 0] =        -0.1014;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n ==  99) {
      nut_data[ 0] =        -0.0585;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0316;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 100) {
      nut_data[ 0] =        -0.0949;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 101) {
      nut_data[ 0] =        -0.0595;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0258;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 102) {
      nut_data[ 0] =         0.0528;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0279;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 103) {
      nut_data[ 0] =        -0.0590;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0252;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 104) {
      nut_data[ 0] =         0.0570;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0244;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 105) {
      nut_data[ 0] =        -0.0502;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0250;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 106) {
      nut_data[ 0] =        -0.0875;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0029;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 107) {
      nut_data[ 0] =        -0.0492;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0275;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 108) {
      nut_data[ 0] =         0.0535;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0228;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 109) {
      nut_data[ 0] =        -0.0467;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0240;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 110) {
      nut_data[ 0] =         0.0591;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0253;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 111) {
      nut_data[ 0] =        -0.0453;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0244;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 112) {
      nut_data[ 0] =         0.0766;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 113) {
      nut_data[ 0] =        -0.0446;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0225;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 114) {
      nut_data[ 0] =        -0.0488;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0207;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 115) {
      nut_data[ 0] =        -0.0468;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0201;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 116) {
      nut_data[ 0] =        -0.0421;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0216;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 117) {
      nut_data[ 0] =         0.0463;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0200;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 118) {
      nut_data[ 0] =        -0.0673;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0014;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 119) {
      nut_data[ 0] =         0.0658;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 120) {
      nut_data[ 0] =        -0.0438;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0188;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 121) {
      nut_data[ 0] =        -0.0390;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0205;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 122) {
      nut_data[ 0] =         0.0639;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =        -0.0019;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 123) {
      nut_data[ 0] =         0.0412;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0176;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 124) {
      nut_data[ 0] =        -0.0361;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0189;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 125) {
      nut_data[ 0] =         0.0360;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0185;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 126) {
      nut_data[ 0] =         0.0588;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0024;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 127) {
      nut_data[ 0] =        -0.0578;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 128) {
      nut_data[ 0] =        -0.0396;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0171;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 129) {
      nut_data[ 0] =         0.0565;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 130) {
      nut_data[ 0] =        -0.0335;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0184;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 131) {
      nut_data[ 0] =         0.0357;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0154;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 132) {
      nut_data[ 0] =         0.0321;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0174;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 133) {
      nut_data[ 0] =        -0.0301;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0162;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 134) {
      nut_data[ 0] =        -0.0334;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0144;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 135) {
      nut_data[ 0] =         0.0493;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 136) {
      nut_data[ 0] =         0.0494;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0019;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 137) {
      nut_data[ 0] =         0.0337;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0143;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 138) {
      nut_data[ 0] =         0.0280;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0144;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 139) {
      nut_data[ 0] =         0.0309;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0134;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 140) {
      nut_data[ 0] =        -0.0263;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0131;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 141) {
      nut_data[ 0] =         0.0253;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0138;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 142) {
      nut_data[ 0] =         0.0245;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0128;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 143) {
      nut_data[ 0] =         0.0416;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0017;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 144) {
      nut_data[ 0] =        -0.0229;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0128;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 145) {
      nut_data[ 0] =         0.0231;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0120;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 146) {
      nut_data[ 0] =        -0.0259;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0109;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 147) {
      nut_data[ 0] =         0.0375;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 148) {
      nut_data[ 0] =         0.0252;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0108;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 149) {
      nut_data[ 0] =        -0.0245;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0104;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 150) {
      nut_data[ 0] =         0.0243;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0104;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 151) {
      nut_data[ 0] =         0.0208;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0112;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 152) {
      nut_data[ 0] =         0.0199;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0102;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 153) {
      nut_data[ 0] =        -0.0208;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0105;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 154) {
      nut_data[ 0] =         0.0335;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0014;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 155) {
      nut_data[ 0] =        -0.0325;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 156) {
      nut_data[ 0] =        -0.0187;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0096;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 157) {
      nut_data[ 0] =         0.0197;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0100;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 158) {
      nut_data[ 0] =        -0.0192;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0094;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 159) {
      nut_data[ 0] =        -0.0188;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0083;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 160) {
      nut_data[ 0] =         0.0276;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 161) {
      nut_data[ 0] =        -0.0286;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 162) {
      nut_data[ 0] =         0.0186;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0079;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 163) {
      nut_data[ 0] =        -0.0219;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0043;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 164) {
      nut_data[ 0] =         0.0276;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 165) {
      nut_data[ 0] =        -0.0153;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0084;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 166) {
      nut_data[ 0] =        -0.0156;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0081;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 167) {
      nut_data[ 0] =        -0.0154;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0078;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 168) {
      nut_data[ 0] =        -0.0174;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0075;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 169) {
      nut_data[ 0] =        -0.0163;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0069;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 170) {
      nut_data[ 0] =        -0.0228;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 171) {
      nut_data[ 0] =         0.0091;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0054;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0004;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 172) {
      nut_data[ 0] =         0.0175;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0075;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 173) {
      nut_data[ 0] =        -0.0159;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0069;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 174) {
      nut_data[ 0] =         0.0141;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0072;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 175) {
      nut_data[ 0] =         0.0147;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0075;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 176) {
      nut_data[ 0] =        -0.0132;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0069;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 177) {
      nut_data[ 0] =         0.0159;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0054;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0028;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0011;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 178) {
      nut_data[ 0] =         0.0213;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 179) {
      nut_data[ 0] =         0.0123;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0064;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 180) {
      nut_data[ 0] =        -0.0118;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0066;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 181) {
      nut_data[ 0] =         0.0144;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0061;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 182) {
      nut_data[ 0] =        -0.0121;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0060;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 183) {
      nut_data[ 0] =        -0.0134;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0056;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 184) {
      nut_data[ 0] =        -0.0105;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0057;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 185) {
      nut_data[ 0] =        -0.0102;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0056;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 186) {
      nut_data[ 0] =         0.0120;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0052;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 187) {
      nut_data[ 0] =         0.0101;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0054;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 188) {
      nut_data[ 0] =        -0.0113;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0059;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 189) {
      nut_data[ 0] =        -0.0106;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0061;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 190) {
      nut_data[ 0] =        -0.0129;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0055;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 191) {
      nut_data[ 0] =        -0.0114;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0057;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 192) {
      nut_data[ 0] =         0.0113;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0049;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 193) {
      nut_data[ 0] =        -0.0102;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0044;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 194) {
      nut_data[ 0] =        -0.0094;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0051;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 195) {
      nut_data[ 0] =        -0.0100;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0056;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 196) {
      nut_data[ 0] =         0.0087;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0047;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 197) {
      nut_data[ 0] =         0.0161;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 198) {
      nut_data[ 0] =         0.0096;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0050;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 199) {
      nut_data[ 0] =         0.0151;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 200) {
      nut_data[ 0] =        -0.0104;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0044;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 201) {
      nut_data[ 0] =        -0.0110;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0048;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 202) {
      nut_data[ 0] =        -0.0100;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0050;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 203) {
      nut_data[ 0] =         0.0092;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0012;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 204) {
      nut_data[ 0] =         0.0082;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0045;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 205) {
      nut_data[ 0] =         0.0082;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0045;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 206) {
      nut_data[ 0] =        -0.0078;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0041;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 207) {
      nut_data[ 0] =        -0.0077;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0043;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 208) {
      nut_data[ 0] =         0.0002;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0054;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 209) {
      nut_data[ 0] =         0.0094;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0040;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 210) {
      nut_data[ 0] =        -0.0093;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0040;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 211) {
      nut_data[ 0] =        -0.0083;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0040;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0010;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 212) {
      nut_data[ 0] =         0.0083;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0036;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 213) {
      nut_data[ 0] =        -0.0091;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0039;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 214) {
      nut_data[ 0] =         0.0128;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 215) {
      nut_data[ 0] =        -0.0079;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0034;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 216) {
      nut_data[ 0] =        -0.0083;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0047;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 217) {
      nut_data[ 0] =         0.0084;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0044;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 218) {
      nut_data[ 0] =         0.0083;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0043;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 219) {
      nut_data[ 0] =         0.0091;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0039;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 220) {
      nut_data[ 0] =        -0.0077;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0039;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 221) {
      nut_data[ 0] =         0.0084;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0043;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 222) {
      nut_data[ 0] =        -0.0092;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0039;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 223) {
      nut_data[ 0] =        -0.0092;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0039;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 224) {
      nut_data[ 0] =        -0.0094;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 225) {
      nut_data[ 0] =         0.0068;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0036;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 226) {
      nut_data[ 0] =        -0.0061;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0032;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 227) {
      nut_data[ 0] =         0.0071;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0031;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 228) {
      nut_data[ 0] =         0.0062;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0034;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 229) {
      nut_data[ 0] =        -0.0063;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0033;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 230) {
      nut_data[ 0] =        -0.0073;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0032;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 231) {
      nut_data[ 0] =         0.0115;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 232) {
      nut_data[ 0] =        -0.0103;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 233) {
      nut_data[ 0] =         0.0063;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0028;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 234) {
      nut_data[ 0] =         0.0074;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0032;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 235) {
      nut_data[ 0] =        -0.0103;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 236) {
      nut_data[ 0] =        -0.0069;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0030;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 237) {
      nut_data[ 0] =         0.0057;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0029;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 238) {
      nut_data[ 0] =         0.0094;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 239) {
      nut_data[ 0] =         0.0064;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0033;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 240) {
      nut_data[ 0] =        -0.0063;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0026;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 241) {
      nut_data[ 0] =        -0.0038;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0020;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 242) {
      nut_data[ 0] =        -0.0043;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0024;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 243) {
      nut_data[ 0] =        -0.0045;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0023;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 244) {
      nut_data[ 0] =         0.0047;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0024;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 245) {
      nut_data[ 0] =        -0.0048;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0025;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 246) {
      nut_data[ 0] =         0.0045;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0026;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 247) {
      nut_data[ 0] =         0.0056;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0025;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 248) {
      nut_data[ 0] =         0.0088;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 249) {
      nut_data[ 0] =        -0.0075;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 250) {
      nut_data[ 0] =         0.0085;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 251) {
      nut_data[ 0] =         0.0049;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0026;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 252) {
      nut_data[ 0] =        -0.0074;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 253) {
      nut_data[ 0] =        -0.0039;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0021;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 254) {
      nut_data[ 0] =         0.0045;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0020;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 255) {
      nut_data[ 0] =         0.0051;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0022;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 256) {
      nut_data[ 0] =        -0.0040;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0021;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 257) {
      nut_data[ 0] =         0.0041;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0021;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 258) {
      nut_data[ 0] =        -0.0042;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0024;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 259) {
      nut_data[ 0] =        -0.0051;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0022;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 260) {
      nut_data[ 0] =        -0.0042;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0022;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 261) {
      nut_data[ 0] =         0.0039;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0021;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 262) {
      nut_data[ 0] =         0.0046;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0018;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 263) {
      nut_data[ 0] =        -0.0053;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0022;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 264) {
      nut_data[ 0] =         0.0082;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 265) {
      nut_data[ 0] =         0.0081;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0001;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 266) {
      nut_data[ 0] =         0.0047;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0019;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 267) {
      nut_data[ 0] =         0.0053;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0023;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 268) {
      nut_data[ 0] =        -0.0045;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0022;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 269) {
      nut_data[ 0] =        -0.0044;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 270) {
      nut_data[ 0] =        -0.0033;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0016;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 271) {
      nut_data[ 0] =        -0.0061;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 272) {
      nut_data[ 0] =         0.0028;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 273) {
      nut_data[ 0] =        -0.0038;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0019;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 274) {
      nut_data[ 0] =        -0.0033;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0021;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 275) {
      nut_data[ 0] =        -0.0060;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 276) {
      nut_data[ 0] =         0.0048;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 277) {
      nut_data[ 0] =         0.0027;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0014;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 278) {
      nut_data[ 0] =         0.0038;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0020;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 279) {
      nut_data[ 0] =         0.0031;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0013;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 280) {
      nut_data[ 0] =        -0.0029;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 281) {
      nut_data[ 0] =         0.0028;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 282) {
      nut_data[ 0] =        -0.0032;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 283) {
      nut_data[ 0] =         0.0045;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 284) {
      nut_data[ 0] =        -0.0044;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0019;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 285) {
      nut_data[ 0] =         0.0028;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 286) {
      nut_data[ 0] =        -0.0051;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 287) {
      nut_data[ 0] =        -0.0036;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0020;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 288) {
      nut_data[ 0] =         0.0044;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0019;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 289) {
      nut_data[ 0] =         0.0026;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0014;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 290) {
      nut_data[ 0] =        -0.0060;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 291) {
      nut_data[ 0] =         0.0035;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0018;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 292) {
      nut_data[ 0] =        -0.0027;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0011;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 293) {
      nut_data[ 0] =         0.0047;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 294) {
      nut_data[ 0] =         0.0036;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 295) {
      nut_data[ 0] =        -0.0036;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0020;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 296) {
      nut_data[ 0] =        -0.0035;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0019;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 297) {
      nut_data[ 0] =        -0.0037;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0019;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 298) {
      nut_data[ 0] =         0.0032;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0016;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 299) {
      nut_data[ 0] =         0.0035;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0014;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 300) {
      nut_data[ 0] =         0.0032;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0013;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 301) {
      nut_data[ 0] =         0.0065;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 302) {
      nut_data[ 0] =         0.0047;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 303) {
      nut_data[ 0] =         0.0032;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0016;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 304) {
      nut_data[ 0] =         0.0037;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0016;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 305) {
      nut_data[ 0] =        -0.0030;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 306) {
      nut_data[ 0] =        -0.0032;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0016;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 307) {
      nut_data[ 0] =        -0.0031;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0013;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 308) {
      nut_data[ 0] =         0.0037;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0016;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 309) {
      nut_data[ 0] =         0.0031;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0013;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 310) {
      nut_data[ 0] =         0.0049;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 311) {
      nut_data[ 0] =         0.0032;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0013;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 312) {
      nut_data[ 0] =         0.0023;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0012;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 313) {
      nut_data[ 0] =        -0.0043;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0018;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 314) {
      nut_data[ 0] =         0.0026;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0011;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 315) {
      nut_data[ 0] =        -0.0032;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0014;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 316) {
      nut_data[ 0] =        -0.0029;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0014;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 317) {
      nut_data[ 0] =        -0.0027;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0012;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 318) {
      nut_data[ 0] =         0.0030;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 319) {
      nut_data[ 0] =        -0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 320) {
      nut_data[ 0] =        -0.0021;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 321) {
      nut_data[ 0] =        -0.0034;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 322) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 323) {
      nut_data[ 0] =        -0.0036;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 324) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 325) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 326) {
      nut_data[ 0] =        -0.0021;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 327) {
      nut_data[ 0] =        -0.0029;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 328) {
      nut_data[ 0] =        -0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 329) {
      nut_data[ 0] =        -0.0020;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 330) {
      nut_data[ 0] =         0.0028;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 331) {
      nut_data[ 0] =         0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 332) {
      nut_data[ 0] =        -0.0022;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0012;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 333) {
      nut_data[ 0] =        -0.0014;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 334) {
      nut_data[ 0] =         0.0024;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0011;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 335) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 336) {
      nut_data[ 0] =         0.0014;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 337) {
      nut_data[ 0] =         0.0024;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 338) {
      nut_data[ 0] =         0.0018;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 339) {
      nut_data[ 0] =        -0.0038;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 340) {
      nut_data[ 0] =        -0.0031;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 341) {
      nut_data[ 0] =        -0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 342) {
      nut_data[ 0] =         0.0029;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 343) {
      nut_data[ 0] =        -0.0018;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 344) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 345) {
      nut_data[ 0] =        -0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 346) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 347) {
      nut_data[ 0] =         0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 348) {
      nut_data[ 0] =         0.0022;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0012;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 349) {
      nut_data[ 0] =         0.0020;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 350) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 351) {
      nut_data[ 0] =        -0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 352) {
      nut_data[ 0] =        -0.0014;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 353) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 354) {
      nut_data[ 0] =         0.0014;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 355) {
      nut_data[ 0] =         0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 356) {
      nut_data[ 0] =        -0.0034;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 357) {
      nut_data[ 0] =        -0.0020;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 358) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 359) {
      nut_data[ 0] =        -0.0018;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 360) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 361) {
      nut_data[ 0] =         0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 362) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 363) {
      nut_data[ 0] =         0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 364) {
      nut_data[ 0] =        -0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 365) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 366) {
      nut_data[ 0] =        -0.0018;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 367) {
      nut_data[ 0] =        -0.0035;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 368) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 369) {
      nut_data[ 0] =        -0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 370) {
      nut_data[ 0] =        -0.0026;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0011;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 371) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 372) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 373) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 374) {
      nut_data[ 0] =        -0.0021;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 375) {
      nut_data[ 0] =        -0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 376) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 377) {
      nut_data[ 0] =        -0.0029;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 378) {
      nut_data[ 0] =        -0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 379) {
      nut_data[ 0] =         0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 380) {
      nut_data[ 0] =         0.0022;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 381) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 382) {
      nut_data[ 0] =        -0.0020;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0011;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 383) {
      nut_data[ 0] =        -0.0020;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 384) {
      nut_data[ 0] =        -0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 385) {
      nut_data[ 0] =         0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 386) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 387) {
      nut_data[ 0] =         0.0014;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 388) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 389) {
      nut_data[ 0] =         0.0025;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 390) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 391) {
      nut_data[ 0] =        -0.0014;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 392) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 393) {
      nut_data[ 0] =        -0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 394) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 395) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 396) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 397) {
      nut_data[ 0] =        -0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 398) {
      nut_data[ 0] =        -0.0022;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 399) {
      nut_data[ 0] =         0.0028;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 400) {
      nut_data[ 0] =         0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 401) {
      nut_data[ 0] =         0.0023;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 402) {
      nut_data[ 0] =         0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 403) {
      nut_data[ 0] =         0.0029;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 404) {
      nut_data[ 0] =        -0.0025;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 405) {
      nut_data[ 0] =         0.0022;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 406) {
      nut_data[ 0] =        -0.0018;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 407) {
      nut_data[ 0] =         0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 408) {
      nut_data[ 0] =        -0.0023;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 409) {
      nut_data[ 0] =         0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 410) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 411) {
      nut_data[ 0] =        -0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 412) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 413) {
      nut_data[ 0] =         0.0021;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 414) {
      nut_data[ 0] =         0.0023;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 415) {
      nut_data[ 0] =        -0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 416) {
      nut_data[ 0] =        -0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 417) {
      nut_data[ 0] =        -0.0022;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 418) {
      nut_data[ 0] =         0.0027;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 419) {
      nut_data[ 0] =         0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 420) {
      nut_data[ 0] =         0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 421) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 422) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 423) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 424) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 425) {
      nut_data[ 0] =         0.0018;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 426) {
      nut_data[ 0] =         0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 427) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 428) {
      nut_data[ 0] =        -0.0023;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 429) {
      nut_data[ 0] =         0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 430) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 431) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 432) {
      nut_data[ 0] =         0.0030;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 433) {
      nut_data[ 0] =         0.0024;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 434) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 435) {
      nut_data[ 0] =        -0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 436) {
      nut_data[ 0] =        -0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 437) {
      nut_data[ 0] =         0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 438) {
      nut_data[ 0] =        -0.0024;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 439) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 440) {
      nut_data[ 0] =        -0.0024;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0011;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 441) {
      nut_data[ 0] =        -0.0023;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 442) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 443) {
      nut_data[ 0] =        -0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 444) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.1988;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.1679;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 445) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0063;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0027;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 446) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 447) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0004;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 448) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 449) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0364;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0176;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 450) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.1044;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0891;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 451) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 452) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 453) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0330;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 454) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 455) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 456) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 457) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 458) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 459) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 460) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 461) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 462) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 463) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 464) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 465) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 466) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 467) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 468) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 469) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 470) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 471) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 472) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 473) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 474) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0012;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0010;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 475) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 476) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 477) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 478) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 479) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 480) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 481) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 482) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 483) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 484) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 485) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 486) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 487) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 488) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 489) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 490) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 491) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 492) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 493) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 494) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 495) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 496) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 497) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 498) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 499) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 500) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 501) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 502) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 503) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 504) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 505) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 506) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 507) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 508) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 509) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 510) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 511) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 512) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 513) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 514) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 515) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 516) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 517) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 518) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 519) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 520) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 521) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 522) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 523) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 524) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 525) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 526) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 527) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 528) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 529) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 530) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 531) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 532) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0013;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0005;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 533) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0030;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0014;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 534) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0162;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0138;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 535) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0075;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 536) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 537) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 538) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 539) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 540) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 541) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 542) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 543) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 544) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 545) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 546) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 547) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 548) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 549) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 550) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 551) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 552) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 553) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 554) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 555) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 556) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 557) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 558) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 559) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 560) {
      nut_data[ 0] =        -0.0001;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0001;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 561) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 562) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0013;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0011;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 563) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 564) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 565) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 566) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 567) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 568) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 569) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 570) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 571) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 572) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 573) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 574) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 575) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 576) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 577) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 578) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 579) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 580) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 581) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 582) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 583) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 584) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 585) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 586) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 587) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 588) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 589) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 590) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 591) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 592) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 593) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 594) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 595) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 596) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0026;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0011;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 597) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0010;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0005;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 598) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 599) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 600) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 601) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 602) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 603) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 604) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 605) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 606) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 607) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 608) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 609) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 610) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 611) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 612) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 613) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 614) {
      nut_data[ 0] =        -0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 615) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 616) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 617) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 618) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 619) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 620) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 621) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 622) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 623) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 624) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 625) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 626) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 627) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 628) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =        -0.0005;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -0.0002;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 629) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 630) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 631) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 632) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 633) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 634) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 635) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 636) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 637) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 638) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 639) {
      nut_data[ 0] =         0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 640) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 641) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 642) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 643) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 644) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 645) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 646) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 647) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 648) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 649) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 650) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 651) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 652) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 653) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 654) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 655) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 656) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 657) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 658) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 659) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 660) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 661) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 662) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 663) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 664) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 665) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 666) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 667) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 668) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 669) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 670) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 671) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 672) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 673) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 674) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 675) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 676) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
    else if (n == 677) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
    }
  }

  else if (*S1 == 'P' && *S2 == 'A') {
    if (n ==   0) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         6.8500;
    }
    else if (n ==   1) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         8.9900;
    }
    else if (n ==   2) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.1100;
    }
    else if (n ==   3) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.1200;
    }
    else if (n ==   4) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         6.7300;
    }
    else if (n ==   5) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.1700;
    }
    else if (n ==   6) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.1300;
    }
    else if (n ==   7) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.1300;
    }
    else if (n ==   8) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.1300;
    }
    else if (n ==   9) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         6.6400;
    }
    else if (n ==  10) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.1400;
    }
    else if (n ==  11) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.2800;
    }
    else if (n ==  12) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.7900;
    }
    else if (n ==  13) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.0500;
    }
    else if (n ==  14) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.2100;
    }
    else if (n ==  15) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.3500;
    }
    else if (n ==  16) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.5300;
    }
    else if (n ==  17) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.6300;
    }
    else if (n ==  18) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        25.1300;
    }
    else if (n ==  19) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.6300;
    }
    else if (n ==  20) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         9.5400;
    }
    else if (n ==  21) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.6500;
    }
    else if (n ==  22) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.0900;
    }
    else if (n ==  23) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.6600;
    }
    else if (n ==  24) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.6600;
    }
    else if (n ==  25) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         8.9100;
    }
    else if (n ==  26) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.6600;
    }
    else if (n ==  27) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        10.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.6600;
    }
    else if (n ==  28) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =         8.9100;
    }
    else if (n ==  29) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.6900;
    }
    else if (n ==  30) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.7500;
    }
    else if (n ==  31) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.7900;
    }
    else if (n ==  32) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        13.9900;
    }
    else if (n ==  33) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        14.1500;
    }
    else if (n ==  34) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        14.3300;
    }
    else if (n ==  35) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -15.3300;
    }
    else if (n ==  36) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -15.5500;
    }
    else if (n ==  37) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        25.7700;
    }
    else if (n ==  38) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        26.3100;
    }
    else if (n ==  39) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.2100;
    }
    else if (n ==  40) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.3000;
    }
    else if (n ==  41) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -27.3300;
    }
    else if (n ==  42) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.3400;
    }
    else if (n ==  43) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.4400;
    }
    else if (n ==  44) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       199.4400;
    }
    else if (n ==  45) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.5500;
    }
    else if (n ==  46) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -27.5500;
    }
    else if (n ==  47) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.5500;
    }
    else if (n ==  48) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -27.5500;
    }
    else if (n ==  49) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.5600;
    }
    else if (n ==  50) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -14.7500;
    }
    else if (n ==  51) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.6100;
    }
    else if (n ==  52) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        27.6700;
    }
    else if (n ==  53) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        28.6900;
    }
    else if (n ==  54) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -28.9200;
    }
    else if (n ==  55) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -29.6000;
    }
    else if (n ==  56) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -31.7400;
    }
    else if (n ==  57) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -205.8300;
    }
    else if (n ==  58) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -31.8100;
    }
    else if (n ==  59) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        31.8100;
    }
    else if (n ==  60) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -31.9700;
    }
    else if (n ==  61) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -32.1000;
    }
    else if (n ==  62) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -33.6400;
    }
    else if (n ==  63) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -34.5700;
    }
    else if (n ==  64) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -35.7000;
    }
    else if (n ==  65) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -37.8500;
    }
    else if (n ==  66) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =   4709189.0800;
    }
    else if (n ==  67) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        60.8800;
    }
    else if (n ==  68) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        60.8800;
    }
    else if (n ==  69) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        60.8800;
    }
    else if (n ==  70) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        60.8800;
    }
    else if (n ==  71) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        60.8800;
    }
    else if (n ==  72) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        60.8800;
    }
    else if (n ==  73) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        60.8800;
    }
    else if (n ==  74) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        60.8800;
    }
    else if (n ==  75) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         6.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        63.4900;
    }
    else if (n ==  76) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         9.0000;
      nut_data[ 7] =        -9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        64.8800;
    }
    else if (n ==  77) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -9.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       -64.8800;
    }
    else if (n ==  78) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -9.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       -64.8800;
    }
    else if (n ==  79) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -9.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       -64.8800;
    }
    else if (n ==  80) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        71.2300;
    }
    else if (n ==  81) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        71.2300;
    }
    else if (n ==  82) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         8.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        72.9900;
    }
    else if (n ==  83) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       -72.9900;
    }
    else if (n ==  84) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        74.9000;
    }
    else if (n ==  85) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        75.6000;
    }
    else if (n ==  86) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        76.9400;
    }
    else if (n ==  87) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        78.3300;
    }
    else if (n ==  88) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        81.1300;
    }
    else if (n ==  89) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        81.1300;
    }
    else if (n ==  90) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        81.1300;
    }
    else if (n ==  91) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        81.1300;
    }
    else if (n ==  92) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         7.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        83.4200;
    }
    else if (n ==  93) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -7.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       -83.4200;
    }
    else if (n ==  94) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        85.9200;
    }
    else if (n ==  95) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        85.9200;
    }
    else if (n ==  96) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        85.9200;
    }
    else if (n ==  97) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        88.4900;
    }
    else if (n ==  98) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         8.0000;
      nut_data[ 7] =        -9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        91.2200;
    }
    else if (n ==  99) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        91.3100;
    }
    else if (n == 100) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        92.7800;
    }
    else if (n == 101) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        92.8900;
    }
    else if (n == 102) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        93.2800;
    }
    else if (n == 103) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        94.2200;
    }
    else if (n == 104) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        94.2200;
    }
    else if (n == 105) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        94.2200;
    }
    else if (n == 106) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        94.3000;
    }
    else if (n == 107) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        95.3300;
    }
    else if (n == 108) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         6.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =        97.3200;
    }
    else if (n == 109) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         6.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =        97.3200;
    }
    else if (n == 110) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       -97.3200;
    }
    else if (n == 111) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        97.4800;
    }
    else if (n == 112) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        99.7200;
    }
    else if (n == 113) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -9.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -100.6300;
    }
    else if (n == 114) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -9.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -100.6300;
    }
    else if (n == 115) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       100.7400;
    }
    else if (n == 116) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       104.2900;
    }
    else if (n == 117) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -3.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       107.2700;
    }
    else if (n == 118) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         7.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       108.1100;
    }
    else if (n == 119) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -5.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       109.3000;
    }
    else if (n == 120) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =        -7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       111.4200;
    }
    else if (n == 121) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =        -7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       111.4200;
    }
    else if (n == 122) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       112.3500;
    }
    else if (n == 123) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       112.3500;
    }
    else if (n == 124) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       112.3500;
    }
    else if (n == 125) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       116.7800;
    }
    else if (n == 126) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       116.7800;
    }
    else if (n == 127) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       116.7800;
    }
    else if (n == 128) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -116.7900;
    }
    else if (n == 129) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       116.9400;
    }
    else if (n == 130) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       118.4200;
    }
    else if (n == 131) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       119.2300;
    }
    else if (n == 132) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        10.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -121.5900;
    }
    else if (n == 133) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        10.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -121.5800;
    }
    else if (n == 134) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        10.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -121.5900;
    }
    else if (n == 135) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       121.7300;
    }
    else if (n == 136) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       121.7500;
    }
    else if (n == 137) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =        16.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       121.7500;
    }
    else if (n == 138) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       121.7700;
    }
    else if (n == 139) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       121.7900;
    }
    else if (n == 140) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        16.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       121.9200;
    }
    else if (n == 141) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       123.0000;
    }
    else if (n == 142) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       123.1400;
    }
    else if (n == 143) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       124.3800;
    }
    else if (n == 144) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       124.5700;
    }
    else if (n == 145) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       125.2700;
    }
    else if (n == 146) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       125.2700;
    }
    else if (n == 147) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       126.9700;
    }
    else if (n == 148) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       126.9800;
    }
    else if (n == 149) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       127.1200;
    }
    else if (n == 150) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       129.0000;
    }
    else if (n == 151) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       129.0000;
    }
    else if (n == 152) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       129.9900;
    }
    else if (n == 153) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         6.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       132.6700;
    }
    else if (n == 154) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -132.6700;
    }
    else if (n == 155) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       132.8700;
    }
    else if (n == 156) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       132.9600;
    }
    else if (n == 157) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       132.9600;
    }
    else if (n == 158) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       132.9800;
    }
    else if (n == 159) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       136.3200;
    }
    else if (n == 160) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       137.1700;
    }
    else if (n == 161) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -9.0000;
      nut_data[ 7] =        12.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -138.9000;
    }
    else if (n == 162) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       139.1100;
    }
    else if (n == 163) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       139.1200;
    }
    else if (n == 164) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       139.1200;
    }
    else if (n == 165) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         5.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -141.6600;
    }
    else if (n == 166) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         1.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       144.2700;
    }
    else if (n == 167) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       145.9800;
    }
    else if (n == 168) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       145.9800;
    }
    else if (n == 169) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       145.9800;
    }
    else if (n == 170) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -145.9800;
    }
    else if (n == 171) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -1.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       147.9700;
    }
    else if (n == 172) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -3.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       151.8700;
    }
    else if (n == 173) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -7.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -153.5600;
    }
    else if (n == 174) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -7.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -153.5600;
    }
    else if (n == 175) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       153.8200;
    }
    else if (n == 176) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       153.8200;
    }
    else if (n == 177) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -5.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       155.9800;
    }
    else if (n == 178) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       160.3200;
    }
    else if (n == 179) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       162.2600;
    }
    else if (n == 180) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       162.2600;
    }
    else if (n == 181) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =        -9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       164.9100;
    }
    else if (n == 182) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       166.7700;
    }
    else if (n == 183) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       168.4200;
    }
    else if (n == 184) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       171.6700;
    }
    else if (n == 185) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       171.6700;
    }
    else if (n == 186) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -171.6800;
    }
    else if (n == 187) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -171.6800;
    }
    else if (n == 188) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       171.7400;
    }
    else if (n == 189) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       172.0100;
    }
    else if (n == 190) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       175.2300;
    }
    else if (n == 191) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       177.0100;
    }
    else if (n == 192) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       177.8000;
    }
    else if (n == 193) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       177.8500;
    }
    else if (n == 194) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       177.8900;
    }
    else if (n == 195) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       178.5900;
    }
    else if (n == 196) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       179.5700;
    }
    else if (n == 197) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =        33.2000;
    }
    else if (n == 198) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -182.2500;
    }
    else if (n == 199) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -182.2500;
    }
    else if (n == 200) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -182.2500;
    }
    else if (n == 201) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.3200;
    }
    else if (n == 202) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.5200;
    }
    else if (n == 203) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.5700;
    }
    else if (n == 204) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.5700;
    }
    else if (n == 205) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.5700;
    }
    else if (n == 206) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =        16.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.6200;
    }
    else if (n == 207) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.6200;
    }
    else if (n == 208) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       182.6200;
    }
    else if (n == 209) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       182.6200;
    }
    else if (n == 210) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       182.6200;
    }
    else if (n == 211) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       182.6300;
    }
    else if (n == 212) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       182.6300;
    }
    else if (n == 213) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.6700;
    }
    else if (n == 214) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.6700;
    }
    else if (n == 215) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.7200;
    }
    else if (n == 216) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       182.9200;
    }
    else if (n == 217) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        15.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       183.0000;
    }
    else if (n == 218) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =        15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       184.9100;
    }
    else if (n == 219) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       185.7700;
    }
    else if (n == 220) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       188.6000;
    }
    else if (n == 221) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       188.6000;
    }
    else if (n == 222) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       189.0400;
    }
    else if (n == 223) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       189.0500;
    }
    else if (n == 224) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       190.6600;
    }
    else if (n == 225) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       190.6700;
    }
    else if (n == 226) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =        13.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       191.0400;
    }
    else if (n == 227) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -3.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       192.4300;
    }
    else if (n == 228) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       194.6300;
    }
    else if (n == 229) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       194.6400;
    }
    else if (n == 230) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       194.6400;
    }
    else if (n == 231) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -194.6400;
    }
    else if (n == 232) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       194.9800;
    }
    else if (n == 233) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       194.9800;
    }
    else if (n == 234) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =        10.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       195.0700;
    }
    else if (n == 235) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       199.4300;
    }
    else if (n == 236) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       199.4400;
    }
    else if (n == 237) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       199.4400;
    }
    else if (n == 238) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       201.8000;
    }
    else if (n == 239) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       201.8100;
    }
    else if (n == 240) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       204.6000;
    }
    else if (n == 241) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       205.4700;
    }
    else if (n == 242) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         6.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       208.3400;
    }
    else if (n == 243) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -208.3500;
    }
    else if (n == 244) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -208.3500;
    }
    else if (n == 245) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       208.8300;
    }
    else if (n == 246) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       208.8400;
    }
    else if (n == 247) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       209.0700;
    }
    else if (n == 248) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       212.1300;
    }
    else if (n == 249) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       219.6700;
    }
    else if (n == 250) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         4.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -219.6800;
    }
    else if (n == 251) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         5.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       220.2300;
    }
    else if (n == 252) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -9.0000;
      nut_data[ 7] =        13.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -224.1300;
    }
    else if (n == 253) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       224.6900;
    }
    else if (n == 254) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       224.7000;
    }
    else if (n == 255) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         3.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       228.9800;
    }
    else if (n == 256) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         5.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -231.4100;
    }
    else if (n == 257) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         1.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       238.4600;
    }
    else if (n == 258) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       243.1600;
    }
    else if (n == 259) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -243.1700;
    }
    else if (n == 260) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -243.1800;
    }
    else if (n == 261) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       243.8300;
    }
    else if (n == 262) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -3.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       259.9800;
    }
    else if (n == 263) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -7.0000;
      nut_data[ 7] =        10.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -264.9400;
    }
    else if (n == 264) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -7.0000;
      nut_data[ 7] =        10.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -264.9500;
    }
    else if (n == 265) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       265.7300;
    }
    else if (n == 266) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       265.7300;
    }
    else if (n == 267) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       265.7300;
    }
    else if (n == 268) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -265.7400;
    }
    else if (n == 269) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -5.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       272.2600;
    }
    else if (n == 270) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -279.9400;
    }
    else if (n == 271) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       285.7600;
    }
    else if (n == 272) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -285.7700;
    }
    else if (n == 273) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       291.5100;
    }
    else if (n == 274) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       291.9500;
    }
    else if (n == 275) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       291.9600;
    }
    else if (n == 276) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -291.9700;
    }
    else if (n == 277) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -291.9700;
    }
    else if (n == 278) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       300.6600;
    }
    else if (n == 279) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -300.6800;
    }
    else if (n == 280) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       305.0600;
    }
    else if (n == 281) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       306.8900;
    }
    else if (n == 282) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       312.5400;
    }
    else if (n == 283) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =        -1.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -317.0600;
    }
    else if (n == 284) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -398.8800;
    }
    else if (n == 285) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       320.9300;
    }
    else if (n == 286) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       323.9200;
    }
    else if (n == 287) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -323.9300;
    }
    else if (n == 288) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -323.9300;
    }
    else if (n == 289) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -323.9400;
    }
    else if (n == 290) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       324.1500;
    }
    else if (n == 291) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       325.1000;
    }
    else if (n == 292) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       325.1100;
    }
    else if (n == 293) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -325.1300;
    }
    else if (n == 294) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       336.8300;
    }
    else if (n == 295) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       336.8500;
    }
    else if (n == 296) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       336.8600;
    }
    else if (n == 297) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -340.1500;
    }
    else if (n == 298) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       342.0100;
    }
    else if (n == 299) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       343.4600;
    }
    else if (n == 300) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       343.4900;
    }
    else if (n == 301) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       353.2400;
    }
    else if (n == 302) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       359.3600;
    }
    else if (n == 303) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        12.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -363.7300;
    }
    else if (n == 304) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        12.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -363.7600;
    }
    else if (n == 305) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       364.8200;
    }
    else if (n == 306) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       365.0200;
    }
    else if (n == 307) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -365.0500;
    }
    else if (n == 308) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       365.2600;
    }
    else if (n == 309) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -365.2700;
    }
    else if (n == 310) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       365.4300;
    }
    else if (n == 311) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       365.4600;
    }
    else if (n == 312) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       365.6700;
    }
    else if (n == 313) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        14.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       366.7600;
    }
    else if (n == 314) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -368.8100;
    }
    else if (n == 315) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -376.7800;
    }
    else if (n == 316) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       378.0900;
    }
    else if (n == 317) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       389.9700;
    }
    else if (n == 318) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       391.8600;
    }
    else if (n == 319) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -392.6200;
    }
    else if (n == 320) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       398.3900;
    }
    else if (n == 321) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       398.8700;
    }
    else if (n == 322) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       398.8800;
    }
    else if (n == 323) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -398.9000;
    }
    else if (n == 324) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       399.3800;
    }
    else if (n == 325) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -3.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       406.6700;
    }
    else if (n == 326) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -412.6600;
    }
    else if (n == 327) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       413.7000;
    }
    else if (n == 328) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       416.6700;
    }
    else if (n == 329) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       416.6900;
    }
    else if (n == 330) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -416.7100;
    }
    else if (n == 331) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -416.7300;
    }
    else if (n == 332) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       418.2700;
    }
    else if (n == 333) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       418.6500;
    }
    else if (n == 334) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       423.7500;
    }
    else if (n == 335) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       439.3300;
    }
    else if (n == 336) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -439.3700;
    }
    else if (n == 337) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       443.9000;
    }
    else if (n == 338) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       450.9900;
    }
    else if (n == 339) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -451.0400;
    }
    else if (n == 340) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       455.0000;
    }
    else if (n == 341) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       469.6800;
    }
    else if (n == 342) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         6.0000;
      nut_data[ 7] =        -9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       484.9800;
    }
    else if (n == 343) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -485.0000;
    }
    else if (n == 344) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -485.0300;
    }
    else if (n == 345) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       487.6400;
    }
    else if (n == 346) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       487.6600;
    }
    else if (n == 347) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       487.6600;
    }
    else if (n == 348) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -487.6900;
    }
    else if (n == 349) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       488.9100;
    }
    else if (n == 350) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -488.9300;
    }
    else if (n == 351) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -488.9600;
    }
    else if (n == 352) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       489.2700;
    }
    else if (n == 353) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -489.3300;
    }
    else if (n == 354) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       501.6700;
    }
    else if (n == 355) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       505.9900;
    }
    else if (n == 356) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -526.8500;
    }
    else if (n == 357) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =        10.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -534.6600;
    }
    else if (n == 358) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =        10.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -534.7200;
    }
    else if (n == 359) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -537.7300;
    }
    else if (n == 360) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       551.1000;
    }
    else if (n == 361) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         4.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -551.1600;
    }
    else if (n == 362) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         5.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       554.6800;
    }
    else if (n == 363) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       583.8900;
    }
    else if (n == 364) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       583.8900;
    }
    else if (n == 365) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       583.9200;
    }
    else if (n == 366) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -583.9600;
    }
    else if (n == 367) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -583.9600;
    }
    else if (n == 368) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -583.9900;
    }
    else if (n == 369) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =        12.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -589.4000;
    }
    else if (n == 370) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         3.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       613.7400;
    }
    else if (n == 371) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =        -3.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -613.8200;
    }
    else if (n == 372) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         5.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -631.4900;
    }
    else if (n == 373) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       638.7900;
    }
    else if (n == 374) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       638.8300;
    }
    else if (n == 375) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         1.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       686.9800;
    }
    else if (n == 376) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -698.4300;
    }
    else if (n == 377) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       727.4700;
    }
    else if (n == 378) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -727.5200;
    }
    else if (n == 379) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -727.5200;
    }
    else if (n == 380) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -727.5200;
    }
    else if (n == 381) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -727.5800;
    }
    else if (n == 382) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       733.4700;
    }
    else if (n == 383) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       733.5300;
    }
    else if (n == 384) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       733.5300;
    }
    else if (n == 385) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =        -1.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       779.9400;
    }
    else if (n == 386) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       814.7100;
    }
    else if (n == 387) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -4.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -814.7800;
    }
    else if (n == 388) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -3.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       901.9800;
    }
    else if (n == 389) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         4.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       925.7300;
    }
    else if (n == 390) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -7.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -964.7000;
    }
    else if (n == 391) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -7.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      -964.7900;
    }
    else if (n == 392) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =       975.1800;
    }
    else if (n == 393) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =       975.2800;
    }
    else if (n == 394) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -975.3800;
    }
    else if (n == 395) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      -975.4800;
    }
    else if (n == 396) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         1.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      -975.4800;
    }
    else if (n == 397) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -5.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1069.3200;
    }
    else if (n == 398) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         5.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -1069.4400;
    }
    else if (n == 399) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         5.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -1069.5600;
    }
    else if (n == 400) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         4.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      1082.9000;
    }
    else if (n == 401) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -1.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1138.7600;
    }
    else if (n == 402) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =        10.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      1152.5200;
    }
    else if (n == 403) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1190.8500;
    }
    else if (n == 404) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -1198.4500;
    }
    else if (n == 405) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1312.8800;
    }
    else if (n == 406) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -1313.0600;
    }
    else if (n == 407) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -1313.2400;
    }
    else if (n == 408) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -1322.7200;
    }
    else if (n == 409) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      1440.5600;
    }
    else if (n == 410) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -1441.0000;
    }
    else if (n == 411) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      1443.7500;
    }
    else if (n == 412) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      1443.9800;
    }
    else if (n == 413) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1443.9700;
    }
    else if (n == 414) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      1443.9800;
    }
    else if (n == 415) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      1454.7100;
    }
    else if (n == 416) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1454.9400;
    }
    else if (n == 417) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -1455.1600;
    }
    else if (n == 418) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -1455.1600;
    }
    else if (n == 419) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -2.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -1455.3900;
    }
    else if (n == 420) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =        10.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      1479.1400;
    }
    else if (n == 421) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =        10.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      1479.3700;
    }
    else if (n == 422) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -1581.4600;
    }
    else if (n == 423) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1616.0300;
    }
    else if (n == 424) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1642.2400;
    }
    else if (n == 425) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -1642.5300;
    }
    else if (n == 426) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1700.1200;
    }
    else if (n == 427) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -1700.4200;
    }
    else if (n == 428) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -1700.7300;
    }
    else if (n == 429) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1851.0900;
    }
    else if (n == 430) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1851.4600;
    }
    else if (n == 431) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      1920.5500;
    }
    else if (n == 432) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      1920.9400;
    }
    else if (n == 433) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      1920.9400;
    }
    else if (n == 434) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      1920.9400;
    }
    else if (n == 435) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -1921.3400;
    }
    else if (n == 436) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2014.6100;
    }
    else if (n == 437) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2023.1200;
    }
    else if (n == 438) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2023.9900;
    }
    else if (n == 439) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2060.8600;
    }
    else if (n == 440) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2061.3100;
    }
    else if (n == 441) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2061.7600;
    }
    else if (n == 442) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      2150.8600;
    }
    else if (n == 443) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2165.3000;
    }
    else if (n == 444) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      2165.3000;
    }
    else if (n == 445) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      2165.8000;
    }
    else if (n == 446) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2165.8000;
    }
    else if (n == 447) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      2165.8000;
    }
    else if (n == 448) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2166.2900;
    }
    else if (n == 449) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         4.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      2179.9300;
    }
    else if (n == 450) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =         9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2268.1500;
    }
    else if (n == 451) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2310.4900;
    }
    else if (n == 452) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2337.2000;
    }
    else if (n == 453) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         1.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2401.5600;
    }
    else if (n == 454) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =        -1.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -2402.7900;
    }
    else if (n == 455) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =       -11.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2411.3600;
    }
    else if (n == 456) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =        11.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -2412.5900;
    }
    else if (n == 457) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2644.7000;
    }
    else if (n == 458) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2645.4400;
    }
    else if (n == 459) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2646.1800;
    }
    else if (n == 460) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2646.9300;
    }
    else if (n == 461) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2862.1500;
    }
    else if (n == 462) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -2863.0200;
    }
    else if (n == 463) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2863.0200;
    }
    else if (n == 464) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -2863.0200;
    }
    else if (n == 465) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -2863.8900;
    }
    else if (n == 466) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2863.8900;
    }
    else if (n == 467) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -2863.8900;
    }
    else if (n == 468) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      2880.2400;
    }
    else if (n == 469) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      2881.1200;
    }
    else if (n == 470) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2881.1200;
    }
    else if (n == 471) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      2881.1200;
    }
    else if (n == 472) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2882.0000;
    }
    else if (n == 473) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2882.8900;
    }
    else if (n == 474) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -2882.8900;
    }
    else if (n == 475) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      2957.3500;
    }
    else if (n == 476) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2957.3500;
    }
    else if (n == 477) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      2957.3500;
    }
    else if (n == 478) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      2958.2800;
    }
    else if (n == 479) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      2958.2800;
    }
    else if (n == 480) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2959.2100;
    }
    else if (n == 481) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -2960.1400;
    }
    else if (n == 482) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -2960.1400;
    }
    else if (n == 483) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      3001.2600;
    }
    else if (n == 484) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3002.2200;
    }
    else if (n == 485) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3118.2700;
    }
    else if (n == 486) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3119.3100;
    }
    else if (n == 487) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3120.3400;
    }
    else if (n == 488) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      3178.3300;
    }
    else if (n == 489) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3287.6100;
    }
    else if (n == 490) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3381.5500;
    }
    else if (n == 491) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3417.0200;
    }
    else if (n == 492) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      3561.6600;
    }
    else if (n == 493) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         3.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      3583.6800;
    }
    else if (n == 494) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         3.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      3585.0400;
    }
    else if (n == 495) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         3.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      3585.0400;
    }
    else if (n == 496) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      3626.7300;
    }
    else if (n == 497) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      3680.3300;
    }
    else if (n == 498) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =        -1.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3716.2400;
    }
    else if (n == 499) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =        11.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -3739.7500;
    }
    else if (n == 500) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4076.3300;
    }
    else if (n == 501) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =       -13.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      4145.7000;
    }
    else if (n == 502) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =        13.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     -4149.3500;
    }
    else if (n == 503) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4164.1800;
    }
    else if (n == 504) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4166.0200;
    }
    else if (n == 505) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4167.8600;
    }
    else if (n == 506) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      4328.6100;
    }
    else if (n == 507) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      4330.6000;
    }
    else if (n == 508) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      4330.6000;
    }
    else if (n == 509) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      4330.6000;
    }
    else if (n == 510) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      4332.5900;
    }
    else if (n == 511) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -4334.5800;
    }
    else if (n == 512) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4334.5800;
    }
    else if (n == 513) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      4391.5700;
    }
    else if (n == 514) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4431.9800;
    }
    else if (n == 515) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4497.2200;
    }
    else if (n == 516) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4528.4900;
    }
    else if (n == 517) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =        15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4655.1200;
    }
    else if (n == 518) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         6.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4895.8600;
    }
    else if (n == 519) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      4943.3100;
    }
    else if (n == 520) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      4945.9100;
    }
    else if (n == 521) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -4948.5100;
    }
    else if (n == 522) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      5002.8200;
    }
    else if (n == 523) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      5237.2300;
    }
    else if (n == 524) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      5240.1500;
    }
    else if (n == 525) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -5243.0700;
    }
    else if (n == 526) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      5373.4700;
    }
    else if (n == 527) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      5376.5400;
    }
    else if (n == 528) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      5376.5400;
    }
    else if (n == 529) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -5382.6900;
    }
    else if (n == 530) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -5413.4000;
    }
    else if (n == 531) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -5692.3600;
    }
    else if (n == 532) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =      5760.4800;
    }
    else if (n == 533) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      5760.4800;
    }
    else if (n == 534) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -5764.0100;
    }
    else if (n == 535) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -5767.5400;
    }
    else if (n == 536) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     -5767.5400;
    }
    else if (n == 537) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        13.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6307.0400;
    }
    else if (n == 538) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        12.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6311.2600;
    }
    else if (n == 539) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      6361.9500;
    }
    else if (n == 540) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =       -18.0000;
      nut_data[ 7] =        16.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6364.5200;
    }
    else if (n == 541) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      6366.2500;
    }
    else if (n == 542) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6404.2400;
    }
    else if (n == 543) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      6494.4100;
    }
    else if (n == 544) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6609.1800;
    }
    else if (n == 545) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6653.3700;
    }
    else if (n == 546) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6658.0800;
    }
    else if (n == 547) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6728.1600;
    }
    else if (n == 548) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        10.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6732.8400;
    }
    else if (n == 549) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      6788.5800;
    }
    else if (n == 550) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =       -10.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6865.2200;
    }
    else if (n == 551) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6870.0800;
    }
    else if (n == 552) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6944.7300;
    }
    else if (n == 553) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6949.8600;
    }
    else if (n == 554) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -6998.7400;
    }
    else if (n == 555) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -7244.2200;
    }
    else if (n == 556) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -7259.0500;
    }
    else if (n == 557) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        18.0000;
      nut_data[ 7] =       -16.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -7295.7300;
    }
    else if (n == 558) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      7360.6500;
    }
    else if (n == 559) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      7361.2200;
    }
    else if (n == 560) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         8.0000;
      nut_data[ 7] =       -12.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -7366.9900;
    }
    else if (n == 561) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         8.0000;
      nut_data[ 7] =       -13.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -7372.7600;
    }
    else if (n == 562) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      7480.8800;
    }
    else if (n == 563) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -7819.7800;
    }
    else if (n == 564) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      8084.7600;
    }
    else if (n == 565) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         1.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =      8197.0300;
    }
    else if (n == 566) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -8277.9100;
    }
    else if (n == 567) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -8285.2000;
    }
    else if (n == 568) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -8437.8500;
    }
    else if (n == 569) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -9135.6800;
    }
    else if (n == 570) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -9.0000;
      nut_data[ 8] =        17.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =      9435.7900;
    }
    else if (n == 571) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =       -17.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     -9454.7300;
    }
    else if (n == 572) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     10180.7400;
    }
    else if (n == 573) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =       -13.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     10624.7400;
    }
    else if (n == 574) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -9.0000;
      nut_data[ 8] =        13.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -10648.7800;
    }
    else if (n == 575) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     10734.7000;
    }
    else if (n == 576) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     10746.9500;
    }
    else if (n == 577) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     10746.9400;
    }
    else if (n == 578) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     10759.2300;
    }
    else if (n == 579) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -10771.5400;
    }
    else if (n == 580) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =    -10771.5300;
    }
    else if (n == 581) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -10857.0100;
    }
    else if (n == 582) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -11119.3300;
    }
    else if (n == 583) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         4.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -11143.6100;
    }
    else if (n == 584) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     11930.1400;
    }
    else if (n == 585) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     11945.2800;
    }
    else if (n == 586) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     11960.4600;
    }
    else if (n == 587) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =       -10.0000;
      nut_data[ 8] =        15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     12565.5000;
    }
    else if (n == 588) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =       -15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -12599.1400;
    }
    else if (n == 589) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        20.0000;
      nut_data[ 7] =       -21.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     12715.0600;
    }
    else if (n == 590) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =       -20.0000;
      nut_data[ 7] =        20.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -12732.2600;
    }
    else if (n == 591) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     12732.5100;
    }
    else if (n == 592) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -13015.7000;
    }
    else if (n == 593) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -13286.2500;
    }
    else if (n == 594) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -13562.8900;
    }
    else if (n == 595) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -13582.4600;
    }
    else if (n == 596) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     13611.1000;
    }
    else if (n == 597) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -13630.8100;
    }
    else if (n == 598) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -13922.1400;
    }
    else if (n == 599) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -5.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -14158.1100;
    }
    else if (n == 600) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     14573.0900;
    }
    else if (n == 601) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =       -15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =       -29.2000;
    }
    else if (n == 602) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =       -15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     14765.9700;
    }
    else if (n == 603) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -9.0000;
      nut_data[ 8] =        15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -14789.1700;
    }
    else if (n == 604) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =        15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =    -14789.1600;
    }
    else if (n == 605) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =        15.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =    -14812.4200;
    }
    else if (n == 606) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         2.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     15294.4000;
    }
    else if (n == 607) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         2.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     15319.2800;
    }
    else if (n == 608) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         2.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     15319.2700;
    }
    else if (n == 609) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -15751.7600;
    }
    else if (n == 610) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -15778.1700;
    }
    else if (n == 611) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -3.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -6.0000;
      nut_data[ 7] =         7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     17462.2000;
    }
    else if (n == 612) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         6.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -17494.6500;
    }
    else if (n == 613) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     18150.9200;
    }
    else if (n == 614) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         5.0000;
      nut_data[ 7] =        -6.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     18185.9900;
    }
    else if (n == 615) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -18430.9600;
    }
    else if (n == 616) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -18467.1200;
    }
    else if (n == 617) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -18503.4200;
    }
    else if (n == 618) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -2.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         9.0000;
      nut_data[ 8] =       -13.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -18801.8200;
    }
    else if (n == 619) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =        -6.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -20462.8200;
    }
    else if (n == 620) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =        -3.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -20870.1400;
    }
    else if (n == 621) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        17.0000;
      nut_data[ 7] =       -16.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -22201.1800;
    }
    else if (n == 622) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     22260.8600;
    }
    else if (n == 623) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =         2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -22313.6300;
    }
    else if (n == 624) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -3.0000;
      nut_data[ 7] =         5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     22862.0700;
    }
    else if (n == 625) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -9.0000;
      nut_data[ 8] =        17.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -24197.4800;
    }
    else if (n == 626) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     25707.3300;
    }
    else if (n == 627) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     25848.5200;
    }
    else if (n == 628) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -26572.5100;
    }
    else if (n == 629) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         2.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =     29904.1900;
    }
    else if (n == 630) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         2.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     29999.4400;
    }
    else if (n == 631) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         2.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     29999.3800;
    }
    else if (n == 632) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         1.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     30688.4800;
    }
    else if (n == 633) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =        -1.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -30788.8600;
    }
    else if (n == 634) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         9.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     32558.5700;
    }
    else if (n == 635) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =        -2.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -34711.7400;
    }
    else if (n == 636) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -34989.3100;
    }
    else if (n == 637) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -35119.8600;
    }
    else if (n == 638) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -2.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     37731.7400;
    }
    else if (n == 639) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     37883.6000;
    }
    else if (n == 640) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     38036.6800;
    }
    else if (n == 641) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =         3.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -38896.8900;
    }
    else if (n == 642) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -42727.3300;
    }
    else if (n == 643) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     52048.3800;
    }
    else if (n == 644) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -74517.0800;
    }
    else if (n == 645) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     76428.1800;
    }
    else if (n == 646) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         8.0000;
      nut_data[ 7] =       -13.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =     86464.2700;
    }
    else if (n == 647) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         8.0000;
      nut_data[ 7] =       -14.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     86463.7900;
    }
    else if (n == 648) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         8.0000;
      nut_data[ 7] =       -13.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     87265.4000;
    }
    else if (n == 649) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        12.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -88082.0100;
    }
    else if (n == 650) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        13.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =    -88081.5200;
    }
    else if (n == 651) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         2.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        11.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -88914.0600;
    }
    else if (n == 652) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -8.0000;
      nut_data[ 7] =        13.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =    -88913.0400;
    }
    else if (n == 653) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =        -2.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    -88997.8100;
    }
    else if (n == 654) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        18.0000;
      nut_data[ 7] =       -17.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     98682.3000;
    }
    else if (n == 655) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         2.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     99100.9400;
    }
    else if (n == 656) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        18.0000;
      nut_data[ 7] =       -16.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =     99727.8500;
    }
    else if (n == 657) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    100155.4200;
    }
    else if (n == 658) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =         1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =        -1.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -107127.7200;
    }
    else if (n == 659) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         3.0000;
      nut_data[ 7] =        -7.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    110464.0100;
    }
    else if (n == 660) {
      nut_data[ 0] =        -2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         2.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         2.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    145247.8300;
    }
    else if (n == 661) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =        10.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    158588.3100;
    }
    else if (n == 662) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    193201.8000;
    }
    else if (n == 663) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         2.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -199810.8000;
    }
    else if (n == 664) {
      nut_data[ 0] =         1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -2.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        19.0000;
      nut_data[ 7] =       -21.0000;
      nut_data[ 8] =         3.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -209413.5500;
    }
    else if (n == 665) {
      nut_data[ 0] =         2.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -7.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -237473.7100;
    }
    else if (n == 666) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =    301927.9600;
    }
    else if (n == 667) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =    311927.5200;
    }
    else if (n == 668) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    311921.2600;
    }
    else if (n == 669) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -322612.1200;
    }
    else if (n == 670) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -1.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -334061.8400;
    }
    else if (n == 671) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =   -334054.6500;
    }
    else if (n == 672) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         2.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =   -346338.7300;
    }
    else if (n == 673) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        -5.0000;
      nut_data[ 7] =         6.0000;
      nut_data[ 8] =         4.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =   -455742.4000;
    }
    else if (n == 674) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         1.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =    562786.5000;
    }
    else if (n == 675) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =    609254.2200;
    }
    else if (n == 676) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -5.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =    609230.3300;
    }
    else if (n == 677) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -651391.3000;
    }
    else if (n == 678) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =        -2.0000;
      nut_data[10] =         6.0000;
      nut_data[11] =        -3.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =   -666209.9000;
    }
    else if (n == 679) {
      nut_data[ 0] =        -1.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =        10.0000;
      nut_data[ 7] =        -3.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -698321.3200;
    }
    else if (n == 680) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         1.0000;
      nut_data[ 3] =        -1.0000;
      nut_data[ 4] =         1.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         3.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =   -699821.5100;
    }
    else if (n == 681) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         4.0000;
      nut_data[ 8] =        -8.0000;
      nut_data[ 9] =         3.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         1.0000;
      nut_data[14] =   -699789.9800;
    }
    else if (n == 682) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -4.0000;
      nut_data[ 8] =         8.0000;
      nut_data[ 9] =        -1.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =   -739551.3200;
    }
    else if (n == 683) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         0.0000;
      nut_data[ 8] =         0.0000;
      nut_data[ 9] =         0.0000;
      nut_data[10] =         0.0000;
      nut_data[11] =        -1.0000;
      nut_data[12] =         2.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =   1170120.0400;
    }
    else if (n == 684) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =       -16.0000;
      nut_data[ 9] =         4.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =   4137408.1100;
    }
    else if (n == 685) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =        -8.0000;
      nut_data[ 8] =        16.0000;
      nut_data[ 9] =        -4.0000;
      nut_data[10] =        -5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         2.0000;
      nut_data[14] =   5464350.6700;
    }
    else if (n == 686) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
      nut_data[ 5] =         0.0000;
      nut_data[ 6] =         0.0000;
      nut_data[ 7] =         8.0000;
      nut_data[ 8] =       -16.0000;
      nut_data[ 9] =         4.0000;
      nut_data[10] =         5.0000;
      nut_data[11] =         0.0000;
      nut_data[12] =         0.0000;
      nut_data[13] =         0.0000;
      nut_data[14] =  34075700.8200;
    }
  }

  else if (*S1 == 'P' && *S2 == 'B') {
    if (n ==   0) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==   1) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==   2) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==   3) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0008;
    }
    else if (n ==   4) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==   5) {
      nut_data[ 0] =        -0.0024;
      nut_data[ 1] =        -0.0012;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0010;
      nut_data[ 4] =         0.0015;
    }
    else if (n ==   6) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==   7) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==   8) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==   9) {
      nut_data[ 0] =         0.0024;
      nut_data[ 1] =        -0.0012;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =        -0.0011;
      nut_data[ 4] =         0.0016;
    }
    else if (n ==  10) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==  11) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  12) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  13) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  14) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  15) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  16) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  17) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  18) {
      nut_data[ 0] =        -0.0021;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0011;
      nut_data[ 4] =         0.0016;
    }
    else if (n ==  19) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  20) {
      nut_data[ 0] =         0.0021;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =        -0.0011;
      nut_data[ 4] =         0.0016;
    }
    else if (n ==  21) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  22) {
      nut_data[ 0] =        -0.0126;
      nut_data[ 1] =        -0.0063;
      nut_data[ 2] =        -0.0027;
      nut_data[ 3] =         0.0055;
      nut_data[ 4] =         0.0083;
    }
    else if (n ==  23) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  24) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  25) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==  26) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0028;
      nut_data[ 2] =         0.0015;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0019;
    }
    else if (n ==  27) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==  28) {
      nut_data[ 0] =         0.0126;
      nut_data[ 1] =        -0.0063;
      nut_data[ 2] =        -0.0027;
      nut_data[ 3] =        -0.0055;
      nut_data[ 4] =         0.0083;
    }
    else if (n ==  29) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  30) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n ==  31) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  32) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  33) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  34) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  35) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  36) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n ==  37) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n ==  38) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  39) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  40) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  41) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  42) {
      nut_data[ 0] =        -0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
    }
    else if (n ==  43) {
      nut_data[ 0] =        -0.0034;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0014;
    }
    else if (n ==  44) {
      nut_data[ 0] =         0.0020;
      nut_data[ 1] =         0.0010;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n ==  45) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  46) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  47) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  48) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==  49) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  50) {
      nut_data[ 0] =         0.0021;
      nut_data[ 1] =         0.0011;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n ==  51) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  52) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n ==  53) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  54) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  55) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  56) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n ==  57) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  58) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  59) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  60) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  61) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  62) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  63) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==  64) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  65) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n ==  66) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  67) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==  68) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  69) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  70) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  71) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==  72) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n ==  73) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  74) {
      nut_data[ 0] =         0.0002;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0006;
    }
    else if (n ==  75) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  76) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  77) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  78) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  79) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n ==  80) {
      nut_data[ 0] =         0.0014;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0006;
    }
    else if (n ==  81) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n ==  82) {
      nut_data[ 0] =         0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
    }
    else if (n ==  83) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0007;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  84) {
      nut_data[ 0] =         0.0017;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n ==  85) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0007;
    }
    else if (n ==  86) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  87) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  88) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n ==  89) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  90) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  91) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n ==  92) {
      nut_data[ 0] =         0.0028;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n ==  93) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0006;
    }
    else if (n ==  94) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  95) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0010;
      nut_data[ 4] =         0.0010;
    }
    else if (n ==  96) {
      nut_data[ 0] =         0.0019;
      nut_data[ 1] =        -0.0023;
      nut_data[ 2] =        -0.0010;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0016;
    }
    else if (n ==  97) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  98) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n ==  99) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 100) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0008;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 101) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 102) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =        -0.0032;
      nut_data[ 2] =        -0.0014;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0019;
    }
    else if (n == 103) {
      nut_data[ 0] =         0.0031;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0013;
      nut_data[ 4] =         0.0018;
    }
    else if (n == 104) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0010;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 105) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 106) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 107) {
      nut_data[ 0] =        -0.0184;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0080;
      nut_data[ 4] =         0.0109;
    }
    else if (n == 108) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 109) {
      nut_data[ 0] =         0.0040;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 110) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0013;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 111) {
      nut_data[ 0] =        -0.0037;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0016;
      nut_data[ 4] =         0.0022;
    }
    else if (n == 112) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 113) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 114) {
      nut_data[ 0] =         0.0034;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0015;
      nut_data[ 4] =         0.0020;
    }
    else if (n == 115) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0017;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 116) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 117) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 118) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0006;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 119) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 120) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 121) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 122) {
      nut_data[ 0] =         0.0370;
      nut_data[ 1] =        -0.0008;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0160;
      nut_data[ 4] =         0.0218;
    }
    else if (n == 123) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =        -0.0015;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 124) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 125) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 126) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 127) {
      nut_data[ 0] =         0.0059;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0024;
    }
    else if (n == 128) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0017;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 129) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 130) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 131) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 132) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 133) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0007;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 134) {
      nut_data[ 0] =         0.0050;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0022;
      nut_data[ 4] =         0.0030;
    }
    else if (n == 135) {
      nut_data[ 0] =        -0.0022;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0010;
      nut_data[ 4] =         0.0015;
    }
    else if (n == 136) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 137) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 138) {
      nut_data[ 0] =        -0.0025;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0011;
      nut_data[ 4] =         0.0015;
    }
    else if (n == 139) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 140) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 141) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 142) {
      nut_data[ 0] =         0.0029;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0013;
      nut_data[ 4] =         0.0017;
    }
    else if (n == 143) {
      nut_data[ 0] =         0.0143;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =        -0.0062;
      nut_data[ 4] =         0.0084;
    }
    else if (n == 144) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 145) {
      nut_data[ 0] =         0.0517;
      nut_data[ 1] =         0.0016;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =        -0.0224;
      nut_data[ 4] =         0.0305;
    }
    else if (n == 146) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0023;
      nut_data[ 4] =         0.0023;
    }
    else if (n == 147) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0114;
      nut_data[ 2] =        -0.0050;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0068;
    }
    else if (n == 148) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
    }
    else if (n == 149) {
      nut_data[ 0] =         0.0030;
      nut_data[ 1] =        -0.0018;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =        -0.0013;
      nut_data[ 4] =         0.0021;
    }
    else if (n == 150) {
      nut_data[ 0] =         0.0067;
      nut_data[ 1] =        -0.0091;
      nut_data[ 2] =        -0.0039;
      nut_data[ 3] =        -0.0029;
      nut_data[ 4] =         0.0066;
    }
    else if (n == 151) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =        -0.0012;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 152) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =        -0.0009;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 153) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 154) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 155) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 156) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 157) {
      nut_data[ 0] =         0.0018;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 158) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 159) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 160) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 161) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0010;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 162) {
      nut_data[ 0] =        -0.0339;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0147;
      nut_data[ 4] =         0.0200;
    }
    else if (n == 163) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =        -0.0023;
      nut_data[ 2] =        -0.0012;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 164) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 165) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 166) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 167) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 168) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0010;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 169) {
      nut_data[ 0] =         0.0083;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0033;
    }
    else if (n == 170) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0024;
      nut_data[ 2] =         0.0013;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0017;
    }
    else if (n == 171) {
      nut_data[ 0] =         0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0008;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 172) {
      nut_data[ 0] =         0.0026;
      nut_data[ 1] =        -0.0014;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =        -0.0011;
      nut_data[ 4] =         0.0017;
    }
    else if (n == 173) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 174) {
      nut_data[ 0] =         0.0074;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0032;
      nut_data[ 4] =         0.0044;
    }
    else if (n == 175) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 176) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 177) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =        -0.0012;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 178) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 179) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =        -0.0173;
      nut_data[ 2] =        -0.0075;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0102;
    }
    else if (n == 180) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0019;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 181) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 182) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 183) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 184) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0010;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 185) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 186) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 187) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 188) {
      nut_data[ 0] =        -0.0021;
      nut_data[ 1] =        -0.0032;
      nut_data[ 2] =        -0.0014;
      nut_data[ 3] =         0.0009;
      nut_data[ 4] =         0.0023;
    }
    else if (n == 189) {
      nut_data[ 0] =        -0.0063;
      nut_data[ 1] =        -0.0016;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =         0.0028;
      nut_data[ 4] =         0.0039;
    }
    else if (n == 190) {
      nut_data[ 0] =         0.0163;
      nut_data[ 1] =        -0.0012;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =        -0.0072;
      nut_data[ 4] =         0.0097;
    }
    else if (n == 191) {
      nut_data[ 0] =        -0.0085;
      nut_data[ 1] =        -0.0070;
      nut_data[ 2] =        -0.0031;
      nut_data[ 3] =         0.0037;
      nut_data[ 4] =         0.0065;
    }
    else if (n == 192) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 193) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 194) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 195) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 196) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0019;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 197) {
      nut_data[ 0] =        -0.0011;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 198) {
      nut_data[ 0] =        -0.0062;
      nut_data[ 1] =        -0.0097;
      nut_data[ 2] =        -0.0042;
      nut_data[ 3] =         0.0027;
      nut_data[ 4] =         0.0068;
    }
    else if (n == 199) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 200) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0015;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 201) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0006;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 202) {
      nut_data[ 0] =         0.0012;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 203) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 204) {
      nut_data[ 0] =        -0.0123;
      nut_data[ 1] =        -0.0416;
      nut_data[ 2] =        -0.0180;
      nut_data[ 3] =         0.0053;
      nut_data[ 4] =         0.0256;
    }
    else if (n == 205) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0019;
      nut_data[ 3] =         0.0006;
      nut_data[ 4] =         0.0020;
    }
    else if (n == 206) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0086;
      nut_data[ 2] =        -0.0019;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0040;
    }
    else if (n == 207) {
      nut_data[ 0] =        -0.0089;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0038;
      nut_data[ 4] =         0.0052;
    }
    else if (n == 208) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 209) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0009;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 210) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =        -0.0032;
      nut_data[ 2] =        -0.0017;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0022;
    }
    else if (n == 211) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 212) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0035;
      nut_data[ 4] =         0.0035;
    }
    else if (n == 213) {
      nut_data[ 0] =         0.0123;
      nut_data[ 1] =        -0.0415;
      nut_data[ 2] =        -0.0180;
      nut_data[ 3] =        -0.0053;
      nut_data[ 4] =         0.0255;
    }
    else if (n == 214) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 215) {
      nut_data[ 0] =        -0.0011;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 216) {
      nut_data[ 0] =         0.0014;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 217) {
      nut_data[ 0] =         0.0061;
      nut_data[ 1] =        -0.0096;
      nut_data[ 2] =        -0.0042;
      nut_data[ 3] =        -0.0027;
      nut_data[ 4] =         0.0068;
    }
    else if (n == 218) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 219) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0020;
      nut_data[ 2] =        -0.0009;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 220) {
      nut_data[ 0] =         0.0080;
      nut_data[ 1] =        -0.0071;
      nut_data[ 2] =        -0.0031;
      nut_data[ 3] =        -0.0035;
      nut_data[ 4] =         0.0063;
    }
    else if (n == 221) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 222) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 223) {
      nut_data[ 0] =         0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 224) {
      nut_data[ 0] =        -0.0154;
      nut_data[ 1] =        -0.0030;
      nut_data[ 2] =        -0.0013;
      nut_data[ 3] =         0.0067;
      nut_data[ 4] =         0.0093;
    }
    else if (n == 225) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0035;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0014;
    }
    else if (n == 226) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 227) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 228) {
      nut_data[ 0] =         0.0054;
      nut_data[ 1] =        -0.0015;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =        -0.0024;
      nut_data[ 4] =         0.0034;
    }
    else if (n == 229) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 230) {
      nut_data[ 0] =         0.0089;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0036;
    }
    else if (n == 231) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0034;
      nut_data[ 2] =         0.0018;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0023;
    }
    else if (n == 232) {
      nut_data[ 0] =         0.0018;
      nut_data[ 1] =        -0.0029;
      nut_data[ 2] =        -0.0013;
      nut_data[ 3] =        -0.0008;
      nut_data[ 4] =         0.0020;
    }
    else if (n == 233) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 234) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 235) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 236) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0021;
      nut_data[ 2] =         0.0011;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0014;
    }
    else if (n == 237) {
      nut_data[ 0] =         0.0393;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0157;
    }
    else if (n == 238) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 239) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 240) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 241) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 242) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 243) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0019;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0013;
    }
    else if (n == 244) {
      nut_data[ 0] =         0.0117;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0051;
      nut_data[ 4] =         0.0069;
    }
    else if (n == 245) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0114;
      nut_data[ 2] =        -0.0049;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0067;
    }
    else if (n == 246) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 247) {
      nut_data[ 0] =         0.0083;
      nut_data[ 1] =         0.0015;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0034;
    }
    else if (n == 248) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 249) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 250) {
      nut_data[ 0] =         0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0007;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 251) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =        -0.0009;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 252) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 253) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0006;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 254) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =        -0.0028;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 255) {
      nut_data[ 0] =        -0.0018;
      nut_data[ 1] =        -0.0010;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0008;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 256) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 257) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0006;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 258) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 259) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 260) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 261) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 262) {
      nut_data[ 0] =         0.0027;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 263) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 264) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0024;
      nut_data[ 2] =        -0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0014;
    }
    else if (n == 265) {
      nut_data[ 0] =         0.0113;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0049;
      nut_data[ 4] =         0.0067;
    }
    else if (n == 266) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 267) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =        -0.0031;
      nut_data[ 2] =        -0.0016;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0021;
    }
    else if (n == 268) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 269) {
      nut_data[ 0] =         0.0016;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 270) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 271) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 272) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 273) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 274) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =        -0.0013;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 275) {
      nut_data[ 0] =        -0.0598;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0239;
    }
    else if (n == 276) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 277) {
      nut_data[ 0] =         0.0012;
      nut_data[ 1] =         0.0055;
      nut_data[ 2] =         0.0029;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0037;
    }
    else if (n == 278) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 279) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 280) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 281) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 282) {
      nut_data[ 0] =        -0.0035;
      nut_data[ 1] =        -0.0048;
      nut_data[ 2] =        -0.0021;
      nut_data[ 3] =         0.0015;
      nut_data[ 4] =         0.0035;
    }
    else if (n == 283) {
      nut_data[ 0] =         0.0001;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
    }
    else if (n == 284) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 285) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 286) {
      nut_data[ 0] =         0.0016;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 287) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 288) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0035;
      nut_data[ 2] =         0.0019;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0024;
    }
    else if (n == 289) {
      nut_data[ 0] =         0.0202;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0087;
      nut_data[ 4] =         0.0119;
    }
    else if (n == 290) {
      nut_data[ 0] =        -0.0019;
      nut_data[ 1] =        -0.0008;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0008;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 291) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0027;
      nut_data[ 2] =        -0.0012;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 292) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 293) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 294) {
      nut_data[ 0] =        -0.0262;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0114;
      nut_data[ 4] =         0.0155;
    }
    else if (n == 295) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 296) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0011;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 297) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 298) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 299) {
      nut_data[ 0] =        -0.0074;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0032;
      nut_data[ 4] =         0.0044;
    }
    else if (n == 300) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 301) {
      nut_data[ 0] =        -0.0014;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0006;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 302) {
      nut_data[ 0] =        -0.0001;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0000;
    }
    else if (n == 303) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 304) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 305) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 306) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 307) {
      nut_data[ 0] =        -0.0019;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 308) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =        -0.0027;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 309) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 310) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 311) {
      nut_data[ 0] =         0.0021;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 312) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =        -0.0030;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0013;
    }
    else if (n == 313) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 314) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 315) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 316) {
      nut_data[ 0] =        -0.0075;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0030;
    }
    else if (n == 317) {
      nut_data[ 0] =        -0.0368;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0147;
    }
    else if (n == 318) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0020;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 319) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 320) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 321) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0007;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 322) {
      nut_data[ 0] =        -0.1223;
      nut_data[ 1] =        -0.0026;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0489;
    }
    else if (n == 323) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 324) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 325) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 326) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 327) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 328) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 329) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0328;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0131;
    }
    else if (n == 330) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 331) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 332) {
      nut_data[ 0] =        -0.0078;
      nut_data[ 1] =         0.0045;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0036;
    }
    else if (n == 333) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 334) {
      nut_data[ 0] =         0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0008;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 335) {
      nut_data[ 0] =        -0.0166;
      nut_data[ 1] =         0.0269;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0126;
    }
    else if (n == 336) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 337) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 338) {
      nut_data[ 0] =        -0.0016;
      nut_data[ 1] =         0.0023;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 339) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 340) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 341) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 342) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 343) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 344) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0045;
      nut_data[ 2] =        -0.0020;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0027;
    }
    else if (n == 345) {
      nut_data[ 0] =        -0.0458;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0198;
      nut_data[ 4] =         0.0270;
    }
    else if (n == 346) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 347) {
      nut_data[ 0] =         0.0014;
      nut_data[ 1] =        -0.0059;
      nut_data[ 2] =        -0.0031;
      nut_data[ 3] =        -0.0008;
      nut_data[ 4] =         0.0040;
    }
    else if (n == 348) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 349) {
      nut_data[ 0] =        -0.0028;
      nut_data[ 1] =         0.0036;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0018;
    }
    else if (n == 350) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 351) {
      nut_data[ 0] =         0.0118;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0052;
      nut_data[ 4] =         0.0070;
    }
    else if (n == 352) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0011;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 353) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 354) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 355) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 356) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 357) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 358) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 359) {
      nut_data[ 0] =        -0.0019;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0010;
      nut_data[ 4] =         0.0013;
    }
    else if (n == 360) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 361) {
      nut_data[ 0] =         0.0030;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =        -0.0013;
      nut_data[ 4] =         0.0018;
    }
    else if (n == 362) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 363) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 364) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =        -0.0032;
      nut_data[ 2] =        -0.0017;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0022;
    }
    else if (n == 365) {
      nut_data[ 0] =         0.1485;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0594;
    }
    else if (n == 366) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0021;
      nut_data[ 2] =         0.0011;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0014;
    }
    else if (n == 367) {
      nut_data[ 0] =         0.0025;
      nut_data[ 1] =         0.0106;
      nut_data[ 2] =         0.0057;
      nut_data[ 3] =        -0.0013;
      nut_data[ 4] =         0.0073;
    }
    else if (n == 368) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 369) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 370) {
      nut_data[ 0] =        -0.0011;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 371) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 372) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 373) {
      nut_data[ 0] =        -0.0028;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0015;
      nut_data[ 4] =         0.0019;
    }
    else if (n == 374) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 375) {
      nut_data[ 0] =         0.0002;
      nut_data[ 1] =         0.0001;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 376) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 377) {
      nut_data[ 0] =        -0.0046;
      nut_data[ 1] =         0.0014;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0019;
    }
    else if (n == 378) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0013;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 379) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0028;
      nut_data[ 2] =         0.0015;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0019;
    }
    else if (n == 380) {
      nut_data[ 0] =        -0.0022;
      nut_data[ 1] =         0.0093;
      nut_data[ 2] =         0.0049;
      nut_data[ 3] =         0.0012;
      nut_data[ 4] =         0.0063;
    }
    else if (n == 381) {
      nut_data[ 0] =         0.0490;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0213;
      nut_data[ 4] =         0.0289;
    }
    else if (n == 382) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0066;
      nut_data[ 2] =         0.0029;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0039;
    }
    else if (n == 383) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 384) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0006;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 385) {
      nut_data[ 0] =        -0.0025;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 386) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 387) {
      nut_data[ 0] =         0.0027;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0014;
      nut_data[ 4] =         0.0018;
    }
    else if (n == 388) {
      nut_data[ 0] =        -0.0068;
      nut_data[ 1] =         0.0039;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0031;
    }
    else if (n == 389) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 390) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 391) {
      nut_data[ 0] =        -0.0016;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0007;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 392) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 393) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 394) {
      nut_data[ 0] =         0.0022;
      nut_data[ 1] =        -0.0087;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0036;
    }
    else if (n == 395) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 396) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 397) {
      nut_data[ 0] =        -0.0020;
      nut_data[ 1] =         0.0034;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 398) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 399) {
      nut_data[ 0] =         0.0017;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0007;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 400) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 401) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 402) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 403) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 404) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0019;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0013;
    }
    else if (n == 405) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0022;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 406) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 407) {
      nut_data[ 0] =         0.0016;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =        -0.0007;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 408) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 409) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =        -0.0010;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 410) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 411) {
      nut_data[ 0] =        -0.0127;
      nut_data[ 1] =         0.0021;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0055;
      nut_data[ 4] =         0.0076;
    }
    else if (n == 412) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0006;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 413) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =        -0.0009;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 414) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0017;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 415) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 416) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0614;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0246;
    }
    else if (n == 417) {
      nut_data[ 0] =         0.0018;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0009;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 418) {
      nut_data[ 0] =         0.0030;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =        -0.0016;
      nut_data[ 4] =         0.0020;
    }
    else if (n == 419) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 420) {
      nut_data[ 0] =         0.0024;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =        -0.0011;
      nut_data[ 4] =         0.0015;
    }
    else if (n == 421) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 422) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 423) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 424) {
      nut_data[ 0] =         0.0050;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0027;
      nut_data[ 4] =         0.0034;
    }
    else if (n == 425) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0016;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 426) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0013;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 427) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 428) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 429) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0032;
      nut_data[ 2] =        -0.0017;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0021;
    }
    else if (n == 430) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 431) {
      nut_data[ 0] =        -0.0024;
      nut_data[ 1] =        -0.0013;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0010;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 432) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 433) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 434) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 435) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0031;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 436) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 437) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 438) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 439) {
      nut_data[ 0] =         0.0117;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0063;
      nut_data[ 4] =         0.0078;
    }
    else if (n == 440) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 441) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 442) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 443) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 444) {
      nut_data[ 0] =        -0.1166;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0505;
      nut_data[ 4] =         0.0687;
    }
    else if (n == 445) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0049;
      nut_data[ 2] =         0.0026;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0033;
    }
    else if (n == 446) {
      nut_data[ 0] =        -0.0027;
      nut_data[ 1] =        -0.0143;
      nut_data[ 2] =        -0.0077;
      nut_data[ 3] =         0.0014;
      nut_data[ 4] =         0.0098;
    }
    else if (n == 447) {
      nut_data[ 0] =         0.0042;
      nut_data[ 1] =         0.0223;
      nut_data[ 2] =         0.0119;
      nut_data[ 3] =        -0.0022;
      nut_data[ 4] =         0.0151;
    }
    else if (n == 448) {
      nut_data[ 0] =        -0.0025;
      nut_data[ 1] =         0.0022;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0013;
    }
    else if (n == 449) {
      nut_data[ 0] =         0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 450) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 451) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 452) {
      nut_data[ 0] =         0.0016;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0009;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 453) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 454) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0010;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 455) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =         0.0006;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 456) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0011;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 457) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 458) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 459) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0015;
      nut_data[ 2] =         0.0008;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 460) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 461) {
      nut_data[ 0] =        -0.0014;
      nut_data[ 1] =        -0.0039;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0017;
    }
    else if (n == 462) {
      nut_data[ 0] =         0.0057;
      nut_data[ 1] =         0.0011;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =        -0.0030;
      nut_data[ 4] =         0.0038;
    }
    else if (n == 463) {
      nut_data[ 0] =         0.0140;
      nut_data[ 1] =         0.0027;
      nut_data[ 2] =         0.0014;
      nut_data[ 3] =        -0.0075;
      nut_data[ 4] =         0.0095;
    }
    else if (n == 464) {
      nut_data[ 0] =         0.0031;
      nut_data[ 1] =         0.0006;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =        -0.0017;
      nut_data[ 4] =         0.0021;
    }
    else if (n == 465) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0007;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 466) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 467) {
      nut_data[ 0] =        -0.0011;
      nut_data[ 1] =        -0.0268;
      nut_data[ 2] =        -0.0116;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0158;
    }
    else if (n == 468) {
      nut_data[ 0] =        -0.0051;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0022;
      nut_data[ 4] =         0.0030;
    }
    else if (n == 469) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 470) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =        -0.0013;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 471) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0005;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 472) {
      nut_data[ 0] =        -0.0086;
      nut_data[ 1] =         0.0153;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0070;
    }
    else if (n == 473) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 474) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 475) {
      nut_data[ 0] =         0.0085;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0037;
      nut_data[ 4] =         0.0050;
    }
    else if (n == 476) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 477) {
      nut_data[ 0] =        -0.2150;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0932;
      nut_data[ 4] =         0.1268;
    }
    else if (n == 478) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =        -0.0049;
      nut_data[ 2] =        -0.0026;
      nut_data[ 3] =        -0.0007;
      nut_data[ 4] =         0.0034;
    }
    else if (n == 479) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0040;
      nut_data[ 2] =         0.0021;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0027;
    }
    else if (n == 480) {
      nut_data[ 0] =        -0.0145;
      nut_data[ 1] =         0.0047;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0061;
    }
    else if (n == 481) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =        -0.0009;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 482) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =        -0.0014;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 483) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 484) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0009;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 485) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 486) {
      nut_data[ 0] =        -0.0037;
      nut_data[ 1] =         0.0035;
      nut_data[ 2] =         0.0019;
      nut_data[ 3] =         0.0020;
      nut_data[ 4] =         0.0034;
    }
    else if (n == 487) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 488) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0010;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 489) {
      nut_data[ 0] =         0.0054;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0029;
      nut_data[ 4] =         0.0036;
    }
    else if (n == 490) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 491) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 492) {
      nut_data[ 0] =        -0.0138;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0055;
    }
    else if (n == 493) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0015;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 494) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =        -0.0004;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 495) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 496) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 497) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 498) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 499) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 500) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 501) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 502) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 503) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 504) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 505) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 506) {
      nut_data[ 0] =         0.0052;
      nut_data[ 1] =         0.0023;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =        -0.0023;
      nut_data[ 4] =         0.0034;
    }
    else if (n == 507) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 508) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =        -0.0047;
      nut_data[ 2] =        -0.0025;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0032;
    }
    else if (n == 509) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0049;
      nut_data[ 2] =         0.0026;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0033;
    }
    else if (n == 510) {
      nut_data[ 0] =        -0.0091;
      nut_data[ 1] =         0.0248;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0106;
    }
    else if (n == 511) {
      nut_data[ 0] =        -0.0013;
      nut_data[ 1] =         0.0052;
      nut_data[ 2] =         0.0028;
      nut_data[ 3] =         0.0007;
      nut_data[ 4] =         0.0036;
    }
    else if (n == 512) {
      nut_data[ 0] =        -0.0050;
      nut_data[ 1] =         0.0194;
      nut_data[ 2] =         0.0103;
      nut_data[ 3] =         0.0027;
      nut_data[ 4] =         0.0133;
    }
    else if (n == 513) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 514) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 515) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 516) {
      nut_data[ 0] =        -0.0053;
      nut_data[ 1] =        -0.0009;
      nut_data[ 2] =        -0.0005;
      nut_data[ 3] =         0.0028;
      nut_data[ 4] =         0.0036;
    }
    else if (n == 517) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 518) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 519) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 520) {
      nut_data[ 0] =         0.0035;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0014;
    }
    else if (n == 521) {
      nut_data[ 0] =        -0.0018;
      nut_data[ 1] =        -0.0436;
      nut_data[ 2] =        -0.0233;
      nut_data[ 3] =         0.0009;
      nut_data[ 4] =         0.0291;
    }
    else if (n == 522) {
      nut_data[ 0] =        -0.0011;
      nut_data[ 1] =        -0.0021;
      nut_data[ 2] =        -0.0011;
      nut_data[ 3] =         0.0006;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 523) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 524) {
      nut_data[ 0] =        -0.0021;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0011;
      nut_data[ 4] =         0.0014;
    }
    else if (n == 525) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 526) {
      nut_data[ 0] =        -0.0133;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0057;
      nut_data[ 4] =         0.0078;
    }
    else if (n == 527) {
      nut_data[ 0] =         0.0051;
      nut_data[ 1] =         0.0114;
      nut_data[ 2] =         0.0061;
      nut_data[ 3] =        -0.0027;
      nut_data[ 4] =         0.0083;
    }
    else if (n == 528) {
      nut_data[ 0] =        -0.0048;
      nut_data[ 1] =        -0.0110;
      nut_data[ 2] =        -0.0059;
      nut_data[ 3] =         0.0026;
      nut_data[ 4] =         0.0080;
    }
    else if (n == 529) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 530) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 531) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0023;
      nut_data[ 2] =         0.0013;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 532) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0047;
      nut_data[ 2] =         0.0025;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0032;
    }
    else if (n == 533) {
      nut_data[ 0] =        -0.0008;
      nut_data[ 1] =        -0.0047;
      nut_data[ 2] =        -0.0025;
      nut_data[ 3] =         0.0004;
      nut_data[ 4] =         0.0032;
    }
    else if (n == 534) {
      nut_data[ 0] =        -0.0449;
      nut_data[ 1] =         0.0430;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0249;
    }
    else if (n == 535) {
      nut_data[ 0] =         0.0273;
      nut_data[ 1] =         0.0080;
      nut_data[ 2] =         0.0043;
      nut_data[ 3] =        -0.0146;
      nut_data[ 4] =         0.0190;
    }
    else if (n == 536) {
      nut_data[ 0] =         0.0023;
      nut_data[ 1] =         0.0007;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =        -0.0013;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 537) {
      nut_data[ 0] =        -0.0040;
      nut_data[ 1] =         0.0057;
      nut_data[ 2] =         0.0030;
      nut_data[ 3] =         0.0021;
      nut_data[ 4] =         0.0046;
    }
    else if (n == 538) {
      nut_data[ 0] =        -0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 539) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 540) {
      nut_data[ 0] =         0.0057;
      nut_data[ 1] =        -0.0028;
      nut_data[ 2] =        -0.0015;
      nut_data[ 3] =        -0.0030;
      nut_data[ 4] =         0.0042;
    }
    else if (n == 541) {
      nut_data[ 0] =        -0.0439;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0176;
    }
    else if (n == 542) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0005;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 543) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 544) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 545) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =        -0.0016;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 546) {
      nut_data[ 0] =        -0.0073;
      nut_data[ 1] =         0.0017;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =         0.0039;
      nut_data[ 4] =         0.0050;
    }
    else if (n == 547) {
      nut_data[ 0] =        -0.0082;
      nut_data[ 1] =         0.0292;
      nut_data[ 2] =         0.0156;
      nut_data[ 3] =         0.0044;
      nut_data[ 4] =         0.0202;
    }
    else if (n == 548) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 549) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 550) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 551) {
      nut_data[ 0] =         0.0084;
      nut_data[ 1] =         0.0298;
      nut_data[ 2] =         0.0159;
      nut_data[ 3] =        -0.0045;
      nut_data[ 4] =         0.0207;
    }
    else if (n == 552) {
      nut_data[ 0] =         0.0076;
      nut_data[ 1] =         0.0017;
      nut_data[ 2] =         0.0009;
      nut_data[ 3] =        -0.0041;
      nut_data[ 4] =         0.0052;
    }
    else if (n == 553) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 554) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 555) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 556) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0014;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 557) {
      nut_data[ 0] =        -0.0068;
      nut_data[ 1] =        -0.0034;
      nut_data[ 2] =        -0.0018;
      nut_data[ 3] =         0.0036;
      nut_data[ 4] =         0.0050;
    }
    else if (n == 558) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 559) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 560) {
      nut_data[ 0] =        -0.0014;
      nut_data[ 1] =         0.0007;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 561) {
      nut_data[ 0] =         0.0046;
      nut_data[ 1] =         0.0066;
      nut_data[ 2] =         0.0035;
      nut_data[ 3] =        -0.0025;
      nut_data[ 4] =         0.0054;
    }
    else if (n == 562) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 563) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 564) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 565) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 566) {
      nut_data[ 0] =         0.0008;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 567) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 568) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0029;
      nut_data[ 2] =         0.0015;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0019;
    }
    else if (n == 569) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 570) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 571) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 572) {
      nut_data[ 0] =        -0.0024;
      nut_data[ 1] =         0.0012;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 573) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 574) {
      nut_data[ 0] =         0.0010;
      nut_data[ 1] =        -0.0022;
      nut_data[ 2] =        -0.0012;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 575) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 576) {
      nut_data[ 0] =         0.0047;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =        -0.0025;
      nut_data[ 4] =         0.0032;
    }
    else if (n == 577) {
      nut_data[ 0] =        -0.0066;
      nut_data[ 1] =        -0.0012;
      nut_data[ 2] =        -0.0006;
      nut_data[ 3] =         0.0035;
      nut_data[ 4] =         0.0045;
    }
    else if (n == 578) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0056;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0023;
    }
    else if (n == 579) {
      nut_data[ 0] =         0.0174;
      nut_data[ 1] =         0.0084;
      nut_data[ 2] =         0.0045;
      nut_data[ 3] =        -0.0093;
      nut_data[ 4] =         0.0129;
    }
    else if (n == 580) {
      nut_data[ 0] =         0.0032;
      nut_data[ 1] =         0.0015;
      nut_data[ 2] =        -0.0008;
      nut_data[ 3] =         0.0017;
      nut_data[ 4] =         0.0024;
    }
    else if (n == 581) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0006;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 582) {
      nut_data[ 0] =        -0.0017;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0009;
      nut_data[ 4] =         0.0012;
    }
    else if (n == 583) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 584) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0006;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 585) {
      nut_data[ 0] =         0.0020;
      nut_data[ 1] =        -0.0070;
      nut_data[ 2] =        -0.0037;
      nut_data[ 3] =        -0.0011;
      nut_data[ 4] =         0.0048;
    }
    else if (n == 586) {
      nut_data[ 0] =        -0.0021;
      nut_data[ 1] =        -0.0078;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0032;
    }
    else if (n == 587) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 588) {
      nut_data[ 0] =         0.0015;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =        -0.0008;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 589) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0008;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 590) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 591) {
      nut_data[ 0] =        -0.0053;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0021;
    }
    else if (n == 592) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 593) {
      nut_data[ 0] =        -0.0015;
      nut_data[ 1] =         0.0022;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 594) {
      nut_data[ 0] =        -0.0349;
      nut_data[ 1] =        -0.0062;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0142;
    }
    else if (n == 595) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0007;
      nut_data[ 2] =         0.0004;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 596) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 597) {
      nut_data[ 0] =         0.0089;
      nut_data[ 1] =        -0.0016;
      nut_data[ 2] =        -0.0009;
      nut_data[ 3] =        -0.0048;
      nut_data[ 4] =         0.0061;
    }
    else if (n == 598) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 599) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 600) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 601) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 602) {
      nut_data[ 0] =         0.0045;
      nut_data[ 1] =        -0.0022;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0020;
    }
    else if (n == 603) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0019;
      nut_data[ 2] =         0.0010;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0013;
    }
    else if (n == 604) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 605) {
      nut_data[ 0] =        -0.0014;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0006;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 606) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 607) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 608) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0008;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 609) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 610) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0003;
      nut_data[ 2] =        -0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 611) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0001;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 612) {
      nut_data[ 0] =         0.0078;
      nut_data[ 1] =        -0.0018;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0032;
    }
    else if (n == 613) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 614) {
      nut_data[ 0] =        -0.0010;
      nut_data[ 1] =         0.0233;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0093;
    }
    else if (n == 615) {
      nut_data[ 0] =        -0.0042;
      nut_data[ 1] =         0.0020;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0019;
    }
    else if (n == 616) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0024;
      nut_data[ 2] =         0.0013;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 617) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 618) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =        -0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 619) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 620) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 621) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 622) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 623) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0008;
      nut_data[ 2] =        -0.0004;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 624) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0001;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 625) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =        -0.0004;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0001;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 626) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =         0.0024;
      nut_data[ 2] =         0.0011;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0016;
    }
    else if (n == 627) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 628) {
      nut_data[ 0] =        -0.0016;
      nut_data[ 1] =         0.0008;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0007;
    }
    else if (n == 629) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 630) {
      nut_data[ 0] =        -0.0009;
      nut_data[ 1] =        -0.0011;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 631) {
      nut_data[ 0] =        -0.0017;
      nut_data[ 1] =        -0.0019;
      nut_data[ 2] =        -0.0010;
      nut_data[ 3] =         0.0009;
      nut_data[ 4] =         0.0017;
    }
    else if (n == 632) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 633) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0004;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 634) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0003;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 635) {
      nut_data[ 0] =         0.0004;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 636) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0131;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0052;
    }
    else if (n == 637) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 638) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0015;
      nut_data[ 2] =         0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 639) {
      nut_data[ 0] =        -0.0460;
      nut_data[ 1] =        -0.0435;
      nut_data[ 2] =        -0.0232;
      nut_data[ 3] =         0.0246;
      nut_data[ 4] =         0.0422;
    }
    else if (n == 640) {
      nut_data[ 0] =         0.0266;
      nut_data[ 1] =        -0.0078;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0111;
    }
    else if (n == 641) {
      nut_data[ 0] =        -0.0006;
      nut_data[ 1] =        -0.0009;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 642) {
      nut_data[ 0] =         0.0013;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0007;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 643) {
      nut_data[ 0] =         0.0015;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0006;
    }
    else if (n == 644) {
      nut_data[ 0] =         0.0006;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 645) {
      nut_data[ 0] =        -0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 646) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =         0.0003;
      nut_data[ 3] =         0.0003;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 647) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =        -0.0012;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 648) {
      nut_data[ 0] =         0.0235;
      nut_data[ 1] =         0.0334;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0163;
    }
    else if (n == 649) {
      nut_data[ 0] =         0.1200;
      nut_data[ 1] =         0.0598;
      nut_data[ 2] =         0.0319;
      nut_data[ 3] =        -0.0641;
      nut_data[ 4] =         0.0895;
    }
    else if (n == 650) {
      nut_data[ 0] =         0.0425;
      nut_data[ 1] =         0.0212;
      nut_data[ 2] =        -0.0133;
      nut_data[ 3] =         0.0269;
      nut_data[ 4] =         0.0355;
    }
    else if (n == 651) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0015;
      nut_data[ 2] =         0.0006;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 652) {
      nut_data[ 0] =        -0.0041;
      nut_data[ 1] =         0.0175;
      nut_data[ 2] =         0.0076;
      nut_data[ 3] =         0.0017;
      nut_data[ 4] =         0.0106;
    }
    else if (n == 653) {
      nut_data[ 0] =         0.0005;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0003;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 654) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0006;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 655) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0008;
      nut_data[ 2] =        -0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 656) {
      nut_data[ 0] =         0.0226;
      nut_data[ 1] =         0.0101;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0099;
    }
    else if (n == 657) {
      nut_data[ 0] =         0.0284;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0151;
      nut_data[ 4] =         0.0189;
    }
    else if (n == 658) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0024;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0010;
    }
    else if (n == 659) {
      nut_data[ 0] =        -0.0007;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 660) {
      nut_data[ 0] =         0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0006;
      nut_data[ 4] =         0.0008;
    }
    else if (n == 661) {
      nut_data[ 0] =         0.0009;
      nut_data[ 1] =        -0.0027;
      nut_data[ 2] =        -0.0014;
      nut_data[ 3] =        -0.0005;
      nut_data[ 4] =         0.0019;
    }
    else if (n == 662) {
      nut_data[ 0] =        -0.0026;
      nut_data[ 1] =        -0.0029;
      nut_data[ 2] =        -0.0016;
      nut_data[ 3] =         0.0014;
      nut_data[ 4] =         0.0026;
    }
    else if (n == 663) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =        -0.0013;
      nut_data[ 2] =        -0.0007;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0009;
    }
    else if (n == 664) {
      nut_data[ 0] =         0.0103;
      nut_data[ 1] =        -0.0060;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0048;
    }
    else if (n == 665) {
      nut_data[ 0] =         0.0026;
      nut_data[ 1] =        -0.0009;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0011;
    }
    else if (n == 666) {
      nut_data[ 0] =         0.0011;
      nut_data[ 1] =        -0.0024;
      nut_data[ 2] =        -0.0011;
      nut_data[ 3] =        -0.0009;
      nut_data[ 4] =         0.0018;
    }
    else if (n == 667) {
      nut_data[ 0] =        -0.1444;
      nut_data[ 1] =         0.2409;
      nut_data[ 2] =        -0.1286;
      nut_data[ 3] =        -0.0771;
      nut_data[ 4] =         0.1874;
    }
    else if (n == 668) {
      nut_data[ 0] =        -0.3084;
      nut_data[ 1] =         0.5123;
      nut_data[ 2] =         0.2735;
      nut_data[ 3] =         0.1647;
      nut_data[ 4] =         0.3989;
    }
    else if (n == 669) {
      nut_data[ 0] =        -0.0491;
      nut_data[ 1] =         0.0128;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0203;
    }
    else if (n == 670) {
      nut_data[ 0] =         0.0031;
      nut_data[ 1] =        -0.0481;
      nut_data[ 2] =        -0.0257;
      nut_data[ 3] =        -0.0017;
      nut_data[ 4] =         0.0322;
    }
    else if (n == 671) {
      nut_data[ 0] =         0.0014;
      nut_data[ 1] =        -0.0218;
      nut_data[ 2] =         0.0117;
      nut_data[ 3] =         0.0008;
      nut_data[ 4] =         0.0146;
    }
    else if (n == 672) {
      nut_data[ 0] =        -0.0012;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0005;
    }
    else if (n == 673) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 674) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0006;
      nut_data[ 2] =         0.0002;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0003;
    }
    else if (n == 675) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 676) {
      nut_data[ 0] =         0.0099;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0053;
      nut_data[ 4] =         0.0066;
    }
    else if (n == 677) {
      nut_data[ 0] =        -0.0462;
      nut_data[ 1] =         0.1604;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0668;
    }
    else if (n == 678) {
      nut_data[ 0] =        -0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0001;
    }
    else if (n == 679) {
      nut_data[ 0] =        -0.0219;
      nut_data[ 1] =         0.0089;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0095;
    }
    else if (n == 680) {
      nut_data[ 0] =        -0.0114;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0061;
      nut_data[ 4] =         0.0076;
    }
    else if (n == 681) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0002;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 682) {
      nut_data[ 0] =         0.0003;
      nut_data[ 1] =        -0.0007;
      nut_data[ 2] =        -0.0003;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0004;
    }
    else if (n == 683) {
      nut_data[ 0] =         0.0000;
      nut_data[ 1] =         0.0005;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0002;
    }
    else if (n == 684) {
      nut_data[ 0] =         0.0125;
      nut_data[ 1] =        -0.0043;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =        -0.0054;
      nut_data[ 4] =         0.0076;
    }
    else if (n == 685) {
      nut_data[ 0] =         0.0056;
      nut_data[ 1] =        -0.0117;
      nut_data[ 2] =        -0.0042;
      nut_data[ 3] =        -0.0040;
      nut_data[ 4] =         0.0078;
    }
    else if (n == 686) {
      nut_data[ 0] =         0.1440;
      nut_data[ 1] =         0.0000;
      nut_data[ 2] =         0.0000;
      nut_data[ 3] =         0.0000;
      nut_data[ 4] =         0.0576;
    }
  }
  return;
}
