#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"

#define ITRF_R   6378136.60000
#define ITRF_e2  6.69439799586e-3

/***************************************************************/
/* Reference Ellipsoid :                                       */
/*  IERS Techinical Notes No. 32                               */
/*        (IERS Convensions 2003)                              */
/*                 ed. by McCarthy and Petit (2003)            */
/*                                                             */
/*      R   : 6363672.6 m                                      */
/*      1/f : 298.25642                                        */
/*      e^2 = 1 - (1-f)^2                                      */
/***************************************************************/

void topocentric_equatorial_rectangular_coordinate(int    *TimUTC,
                                                   double  UT1_UTC,
                                                   double lambda,
                                                   double fai,
                                                   double h,
                                                   double *earth)
{
  double s[3];

  J_system_geocentric_equatorial_rectangular_coordinate(lambda, fai, h, s);
  earth_rotation(TimUTC, UT1_UTC,  s,   earth);

  return;
}



void J_system_geocentric_equatorial_rectangular_coordinate
                                                   (double lambda,
                                                    double fai,
                                                    double h,
                                                    double *s)
{
  double N;

/** J system geocentric equatorial rectangular coordinate **/

  N = (double)ITRF_R / sqrt(1.0 - (double)ITRF_e2*sin(fai)*sin(fai));
  s[0] = (N + h) * cos (fai) * cos (lambda);
  s[1] = (N + h) * cos (fai) * sin (lambda);
  s[2] = (N * (1.0 - (double)ITRF_e2) + h) * sin (fai);

  return;
}


void J_system_geocentric_equatorial_rectangular_coordinate2llh
                                           (double *lambda,
                                            double *fai,
                                            double *h,
                                            double *s)
{
  double R, tan_fai0, tan_fai1;
  double N, r;

/** J system geocentric equatorial rectangular coordinate **/

  *lambda = atan2(s[1], s[0]);
  R = sqrt(s[0]*s[0] + s[1]*s[1]);

  tan_fai0 = s[2] / R;
  while (1) {
    tan_fai1 = s[2] / R + (double)ITRF_R * (double)ITRF_e2 / R
                  * tan_fai0 / sqrt(1.0 + (1.0-ITRF_e2)*tan_fai0*tan_fai0);
    if (fabs((tan_fai1 - tan_fai0) / tan_fai1) < 1.0e-10) {
      break;
    } else {
      tan_fai0 = tan_fai1;
    }
  }

  *fai = atan(tan_fai1);
  N = (double)ITRF_R / sqrt(1.0 - ITRF_e2*sin(*fai)*sin(*fai));
  *h = sqrt(R*R + pow((s[2] + N*ITRF_e2*sin(*fai)), 2.0)) - N;

  return;
}
