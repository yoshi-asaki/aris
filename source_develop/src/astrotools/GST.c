#include <math.h>
#include <stdio.h>
#include "astrotools.h"


/***************************************************************/
/*                                                             */
/* The IAU 2000 Theory of Precession/Nutation: Precession      */
/*  IERS Techinical Notes No. 32                               */
/*        (IERS Convensions 2003)                              */
/*                 ed. by McCarthy and Petit (2003)            */
/*                                                             */
/***************************************************************/

double GST(int *TimUTC, double dT, double UT1_UTC)
{
  double t, theta;
  double dpi = 3.141592653589793238462643;
  double gst_offset = 0.7790572732640;
  double gst_factor = 1.00273781191135448;

  t = MJD(TimUTC[0], TimUTC[1], TimUTC[2],
          TimUTC[3], TimUTC[4], TimUTC[5], 0.0) - 51544.5;
  t *= gst_factor;
  t -= (double)((int)t);
  t += (dT + UT1_UTC) / 86400.0 * gst_factor;

/****
     Tu = Julian day in UT1 - 2451545.0 
        = MJD in UTC + 2400000.5 + UT1_UTC - 2451545.0 
        = MJD in UTC + UT1_UTC - 51544.5;
     theta = 2.0 * dpi * (0.7790572732640 + 1.00273781191135448 * Tu)
****/

  theta = 2.0 * dpi * (gst_offset + t);

  return (theta);
}
