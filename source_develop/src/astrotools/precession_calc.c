#include <stdio.h>
#include <math.h>
#include "astrotools.h"

void precession_calc(int    *TimUTC,  double dT,
                     double *tueta,   double *theta,   double *z     )
/***************************************************************/
/*                                                             */
/* The IAU 2000 Theory of Precession/Nutation: Precession      */
/*  IERS Techinical Notes No. 32                               */
/*        (IERS Convensions 2003)                              */
/*                 ed. by McCarthy and Petit (2003)            */
/*                                                             */
/* Since, due to the theoritical bases, the original develop-  */
/* ment of the precession quantities as function of time can   */
/* be considered as being expressed in TDB, TT is used in the  */
/* following expressions in place of TDB. The largest term in  */
/* the difference TDB-TT being 1.7 ms X sin(l'), the resulting */
/* error in the precession quantity psi_A is periodic, with    */
/* an annual period and an amplitude of 2.7" X 10e-9, which    */
/* is significantly under the required microarcsecond          */
/* accuracy.                                                   */
/*                                                             */
/***************************************************************/

{
  double Delta_TT, Delta_UTC_TAI;
  double t1, t2, t3, t4, t5;
  double dpi = 3.141592653589793238462643;
  double DPI;

  DPI = dpi / 180.0 / 3600.0;

  Delta_UTC_TAI = (double)UTC_minus_TAI(TimUTC) + 32.0;
  /* -32 (sec) : UTC - TAI at J2000  */

  Delta_TT   = MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                   TimUTC[3], TimUTC[4], TimUTC[5], dT)
             - MJD(     2000,         1,         1,
                          12,         0,         0, 0.0)
             + (double)Delta_UTC_TAI / 86400.0;
  /* 32.184 sec to convert from TAI to TT is canceled. */

  t1 = Delta_TT / 36525.0;
  t2 = t1 * t1;
  t3 = t2 * t1;
  t4 = t3 * t1;
  t5 = t4 * t1;

  *tueta =  2.5976176 + 2306.0809506*t1  + 0.3019015*t2   + 0.0179663*t3
                      -    0.0000327*t4  - 0.0000002*t5;
  *theta =              2004.1917476*t1  - 0.4269353*t2   - 0.0418251*t3
                      -    0.0000601*t4  - 0.0000001*t5;
  *z     = -2.5976176 + 2306.0803226*t1  + 1.0947790*t2   + 0.0182273*t3
                      +    0.0000470*t4  - 0.0000003*t5;

  *tueta *= DPI;
  *z     *= DPI;
  *theta *= DPI;

  return;
}
