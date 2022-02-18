#include <stdio.h>
#include <math.h>
#include "astrotools.h"


int  UTC_minus_TAI(int    *TimUTC)
{
  int    utc_tai;
  double mjd;

  mjd = MJD(TimUTC[0], TimUTC[1], TimUTC[2],
            TimUTC[3], TimUTC[4], TimUTC[5], 0.0);

  if (mjd <  MJD(1972,  1,  1,  0,  0,  0, 0.0)) {
    printf("WARNING: ");
    printf("UTC_minus_TAI: TAI had not been difined at that time. \n");
    printf("WARNING: ");
    printf("Instead, UTC-TAI is set to -10 sec. \n");
    utc_tai = -10;
  } else if (mjd >= MJD(1972,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1972,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -10;
  } else if (mjd >= MJD(1972,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1973,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -11;
  } else if (mjd >= MJD(1973,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1974,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -12;
  } else if (mjd >= MJD(1974,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1975,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -13;
  } else if (mjd >= MJD(1975,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1976,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -14;
  } else if (mjd >= MJD(1976,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1977,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -15;
  } else if (mjd >= MJD(1977,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1978,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -16;
  } else if (mjd >= MJD(1978,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1979,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -17;
  } else if (mjd >= MJD(1979,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1980,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -18;
  } else if (mjd >= MJD(1980,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1981,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -19;
  } else if (mjd >= MJD(1981,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1982,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -20;
  } else if (mjd >= MJD(1982,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1983,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -21;
  } else if (mjd >= MJD(1983,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1985,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -22;
  } else if (mjd >= MJD(1985,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1988,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -23;
  } else if (mjd >= MJD(1988,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1990,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -24;
  } else if (mjd >= MJD(1990,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1991,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -25;
  } else if (mjd >= MJD(1991,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1992,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -26;
  } else if (mjd >= MJD(1992,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1993,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -27;
  } else if (mjd >= MJD(1993,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1994,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -28;
  } else if (mjd >= MJD(1994,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1996,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -29;
  } else if (mjd >= MJD(1996,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1997,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -30;
  } else if (mjd >= MJD(1997,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(1999,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -31;
  } else if (mjd >= MJD(1999,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(2006,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -32;
  } else if (mjd >= MJD(2006,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(2009,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -33;
  } else if (mjd >= MJD(2009,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(2012,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -34;
  } else if (mjd >= MJD(2012,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(2015,  7,  1,  0,  0,  0, 0.0)) {
    utc_tai = -35;
  } else if (mjd >= MJD(2015,  7,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(2017,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -36;
  } else if (mjd >= MJD(2017,  1,  1,  0,  0,  0, 0.0) &&
             mjd <  MJD(9999,  1,  1,  0,  0,  0, 0.0)) {
    utc_tai = -37;
  }

  return utc_tai;
}
