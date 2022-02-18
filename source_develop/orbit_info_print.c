#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


void orbit_info_print(FILE  *log_fp,
               int SRT_NUM,
               struct srt_orbit_parameter *srt,
               int *TimUTC, double UT1_UTC)
{
  int    iant;
  int    timUTC[6], ut1_utc;
  float  DPI;

/*
--------------------------------------------------
*/

  DPI = (float)(dpi / 180.0);

/*
------------------------------------------
*/

  for (iant=0; iant<SRT_NUM; iant++) {
    fprintf(log_fp, "ORBIT : inclination  [deg]      : %f\n",
            srt[iant].inclination / DPI);
    fprintf(log_fp, "ORBIT : Omega        [deg]      : %f\n",
            srt[iant].Omega       / DPI);
    fprintf(log_fp, "ORBIT : omega        [deg]      : %f\n",
            srt[iant].omega       / DPI);
    fprintf(log_fp, "ORBIT : ecicen                  : %f\n", srt[iant].e);
    fprintf(log_fp, "ORBIT : apogee       [m]        : %f\n", srt[iant].apogee);
    MJD2date(srt[iant].t0, &timUTC[0], &timUTC[1], &timUTC[2],
                           &timUTC[3], &timUTC[4], &timUTC[5]);
    fprintf(log_fp,
            "ORBIT : time of perigee passage : %4d.%2d.%2d %2d:%2d:%2d\n",
            timUTC[0], TimUTC[1], timUTC[2], timUTC[3], timUTC[4], timUTC[5]);
  }

/*
----
*/

  return;
}
