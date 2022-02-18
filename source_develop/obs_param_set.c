#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>



int  obs_param_set(_Bool  *ERROR_FLAG,
                   int    SRT_NUM,
                   struct char_srt_info       *ch_srt,
                   struct srt_orbit_parameter *srt,
                   char   *ch_grt_el_lim,
                   double *grt_elevation_limit,
                   double sep_angle_limit_from_earth_limb,
                   struct char_obs_time *ch_obs_t,
                   int    *TimUTC, double *UT1_UTC,
                   double *obs_duration,
                   int    SRCPROC_MODE,
                   struct char_src_info    *ch_src,
                   struct pair_src_info    *pair_src,
                   struct source_parameter *src,
                   struct source_parameter *sun)
{
  int     i, iant;
  int     SEP_MODE, POS_MODE;
  int     timUTC[6];
  double  inner_sun_angle;
  double  f1, f2;
  double  DPI;

/*
-----------------------------------------------------
*/

  DPI = dpi / 180.0;

/*
-----------------------------------------------------
*/

  sscanf(ch_obs_t->start_t[0], "%4d", &TimUTC[0]);
  sscanf(ch_obs_t->start_t[1], "%2d", &TimUTC[1]);
  sscanf(ch_obs_t->start_t[2], "%2d", &TimUTC[2]);
  sscanf(ch_obs_t->start_t[3], "%2d", &TimUTC[3]);
  sscanf(ch_obs_t->start_t[4], "%2d", &TimUTC[4]);
  sscanf(ch_obs_t->start_t[5], "%2d", &TimUTC[5]);
  *UT1_UTC       = (double)   0.0;
  sscanf(ch_obs_t->obsd, "%lf", obs_duration);
  *obs_duration *= (double)3600.0;
  ch_time_set(TimUTC, ch_obs_t);

/*
============================================================
*/

  sscanf(ch_grt_el_lim, "%lf", grt_elevation_limit);
  *grt_elevation_limit *= DPI;

/*
============================================================
*/

  in__src_proc(SRCPROC_MODE, &SEP_MODE, &POS_MODE);
  source_position(src, pair_src, ch_src, SEP_MODE, POS_MODE);

/*
============================================================
*/

  inner_sun_angle
      = sun_angle_check(TimUTC, *obs_duration, &src[0], sun) / DPI;
  if (SRT_NUM >= 1) {
    if (inner_sun_angle < 50.0) {
      printf("ERROR: Inner Sun Angle is less than 50 deg. (%5.1lf [deg])\n",
             inner_sun_angle);
      return -1;
    }
  }

/*
============================================================
*/

  for (iant=0; iant<SRT_NUM; iant++) {
    sscanf(ch_srt[iant].apo, "%lf", &srt[iant].apogee);
    srt[iant].apogee *= 1.0e3;
    sscanf(ch_srt[iant].per, "%lf", &srt[iant].perigee);
    srt[iant].perigee *= 1.0e3;
    sscanf(ch_srt[iant].inc, "%lf", &srt[iant].inclination);
    srt[iant].inclination *= DPI;
    sscanf(ch_srt[iant].OMG, "%lf", &srt[iant].Omega);
    srt[iant].Omega       *= DPI;
    sscanf(ch_srt[iant].omg, "%lf", &srt[iant].omega);
    srt[iant].omega       *= DPI;
    sscanf(ch_srt[iant].t_0, "%4d%2d%2d%2d%2d%2d",
           &timUTC[0], &timUTC[1], &timUTC[2],
           &timUTC[3], &timUTC[4], &timUTC[5]);
    srt[iant].t0 = MJD(timUTC[0], timUTC[1], timUTC[2],
                       timUTC[3], timUTC[4], timUTC[5], 0.0);
    sscanf(ch_srt[iant].d_OMG, "%lf", &srt[iant].d_Omega);
    srt[iant].d_Omega     *= (DPI / 365.25);
    sscanf(ch_srt[iant].d_omg, "%lf", &srt[iant].d_omega);
    srt[iant].d_omega     *= (DPI / 365.25);

    srt[iant].a             = earth_radius
                          + (srt[iant].apogee + srt[iant].perigee) / 2.0;
    srt[iant].e             = 1.0
                        - (earth_radius + srt[iant].perigee) / srt[iant].a;
    srt[iant].ef            = sqrt(1.0 - pow(srt[iant].e, 2.0));
    srt[iant].n             = sqrt(GM / pow(srt[iant].a, 3.0));
    f1                      = earth_radius * earth_radius * sqrt(GM)
                            * 86400.0 / 360.0
                            / pow(srt[iant].ef,     2.0)
                            / pow(srt[iant].a,  7.0/2.0);
    f2                      = cos(srt[iant].inclination);

/**** ASTRO-G ****/
/****
    srt[iant].d_Omega       =  0.0;
    srt[iant].d_omega       =  0.0;
    srt[iant].d_Omega       = -0.584 * f1 * f2;
    srt[iant].d_omega       =  0.292 * f1 * (5.0 * f2 * f2 - 1.0);
    srt[iant].d_Omega       = -167.00 / 365.25 * DPI;
    srt[iant].d_omega       =  261.00 / 365.25 * DPI;
****/

    srt[iant].ODDA          = 0.0;
    srt[iant].initial_phase = 2.0 * dpi * random_val0();

    srt[iant].uplink_freq   = 40.0e9;
    srt[iant].downlink_freq = 37.0e9;
  }

  if (SRT_NUM == 0) {
    ERROR_FLAG[SRTAER] = false;
  }

  return (1);
}
