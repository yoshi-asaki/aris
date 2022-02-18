#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


_Bool tracking_status(_Bool  TRK_STATUS_CHK,  float  *WEIGHT_FLAG_DATA,
              int *ntrk,
              int *TimUTC, double UT1_UTC,
              struct srt_orbit_parameter      srt,
              struct source_parameter         src,
              struct source_parameter         sun,
              double OBS_T, double OBS_p,
              int    TRK_NUM,
              double *ant_XYZ,     double *err_XYZ,  double *V_xyz_CIP,
              struct antenna_parameter    *trk_pos,
              struct srt_data_link        *srt_link,
              double *AZ, double *EL, double *dAZdt, double *dELdt,
              double *srt_AZ, double *srt_EL,
              double sep_angle_limit,
              double *ptrk, double *pdis,
              double *init_l)
{
  int    itrk, I;
  int    trk_status[ANTMAX], trk_priority[ANTMAX];
  double L, R;
  double Oe_rdc_CIP[3];
  double srt_dAZdt, srt_dELdt;
  _Bool  TRK_SWT;

/*
--------
*/

  *ntrk   = -1;

/*
-- SRT Orbit Status Check      --
*/

  spacecraft_position(srt, TimUTC, UT1_UTC, (double)1,
                      ant_XYZ, err_XYZ, Oe_rdc_CIP, V_xyz_CIP, init_l);
  if (TRK_STATUS_CHK == false) {
    return true;
  }
/****SPECIAL TEST  2014/06/01****/
/****
  return true;
****/
/****SPECIAL TEST  2014/06/01****/

/*
----
*/

  if (earth_eclipse(sun.s, ant_XYZ, &R, &L, sep_angle_limit) == NIGHT) {
    *WEIGHT_FLAG_DATA += SRT_EARTH_ECLIPS_SUN;
  }

  if (earth_eclipse(src.s, ant_XYZ, &R, &L, sep_angle_limit) == NIGHT) {
    *WEIGHT_FLAG_DATA += SRT_EARTH_ECLIPS_SRC;
  }

  if (TRK_NUM == 0) {
    return true;
  }

  for (itrk=0; itrk<TRK_NUM; itrk++) {
    trk_priority[itrk] = trk_pos[itrk].priority % ON_TRK;
    trk_status[itrk]   = trk_pos[itrk].priority / ON_TRK;
  }

/*
-- Tracking Continuation Check --
*/

  TRK_SWT = false;
  for (I=1; I<=TRK_NUM; I++) {
    for (itrk=0; itrk<TRK_NUM; itrk++) {
      if (trk_priority[itrk] == I && trk_status[itrk] == 1) {
        if (tracking_condition(TimUTC, UT1_UTC, trk_pos[itrk],
                    OBS_T, OBS_p, false, ant_XYZ, V_xyz_CIP,
                    AZ, EL, dAZdt, dELdt,
                    srt_link, src, sun,
                    srt_AZ, srt_EL, &srt_dAZdt, &srt_dELdt, ptrk, pdis,
                    srt.BODY_X_SUN) == 0) {
          *ntrk   = itrk;
          TRK_SWT = true;
          break;
        }
      }
    }

    if (TRK_SWT == true) {
      for (itrk=0; itrk<TRK_NUM; itrk++) {
        if (trk_priority[itrk] == I && trk_status[itrk] != 1) {
          if (tracking_condition(TimUTC, UT1_UTC, trk_pos[itrk],
                      OBS_T, OBS_p, false, ant_XYZ, V_xyz_CIP,
                      AZ, EL, dAZdt, dELdt,
                      srt_link, src, sun,
                      srt_AZ, srt_EL, &srt_dAZdt, &srt_dELdt, ptrk, pdis,
                      srt.BODY_X_SUN) == 0) {
            *ntrk   = itrk;
            TRK_SWT = true;
            break;
          }
        }
      }
    }

    if (TRK_SWT == true) {
      break;
    }
  }

/*
--------
*/

  if (TRK_SWT == true) {
    for (itrk=0; itrk<TRK_NUM; itrk++) {
      trk_pos[itrk].priority = trk_priority[itrk];
      trk_pos[itrk].UFL      = true;
    }
    trk_pos[*ntrk].priority += ON_TRK;
  } else if (TRK_SWT == false) {
    *WEIGHT_FLAG_DATA += (float)SRT_TRACKING_CONDITION_LIMIT;
  }

  return TRK_SWT;
}
