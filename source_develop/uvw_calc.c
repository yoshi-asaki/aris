#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __DEBUG__
****/

int  uvw_calc(float  src_flag,
              int ANT_NUM, int GRT_NUM,
              int nobs,    int iobs,
              int *TimUTC, double UT1_UTC,
              struct antenna_parameter        *ant_prm,
              struct source_parameter         src,
              struct source_parameter         sun,
              struct phase_screen_parameter   *wvc,
              struct phase_screen_parameter   *dry,
              struct phase_screen_parameter   *ion,
              double *Cw,            double *Cd,
              double Ci[][86400],
              double CI,
              double *wvc_ds,      
              double *dry_ds,
              double *ion_ds,
              double *fqs_ds,
              _Bool  ERROR_PROC_SWT, _Bool  *ERROR_FLAG,
              double W[][3],         double dUT1,
              struct st_observable  *int_obs,
              double OBS_T, double OBS_p,
              struct atmospheric_zenith_error *dz,
              _Bool  SRT_STATUS_CHK,
              struct srt_orbit_parameter      *srt,
              double *init_l,
              int    TRK_NUM,
              struct antenna_parameter        *trk_pos,
              struct srt_data_link  *srt_link,
              double sep_angle_limit, struct TID TID)
{
  int    i, I, iant, jant, nant, itrk, ntrk, nsec;
  float  WEIGHT_FLAG_DATA;
  double ant_XYZ[3], err_XYZ[3];
  double V_xyz[3];
  double  delay,    rate;
  double ddelay,   drate;
  double AZ, EL, dAZdt, dELdt;
  double srt_AZ, srt_EL;
  double s[3], s_e[3];
  double p[3], ptrk[3];
  struct source_parameter SRC, SUN;
  double earth_center[3];
  double ehta, arc;
  double ion_const, ids_const;
  double delta_z, rho;
  double vtmp[3];
  double z, Z, r;
  double gst, gst_e;     /* Greenwich Siderial Time */
  double S[3];
  double hnmf, wnmf;

/*
===================================================================
*/

  ion_const = TEC_CONST / speed_of_light;
  ids_const = CI * ion_const;

  gst       = GST(TimUTC, 0.0,  UT1_UTC);
  if (ERROR_FLAG[EOPERR] == true) {
    gst_e   = GST(TimUTC, dUT1, UT1_UTC);
  } else {
    gst_e   = gst;
  }

  SRC = src;
  drotate(SRC.s,   -gst,   "z");
  SRC.RA   = atan2(SRC.s[1], SRC.s[0]);
  SRC.DC   = atan2(SRC.s[2], vlen2(SRC.s));

  drotate(SRC.s_e, -gst_e, "z");
  SRC.RA_e = atan2(SRC.s_e[1], SRC.s_e[0]);
  SRC.DC_e = atan2(SRC.s_e[2], vlen2(SRC.s_e));

/*
===================================================================
*/

  for (i=0; i<3; i++) {
    earth_center[i] = 0.0;
  }

  nant = 0;
  for (iant=0; iant<ANT_NUM; iant++) {
    I = iant * nobs + iobs;
    int_obs[I].ion_delay = 0.0;
    int_obs[I].grp_delay = 0.0;
    WEIGHT_FLAG_DATA     = 0.0;

/*
-------- GRT --------
*/

    if (iant < GRT_NUM) {
      azel_position(TimUTC, UT1_UTC,
                    ant_prm[iant].LLH[0],
                    ant_prm[iant].LLH[1],
                    ant_prm[iant].LLH[2],
                    OBS_T, OBS_p, false, src.RA, src.DC,
                    &(int_obs[I].az), &(int_obs[I].el),
                    &dAZdt, &dELdt, 0.0, vtmp);
#ifdef __DEBUG__
      printf("__DEBUG__  %d   (%lf  %lf  %lf)    (%lf  %lf  %lf)    %lf  %lf\n",
                                          iant,
                                          ant_prm[iant].LLH[0],
                                          ant_prm[iant].LLH[1],
                                          ant_prm[iant].LLH[2],
                                          ant_prm[iant].XYZ[0],
                                          ant_prm[iant].XYZ[1],
                                          ant_prm[iant].XYZ[2],
                                          int_obs[I].az*180.0/dpi, int_obs[I].el*180.0/dpi);
#endif

/**** SPECIAL TEST 2014/06/01 ****/
/****
      if (int_obs[I].el < ant_prm[iant].ELLIM) {
        WEIGHT_FLAG_DATA += GRT_ELEVATION_LIMIT;
      }
****/
/****
****/
/**** SPECIAL TEST 2014/06/01 ****/

/*
-------- SRT --------
*/

    } else if (iant >= GRT_NUM) {
      jant = iant - GRT_NUM;

/*
-- Spacecraft Position & Tracking Status --
*/

      tracking_status(SRT_STATUS_CHK, &WEIGHT_FLAG_DATA, &ntrk,
                      TimUTC, UT1_UTC, srt[jant], src, sun,
                      OBS_T,  OBS_p, TRK_NUM, ant_XYZ, err_XYZ, V_xyz,
                      trk_pos, srt_link,
                      &AZ, &EL, &dAZdt, &dELdt,
                      &srt_AZ, &srt_EL,
                      sep_angle_limit, ptrk, p, init_l+jant);
      ant_prm[iant].XYZ[0] = ant_XYZ[0];
      ant_prm[iant].XYZ[1] = ant_XYZ[1];
      ant_prm[iant].XYZ[2] = ant_XYZ[2];
      drotate(ant_prm[iant].XYZ, -gst, "z");
      int_obs[I].az = 0.0;
      int_obs[I].el = 0.5 * dpi;

#ifdef __DEBUG__
      printf("__DEBUG__  %lf, %lf, %lf, %lf, %lf, %lf\n",
              ant_XYZ[0], ant_XYZ[1], ant_XYZ[2],
              V_xyz[0],   V_xyz[1],   V_xyz[2]);
#endif /* __DEBUG__ */

    }

    weight_flag_add(src_flag, &(int_obs[I].wt), WEIGHT_FLAG_DATA);
    if (weight_flag_chk(src_flag, &(int_obs[I].wt),
                                     GRT_ELEVATION_LIMIT) == false) {

      position2uvw(&(int_obs[I].u), &(int_obs[I].v), &(int_obs[I].w),
                   &delay, &rate, ant_prm[iant].XYZ, earth_center,
                   SRC.RA, SRC.DC, SRC.s);
/****
      printf("__DEGBU_UVW_CALC_(1)__ %lf  %lf  %lf\n", int_obs[I].u, int_obs[I].v, int_obs[I].w);
      printf("__DEGBU_UVW_CALC_(2)__ %lf  %lf  %lf\n", ant_prm[iant].XYZ[0],
                                                       ant_prm[iant].XYZ[1],
                                                       ant_prm[iant].XYZ[2]);
****/
      if (int_obs[I].wt == src_flag) {
        int_obs[I].wt += 1.0;
      }

/*
===========================================================
*/

      if (ERROR_PROC_SWT == true) {


/*
--------------------------
*/

        if (ERROR_FLAG[IONTRB] == true) {
          sun_position(TimUTC, UT1_UTC, &SUN.RA, &SUN.DC);
          SUN.s[0] = cos(SUN.DC) * cos(SUN.RA);
          SUN.s[1] = cos(SUN.DC) * sin(SUN.RA);
          SUN.s[2] = sin(SUN.DC);
          drotate(SUN.s,   -gst,   "z");
        }

/*
--------------------------
*/

        if (ERROR_FLAG[APOSER] == true || ERROR_FLAG[RPOSER] == true || 
            ERROR_FLAG[EOPERR] == true) {
          if (iant < GRT_NUM) {
            err_XYZ[0]  = ant_prm[iant].ERR[0];
            err_XYZ[1]  = ant_prm[iant].ERR[1];
            err_XYZ[2]  = ant_prm[iant].ERR[2];
          } else if (iant >= GRT_NUM) { 
            drotate(err_XYZ, -gst, "z");
            err_XYZ[0] += ant_prm[iant].XYZ[0];
            err_XYZ[1] += ant_prm[iant].XYZ[1];
            err_XYZ[2] += ant_prm[iant].XYZ[2];
            dvector_calc(W, err_XYZ);
          }
          position2uvw(&(int_obs[I].u), &(int_obs[I].v), &(int_obs[I].w),
                       &ddelay, &drate, err_XYZ, earth_center,
                       SRC.RA_e, SRC.DC_e, SRC.s_e);
          ddelay -= delay;
          drate  -= rate;
          int_obs[I].grp_delay += ddelay;
        }

/*
--------------------------
*/

        if (ERROR_FLAG[FQSERR] == true) {
          int_obs[I].grp_delay += fqs_ds[I];
        }

/*
--------------------------
*/

        if (iant < GRT_NUM) {

/*
--------------------------
*/

          if (ERROR_FLAG[TWVTRB] == true) {
            spherical_geometry(wvc[iant].H_d, int_obs[I].el, &z,
                               &rho, &EL, &Z, &delta_z);
            ddelay = Cw[iant] * sqrt(airmass(EL)) * *(wvc_ds + I);
            int_obs[I].grp_delay += ddelay;
          }

/*
--------------------------
*/

          if (ERROR_FLAG[DRYTRB] == true) {
            spherical_geometry(dry[iant].H_d, int_obs[I].el, &z,
                               &rho, &EL, &Z, &delta_z);
            ddelay = Cd[iant] * sqrt(airmass(EL)) * *(dry_ds + I);
            int_obs[I].grp_delay += ddelay;
          }

/*
--------------------------
*/

          if (ERROR_FLAG[IONTRB] == true) {
            spherical_geometry(ion[iant].H_d, int_obs[I].el, &z,
                               &rho, &EL, &Z, &delta_z);

            r = sin(delta_z) / sin(z) * (ion[iant].H_d + earth_radius);
            S[0] = r * cos(int_obs[I].az) * cos(int_obs[I].el);
            S[1] = r * sin(int_obs[I].az) * cos(int_obs[I].el);
            S[2] = r * sin(int_obs[I].el);

            drotate(S, ant_prm[iant].LLH[1]-0.5*dpi, "x");
            drotate(S, ant_prm[iant].LLH[0]-0.5*dpi, "z");
            S[0] += ant_prm[iant].XYZ[0];
            S[1] += ant_prm[iant].XYZ[1];
            S[2] += ant_prm[iant].XYZ[2];

            nsec = (int)(86400.0 * (atan2(
                   (SUN.s[0]*S[1] - SUN.s[1]*S[0]),
                   (SUN.s[0]*S[0] + SUN.s[1]*S[1])) / 2.0 / dpi + 0.5));
            while (nsec < 0 || nsec >= 86400) {
              if (nsec < 0) {
                nsec += 86400;
              } else if (nsec >= 86400) {
                nsec -= 86400;
              }
            }

            if (ant_prm[iant].LLH[1] >= 0.0) {
              ddelay
                 = Ci[0][nsec] * sqrt(airmass(EL)) * ids_const * *(ion_ds + I);
            } else {
              ddelay
                 = Ci[1][nsec] * sqrt(airmass(EL)) * ids_const * *(ion_ds + I);
            }
            int_obs[I].ion_delay -= ddelay;
          }

/*
--------------------------
*/

          if (ERROR_FLAG[TDSECZ] == true) {
            nmf20(&hnmf, &wnmf,           ant_prm[iant].nmfh,
                  ant_prm[iant].nmfh_hcr, ant_prm[iant].nmfw,
                  ant_prm[iant].LLH[2],   int_obs[I].el);
            ddelay = dz[iant].trp * wnmf;
            int_obs[I].grp_delay += ddelay;
          }

/*
--------------------------
*/

          if (ERROR_FLAG[IDSECZ] == true) {
            spherical_geometry(ion[iant].H_s, int_obs[I].el, &z,
                               &rho, &EL, &Z, &delta_z);
            ddelay = ion_const * dz[iant].tec / sin(EL);;
            int_obs[I].ion_delay -= ddelay;
          }

/********
          if (GLOBAL_SWT[1] == false) {
            int_obs[I].grp_delay = 0.0;
            int_obs[I].ion_delay = 0.0;
          }
*********/

/*
--------------------------
*/

          if (ERROR_FLAG[AMPERR] == true) {
/********
            if (ERROR_FLAG[TWVTRB] == true) {
              spherical_geometry(wvc[iant].H_d, int_obs[I].el, &z,
                                 &rho, &EL, &Z, &delta_z);
              int_obs[I].amp_error +=
                     1.00 * sqrt(airmass(EL)) * *(wvc_ds + I);
            }

            if (ERROR_FLAG[TDSECZ] == true) {
              nmf20(&hnmf, &wnmf,           ant_prm[iant].nmfh,
                    ant_prm[iant].nmfh_hcr, ant_prm[iant].nmfw,
                    ant_prm[iant].LLH[2],   int_obs[I].el);
              int_obs[I].amp_error += ant_prm[iant].d_gain * wnmf;
            } else {
              int_obs[I].amp_error += ant_prm[iant].d_gain;
            }
********/
            int_obs[I].amp_error += ant_prm[iant].d_gain;
          }

/*
--------------------------
*/

        } else if (iant >= GRT_NUM) {
          jant = iant - GRT_NUM;
          srt[jant].uplink_freq   = 40.0e9;
          srt[jant].downlink_freq = 37.0e9;

/*
--------------------------
*/

#ifdef __AAAA__
          if (ERROR_FLAG[IONTRB] == true && ntrk != -1) {
            spherical_geometry(ion[iant].H_d, int_obs[I].el, &z,
                               &rho, &EL, &Z, &delta_z);

            r = sin(delta_z) / sin(z) * (ion[iant].H_d + earth_radius);
            S[0] = r * cos(int_obs[I].az) * cos(int_obs[I].el);
            S[1] = r * sin(int_obs[I].az) * cos(int_obs[I].el);
            S[2] = r * sin(int_obs[I].el);

            drotate(S, trk_pos[ntrk].LLH[1]-0.5*dpi, "x");
            drotate(S, trk_pos[ntrk].LLH[0]-0.5*dpi, "z");
            S[0] += trk_pos[iant].XYZ[0];
            S[1] += trk_pos[iant].XYZ[1];
            S[2] += trk_pos[iant].XYZ[2];

            arc = r * (0.5 * dpi - atan2(S[2], vlen2(S)));
            ehta = 2.0 * dpi * (TID.v[1] * (double)iobs - arc) / TID.lambda;
            *(ion_ds + I) += TID.amp * sin(ehta);

            nsec = (int)(86400.0 * (atan2(
                   (SUN.s[0]*S[1] - SUN.s[1]*S[0]),
                   (SUN.s[0]*S[0] + SUN.s[1]*S[1])) / 2.0 / dpi + 0.5));
            while (nsec < 0 || nsec >= 86400) {
              if (nsec < 0) {
                nsec += 86400;
              } else if (nsec >= 86400) {
                nsec -= 86400;
              }
            }

            if (trk_pos[ntrk].LLH[1] >= 0.0) {
              ddelay
                 = Ci[0][nsec] * sqrt(airmass(EL)) * ids_const * *(ion_ds + I);
            } else {
              ddelay
                 = Ci[1][nsec] * sqrt(airmass(EL)) * ids_const * *(ion_ds + I);
            }
/****
            int_obs[I].local_phs -= 2.0 * dpi * ddelay
                 * ((1.0 / srt[jant].uplink_freq)
                  - (1.0 / srt[jant].downlink_freq));
****/
            int_obs[I].local_phs -= 0.5 * ddelay
                  / pow(srt[jant].uplink_freq, 2.0) * 
              (pow(srt[jant].downlink_freq/srt[jant].uplink_freq, 2.0) - 1.0)
             / pow(srt[jant].downlink_freq/srt[jant].uplink_freq, 2.0);
          }
#endif

/*
--------------------------
*/

          if (ERROR_FLAG[IDSECZ] == true) {
            spherical_geometry(ion[iant].H_s, int_obs[I].el, &z,
                               &rho, &EL, &Z, &delta_z);
            ddelay = ion_const * dz[iant].tec / sin(EL);;
            int_obs[I].local_phs -= 0.5 * ddelay
                  / pow(srt[jant].uplink_freq, 2.0) * 
             (pow(srt[jant].downlink_freq/srt[jant].uplink_freq, 2.0) - 1.0)
             / pow(srt[jant].downlink_freq/srt[jant].uplink_freq, 2.0);
          }

/*
--------------------------
*/

          if (ERROR_FLAG[AMPERR] == true) {
            int_obs[I].amp_error += ant_prm[iant].d_gain;
          }

/*
--------------------------
*/

        }

/*
--------------------------
*/

      }
      nant++;
    }

/*
--------------------------
*/

/********
    if (iant >= GRT_NUM) { 
      ant_prm[iant].XYZ[0] = 0.0;
      ant_prm[iant].XYZ[1] = 0.0;
      ant_prm[iant].XYZ[2] = 0.0;
      ant_prm[iant].ERR[0] = 0.0;
      ant_prm[iant].ERR[1] = 0.0;
      ant_prm[iant].ERR[2] = 0.0;
    }
********/
  }
  return (nant);
}
