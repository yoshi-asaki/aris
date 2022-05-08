#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __DEBUG__
#define __ANT_DEBUG__
****/

int  antenna_visibility(
                  int  ANT_NUM, int  GRT_NUM,  int  SRT_NUM,
                  struct antenna_parameter *ant_prm,
                  char   antenna_code[][10],
                  struct source_parameter  *src,
                  struct source_parameter  sun,
                  struct phase_screen_parameter *wvc,
                  struct phase_screen_parameter *dry,
                  struct phase_screen_parameter *ion,
                  double *Cw, double *Cd,
                  double Ci[][86400],
                  double *fqs_ds,
                  float  *src_flag,
                  struct srt_orbit_parameter *srt,
                  struct antenna_parameter  *trk_pos,
                  double sep_angle_limit_from_earth_limb,
                  struct comment_param   *cmnt,
                  char   comment[][NCOMLEN],
                  _Bool  TV_SWT, int *pgid)
{
  int    i, j, k, n, I, ns, iant;
  int    SRCPROC_MODE;
  int    sobs, eobs, iobs, nobs, IOBS, NOBS;
  int    grt_num;
  int    TRK_NUM=0, itrk, DOY;
  int    TimUTC[6], timUTC[6];
  int    trk_num, priority[TRKMAX], trk_priority[TRKMAX];
  char   trk_name[TRKMAX][10], c_dummy = 0;
  char   antenna_list_file_name[500];
  FILE   *fp;
  double UT1_UTC,   ut1_utc;
  double obs_duration, obs_dh, f;
  double mjd, uv_max[SRC_NUM];
  double grt_elevation_limit;
  double dUT1=0.0, W[3][3];
  double OBS_T=290.0, OBS_p=1014.0;
  double CI=0.0;
  struct TID  MSTID;

  struct st_observable          *int_obs[SRC_NUM];
  double *wvc_ds[SRC_NUM], *dry_ds[SRC_NUM], *ion_ds[SRC_NUM];
  struct atmospheric_zenith_error *dz=NULL;

  float  ftmp1,   ftmp2;
  float  pgx[10], pgy[10];
  double inner_sun_angle;

  int    d_t = 60;
  float  pgxmin, pgxmax, pgymin, pgymax;
  double wave_length[2];
  double DPI;
  double init_l_p[SRTMAX], init_l_m[SRTMAX];

  float  bttn_box[20][4];
  float  cursor_pos[2];
  int    idum;
  _Bool  PROC_SWT, ERROR_FLAG[ERROR_NUM];
  char   string[100];
  struct array_parameter array;
  struct char_src_info  ch_src;
  struct pair_src_info  pair_src;
  struct char_srt_info  ch_srt[SRTMAX];
  struct char_obs_time  ch_obs_t;
  char   ch_grt_el_lim[20];
  int    wave_id[SRC_NUM_P1];
  struct srt_data_link  *srt_link=NULL;
  char   tv_time[30];
  float  *gcode_pos;
  float  scode_pos[2][SRTMAX];
  double fov_st;

/*
===================================================================
*/

  DPI = dpi / 180.0;

/*
-------------------------------------------
*/

  if (SRT_NUM >= 1) {
    if ((srt_link = (struct srt_data_link *)
                        calloc(1, sizeof(struct srt_data_link))) == NULL) {
      printf("ARIS: ANTENNA_VISIBILITY: ERROR: memory alloc for srt_link.\n");
      return (0);
    }
    get_srt_link(srt_link, &fov_st);
  }

/*
-------------------------------------------
*/

  wave_length[0] = 1.0;
  wave_length[1] = 1.0;

  bttn_box[0][0] = 0.70;
  bttn_box[0][1] = 0.90;
  bttn_box[0][2] = 0.18;
  bttn_box[0][3] = 0.21;

  bttn_box[1][0] = 0.70;
  bttn_box[1][1] = 0.90;
  bttn_box[1][2] = 0.14;
  bttn_box[1][3] = 0.17;

  bttn_box[2][0] = 0.70;
  bttn_box[2][1] = 0.90;
  bttn_box[2][2] = 0.10;
  bttn_box[2][3] = 0.13;

  bttn_box[3][0] = 0.65;
  bttn_box[3][1] = 0.90;
  bttn_box[3][2] = 0.28;
  bttn_box[3][3] = 0.31;

  bttn_box[4][0] = 0.65;
  bttn_box[4][1] = 0.90;
  bttn_box[4][2] = 0.24;
  bttn_box[4][3] = 0.27;

/*
-------------------
*/

  PROC_SWT = true;
  while (1) {
    for (i=0; i<ERROR_NUM; i++) {
      ERROR_FLAG[i] = false;
    }
    for (ns=0; ns<SRC_NUM; ns++) {
      wvc_ds[ns] = NULL;
      dry_ds[ns] = NULL;
      ion_ds[ns] = NULL;
    }
    obs_param_file_io(ERROR_FLAG, antenna_list_file_name,
                      &array, &ANT_NUM, &GRT_NUM, &SRT_NUM,
                      srt, &grt_elevation_limit,
                      sep_angle_limit_from_earth_limb,
                      TimUTC, &UT1_UTC, &obs_duration,
                      &SRCPROC_MODE, src, &sun, antenna_code,
                      ch_grt_el_lim, &pair_src, &ch_src,
                      ch_srt, &ch_obs_t, __READ__);
    for (iant=0; iant<GRT_NUM; iant++) {
      if (strncmp(antenna_code[iant], (ant_prm+iant)->IDC, 
                  strlen(antenna_code[iant])) == 0) {
        ant_prm[iant].UFL = true;
      }
    }
    antenna_selection(&ANT_NUM, &GRT_NUM, &SRT_NUM,
                      wave_id, grt_elevation_limit, ant_prm,
                      antenna_code, "aris_input/antenna.prm", true);
    if ((gcode_pos = (float *)calloc(GRT_NUM, sizeof(float))) == NULL) {
      printf("ERROR: antenna_visibility: Memory allocation error.\n");
      return (0);
    }

#ifdef __ANT_DEBUG__
    for (iant=0; iant<ANT_NUM; iant++) {
      printf("__ANT_DEBUG__  %d  %s  %s  %d  %lf %lf %lf\n", iant+1,
            antenna_code[iant], ant_prm[iant].IDC, ant_prm[iant].UFL,
            ant_prm[iant].XYZ[0], ant_prm[iant].XYZ[1], ant_prm[iant].XYZ[2]);
    }
    printf("__ANT_DEBUG__   %lf  %lf  %lf\n",
             src[0].s[0], src[0].s[1], src[0].s[2]);
    printf("__ANT_DEBUG__   %lf  %lf  %lf\n",
             src[1].s[0], src[1].s[1], src[1].s[2]);
    printf("__ANT_DEBUG__ HIT RETURN KEY : ");
    getchar();
#endif /* __ANT_DEBUG__ */

    for (iant=0; iant<GRT_NUM; iant++) {
      gcode_pos[iant] = 0.02 + (float)iant / (float)ANT_NUM;
    }

    if (SRT_NUM >= 1) {
      for (iant=GRT_NUM; iant<ANT_NUM; iant++) {
        scode_pos[0][iant-GRT_NUM] = 0.02 + (float)iant / (float)ANT_NUM;
        scode_pos[1][iant-GRT_NUM] = 0.01 + (float)iant / (float)ANT_NUM;
      }

      trk_num = 0;
      for (itrk=0; itrk<TRKMAX; itrk++) {
        trk_priority[itrk] = 0;
        priority[itrk]     = 0;
      }

      if ((fp = fopen("aris_input/sim.prm", "r")) != NULL) {
        while (1) {
          if (fgets(string, sizeof(string), fp) == NULL) {
            break;
          } else {
            string[strlen(string)-1] = '\0';
          }

          if (strncmp(string, "TRACKING NETWORK    ", 20) == 0) {
            trk_num = tracking_station_name_read(trk_num,  string+20,
                                                 trk_name, priority);
          }
        }
      }
      TRK_NUM = tracking_init(trk_num, priority, trk_priority,
                              trk_name, trk_pos);
    }

/*
--------
*/

/****
    obs_duration *= 3600.0;
****/
    nobs = (int)lrint(obs_duration);
    grt_elevation_limit *= DPI;

    if (pgid[0] < 0) {
      pgid[0] = cpgopen("/xs");
    }
    cpgslct(pgid[0]);
    cpgpap(1.5*pgpap_prm, 1.0);

    cpgbbuf();

    if (SRT_NUM >= 1) {
      cpgsch(1.0);
      cpgsvp(0.00, 1.00, 0.00, 1.00);
      cpgswin(0.00, 1.00, 0.00, 1.00);
      tracking_button_disp(5, 0.960, 0.03, bttn_box,
                           TRK_NUM, trk_priority, trk_name, trk_pos);
    }

    pgxmin = (float)(3600*TimUTC[3] + 60*TimUTC[4] + TimUTC[5]       );
    pgxmax = (float)(3600*TimUTC[3] + 60*TimUTC[4] + TimUTC[5] + nobs);
    pgymin = 0.0;
    pgymax = 1.0;

    cpgsvp(0.15, 0.85, 0.50, 0.87);
    cpgsch(0.8);
    cpgswin(pgxmin, pgxmax, pgymin, pgymax);
    cpgtbox("BCSTNZHI", 0, 0, "BC", 0, 0);

    cpgsch(0.65);
    cpgsci(1);
    for (iant=0; iant<GRT_NUM; iant++) {
      cpgptxt(pgxmin, gcode_pos[iant], 0.0, 1.1, ant_prm[iant].IDC);
    }
    for (iant=GRT_NUM; iant<ANT_NUM; iant++) {
      sprintf(string, "%s[+X]", ant_prm[iant].IDC);
      cpgptxt(pgxmin, scode_pos[0][iant-GRT_NUM], 0.0,  1.1, string);
      sprintf(string, "%s[-X]", ant_prm[iant].IDC);
      cpgsci(2);
      cpgptxt(pgxmax, scode_pos[1][iant-GRT_NUM], 0.0, -0.1, string);
      cpgsci(1);
    }
    cpgsch(1.0);

    if (SRT_NUM > 0) {
      inner_sun_angle
        = sun_angle_check(TimUTC, obs_duration, &src[0], &sun) / DPI;
      if (inner_sun_angle < 50.0) {
        sprintf(string,
                "ERROR: Inner Sun Angle is less than 50 deg. (%5.1lf [deg])",
                inner_sun_angle);
        printf("%s\n", string); 
/********
        if (TV_SWT == true) {
          comment_disp(cmnt, comment, string, true);
        } else {
          printf("%s\n", string); 
        }
********/
      }
    }

/*
--------
*/

    for (ns=0; ns<SRC_NUM; ns++) {
      if ((int_obs[ns] = (struct st_observable *)
                     calloc(ANT_NUM, sizeof(struct st_observable))) == NULL) {
        printf("ARIS: ERROR: memory allocation error: int_obs(1).\n");
        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        return -1;
      }
    }

    IOBS = 0;
    NOBS = 1;
    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;
    for (iant=0; iant<SRT_NUM; iant++) {
      init_l_p[iant] = dpi;
      init_l_m[iant] = dpi;
    }

    for (iobs=0; iobs<nobs; iobs+=d_t) {
      timUTC[5] = TimUTC[5] + iobs;
      ftmp1 = (float)(3600*TimUTC[3] + 60*TimUTC[4] + TimUTC[5] + iobs);

      if (iobs % 86400 < d_t) {
        if (SRT_NUM > 0) {
          inner_sun_angle
            = sun_angle_check(timUTC, obs_duration, &src[0], &sun) / DPI;
          if (inner_sun_angle < 50.0) {
            sprintf(string,
                "ERROR: Inner Sun Angle is less than 50 deg. (%5.1lf [deg])",
                inner_sun_angle);
            printf("%s\n", string); 
/********
            if (TV_SWT == true) {
              comment_disp(cmnt, comment, string, true);
            } else {
              printf("%s\n", string); 
            }
********/
          }
        }
      }

      ns = 0;
      for (iant=0; iant<ANT_NUM; iant++) {
        int_obs[ns][iant].wt        = src_flag[ns];
        int_obs[ns][iant].amp_error = 1.0;
      }
      for (iant=0; iant<SRT_NUM; iant++) {
        srt[iant].BODY_X_SUN = +1;
      }
      uvw_calc(100.0 * (float)(ns + 1),
               ANT_NUM, GRT_NUM,  NOBS, IOBS, timUTC, ut1_utc,
               ant_prm, src[0], sun,
               wvc, dry, ion, Cw, Cd, Ci, CI,
               wvc_ds[0], dry_ds[0], ion_ds[0], fqs_ds,
               false,  ERROR_FLAG, W, dUT1, int_obs[ns],
               OBS_T, OBS_p, dz, true, srt, init_l_p,
               TRK_NUM, trk_pos, srt_link, sep_angle_limit_from_earth_limb,
               MSTID);
      cpgsci(1);
      for (iant=0; iant<GRT_NUM; iant++) {
        if (int_obs[ns][iant].wt > 100.0) {
          cpgpt(1, &ftmp1, &gcode_pos[iant], 1);
        }
      }
      if (SRT_NUM >= 1) {
        cpgsci(1);
        for (iant=GRT_NUM; iant<ANT_NUM; iant++) {
          if (int_obs[ns][iant].wt > 100.0) {
            cpgpt(1, &ftmp1, &scode_pos[0][iant-GRT_NUM], 1);
          }
        }
        cpgsci(3);
        if (iobs % (4 * d_t) == 0) {
          for (iant=GRT_NUM; iant<ANT_NUM; iant++) {
            ftmp2 = vlen3(ant_prm[iant].XYZ)
                      / (srt[iant-GRT_NUM].apogee + earth_radius);
            cpgpt(1, &ftmp1, &ftmp2, 1);
          }
        }

        ns = 0;
        for (iant=0; iant<ANT_NUM; iant++) {
          int_obs[ns][iant].wt        = src_flag[ns];
          int_obs[ns][iant].amp_error = 1.0;
        }
        for (iant=0; iant<SRT_NUM; iant++) {
          srt[iant].BODY_X_SUN = -1;
        }
        uvw_calc(100.0 * (float)(ns + 1),
                 ANT_NUM, GRT_NUM,  NOBS, IOBS, timUTC, ut1_utc,
                 ant_prm, src[0], sun,
                 wvc, dry, ion, Cw, Cd, Ci, CI,
                 wvc_ds[0], dry_ds[0], ion_ds[0], fqs_ds,
                 false,  ERROR_FLAG, W, dUT1, int_obs[ns],
                 OBS_T, OBS_p, dz, true, srt, init_l_m,
                 TRK_NUM, trk_pos, srt_link, sep_angle_limit_from_earth_limb,
                 MSTID);
        cpgsci(2);
        for (iant=GRT_NUM; iant<ANT_NUM; iant++) {
          if (int_obs[ns][iant].wt > 100.0) {
            cpgpt(1, &ftmp1, &scode_pos[1][iant-GRT_NUM], 1);
          }
        }
        cpgsci(1);
      }
    }

    free (gcode_pos);
    if (SRT_NUM >= 1) {
      free (int_obs[0]);
      free (int_obs[1]);
    }

    cpgebuf();

/*
----------------------------------------
*/

    cpgsci(7);
    pgxregion(pgxmin, pgxmax, pgx, pgy);
    pgy[0] = pgymin + 0.005 * (pgymax - pgymin);
    pgy[1] = pgymax - 0.005 * (pgymax - pgymin);
    cpgsfs(4);
    cpgrect(pgx[0], pgx[1], pgy[0], pgy[1]);
    cpgsfs(1);

    sobs = (int)pgx[0];
    eobs = (int)pgx[1];
    obs_duration = lrint(eobs - sobs) + 1;
    NOBS = (int)(obs_duration / d_t);
    obs_dh = obs_duration / 3600.0;

    n = (int)log10(obs_dh);
    i = 0;
    for (j=n; j>=0; j--) {
      f = pow(10.0, (double)j);
      k = (int)floor(obs_dh / f);
      ch_obs_t.obsd[i] = '0' + k;
      obs_dh -= (double)k * f;
      i++;
    }
    ch_obs_t.obsd[i++] = '.';
    for (j=0; j<2; j++) {
      obs_dh *= 10.0;
      ch_obs_t.obsd[i++] = '0' + (int)obs_dh % 10;
    }
    ch_obs_t.obsd[i] = '\0';

    I = NOBS * ANT_NUM;
    for (ns=0; ns<SRC_NUM; ns++) {
      if ((int_obs[ns] = (struct st_observable *)
                     calloc(I, sizeof(struct st_observable))) == NULL) {
        printf("ARIS: ERROR: memory allocation error: int_obs(2).\n");
        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        return -1;
      }
    }

    for (ns=0; ns<SRC_NUM; ns++) {
      for (i=0; i<I; i++) {
        int_obs[ns][i].wt        = src_flag[ns];
        int_obs[ns][i].amp_error = 1.0;
      }
    }

    mjd = MJD(timUTC[0], timUTC[1], timUTC[2], 0, 0, sobs, 0.0);
    MJD2date(mjd, &timUTC[0], &timUTC[1], &timUTC[2],
                  &timUTC[3], &timUTC[4], &timUTC[5]);
    ch_time_set(timUTC, &ch_obs_t);
    for (iant=0; iant<SRT_NUM; iant++) {
      init_l_p[iant] = dpi;
      init_l_m[iant] = dpi;
    }

    for (iobs=0; iobs<NOBS; iobs++) {
      timUTC[5] = d_t * iobs;
      if (SRT_NUM > 0) {
        inner_sun_angle
          = sun_angle_check(timUTC, obs_duration, &src[0], &sun) / DPI;
        if (inner_sun_angle < 50.0) {
          sprintf(string,
                "ERROR: Inner Sun Angle is less than 50 deg. (%5.1lf [deg])",
                inner_sun_angle);
          printf("%s\n", string); 
/********
          if (TV_SWT == true) {
            comment_disp(cmnt, comment, string, true);
          } else {
            printf("%s\n", string); 
          }
********/
        }
      }

      ns = 0;
      for (iant=0; iant<SRT_NUM; iant++) {
        srt[iant].BODY_X_SUN = +1;
      }
      uvw_calc(100.0 * (float)(ns + 1),
               ANT_NUM, GRT_NUM,  NOBS, iobs, timUTC, ut1_utc,
               ant_prm, src[0], sun,
               wvc, dry, ion, Cw, Cd, Ci, CI,
               wvc_ds[0], dry_ds[0], ion_ds[0], fqs_ds,
               false,  ERROR_FLAG, W, dUT1, int_obs[ns],
               OBS_T, OBS_p, dz, true, srt, init_l_p,
               TRK_NUM, trk_pos, srt_link,
               sep_angle_limit_from_earth_limb, MSTID);
      ns = 1;
      for (iant=0; iant<SRT_NUM; iant++) {
        srt[iant].BODY_X_SUN = -1;
      }
      uvw_calc(100.0 * (float)(ns + 1),
               ANT_NUM, GRT_NUM,  NOBS, iobs, timUTC, ut1_utc,
               ant_prm, src[0], sun,
               wvc, dry, ion, Cw, Cd, Ci, CI,
               wvc_ds[0], dry_ds[0], ion_ds[0], fqs_ds,
               false,  ERROR_FLAG, W, dUT1, int_obs[ns],
               OBS_T, OBS_p, dz, true, srt, init_l_m,
               TRK_NUM, trk_pos, srt_link,
               sep_angle_limit_from_earth_limb, MSTID);
    }

/*
----------------------
*/

    cpgsci(1);
    cpgsls(1);
    cpgsch(0.65);
    cpgsvp(0.10, 0.40, 0.10, 0.40);
    for (ns=0; ns<SRC_NUM; ns++) {
      uv_max[ns] = 0.0;
    }
    uv_display(NOBS, uv_max, ANT_NUM, 0, ANT_NUM, 0, ANT_NUM, 
               timUTC, ut1_utc, ant_prm, src, src_flag,
               int_obs, wave_length, 0, false, 0.65, &c_dummy, F__STRUCTURE);
    if (SRT_NUM >= 1) {
      free (int_obs[0]);
      free (int_obs[1]);
    }

    cpgsvp(0.00, 1.00, 0.00, 1.00);
    pgxmin = 0.0;
    pgxmax = 1.0;
    pgymin = 0.0;
    pgymax = 1.0;
    cpgsch(0.65);
    cpgswin(pgxmin, pgxmax, pgymin, pgymax);

    cpgsci(1);
    sprintf(tv_time, "%4s/%2s/%2s   %2s:%2s:%2s",
            ch_obs_t.start_t[0], ch_obs_t.start_t[1], ch_obs_t.start_t[2],
            ch_obs_t.start_t[3], ch_obs_t.start_t[4], ch_obs_t.start_t[5]);
    tv_button_disp(bttn_box[3], tv_time);
    sprintf(tv_time, "%s hours", ch_obs_t.obsd);
    tv_button_disp(bttn_box[4], tv_time);

    off_button(&idum, "Continue\0",  bttn_box[0]);
    off_button(&idum, "Main Menu\0", bttn_box[1]);
    off_button(&idum, "Exit\0",      bttn_box[2]);

    while (1) {
      cpgcurs(cursor_pos, cursor_pos+1, string);
      if (_button_chk(cursor_pos, bttn_box[0]) == true) {
        on_button(&idum, "Continue\0", bttn_box[0]);
        PROC_SWT = true;
        break;
      }
      if (_button_chk(cursor_pos, bttn_box[1]) == true) {
        on_button(&idum, "Main Menu\0",    bttn_box[1]);
        PROC_SWT = false;
        break;
      }
      if (_button_chk(cursor_pos, bttn_box[2]) == true) {
        on_button(&idum, "Exit\0",     bttn_box[2]);
        PROC_SWT = false;
        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        return 0;
      }
    }

    if (PROC_SWT == false) {
      break;
    }
  }

/*
------------------------------
*/

  sscanf(ch_obs_t.start_t[0], "%d", TimUTC);
  sscanf(ch_obs_t.start_t[1], "%d", TimUTC+1);
  sscanf(ch_obs_t.start_t[2], "%d", TimUTC+2);
  sscanf(ch_obs_t.start_t[3], "%d", TimUTC+3);
  sscanf(ch_obs_t.start_t[4], "%d", TimUTC+4);
  sscanf(ch_obs_t.start_t[5], "%d", TimUTC+5);

/*
------------------------------
*/

  grt_num = 0;
  for (iant=0; iant<GRT_NUM; iant++) {
    if (ant_prm[iant].UFL == true) {
      sprintf(antenna_code[grt_num], "%s", (ant_prm+iant)->IDC);
      grt_num++;
    }
  }
  GRT_NUM = grt_num;
  obs_param_file_io(ERROR_FLAG, antenna_list_file_name,
                    &array, &ANT_NUM, &GRT_NUM, &SRT_NUM,
                    srt, &grt_elevation_limit,
                    sep_angle_limit_from_earth_limb,
                    TimUTC, &UT1_UTC, &obs_duration,
                    &SRCPROC_MODE, src, &sun, antenna_code,
                    ch_grt_el_lim, &pair_src, &ch_src,
                    ch_srt, &ch_obs_t, _WRITE__);
  if (SRT_NUM >= 1) {
    free (srt_link);
  }

  return 1;
}
