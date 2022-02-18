#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __FILE_OUT__
****/

void link_TRK_sim(int SRT_NUM, int  nobs,
               struct srt_orbit_parameter *srt,
               struct antenna_parameter *trk_pos,
               struct source_parameter src,
               struct source_parameter sun,
               struct srt_data_link *srt_link,
               int *TimUTC, double UT1_UTC,
               double OBS_T, double OBS_p,
               double sep_angle_limit,
               int  TRK_NUM, int  DISP_SWT)
{
  int    i, I, J, K, II, iant;
  int    iobs, NOBS;
  int    itrk, ERR_FLG;
  int    ntrk;
  _Bool  TRK_SWT;
  float  *pgx, *pg1, *pg2, *pg3;
  float  xmin, xmax, ymin, ymax;
  _Bool  *status[3];
  double P[SRTMAX][3], E[3], V_xyz[SRTMAX][3], Oe[3];
  double ptrk[3], p[3];
  int    timUTC[6];
  double ut1_utc;
  double AZ, EL, dAZdt, dELdt;
  double srt_AZ, srt_EL, srt_dAZdt, srt_dELdt;
  float  DPI;
  double R, L;
  float  WEIGHT_FLAG_DATA;
  double q[4];
  double init_l[SRTMAX];

#ifdef __FILE_OUT__
  FILE   *log_fp;
#endif

/*
--------------------------------------------------
*/

  DPI = (float)(dpi / 180.0);

  for (i=0; i<6; i++) {
    timUTC[i] = TimUTC[i];
  }
  ut1_utc = UT1_UTC;

  NOBS = nobs / TIME_STEP;

  pgx       = calloc(NOBS, sizeof(float));
  pg1       = calloc(NOBS, sizeof(float));
  pg2       = calloc(NOBS, sizeof(float));
  pg3       = calloc(NOBS, sizeof(float));
  status[0] = calloc(NOBS, sizeof(_Bool));
  status[1] = calloc(NOBS, sizeof(_Bool));
  status[2] = calloc(NOBS, sizeof(_Bool));

  cpgpap(1.0*pgpap_prm, 1.0);
  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);

/*
------------------------------------------
*/

  xmin = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               TimUTC[3], TimUTC[4], TimUTC[5], UT1_UTC)
                         - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,         0,         0,         UT1_UTC));
  xmax = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               TimUTC[3], TimUTC[4], TimUTC[5]+nobs, UT1_UTC)
                         - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,         0,         0,         UT1_UTC));
  cpgsch(0.5);
  cpgsci(1);

  II = 0;
  K = 0;
  for (itrk=0; itrk<TRK_NUM; itrk++) {
    for (i=0; i<3; i++) {
      for (iobs=0; iobs<NOBS; iobs++) {
        status[i][iobs] = true;
      }
    }
    J = 0;
    for (iant=0; iant<SRT_NUM; iant++) {
      init_l[iant] = dpi;
      for (iobs=0; iobs<NOBS; iobs++) {
        timUTC[5] = TimUTC[5] + iobs * TIME_STEP;
        pgx[iobs] = 86400.0 * (float)(MJD(timUTC[0], timUTC[1], timUTC[2],
                                          timUTC[3], timUTC[4], timUTC[5],
                                          ut1_utc)
                                    - MJD(timUTC[0], timUTC[1], timUTC[2],
                                                  0,         0,         0,
                                          ut1_utc));
        spacecraft_position(srt[iant], timUTC, ut1_utc, (double)TIME_STEP,
                            P[iant], E, Oe, V_xyz[iant], init_l+iant);
        ERR_FLG = tracking_condition(timUTC, ut1_utc, trk_pos[itrk],
                                     OBS_T, OBS_p, false, P[iant], V_xyz[iant],
                                     &AZ, &EL, &dAZdt, &dELdt,
                                     srt_link, src, sun,
                                     &srt_AZ, &srt_EL,
                                     &srt_dAZdt, &srt_dELdt,
                                     ptrk, p, srt[iant].BODY_X_SUN);
        if (ERR_FLG / 4 == 1) {
          status[0][iobs] = false;
          ERR_FLG %= 4;
        }
        if (ERR_FLG / 2 == 1) {
          status[1][iobs] = false;
          ERR_FLG %= 2;
        }
        if (ERR_FLG == 1) {
          status[2][iobs] = false;
        }

        pg3[iobs] = (float)vlen3(p) * 1.0e-7;
        if (DISP_SWT == TRACKING) {
          pg1[iobs] = (float)(    AZ / DPI);
          pg2[iobs] = (float)(    EL / DPI);
        } else if (DISP_SWT == ONBOARD) {
          pg1[iobs] = (float)(srt_AZ / DPI);
          pg2[iobs] = (float)(srt_EL / DPI);
        }
      }
    }

/*
---- Distance ----
*/

    cpgsvp(0.05+0.50*(float)K, 0.45+0.50*(float)K,
           0.92-0.05*(float)J-0.17*(float)II,
           0.97-0.05*(float)J-0.17*(float)II);
    ymin = 0.01;
    ymax = 4.10;
    cpgswin(xmin, xmax, ymin, ymax);
    cpgtbox("BCSTZHI", 0, 0, "BCNTS", 4.0, 8);
    for (iobs=0; iobs<NOBS; iobs++) {
      cpgsci(1);
      if (DISP_SWT == TRACKING && status[0][iobs] == true &&
         status[1][iobs] == true) {
        cpgsci(2);
      } else if (DISP_SWT == ONBOARD && status[2][iobs] == true) {
        cpgsci(2);
      }
      cpgpt(1, pgx+iobs, pg3+iobs, 1);
    }
    cpgsci(1);
    cpgtext(0.90*xmin+0.10*xmax, 0.30*ymin+0.70*ymax, trk_pos[itrk].IDC);
    cpgtext(0.70*xmin+0.30*xmax, 0.30*ymin+0.70*ymax, "(DIST)");
    J++;

/*
---- AZ ----
*/

    cpgsvp(0.05+0.50*(float)K, 0.45+0.50*(float)K,
           0.92-0.05*(float)J-0.17*(float)II,
           0.97-0.05*(float)J-0.17*(float)II);
    ymin = -179.9;
    ymax =  180.0;
    cpgswin(xmin, xmax, ymin, ymax);
    cpgtbox("BCSTZHI", 0, 0, "BCNTS", 180, 2);
    for (iobs=0; iobs<NOBS; iobs++) {
      cpgsci(1);
      if (DISP_SWT == TRACKING && status[0][iobs] == true &&
          status[1][iobs] == true) {
        cpgsci(2);
      } else if (DISP_SWT == ONBOARD && status[2][iobs] == true) {
        cpgsci(2);
      }
      cpgpt(1, pgx+iobs, pg1+iobs, 1);
    }
    cpgsci(1);
    cpgtext(0.70*xmin+0.30*xmax, 0.30*ymin+0.70*ymax, "(AZ)");
    J++;

/*
---- EL ----
*/

    cpgsvp(0.05+0.50*(float)K, 0.45+0.50*(float)K,
           0.92-0.05*(float)J-0.17*(float)II,
           0.97-0.05*(float)J-0.17*(float)II);
    if (DISP_SWT == TRACKING) {
      ymin = 0.0;
    } else if (DISP_SWT == ONBOARD) {
      ymin = -30.0;
    }
    ymax = 90.0;
    cpgswin(xmin, xmax, ymin, ymax);
    if (II == 3 || itrk == TRK_NUM-1) {
      cpgtbox("BCSTNZHI", 0, 0, "BCNTS", 60, 3);
    } else {
      cpgtbox("BCSTZHI", 0, 0, "BCNTS", 60, 3);
    }
    for (iobs=0; iobs<NOBS; iobs++) {
      cpgsci(1);
      if (DISP_SWT == TRACKING && status[0][iobs] == true &&
          status[1][iobs] == true) {
        cpgsci(2);
      } else if (DISP_SWT == ONBOARD && status[2][iobs] == true) {
        cpgsci(2);
      }
      cpgpt(1, pgx+iobs, pg2+iobs, 1);
    }
    cpgsci(1);
    cpgtext(0.70*xmin+0.30*xmax, 0.30*ymin+0.70*ymax, "(EL)");
    J++;

/*
------------
*/

    II++;
    if (II == 4) {
      II = 0;
      K = 1;
    }
  }

/*
----------------------------------------
*/

#ifdef __FILE_OUT__

  if ((log_fp = fopen("aris_log/link_TRK_status.dat", "w")) == NULL) {
    printf("link_TRK_sim: warning; ./aris_log/link_TRK_status.dat ");
    printf("cannot open.\n");
  } else {
    fprintf(log_fp, "----HEADER----\n");
    fprintf(log_fp, "OBSERVATION EPOCH  : %4d.%2d.%2d %2d:%2d:%2d\n",
            TimUTC[0], TimUTC[1], TimUTC[2], TimUTC[3], TimUTC[4], TimUTC[5]);
    for (iant=0; iant<SSRT_NUM; iant++) {
      attitude_Q(srt[iant].BODY_X_SUN, src, sun, q);
      fprintf(log_fp, "SATELLITE[%d] ATTITUDE: %lf  %lf  %lf  %lf\n",
            iant+1, q[0], q[1], q[2], q[3]);
    }
    orbit_info_print(log_fp, SRT_NUM, srt, TimUTC, UT1_UTC);

/*
----
*/

    fprintf(log_fp, "----DATA------\n");
    fprintf(log_fp, "  UTC   D/N      SATELLITE POSITON      ");
    for (itrk=0; itrk<TRK_NUM; itrk++) {
      fprintf(log_fp, "%8s            ", trk_pos[itrk].IDC);
    }
    fprintf(log_fp, "\n");

    fprintf(log_fp, "              X        Y        Z     ");
    for (itrk=0; itrk<TRK_NUM; itrk++) {
      fprintf(log_fp, "  DIST     AZ     EL");
    }
    fprintf(log_fp, "\n");

    fprintf(log_fp, "   [s]        [km]     [km]     [km]  ");
    for (itrk=0; itrk<TRK_NUM; itrk++) {
      fprintf(log_fp, "  [km]  [deg]  [deg]");
    }
    fprintf(log_fp, "\n");
    fprintf(log_fp, "--------------\n");

    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;

    for (iobs=0; iobs<NOBS; iobs++) {
      timUTC[5] = TimUTC[5] + iobs * TIME_STEP;
      fprintf(log_fp, "%6d  ",
          (int)lrint(
            86400.0 * (float)(MJD(timUTC[0], timUTC[1], timUTC[2],
                                  timUTC[3], timUTC[4], timUTC[5], ut1_utc)
                            - MJD(timUTC[0], timUTC[1], timUTC[2],
                                  0,         0,         0,         ut1_utc))));
      for (iant=0; iant<SRT_NUM; iant++) {
        spacecraft_position(srt[iant], timUTC, ut1_utc, (double)TIME_STEP,
                            P[iant], E, Oe, V_xyz[iant], init_l+iant);
      }
      for (itrk=0; itrk<TRK_NUM; itrk++) {
        for (iant=0; iant<SRT_NUM; iant++) {
          ERR_FLG = tracking_condition(timUTC, ut1_utc, trk_pos[itrk],
                                       OBS_T, OBS_p, false, P[iant], V_xyz[iant],
                                       &AZ, &EL, &dAZdt, &dELdt,
                                       srt_link, src, sun,
                                       &srt_AZ, &srt_EL,
                                       &srt_dAZdt, &srt_dELdt,
                                       ptrk, p, srt[iant].BODY_X_SUN);
          if (itrk == 0) {
            fprintf(log_fp, "%1d  %8.1f %8.1f %8.1f  ",
                    earth_eclipse(sun.s, P[iant], &R, &L, sep_angle_limit),
                    P[iant][0]*1.0e-3, P[iant][1]*1.0e-3, P[iant][2]*1.0e-3);
          }

          fprintf(log_fp, "%5.0f ", (float)vlen3(p) * 1.0e-3);
          if (DISP_SWT == TRACKING) {
            fprintf(log_fp, "%6.1f %6.1f ",
                  (float)(    AZ / DPI), (float)(    EL / DPI));
          } else if (DISP_SWT == ONBOARD) {
            fprintf(log_fp, "%6.1f %6.1f ",
                  (float)(srt_AZ / DPI), (float)(srt_EL / DPI));
          }
        }
      }
      fprintf(log_fp, "\n");
    }

    fclose (log_fp);
  }

#endif /* __FILE_OUT__ */

/*
----------------------------------------
*/

#ifdef __FILE_OUT__

  if ((log_fp = fopen("aris_log/link_info.dat", "w")) == NULL) {
    printf("link_TRK_sim: warning; ./aris_log/link_info.dat cannot open.\n");
  } else {
    fprintf(log_fp, "----HEADER----\n");
    fprintf(log_fp, "OBSERVATION EPOCH  : %4d.%2d.%2d %2d:%2d:%2d\n",
            TimUTC[0], TimUTC[1], TimUTC[2], TimUTC[3], TimUTC[4], TimUTC[5]);
    for (iant=0; iant<SRT_NUM; iant++) {
      attitude_Q(srt[iant].BODY_X_SUN, src, sun, q);
      fprintf(log_fp, "SATELLITE[%d] ATTITUDE: %lf  %lf  %lf  %lf\n",
              iant+1, q[0], q[1], q[2], q[3]);
    }
    orbit_info_print(log_fp, SRT_NUM, srt, TimUTC, UT1_UTC);

/*
----
*/

    fprintf(log_fp, "----DATA------\n");
    fprintf(log_fp, "  UTC   D/N      SATELLITE POSITON      ");
    fprintf(log_fp, "\n");
  
    fprintf(log_fp, "              X        Y        Z     ");
    fprintf(log_fp, "  DIST             AZ     EL  srtAZ  srtEL");
    fprintf(log_fp, "\n");

    fprintf(log_fp, "   [s]        [km]     [km]     [km]  ");
    fprintf(log_fp, "  [km]          [deg]  [deg]  [deg]  [deg]");
    fprintf(log_fp, "\n");
    fprintf(log_fp, "--------------\n");

    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;

    for (iobs=0; iobs<NOBS; iobs++) {
      timUTC[5] = TimUTC[5] + iobs * TIME_STEP;
      fprintf(log_fp, "%6d  ",
          (int)lrint(
            86400.0 * (float)(MJD(timUTC[0], timUTC[1], timUTC[2],
                                  timUTC[3], timUTC[4], timUTC[5], ut1_utc)
                            - MJD(timUTC[0], timUTC[1], timUTC[2],
                                  0,         0,         0,         ut1_utc))));
      for (iant=0; iant<SRT_NUM; iant++) {
        TRK_SWT = tracking_status(true, &WEIGHT_FLAG_DATA, &ntrk,
                        timUTC, ut1_utc, srt[iant], src, sun,
                        OBS_T,  OBS_p, TRK_NUM, P[iant], E, V_xyz[iant],
                        trk_pos, srt_link,
                        &AZ, &EL, &dAZdt,  &dELdt,
                        &srt_AZ, &srt_EL,
                        sep_angle_limit, ptrk, p, init_l+iant);

        fprintf(log_fp, "%1d  %8.1f %8.1f %8.1f ",
                earth_eclipse(sun.s, P[iant], &R, &L, sep_angle_limit),
                P[iant][0]*1.0e-3, P[iant][1]*1.0e-3, P[iant][2]*1.0e-3);
        fprintf(log_fp, "%5.0f ", (float)vlen3(p) * 1.0e-3);
        fprintf(log_fp, "%8s ", trk_pos[ntrk].IDC);

        if (TRK_SWT == true) {
          fprintf(log_fp, "%6.1f %6.1f ",
                (float)(    AZ / DPI), (float)(    EL / DPI));
          fprintf(log_fp, "%6.1f %6.1f",
                (float)(srt_AZ / DPI), (float)(srt_EL / DPI));
        } else {
          fprintf(log_fp, "%6.1f %6.1f ", 0.0, 0.0);
          fprintf(log_fp, "%6.1f %6.1f",  0.0, 0.0);
        }
        fprintf(log_fp, "\n");
      }
    }
    fclose (log_fp);
  }

#endif /* __FILE_OUT__ */

/*
----------------------------------------
*/

  free (pgx);
  free (pg1);
  free (pg2);
  free (pg3);
  free (status[0]);
  free (status[1]);
  free (status[2]);

  return;
}
