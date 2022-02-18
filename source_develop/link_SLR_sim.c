#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __ELEVATION_DEF__
****/
#ifdef __ELEVATION_DEF__
  #define __LOW_EL_DATA_PLOT__
#endif
#ifndef __ELEVATION_DEF__
  #define __DAY_NIGHT_DEF__
#endif


/****
#define __FILE_OUT__
****/

void link_SLR_sim( int SRT_NUM, int  nobs,
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
  int    i, j, I, J, iant;
  int    NOBS;
  int    itrk, ntrk, ERR_FLG;
  float  *pgx, *pg1, *pg2;
  float  xmin, xmax, ymin, ymax;
  double P[3], E[3], V_xyz[3], Oe[3];
  double ptrk[3], p[3], pptrk[3], pp[3];
  int    timUTC[6];
  double ut1_utc;
  double AZ, EL, dAZdt, dELdt;
  double DPI;
  double srt_AZ, srt_EL, srt_dAZdt, srt_dELdt;
  double R, L;
  double D;
  int    SLR_NUM, islr;
  struct antenna_parameter slr_pos[ANTMAX];
  float  WEIGHT_FLAG_DATA;
  _Bool  TRK_SWT;

  int    ncnt[3];
  float  *slr_prm[6][3];
  float  nslr[3][30], dist[30];

  double alpha, delta;
  double AZ_sun, EL_sun, dAZdt_sun, dELdt_sun;
  double vtmp[3];
  double init_l[SRTMAX];


#ifdef __FILE_OUT__
  FILE   *fp;
#endif

/*
--------------------------------------------------
*/

  DPI = dpi / 180.0;

  for (i=0; i<6; i++) {
    timUTC[i] = TimUTC[i];
  }
  ut1_utc = UT1_UTC;

  NOBS = nobs / TIME_STEP;
  SLR_NUM = slr_config(slr_pos);

  cpgpap(1.2*pgpap_prm, 1.2);
  cpgsch(0.5);
  cpgsci(1);

/*
------------------------------------------
*/

  if (DISP_SWT == 0) {

    pgx = (float *)calloc(NOBS, sizeof(float));
    pg1 = (float *)calloc(NOBS, sizeof(float));
    pg2 = (float *)calloc(NOBS, sizeof(float));

    xmin = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                                 TimUTC[3], TimUTC[4], TimUTC[5], UT1_UTC)
                           - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,         0,         0,         UT1_UTC));
    xmax = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                                 TimUTC[3], TimUTC[4], TimUTC[5]+nobs, UT1_UTC)
                           - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                                 0,         0,         0,         UT1_UTC));

/*
----------------------------------------
*/

    cpgsci(1);
    cpgsvp(0.10, 0.90, 0.80, 0.90);
    ymin = 0.00;
    ymax = 4.00;
    cpgswin(xmin, xmax, ymin, ymax);
    cpgtbox("BCSNTZHI", 0, 0, "BCNTS", 0, 0);
    cpglab("", "Altitude [km]/10000", "");

    for (iant=0; iant<SRT_NUM; iant++) {
      init_l[iant] = dpi;
      for (i=0; i<NOBS; i++) {
        timUTC[5] = TimUTC[5] + i * TIME_STEP;
        pgx[i] = 86400.0
                   * (float)(MJD(timUTC[0], timUTC[1], timUTC[2],
                                 timUTC[3], timUTC[4], timUTC[5], ut1_utc)
                           - MJD(timUTC[0], timUTC[1], timUTC[2],
                                 0,         0,         0,         ut1_utc));
        spacecraft_position(srt[iant], timUTC, ut1_utc, (double)TIME_STEP,
                            P, E, Oe, V_xyz, init_l+iant);
        pg1[i] = (float)vlen3(P) * 1.0e-7;
      }
      cpgpt(NOBS, pgx, pg1, 1);
    }

/*
----------------------------------------
*/

    cpgsci(1);
    cpgsvp(0.10, 0.90, 0.65, 0.75);
    ymin = 0.00;
    ymax = 90.00;
    cpgswin(xmin, xmax, ymin, ymax);
    cpgtbox("BCSNTZHI", 0, 0, "BCNTS", 0, 0);
    cpglab("", "Ka Trk EL [deg]", "");

    for (iant=0; iant<SRT_NUM; iant++) {
      init_l[iant] = dpi;
    }
    for (i=0; i<NOBS; i++) {
      timUTC[5] = TimUTC[5] + i * TIME_STEP;
      for (iant=0; iant<SRT_NUM; iant++) {
        TRK_SWT = tracking_status(true, &WEIGHT_FLAG_DATA, &ntrk,
                        timUTC, ut1_utc, srt[iant], src, sun,
                        OBS_T,  OBS_p, TRK_NUM, P, E, V_xyz,
                        trk_pos, srt_link,
                        &AZ, &EL, &dAZdt,  &dELdt,
                        &srt_AZ, &srt_EL,
                        sep_angle_limit, ptrk, p, init_l+iant);

        if (TRK_SWT == true) {
          pg1[i] = (float)( EL / DPI);
          cpgsci(ntrk+1);
          cpgpt(1, pgx+i, pg1+i, ntrk+2);
          cpgsci(1);
        }
      }
    }

/*
----------------------------------------
*/

    cpgsci(1);
    cpgsvp(0.10, 0.90, 0.50, 0.60);
    ymin = 0.00;
    ymax = 4.00;
    cpgswin(xmin, xmax, ymin, ymax);
    cpgtbox("BCSNTZHI", 0, 0, "BCNTS", 0, 0);
    cpglab("", "SLR Dis [km]/10000", "");

    for (i=0; i<NOBS; i++) {
      timUTC[5] = TimUTC[5] + i * TIME_STEP;
      init_l[iant] = dpi;
      for (iant=0; iant<SRT_NUM; iant++) {
        TRK_SWT = tracking_status(true, &WEIGHT_FLAG_DATA, &ntrk,
                        timUTC, ut1_utc, srt[iant], src, sun,
                        OBS_T,  OBS_p, TRK_NUM, P, E, V_xyz,
                        trk_pos, srt_link,
                        &AZ, &EL, &dAZdt,  &dELdt,
                        &srt_AZ, &srt_EL,
                        sep_angle_limit, ptrk, p, init_l+iant);

        if (TRK_SWT == true) {
          for (islr=0; islr<SLR_NUM; islr++) {
            ERR_FLG = tracking_condition(timUTC, ut1_utc, slr_pos[islr],
                                         OBS_T, OBS_p, false, P, V_xyz,
                                         &AZ, &EL, &dAZdt, &dELdt,
                                         srt_link, src, sun,
                                         &srt_AZ, &srt_EL,
                                         &srt_dAZdt, &srt_dELdt,
                                         pptrk, pp, srt[iant].BODY_X_SUN);
            pg1[i] = (float)vlen3(pp) * 1.0e-7;
            if (EL > 0.0) {
              cpgsci(islr%8+1);
              cpgpt(1, pgx+i, pg1+i, islr+2);
              cpgsci(1);
            }
          }
        }
      }
    }

/*
----------------------------------------
*/

    cpgsci(1);
    cpgsvp(0.10, 0.90, 0.35, 0.45);
    ymin = 0.00;
    ymax = 90.00;
    cpgswin(xmin, xmax, ymin, ymax);
    cpgtbox("BCSNTZHI", 0, 0, "BCNTS", 0, 0);
    cpglab("", "SLR EL [deg]", "");

    for (iant=0; iant<SRT_NUM; iant++) {
      init_l[iant] = dpi;
    }
    for (i=0; i<NOBS; i++) {
      timUTC[5] = TimUTC[5] + i * TIME_STEP;
      for (iant=0; iant<SRT_NUM; iant++) {
        TRK_SWT = tracking_status(true, &WEIGHT_FLAG_DATA, &ntrk,
                        timUTC, ut1_utc, srt[iant], src, sun,
                        OBS_T,  OBS_p, TRK_NUM, P, E, V_xyz,
                        trk_pos, srt_link,
                        &AZ, &EL, &dAZdt,  &dELdt,
                        &srt_AZ, &srt_EL,
                        sep_angle_limit, ptrk, p, init_l+iant);

        if (TRK_SWT == true) {
          for (islr=0; islr<SLR_NUM; islr++) {
            ERR_FLG = tracking_condition(timUTC, ut1_utc, slr_pos[islr],
                                         OBS_T, OBS_p, false, P, V_xyz,
                                         &AZ, &EL, &dAZdt, &dELdt,
                                         srt_link, src, sun,
                                         &srt_AZ, &srt_EL,
                                         &srt_dAZdt, &srt_dELdt,
                                         pptrk, pp, srt[iant].BODY_X_SUN);

            pg1[i] = (float)(    EL / DPI);
            if (EL > 0.0) {
              cpgsci(islr%8+1);
              cpgpt(1, pgx+i, pg1+i, islr+2);
              cpgsci(1);
            }
          }
        }
      }
    }

/*
----------------------------------------
*/

    cpgsci(1);
    cpgsvp(0.10, 0.90, 0.20, 0.30);
    ymin =  -5.00;
    ymax =  90.00;
    cpgswin(xmin, xmax, ymin, ymax);
    cpgtbox("BCSNTZHI", 0, 0, "BCNTS", 0, 0);
    cpglab("", "SLR acos [deg]", "");

    for (iant=0; iant<SRT_NUM; iant++) {
      init_l[iant] = dpi;
    }
    for (i=0; i<NOBS; i++) {
      timUTC[5] = TimUTC[5] + i * TIME_STEP;
      for (iant=0; iant<SRT_NUM; iant++) {
        TRK_SWT = tracking_status(true, &WEIGHT_FLAG_DATA, &ntrk,
                        timUTC, ut1_utc, srt[iant], src, sun,
                        OBS_T,  OBS_p, TRK_NUM, P, E, V_xyz,
                        trk_pos, srt_link,
                        &AZ, &EL, &dAZdt,  &dELdt,
                        &srt_AZ, &srt_EL,
                        sep_angle_limit, ptrk, p, init_l+iant);

        if (TRK_SWT == true) {
          for (islr=0; islr<SLR_NUM; islr++) {
            ERR_FLG = tracking_condition(timUTC, ut1_utc, slr_pos[islr],
                                         OBS_T, OBS_p, false, P, V_xyz,
                                         &AZ, &EL, &dAZdt, &dELdt,
                                         srt_link, src, sun,
                                         &srt_AZ, &srt_EL,
                                         &srt_dAZdt, &srt_dELdt,
                                         pptrk, pp, srt[iant].BODY_X_SUN);
            pg1[i] = acos((p[0]*pp[0] + p[1]*pp[1] + p[2]*pp[2])
                                 / vlen3(p) / vlen3(pp)) / DPI;
            if (EL > 0.0) {
              cpgsci(islr%8+1);
              cpgpt(1, pgx+i, pg1+i, islr+2);
              cpgsci(1);
            }
          }
        }
      }
    }
    free (pgx);
    free (pg1);
    free (pg2);
  }

/*
----------------------------------------
*/

  if (DISP_SWT == 1) {

#ifdef __FILE_OUT__
    fp = fopen("slr_trk.dat", "w");
    fprintf(fp, "----HEADER----\n");
    fprintf(fp, "OBSERVATION EPOCH  : %4d.%2d.%2d %2d:%2d:%2d\n",
            TimUTC[0], TimUTC[1], TimUTC[2], TimUTC[3], TimUTC[4], TimUTC[5]);
    orbit_info_print(fp, SRT_NUM, srt, TimUTC, UT1_UTC);
    fprintf(fp, "----\n");
    fprintf(fp, "%2d  : THE NUMBER OF TRACKING STATIONS\n", TRK_NUM);
    fprintf(fp, "KATRK : ID ST_NAME       X[m]          Y[m]          Z[m]");
    fprintf(fp, "     LO[deg] LA[deg]    H[m]\n");
    for (i=0; i<TRK_NUM; i++) {
      fprintf(fp,
              "KATRK : %2d %8s %13.4lf %13.4lf %13.4lf %7.2lf %7.2lf %9.4lf\n",
                  i+1, trk_pos[i].IDC,
                  trk_pos[i].XYZ[0],
                  trk_pos[i].XYZ[1],
                  trk_pos[i].XYZ[2],
                  -1.0 * trk_pos[i].LLH[0] / DPI,
                  trk_pos[i].LLH[1] / DPI,
                  trk_pos[i].LLH[2]);
    }
    fprintf(fp, "----\n");
    fprintf(fp, "%2d  : THE NUMBER OF SLR STATIONS\n", SLR_NUM);
    fprintf(fp, "SLRST : ID ST_NAME       X[m]          Y[m]          Z[m]");
    fprintf(fp, "     LO[deg] LA[deg]    H[m]\n");
    for (i=0; i<SLR_NUM; i++) {
      fprintf(fp,
              "SLRST : %2d %8s %13.4lf %13.4lf %13.4lf %7.2lf %7.2lf %9.4lf\n",
                  i+1, slr_pos[i].IDC,
                  slr_pos[i].XYZ[0],
                  slr_pos[i].XYZ[1],
                  slr_pos[i].XYZ[2],
                  -1.0 * slr_pos[i].LLH[0] / DPI,
                  slr_pos[i].LLH[1] / DPI,
                  slr_pos[i].LLH[2]);
    }
    fprintf(fp, "----DATA----\n");

#endif

    for (I=0; I<6; I++) {
      for (J=0; J<3; J++) {
        slr_prm[I][J] = (float *)calloc(SLR_NUM*NOBS, sizeof(float));
      }
    }

    for (I=0; I<3; I++) {
      ncnt[I] = 0;
    }

/*
----------------------------------------
*/

    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;
    for (iant=0; iant<SRT_NUM; iant++) {
      init_l[iant] = dpi;
    }

    for (i=0; i<NOBS; i++) {
      timUTC[5] = TimUTC[5] + i * TIME_STEP;
      for (iant=0; iant<SRT_NUM; iant++) {
        TRK_SWT = tracking_status(true, &WEIGHT_FLAG_DATA, &ntrk,
                        timUTC, ut1_utc, srt[iant], src, sun,
                        OBS_T,  OBS_p, TRK_NUM, P, E, V_xyz,
                        trk_pos, srt_link,
                        &AZ, &EL, &dAZdt,  &dELdt,
                        &srt_AZ, &srt_EL,
                        sep_angle_limit, ptrk, p, init_l+iant);
        D = vlen3(P);
        if (TRK_SWT == true) {
#ifdef __FILE_OUT__
          fprintf(fp, "%5d,%8.2f,%8s,",
                  3600*timUTC[3] + 60*timUTC[4] + timUTC[5],
                  D * 1.0e-3, trk_pos[ntrk].IDC);
#endif
          for (islr=0; islr<SLR_NUM; islr++) {
            ERR_FLG = tracking_condition(timUTC, ut1_utc, slr_pos[islr],
                                         OBS_T, OBS_p, false, P, V_xyz,
                                         &AZ, &EL, &dAZdt, &dELdt,
                                         srt_link, src, sun,
                                         &srt_AZ, &srt_EL,
                                         &srt_dAZdt, &srt_dELdt,
                                         pptrk, pp, srt[iant].BODY_X_SUN);
            xyz2radec_rad(sun.s, &alpha, &delta);
            azel_position(timUTC, ut1_utc,
                          slr_pos[islr].LLH[0],
                          slr_pos[islr].LLH[1],
                          slr_pos[islr].LLH[2],
                          OBS_T, OBS_p, false, alpha, delta,
                          &AZ_sun, &EL_sun, &dAZdt_sun, &dELdt_sun,
                          0.0, vtmp);

            if (EL > 0.0) {
              if (D <= 15.0e6) {
                I = 0;
              } else if (D > 15.0e6 && D <= 25.0e6) {
                I = 1;
              } else if (D > 25.0e6) {
                I = 2;
              }

              slr_prm[0][I][ncnt[I]] = (float)(EL / DPI);
              slr_prm[1][I][ncnt[I]] = (float)vlen3(pp);
              slr_prm[2][I][ncnt[I]] = (float)sqrt(pow(ptrk[0]-pptrk[0], 2.0)
                                                 + pow(ptrk[1]-pptrk[1], 2.0)
                                                 + pow(ptrk[2]-pptrk[2], 2.0));
              slr_prm[3][I][ncnt[I]] = (float)islr + 0.5;
              slr_prm[4][I][ncnt[I]] = acos(
                                       (p[0]*pp[0] + p[1]*pp[1] + p[2]*pp[2])
                                       / vlen3(p) / vlen3(pp)) / DPI;
              slr_prm[5][I][ncnt[I]] = (float)(EL_sun / DPI);

#ifdef __FILE_OUT__
              if (islr == SLR_NUM - 1) {
                fprintf(fp, "%6.2f,%6.2f,%8.2f,%6.2f\n",
                        slr_prm[0][I][ncnt[I]],
                        slr_prm[4][I][ncnt[I]],
                        slr_prm[1][I][ncnt[I]]*1.0e-3,
                        slr_prm[5][I][ncnt[I]]);
              } else {
                fprintf(fp, "%6.2f,%6.2f,%8.2f,%6.2f,",
                        slr_prm[0][I][ncnt[I]],
                        slr_prm[4][I][ncnt[I]],
                        slr_prm[1][I][ncnt[I]]*1.0e-3,
                        slr_prm[5][I][ncnt[I]]);
              }
#endif

              slr_prm[1][I][ncnt[I]] *= 1.0e-6;
              slr_prm[2][I][ncnt[I]] *= 1.0e-6;

              ncnt[I]++;
            } else {
#ifdef __FILE_OUT__
              if (islr == SLR_NUM - 1) {
                fprintf(fp, "%6.2f,%6.2f,%8.2f,%6.2f\n", 0.0, 0.0, 0.0, 0.0);
              } else {
                fprintf(fp, "%6.2f,%6.2f,%8.2f,%6.2f,",  0.0, 0.0, 0.0, 0.0);
              }
#endif
            }
          }
        }
      }
    }

#ifdef __FILE_OUT__
    fclose (fp);
#endif

/*
----------------------------------------
*/

/**** Statistics (Duration) ****/

    for (I=0; I<3; I++) {
      cpgsci(1);
      cpgsvp(0.10+0.28*(float)I, 0.33+0.28*(float)I, 0.80, 0.95);
      cpgswin(0.0, 13.0, 0.0, 400.0*(float)TIME_STEP/60.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      if (I == 0) {
        cpglab("", "Duration [min]", "R <= 15,000 [km]");
      } else if (I == 1) {
        cpglab("", "", "15,000 < R <= 25,000 [km]");
      } else {
        cpglab("", "", "R > 25,000 [km]");
      }
      for (j=0; j<22; j++) {
        nslr[0][j] = 0.0;
        nslr[1][j] = 0.0;
        nslr[2][j] = 0.0;
        dist[j]    = 0.0;
      }
      for (i=0; i<ncnt[I]; i++) {
        for (j=0; j<22; j++) {
          if (slr_prm[3][I][i] > (float)j && slr_prm[3][I][i] < (float)(j+1)) {
            if (nslr[0][j] == 0.0 && nslr[1][j] == 0) {
              dist[j] = slr_prm[2][I][i];
            }
#ifdef __ELEVATION_DEF__
            if (slr_prm[0][I][i] >= 25.0) {
              nslr[0][j] += 1.0;
            } else if (slr_prm[0][I][i] >= 20.0 && slr_prm[0][I][i] < 25.0) {
              nslr[1][j] += 1.0;
            } else if (slr_prm[0][I][i] < 20.0) {
              nslr[2][j] += 1.0;
            }
#elif defined __DAY_NIGHT_DEF__
            if (slr_prm[0][I][i] >= 20.0) {
              if (slr_prm[5][I][i] >= -5.0) {
                nslr[0][j] += 1.0;
              } else if (slr_prm[5][I][i] < -5.0) {
                nslr[1][j] += 1.0;
              }
            }
#endif
          }
        }
      }
      i = 0;
      for (j=0; j<22; j++) {
        if (nslr[0][j] != 0.0 || nslr[1][j] != 0.0 || nslr[2][j] != 0.0) {
          dist[i] = dist[j];
          nslr[0][i] = nslr[0][j];
          nslr[1][i] = nslr[1][j];
          nslr[2][i] = nslr[2][j];
          i++;
        }
      }

      cpgsfs(2);
      cpgsci(1);
      for (j=0; j<i; j++) {
        cpgrect(dist[j]-0.05, dist[j]+0.05,
                0.0, nslr[0][j]*(float)TIME_STEP/60.0);
      }
      cpgsci(2);
      for (j=0; j<i; j++) {
        cpgrect(dist[j]-0.05, dist[j]+0.05,
                nslr[0][j]*(float)TIME_STEP/60.0,
               (nslr[0][j]+nslr[1][j])*(float)TIME_STEP/60.0);
      }
#ifdef __ELEVATION_DEF__
  #ifdef __LOW_EL_DATA_PLOT__
      cpgsci(7);
      for (j=0; j<i; j++) {
        cpgrect(dist[j]-0.05, dist[j]+0.05,
               (nslr[0][j]+nslr[1][j])*(float)TIME_STEP/60.0,
               (nslr[0][j]+nslr[1][j]+nslr[2][j])*(float)TIME_STEP/60.0);
      }
  #endif
#endif
      cpgsfs(1);
      cpgsci(1);
    }

/*
--------
*/

/**** Range ****/

    for (I=0; I<3; I++) {
      cpgsci(1);
      cpgsvp(0.10+0.28*(float)I, 0.33+0.28*(float)I, 0.62, 0.77);
      cpgswin(0.0, 13.0, 0.0, 35.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      if (I == 0) {
        cpglab("", "Range [x10\\u6\\d m]", "");
      }
      for (i=0; i<ncnt[I]; i++) {
#ifdef __ELEVATION_DEF__
        if (slr_prm[0][I][i] >= 25.0) {
          cpgsci(1);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[1][I]+i, 1);
        } else if (slr_prm[0][I][i] >= 20.0 && slr_prm[0][I][i] < 25.0) {
          cpgsci(2);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[1][I]+i, 1);
  #ifdef __LOW_EL_DATA_PLOT__
        } else if (slr_prm[0][I][i] < 20.0) {
          cpgsci(2);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[1][I]+i, 1);
  #endif
        }
#elif defined __DAY_NIGHT_DEF__
        if (slr_prm[0][I][i] >= 20.0) {
          if (slr_prm[5][I][i] >= -5.0) {
            cpgsci(1);
            cpgpt(1, slr_prm[2][I]+i, slr_prm[1][I]+i, 1);
          } else if (slr_prm[5][I][i] < -5.0) {
            cpgsci(2);
            cpgpt(1, slr_prm[2][I]+i, slr_prm[1][I]+i, 1);
          }
        }
#endif
      }
      cpgsci(1);
    }

/*
--------
*/

/**** Elevation ****/

    for (I=0; I<3; I++) {
      cpgsci(1);
      cpgsvp(0.10+0.28*(float)I, 0.33+0.28*(float)I, 0.44, 0.59);
      cpgswin(0.0, 13.0, 0.0, 90.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      if (I == 0) {
        cpglab("", "Elevation [deg]", "");
      }
      for (i=0; i<ncnt[I]; i++) {
#ifdef __ELEVATION_DEF__
        if (slr_prm[0][I][i] >= 25.0) {
          cpgsci(1);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[0][I]+i, 1);
        } else if (slr_prm[0][I][i] >= 20.0 && slr_prm[0][I][i] < 25.0) {
          cpgsci(2);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[0][I]+i, 1);
  #ifdef __LOW_EL_DATA_PLOT__
        } else if (slr_prm[0][I][i] < 20.0) {
          cpgsci(7);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[0][I]+i, 1);
  #endif
        }
#elif defined __DAY_NIGHT_DEF__
        if (slr_prm[0][I][i] >= 20.0) {
          if (slr_prm[5][I][i] >= -5.0) {
            cpgsci(1);
            cpgpt(1, slr_prm[2][I]+i, slr_prm[0][I]+i, 1);
          } else if (slr_prm[5][I][i] < -5.0) {
            cpgsci(2);
            cpgpt(1, slr_prm[2][I]+i, slr_prm[0][I]+i, 1);
          }
        }
#endif
      }
    }

/*
--------
*/

/**** Incident Angle ****/

    for (I=0; I<3; I++) {
      cpgsci(1);
      cpgsvp(0.10+0.28*(float)I, 0.33+0.28*(float)I, 0.26, 0.41);
      cpgswin(0.0, 13.0, 0.0, 90.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      if (I == 0) {
        cpglab("", "Separation [deg]", "");
      } else if (I == 1) {
        cpglab("Distance from tracking stations [x10\\u6\\d m]", "", "");
      }
      for (i=0; i<ncnt[I]; i++) {
#ifdef __ELEVATION_DEF__
        if (slr_prm[0][I][i] >= 25.0) {
          cpgsci(1);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[4][I]+i, 1);
        } else if (slr_prm[0][I][i] >= 20.0 && slr_prm[0][I][i] < 25.0) {
          cpgsci(2);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[4][I]+i, 1);
  #ifdef __LOW_EL_DATA_PLOT__
        } else if (slr_prm[0][I][i] < 20.0) {
          cpgsci(7);
          cpgpt(1, slr_prm[2][I]+i, slr_prm[4][I]+i, 1);
  #endif
        }
#elif defined __DAY_NIGHT_DEF__
        if (slr_prm[0][I][i] >= 20.0) {
          if (slr_prm[5][I][i] >= -5.0) {
            cpgsci(1);
            cpgpt(1, slr_prm[2][I]+i, slr_prm[4][I]+i, 1);
          } else if (slr_prm[5][I][i] < -5.0) {
            cpgsci(2);
            cpgpt(1, slr_prm[2][I]+i, slr_prm[4][I]+i, 1);
          }
        }
#endif
      }
      cpgsci(1);
    }

/*
--------
*/

/**** Incident Angle (VS. Range) ****/

    for (I=0; I<3; I++) {
      cpgsci(1);
      cpgsvp(0.10+0.28*(float)I, 0.33+0.28*(float)I, 0.04, 0.19);
      cpgswin(0.0, 35.0, 0.0, 90.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      if (I == 0) {
        cpglab("",     "Separation [deg]", "");
      } else if (I == 1) {
        cpglab("Range [x10\\u6\\d m]", "", "");
      } else {
        cpglab("", "",                     "");
      }

      for (i=0; i<ncnt[I]; i++) {
#ifdef __ELEVATION_DEF__
        if (slr_prm[0][I][i] >= 25.0) {
          cpgsci(1);
          cpgpt(1, slr_prm[1][I]+i, slr_prm[4][I]+i, 1);
        } else if (slr_prm[0][I][i] >= 20.0 && slr_prm[0][I][i] < 25.0) {
          cpgsci(2);
          cpgpt(1, slr_prm[1][I]+i, slr_prm[4][I]+i, 1);
  #ifdef __LOW_EL_DATA_PLOT__
        } else if (slr_prm[0][I][i] < 20.0) {
          cpgsci(7);
          cpgpt(1, slr_prm[1][I]+i, slr_prm[4][I]+i, 1);
  #endif
        }
#elif defined __DAY_NIGHT_DEF__
        if (slr_prm[0][I][i] >= 20.0) {
          if (slr_prm[5][I][i] >= -5.0) {
            cpgsci(1);
            cpgpt(1, slr_prm[1][I]+i, slr_prm[4][I]+i, 1);
          } else if (slr_prm[5][I][i] < -5.0) {
            cpgsci(2);
            cpgpt(1, slr_prm[1][I]+i, slr_prm[4][I]+i, 1);
          }
        }
#endif
      }
      cpgsci(1);
    }

/*
----------------------------------------
*/

    for (I=0; I<5; I++) {
      for (J=0; J<3; J++) {
        free (slr_prm[I][J]);
      }
    }
  }

  return;
}
