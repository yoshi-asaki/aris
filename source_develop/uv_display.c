#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <aris.h>


int     uv_display(
              int nobs,      double *uv_max,
              int ANT_NUM,
              int BGN_ANT_I, int END_ANT_I,
              int BGN_ANT_J, int END_ANT_J,
              int *TimUTC,   double UT1_UTC,
              struct antenna_parameter *ant_prm,
              struct source_parameter  *src,
              float *src_flag, struct st_observable *int_obs[],
              double *wave_length, int  nswt,
              _Bool  WAVE_NORM,   float pgsch_prm,
              char   *ascii_out,  int    DISP_SWT)
{
  int    iant, jant, iobs;
  int    i, j, I, J, k, l, m, n, ns;
  int    n_blank;
  float  *pgx, *pgy;
  float  scaling_factor;
  double uv_length, wl;
  char   xlab[30], ylab[30];
  char   RA[20], DC[20];
  int    FACTOR;
  FILE   *log_fp;

/*
---------------------------------
*/

  if (*uv_max <= 0.0) {
    for (iant=0; iant<ANT_NUM; iant++) {
      for (jant=iant+1; jant<ANT_NUM; jant++) {
        if (iant >= BGN_ANT_I && iant < END_ANT_I &&
            jant >= BGN_ANT_J && jant < END_ANT_J) {
          for (iobs=0; iobs<nobs; iobs++) {
            I = iant * nobs + iobs;
            J = jant * nobs + iobs;
            for (ns=0; ns<SRC_NUM; ns++) {
              if (int_obs[ns][I].wt >= src_flag[ns] &&
                  int_obs[ns][J].wt >= src_flag[ns]) {
                uv_length =
                     pow(diff(int_obs[ns][I].u, int_obs[ns][J].u), 2.0) +
                     pow(diff(int_obs[ns][I].v, int_obs[ns][J].v), 2.0);
                if (uv_length > *uv_max) {
                  *uv_max = uv_length;
                }
              }
            }
          }
        }
      }
    }
    *uv_max = sqrt(*uv_max);
  }
  if (DISP_SWT == NO_STRUCTURE) {
    return (__GO__);
  }

/*
---------------------------------
*/

  if (ascii_out[0] == '!') {
    if ((log_fp = fopen(ascii_out+1, "w")) == NULL) {
      printf("Warning; ORBIT_DISP: ./aris_log/uvw.log cannot be made.\n");
    } else {
      fprintf(log_fp, "########\n");
      fprintf(log_fp, "# UVW\n");
      fprintf(log_fp, "# FORMAT\n");
      fprintf(log_fp, "# I4,I4,I8,D14.2,D15.3,D15.3,D15.3,D15.3,D15.3,D15.3\n");
      fprintf(log_fp, "#    (ANT_ID1, ANT_ID2, DATA_ID, ");
      fprintf(log_fp, "U(SOURCE11), V(SOURCE1), W(SOURCE1), U(SOURCE2), V(SOURCE2), W(SOURCE2)\n");
      fprintf(log_fp, "#\n");
      fprintf(log_fp, "# Wave length for source 1: %lf [m]\n", wave_length[0]);
      fprintf(log_fp, "# Wave length for source 2: %lf [m]\n", wave_length[1]);
      fprintf(log_fp, "#\n");
      output_star_position(RA, DC, src->s);
      fprintf(log_fp, "# Source position 1       : %s, %s\n", RA, DC);
      output_star_position(RA, DC, (src+1)->s);
      fprintf(log_fp, "# Source position 2       : %s, %s\n", RA, DC);
      fprintf(log_fp, "#\n");
      fprintf(log_fp, "# Observation start time  : %4d/%2d/%2d %2d:%2d:%2d\n",
                         TimUTC[0], TimUTC[1], TimUTC[2],
                         TimUTC[3], TimUTC[4], TimUTC[5]);
      fprintf(log_fp, "#\n");
      for (iant=0; iant<ANT_NUM; iant++) {
        fprintf(log_fp, "# ANT%2d : %s\n", iant+1, ant_prm[iant].IDC);
      }
      fprintf(log_fp, "#\n");
      fprintf(log_fp, "# The number of fringe data per baseline : %d\n", nobs);
      fprintf(log_fp, "# A1,  A2, OBS NUM,             ");
      fprintf(log_fp, "U1,             V1,             W1,             U2,             V2,             W2\n");
      fprintf(log_fp, "#\n");
      for (iant=0; iant<ANT_NUM; iant++) {
        for (jant=iant+1; jant<ANT_NUM; jant++) {
          if (iant >= BGN_ANT_I && iant < END_ANT_I &&
              jant >= BGN_ANT_J && jant < END_ANT_J) {
            for (iobs=0; iobs<nobs; iobs++) {
              I = iant * nobs + iobs;
              J = jant * nobs + iobs;
              fprintf(log_fp, "%4d,%4d,%8d,", iant+1, jant+1, iobs);
              for (ns=0; ns<SRC_NUM; ns++) {
/****
                if (int_obs[ns][I].wt >= src_flag[ns] &&
                    int_obs[ns][J].wt >= src_flag[ns]) {
****/
                  if (ns < SRC_NUM - 1) {
                    fprintf(log_fp, "%15.3lf,%15.3lf,%15.3lf,",
                    diff(int_obs[ns][I].u, int_obs[ns][J].u)/wave_length[ns],
                    diff(int_obs[ns][I].v, int_obs[ns][J].v)/wave_length[ns],
                    diff(int_obs[ns][I].w, int_obs[ns][J].w)/wave_length[ns]);
                  } else {
                    fprintf(log_fp, "%15.3lf,%15.3lf,%15.3lf\n",
                    diff(int_obs[ns][I].u, int_obs[ns][J].u)/wave_length[ns],
                    diff(int_obs[ns][I].v, int_obs[ns][J].v)/wave_length[ns],
                    diff(int_obs[ns][I].w, int_obs[ns][J].w)/wave_length[ns]);
                  }
/****
                }
****/
              }
            }
          }
        }
      }
      fclose (log_fp);
    }
    return (__GO__);
  }

/*
---------------------------------
*/

  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);
  cpgbbuf();

  if ((pgx = calloc(nobs, sizeof(float))) == NULL) {
    printf("ERROR: UV_DISPLAY: calloc fail for pgx.\n");
    return (__NG__);
  }
  if ((pgy = calloc(nobs, sizeof(float))) == NULL) {
    printf("ERROR: UV_DISPLAY: calloc fail for pgy.\n");
    free (pgx);
    return (__NG__);
  }

  cpgsch(pgsch_prm);
  if (WAVE_NORM == true) {
    if (wave_length[0] <= wave_length[1]) {
      wl = wave_length[0];
    } else {
      wl = wave_length[1];
    }

    scaling_factor = 1.0;
    sprintf(xlab, "\\fiu\\fn [\\gl\\fn]");
    sprintf(ylab, "\\fiv\\fn [\\gl\\fn]");

    for (i=1; i<=3; i++) {
      FACTOR = (int)(*uv_max/wl * pow(10.0, (double)(-3*i)));
      if (FACTOR >= 1 && FACTOR < 1000) {
        if (i == 1) {
          scaling_factor = (float)pow(10.0, (double)(-3*i));
          sprintf(xlab, "\\fiu\\fn [k\\gl\\fn]");
          sprintf(ylab, "\\fiv\\fn [k\\gl\\fn]");
          break;
        } else if (i == 2) {
          scaling_factor = (float)pow(10.0, (double)(-3*i));
          sprintf(xlab, "\\fiu\\fn [M\\gl\\fn]");
          sprintf(ylab, "\\fiv\\fn [M\\gl\\fn]");
          break;
        } else if (i == 3) {
          scaling_factor = (float)pow(10.0, (double)(-3*i));
          sprintf(xlab, "\\fiu\\fn [G\\gl\\fn]");
          sprintf(ylab, "\\fiv\\fn [G\\gl\\fn]");
          break;
        }
      }
    }
/**
    cpgsvp(0.15, 0.90, 0.15, 0.90);
**/
    cpgswin(-1.2* *uv_max/wl*scaling_factor, 1.2* *uv_max/wl*scaling_factor,
            -1.2* *uv_max/wl*scaling_factor, 1.2* *uv_max/wl*scaling_factor);
    cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
    cpglab(xlab, ylab, 
           "(\\fiu\\fn, \\fiv\\fn) plot (normalized by wave length)");
  } else {
    scaling_factor = 1.0;
    sprintf(xlab, "\\fiu\\fn [m]");
    sprintf(ylab, "\\fiv\\fn [m]");

    FACTOR = (int)(*uv_max * pow(10.0, -3.0));
    if (FACTOR >= 1 && FACTOR < 10) {
      scaling_factor = (float)pow(10.0, -3);
      sprintf(xlab, "\\fiu\\fn [km]");
      sprintf(ylab, "\\fiv\\fn [km]");
    }
/**
    cpgsvp(0.15, 0.90, 0.15, 0.90);
**/
    cpgswin(-1.2* *uv_max*scaling_factor, 1.2* *uv_max*scaling_factor,
            -1.2* *uv_max*scaling_factor, 1.2* *uv_max*scaling_factor);
    cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
    cpglab(xlab, ylab, "(\\fiu\\fn, \\fiv\\fn) plot\0");
  }

/*
--------
*/

  for (ns=0; ns<SRC_NUM; ns++) {
    cpgsci(ns+1);
    cpgsls(3*ns+1);

    for (iant=0; iant<ANT_NUM; iant++) {
      for (jant=iant+1; jant<ANT_NUM; jant++) {
        if (iant >= BGN_ANT_I && iant < END_ANT_I &&
            jant >= BGN_ANT_J && jant < END_ANT_J) {
          if (DISP_SWT == F__STRUCTURE) {
            iobs = 0;
            while (1) {
              l = 0;
              I = iant * nobs + iobs;
              J = jant * nobs + iobs;
              for (k=iobs; k<nobs; k++) {
                if (int_obs[ns][I].wt >= src_flag[ns] &&
                    int_obs[ns][J].wt >= src_flag[ns]) {
                  pgx[l] = scaling_factor
                               * diff(int_obs[ns][I].u, int_obs[ns][J].u);
                  pgy[l] = scaling_factor
                               * diff(int_obs[ns][I].v, int_obs[ns][J].v);
                  l++;
                } else {
                  iobs = k + 1;
                  break;
                }
                I++;
                J++;
              }
              if (WAVE_NORM == true) {
                for (m=0; m<l; m++) {
                  pgx[m] /= wave_length[ns];
                  pgy[m] /= wave_length[ns];
                }
              }
              cpgline(l, pgx, pgy);
              for (m=0; m<l; m++) {
                pgx[m] *= -1.0;
                pgy[m] *= -1.0;
              }
              cpgline(l, pgx, pgy);

              if (k >= nobs - 1) {
                break;
              }
            }

          } else if (DISP_SWT == C__STRUCTURE) {
            if (nswt != 0) {
              iobs = 0;
              while (1) {
                l = 0;
                m = 0;
                pgx[l] = 0.0;
                pgy[l] = 0.0;
                n_blank = 0;
                I = iant * nobs + iobs;
                J = jant * nobs + iobs;
                for (k=iobs; k<nobs; k++) {
                  if (int_obs[ns][I].wt >= src_flag[ns] &&
                      int_obs[ns][J].wt >= src_flag[ns]) {
                    pgx[l] += scaling_factor
                                 * diff(int_obs[ns][I].u, int_obs[ns][J].u);
                    pgy[l] += scaling_factor
                                 * diff(int_obs[ns][I].v, int_obs[ns][J].v);
                    m++;
                    n_blank = 0;
                  } else {
                    if (n_blank == 0) {
                      if (m != 0) {
                        pgx[l] /= (float)m;
                        pgy[l] /= (float)m;
                        m = 0;
                        l++;
                        pgx[l] = 0.0;
                        pgy[l] = 0.0;
                      }
                      n_blank++;
                    } else {
                      n_blank++;
                    }
                  }

                  if (n_blank > nswt) {
                    iobs = k + 1;
                    break;
                  }
                  I++;
                  J++;
                }

                if (WAVE_NORM == true) {
                  for (m=0; m<l; m++) {
                    pgx[m] /= wave_length[ns];
                    pgy[m] /= wave_length[ns];
                  }
                }
                cpgline(l, pgx, pgy);
                for (m=0; m<l; m++) {
                  pgx[m] *= -1.0;
                  pgy[m] *= -1.0;
                }
                cpgline(l, pgx, pgy);

                if (k >= nobs - 1) {
                  break;
                }
              }
            } else if (nswt == 0) {
              nswt = 60;
              iobs = 0;
              n    = 0;
              while (1) {
                for (k=iobs; k<nobs/nswt; k++) {
                  pgx[n] = 0.0;
                  pgy[n] = 0.0;
                  m = 0;
                  I = iant * nobs + k * nswt;
                  J = jant * nobs + k * nswt;
                  for (i=0; i<nswt; i++) {
                    if (int_obs[ns][I].wt >= src_flag[ns] &&
                        int_obs[ns][J].wt >= src_flag[ns]) {
                      pgx[n] += scaling_factor
                                 * diff(int_obs[ns][I].u, int_obs[ns][J].u);
                      pgy[n] += scaling_factor
                                 * diff(int_obs[ns][I].v, int_obs[ns][J].v);
                      m++;
                    }
                    I++;
                    J++;
                  }
                  if (m != 0) {
                    pgx[n] /= (float)m;
                    pgy[n] /= (float)m;
                    if (WAVE_NORM == true) {
                      pgx[n] /= wave_length[ns];
                      pgy[n] /= wave_length[ns];
                    }
                    n++;
                  } else if (m == 0) {
                    cpgline(n, pgx, pgy);
                    for (l=0; l<n; l++) {
                      pgx[l] *= -1.0;
                      pgy[l] *= -1.0;
                    }
                    cpgline(n, pgx, pgy);
                    n = 0;
                  }
                  iobs += nswt;
                }
                if (iobs >= nobs / nswt) {
                  cpgline(n, pgx, pgy);
                  for (l=0; l<n; l++) {
                    pgx[l] *= -1.0;
                    pgy[l] *= -1.0;
                  }
                  cpgline(n, pgx, pgy);
                  break;
                }
              }
              nswt = 0;
            }
          }
        }
      }
    }
  }
  cpgsch(1.2);
  cpgsls(1);
  cpgebuf();

/*
--------
*/

  free (pgx);
  free (pgy);

  return (__GO__);
}
