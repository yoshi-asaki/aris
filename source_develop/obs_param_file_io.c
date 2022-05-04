#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


int  obs_param_file_io(
                       _Bool *ERROR_FLAG,  char *antenna_list_file_name,
                       int  *ARRAY_TYPE,
                       int  *ARRAY_ID,
                       int  *ANT_NUM,  int  *GRT_NUM,  int  *SRT_NUM,
                       struct srt_orbit_parameter *srt,
                       double *grt_elevation_limit,
                       double sep_angle_limit_from_earth_limb,
                       int  *TimUTC, double *UT1_UTC, double *obs_duration,
                       int  *SRCPROC_MODE,
                       struct source_parameter   *src,
                       struct source_parameter   *sun,
                       char    ant_code[][10],
                       char    *ch_grt_el_lim,
                       struct  pair_src_info     *pair_src,
                       struct  char_src_info     *ch_src,
                       struct  char_srt_info     *ch_srt,
                       struct  char_obs_time     *ch_obs_t,
                       int    io_swt)
{
  FILE    *fp;
  int     i, idum, iant;
  char    string[100], cdum[100];

/*
==================================================
*/

  if (io_swt == 0) {
    for (i=0; i<ERROR_NUM; i++) {
      ERROR_FLAG[i] = false;
    }
    sprintf(ch_obs_t->start_t[0], "2015");
    sprintf(ch_obs_t->start_t[1], "01");
    sprintf(ch_obs_t->start_t[2], "01");
    sprintf(ch_obs_t->start_t[3], "00");
    sprintf(ch_obs_t->start_t[4], "00");
    sprintf(ch_obs_t->start_t[5], "00");
    sprintf(ch_obs_t->obsd, "8.0");
    *obs_duration = 8.0;
    *ARRAY_ID = 0;
    *SRT_NUM  = 0;
    *GRT_NUM  = 0;

    ch_src->tgt_ra[0] = '\0';
    ch_src->tgt_dc[0] = '\0';
    ch_src->ref_ra[0] = '\0';
    ch_src->ref_dc[0] = '\0';
    ch_src->mid_ra[0] = '\0';
    ch_src->mid_dc[0] = '\0';
    pair_src->dlt_ra  = (double)0.0;
    pair_src->dlt_dc  = (double)0.0;
    sprintf(ch_src->dlt_ra, "%lf", pair_src->dlt_ra);
    sprintf(ch_src->dlt_dc, "%lf", pair_src->dlt_dc);
    sprintf(ch_grt_el_lim, "20.0");

    if ((fp = fopen("aris_input/input.prm", "r")) != NULL) {
      while (1) {
        if (fgets(string, sizeof(string), fp) == NULL) {
          break;
        } else {
          string[strlen(string)-1] = '\0';
        }
        if (strncmp(string, "ERROR FLAG          ", 20) == 0) {
          for (i=0; i<ERROR_NUM; i++) {
            sscanf(string+20+i, "%1d", &idum);
            if (idum == 1) {
              ERROR_FLAG[i] = true;
            } else {
              ERROR_FLAG[i] = false;
            }
          }
        } else if (strncmp(string, "OBS DATE            ", 20) == 0) {
          char_ncopy(ch_obs_t->start_t[0], string+20, 4);
          char_ncopy(ch_obs_t->start_t[1], string+24, 2);
          char_ncopy(ch_obs_t->start_t[2], string+26, 2);
        } else if (strncmp(string, "START TIME          ", 20) == 0) {
          char_ncopy(ch_obs_t->start_t[3], string+20, 2);
          char_ncopy(ch_obs_t->start_t[4], string+22, 2);
          char_ncopy(ch_obs_t->start_t[5], string+24, 2);
        } else if (strncmp(string, "OBS DURATION        ", 20) == 0) {
          if (strlen(string+20) == 0) {
            sprintf(ch_obs_t->obsd, "0");
            *obs_duration = 0.0;
          } else {
            sscanf(string+20, "%lf", obs_duration);
            char_copy(ch_obs_t->obsd, string+20);
          }
        } else if (strncmp(string, "GRT MINIMUM EL      ", 20) == 0) {
          char_copy(ch_grt_el_lim, string+20);
          sscanf(string+20, "%lf", grt_elevation_limit);
        } else if (strncmp(string, "ANTENNA LIST FILE   ", 20) == 0) {
          char_copy(antenna_list_file_name, string+20);
        }
      }
      fclose (fp);
    }

/******************
    *ANT_NUM = array_config(*ARRAY_ID, wave_id, *SRT_NUM, GRT_NUM,
                        ant_prm,   "",      OFF,
                        antenna_list_file_name,     true);
    for (iant=0; iant<*ANT_NUM; iant++) {
      ant_prm[iant].UFL = OFF;
    }
*******************/


    if ((fp = fopen("aris_input/input.prm", "r")) != NULL) {
      while (1) {
        if (fgets(string, sizeof(string), fp) == NULL) {
          break;
        } else {
          string[strlen(string)-1] = '\0';
        }
        if (strncmp(string, "ARRAY               ", 20) == 0) {
          sscanf(string+20, "%d", ARRAY_ID);
        } else if (strncmp(string, "STATION             ", 20) == 0) {
          sprintf(ant_code[*GRT_NUM], "%s", string+20);
          (*GRT_NUM)++;
        } else if (strncmp(string, "SRT NUMBER          ", 20) == 0) {
          sscanf(string+20, "%d", SRT_NUM);
          if (strlen(string+20) == 0) {
            *SRT_NUM == 0;
          }
          if (*SRT_NUM > SRTMAX) {
            *SRT_NUM = SRTMAX;
          }
        } else if (strncmp(string, "SRT APOGEE          ", 20) == 0) {
          sscanf(string+20, "%d %s", &idum, cdum);
          if (idum >= 1 && idum <= SRTMAX) {
            i = idum - 1;
            if (strlen(string+22) == 0) {
              sprintf(ch_srt[i].apo, "0");
            } else {
              char_copy(ch_srt[i].apo, cdum);
            }
          }
        } else if (strncmp(string, "SRT PERIGEE         ", 20) == 0) {
          sscanf(string+20, "%d %s", &idum, cdum);
          if (idum >= 1 && idum <= SRTMAX) {
            i = idum - 1;
            if (strlen(string+22) == 0) {
              sprintf(ch_srt[i].per, "0");
            } else {
              char_copy(ch_srt[i].per, cdum);
            }
          }
        } else if (strncmp(string, "SRT INCLINATION     ", 20) == 0) {
          sscanf(string+20, "%d %s", &idum, cdum);
          if (idum >= 1 && idum <= SRTMAX) {
            i = idum - 1;
            if (strlen(string+22) == 0) {
              sprintf(ch_srt[i].inc, "0");
            } else {
              char_copy(ch_srt[i].inc, cdum);
            }
          }
        } else if (strncmp(string, "SRT OMEGA           ", 20) == 0) {
          sscanf(string+20, "%d %s", &idum, cdum);
          if (idum >= 1 && idum <= SRTMAX) {
            i = idum - 1;
            if (strlen(string+22) == 0) {
              sprintf(ch_srt[i].OMG, "0");
            } else {
              char_copy(ch_srt[i].OMG, cdum);
            }
          }
        } else if (strncmp(string, "SRT SMALL OMEGA     ", 20) == 0) {
          sscanf(string+20, "%d %s", &idum, cdum);
          if (idum >= 1 && idum <= SRTMAX) {
            i = idum - 1;
            if (strlen(string+22) == 0) {
              sprintf(ch_srt[i].omg, "0");
            } else {
              char_copy(ch_srt[i].omg, cdum);
            }
          }
        } else if (strncmp(string, "SRT T0              ", 20) == 0) {
          sscanf(string+20, "%d %s", &idum, cdum);
          if (idum >= 1 && idum <= SRTMAX) {
            i = idum - 1;
            if (strlen(string+22) == 0) {
              sprintf(ch_srt[i].t_0, "0");
            } else {
              char_copy(ch_srt[i].t_0, cdum);
            }
          }
        } else if (strncmp(string, "SRT d(OMEGA)/dt     ", 20) == 0) {
          sscanf(string+20, "%d %s", &idum, cdum);
          if (idum >= 1 && idum <= SRTMAX) {
            i = idum - 1;
            if (strlen(string+22) == 0) {
              sprintf(ch_srt[i].d_OMG, "0");
            } else {
              char_copy(ch_srt[i].d_OMG, cdum);
            }
          }
        } else if (strncmp(string, "SRT d(omega)/dt     ", 20) == 0) {
          sscanf(string+20, "%d %s", &idum, cdum);
          if (idum >= 1 && idum <= SRTMAX) {
            i = idum - 1;
            if (strlen(string+22) == 0) {
              sprintf(ch_srt[i].d_omg, "0");
            } else {
              char_copy(ch_srt[i].d_omg, cdum);
            }
          }
        } else {
          if (strncmp(string, "SOURCE POS CAL TYPE ", 20) == 0) {
            sscanf(string+20, "%d", SRCPROC_MODE);
          }
          if (       strncmp(string, "SOURCE TGT RA       ", 20) == 0) {
            char_copy(ch_src->tgt_ra, string+20);
          } else if (strncmp(string, "SOURCE TGT DC       ", 20) == 0) {
            char_copy(ch_src->tgt_dc, string+20);
          } else if (strncmp(string, "SOURCE REF RA       ", 20) == 0) {
            char_copy(ch_src->ref_ra, string+20);
          } else if (strncmp(string, "SOURCE REF DC       ", 20) == 0) {
            char_copy(ch_src->ref_dc, string+20);
          }

          if (       strncmp(string, "SOURCE MID RA       ", 20) == 0) {
            char_copy(ch_src->mid_ra, string+20);
          } else if (strncmp(string, "SOURCE MID DC       ", 20) == 0) {
            char_copy(ch_src->mid_dc, string+20);
          } else if (strncmp(string, "SOURCE DELTA RA     ", 20) == 0) {
            if (strlen(string+20) == 0) {
              sprintf(ch_src->dlt_ra, "0");
              pair_src->dlt_ra = (double)0.0;
            } else {
              sscanf(string+20, "%lf", &pair_src->dlt_ra);
              char_copy(ch_src->dlt_ra, string+20);
            }
          } else if (strncmp(string, "SOURCE DELTA DC     ", 20) == 0) {
            if (strlen(string+20) == 0) {
              sprintf(ch_src->dlt_dc, "0");
              pair_src->dlt_dc = (double)0.0;
            } else {
              sscanf(string+20, "%lf", &pair_src->dlt_dc);
              char_copy(ch_src->dlt_dc, string+20);
            }
          }

          if (       strncmp(string, "SOURCE SEPARATION   ", 20) == 0) {
            if (strlen(string+20) == 0) {
              sprintf(ch_src->sepang, "0");
              pair_src->sepang = (double)0.0;
            } else {
              sscanf(string+20, "%lf", &pair_src->sepang);
              char_copy(ch_src->sepang, string+20);
            }
          } else if (strncmp(string, "POSITION ANGLE      ", 20) == 0) {
            if (strlen(string+20) == 0) {
              sprintf(ch_src->posang, "0");
              pair_src->posang = (double)0.0;
            } else {
              sscanf(string+20, "%lf", &pair_src->posang);
              char_copy(ch_src->posang, string+20);
            }
          }
        }
      }
      *ANT_NUM = *GRT_NUM + *SRT_NUM;
      fclose (fp);
    }

    obs_param_set(ERROR_FLAG, *SRT_NUM,
                  ch_srt, srt,
                  ch_grt_el_lim, grt_elevation_limit,
                  sep_angle_limit_from_earth_limb,
                  ch_obs_t, TimUTC, UT1_UTC, obs_duration, *SRCPROC_MODE,
                  ch_src, pair_src, src, sun);

/*
==================================================
*/

  } else if (io_swt == 1) {
    number_char_cut(ch_obs_t->obsd);
    number_char_cut(ch_src->sepang);
    number_char_cut(ch_src->posang);
    number_char_cut(ch_grt_el_lim);
    for (i=0; i<SRTMAX; i++) {
      number_char_cut(ch_srt[i].apo);
      number_char_cut(ch_srt[i].per);
      number_char_cut(ch_srt[i].inc);
      number_char_cut(ch_srt[i].OMG);
      number_char_cut(ch_srt[i].omg);
      number_char_cut(ch_srt[i].d_OMG);
      number_char_cut(ch_srt[i].d_omg);
    }

    if ((fp = fopen("aris_input/input.prm", "w")) != NULL) {
      fprintf(fp, "ERROR FLAG          ");
      for (i=0; i<ERROR_NUM; i++) {
        fprintf(fp, "%1d", ERROR_FLAG[i]);
      }
      fprintf(fp, "\n");
      fprintf(fp, "SOURCE POS CAL TYPE %d\n",  *SRCPROC_MODE);
      fprintf(fp, "SOURCE TGT RA       %s\n",  ch_src->tgt_ra);
      fprintf(fp, "SOURCE TGT DC       %s\n",  ch_src->tgt_dc);
      fprintf(fp, "SOURCE REF RA       %s\n",  ch_src->ref_ra);
      fprintf(fp, "SOURCE REF DC       %s\n",  ch_src->ref_dc);
      fprintf(fp, "SOURCE MID RA       %s\n",  ch_src->mid_ra);
      fprintf(fp, "SOURCE MID DC       %s\n",  ch_src->mid_dc);
      fprintf(fp, "SOURCE DELTA RA     %s\n",  ch_src->dlt_ra);
      fprintf(fp, "SOURCE DELTA DC     %s\n",  ch_src->dlt_dc);
      fprintf(fp, "SOURCE SEPARATION   %s\n",  ch_src->sepang);
      fprintf(fp, "POSITION ANGLE      %s\n",  ch_src->posang);
      fprintf(fp, "OBS DATE            %s%s%s\n",
          ch_obs_t->start_t[0], ch_obs_t->start_t[1], ch_obs_t->start_t[2]);
      fprintf(fp, "START TIME          %s%s%s\n",
          ch_obs_t->start_t[3], ch_obs_t->start_t[4], ch_obs_t->start_t[5]);
      fprintf(fp, "OBS DURATION        %s\n",  ch_obs_t->obsd);
      fprintf(fp, "GRT MINIMUM EL      %s\n",  ch_grt_el_lim);
      fprintf(fp, "ANTENNA LIST FILE   %s\n",  antenna_list_file_name);
      fprintf(fp, "ARRAY               %d\n",  *ARRAY_ID);
      for (i=0; i<*GRT_NUM; i++) {
        fprintf(fp, "STATION             %s\n",  ant_code[i]);
      }
      fprintf(fp, "SRT NUMBER          %d\n",  *SRT_NUM);
      for (i=0; i<SRTMAX; i++) {
        fprintf(fp, "SRT APOGEE          %1d %s\n",   i+1, ch_srt[i].apo);
        fprintf(fp, "SRT PERIGEE         %1d %s\n",   i+1, ch_srt[i].per);
        fprintf(fp, "SRT INCLINATION     %1d %s\n",   i+1, ch_srt[i].inc);
        fprintf(fp, "SRT OMEGA           %1d %s\n",   i+1, ch_srt[i].OMG);
        fprintf(fp, "SRT SMALL OMEGA     %1d %s\n",   i+1, ch_srt[i].omg);
        fprintf(fp, "SRT T0              %1d %s\n",   i+1, ch_srt[i].t_0);
        fprintf(fp, "SRT d(OMEGA)/dt     %1d %s\n",   i+1, ch_srt[i].d_OMG);
        fprintf(fp, "SRT d(omega)/dt     %1d %s\n",   i+1, ch_srt[i].d_omg);
      }
      fclose (fp);
    } else {
      printf("CAUTION: Input parameters cannot be saved. ");
      printf("Make directory \"./aris_input/\".\n");
    }
  }

/*
==================================================
*/

  return (1);
}
