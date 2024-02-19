#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

#define NCOM      24
#define MENU_NUM  23

/****
#define __UV_SCALE_FIX__
****/
#ifdef __UV_SCALE_FIX__
  #define UV_MAX 3.0e7
#endif

/****
#define __FITS_SAVE_DEBUG__
****/

int  menu_config(int ANT_NUM, int GRT_NUM,  int SRT_NUM,  int TRK_NUM,
                 int BGN_ANT_I, int END_ANT_I,
                 int BGN_ANT_J, int END_ANT_J,
                 int nobs,
                 int  *TimUTC, double UT1_UTC, _Bool *ERROR_FLAG,
                 int nswt, int nont,
                 struct antenna_parameter *ant_prm,
                 struct antenna_error_parameter *ant_err,
                 struct atmospheric_zenith_error *dz,
                 struct st_observable *int_obs[],
                 double elevation_limit,
                 double *wave_length, double *nu,
                 struct EOP_data EOP,
                 int nbase,
                 float  *src_flag,
                 struct fringe **frng,
                 float  **fringe_weight,
                 double inttim,
                 int    nfrq,
                 double band_width,
                 struct baseline_uvw *bluvw[],
                 struct source_parameter *src,
                 struct source_parameter sun,
                 double *pix_uvl, double *pix_mas,
                 struct srt_orbit_parameter *srt,
                 struct antenna_parameter   *trk_pos,
                 struct srt_data_link       *srt_link,
                 double sep_angle_limit_from_earth_limb,
                 double OBS_T,   double  OBS_p,
                 double **wvc_ds, double **ion_ds,
                 _Bool  TV_SWT, float *cursor_pos,
                 int cpgid1, int cpgid2 )
{
  int    i, j, I, J, ID, iobs, idum;
  int    NX, NY;
  int    proc_mode;
  int    iant, ns, nseq;
  int    timUTC[6];
  double ut1_utc;
  double *tim, obs_start_time_utc;
  double uv_max;
  char   fits_fname[3][40];
  _Bool  fits_save_flag[3];
  char   strtmp[2] = {0, 0};

  int    m_block;
  float  **bttn_box=NULL;
  float  pitch=0.03;

  int    menu_num, MENU_ID, MENU_SECTION;

  int    output_dev_num, PGDEV_ID, PGDEV_SECTION;
  char   pgdev[10][20], ascii_out[100];
  float  sel_pos = 0.220;

  char   string[NCOMLEN];
  char   comstr[NCOM][NCOMLEN];
  struct comment_param cmnt;

  int    FRNG_NUM;
  float  rms_phase[2], coherence[2];

  FILE   *log_fp;

  struct proc_menu {
    char   name[80];
    int    ID;
    _Bool  flag;
  } proc_menu[MENU_NUM];

/*
------------------------------------------
*/

  MENU_ID  = 0;
  FRNG_NUM = 0;

/*
------------------------------------------
*/

  menu_num = MENU_NUM - 1;

  ID =  0;
  sprintf((proc_menu+ID)->name,
                 "ON source duration");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID =  1;
  sprintf((proc_menu+ID)->name,
                 "GRT error quantities");
  proc_menu[ID].ID = ID;
  if (GRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID =  2;
  sprintf((proc_menu+ID)->name,
                 "Az-El plot");
  proc_menu[ID].ID = ID;
  if (GRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID =  3;
  sprintf((proc_menu+ID)->name,
                 "(u, v) Plane");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID =  4;
  sprintf((proc_menu+ID)->name,
                 "Fringe Phase & Amplitude");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID =  5;
  sprintf((proc_menu+ID)->name,
                 "FITS-IDI save");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID =  6;
  sprintf((proc_menu+ID)->name,
                 "Simple FFT Imaging");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID =  7;
  sprintf((proc_menu+ID)->name,
                 "Fringe Fitting");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }
  proc_menu[ID].flag = false;

  ID =  8;
  sprintf((proc_menu+ID)->name,
                 "Orbit of Spacecrafts");
  proc_menu[ID].ID = ID;
/****** xxxxxxxxxxxxxxxx
  if (SRT_NUM >= 1) {
********/
  if (SRT_NUM >= 0) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID =  9;
  sprintf((proc_menu+ID)->name,
                 "Ground Tracking AZ-EL Conditions");
  proc_menu[ID].ID = ID;
  if (SRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 10;
  sprintf((proc_menu+ID)->name,
                 "On-board Link AZ-EL Conditions");
  proc_menu[ID].ID = ID;
  if (SRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 11;
  sprintf((proc_menu+ID)->name,
                 "Ka-link Schedule");
  proc_menu[ID].ID = ID;
  if (SRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 12;
  sprintf((proc_menu+ID)->name,
                 "SLR/Link Conditions (time line)");
  proc_menu[ID].ID = ID;
  if (SRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 13;
  sprintf((proc_menu+ID)->name,
                 "SLR/Link Conditions (statistics)");
  proc_menu[ID].ID = ID;
  if (SRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 14;
  sprintf((proc_menu+ID)->name,
                 "GPS/Link Conditions (time line) (NOT IMPLEMENTED)");
  proc_menu[ID].ID = ID;
  if (SRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }
  proc_menu[ID].flag = false;

  ID = 15;
  sprintf((proc_menu+ID)->name,
                 "GPS/Link Conditions (statistics) (NOT IMPLEMENTED)");
  proc_menu[ID].ID = ID;
  if (SRT_NUM >= 1) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }
  proc_menu[ID].flag = false;

  ID = 16;
  sprintf((proc_menu+ID)->name,
                 "Path Length Error time series");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2 &&
      (ERROR_FLAG[APOSER] == true || ERROR_FLAG[EOPERR] == true ||
       ERROR_FLAG[TDSECZ] == true || ERROR_FLAG[IDSECZ] == true ||
       ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true ||
       ERROR_FLAG[IONTRB] == true ||
       ERROR_FLAG[FQSERR] == true || ERROR_FLAG[LOPOFS] == true)) {
    proc_menu[ID].flag =  true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 17;
  sprintf((proc_menu+ID)->name,
                 "Path Length Error power spectrum (target)");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2 &&
      (ERROR_FLAG[APOSER] == true || ERROR_FLAG[EOPERR] == true ||
       ERROR_FLAG[TDSECZ] == true || ERROR_FLAG[IDSECZ] == true ||
       ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true ||
       ERROR_FLAG[IONTRB] == true ||
       ERROR_FLAG[FQSERR] == true || ERROR_FLAG[LOPOFS] == true)) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 18;
  sprintf((proc_menu+ID)->name,
                 "Path Length Error Allan standard deviation (target)");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2 &&
      (ERROR_FLAG[APOSER] == true || ERROR_FLAG[EOPERR] == true ||
       ERROR_FLAG[TDSECZ] == true || ERROR_FLAG[IDSECZ] == true ||
       ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true ||
       ERROR_FLAG[IONTRB] == true ||
       ERROR_FLAG[FQSERR] == true || ERROR_FLAG[LOPOFS] == true)) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 19;
  sprintf((proc_menu+ID)->name,
                 "Coherence (target)");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2 &&
      (ERROR_FLAG[APOSER] == true || ERROR_FLAG[EOPERR] == true ||
       ERROR_FLAG[TDSECZ] == true || ERROR_FLAG[IDSECZ] == true ||
       ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true ||
       ERROR_FLAG[IONTRB] == true ||
       ERROR_FLAG[FQSERR] == true || ERROR_FLAG[LOPOFS] == true)) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 20;
  sprintf((proc_menu+ID)->name,
                 "Phase Spatial Structure Function (target)");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2 &&
      (ERROR_FLAG[APOSER] == true || ERROR_FLAG[EOPERR] == true ||
       ERROR_FLAG[TDSECZ] == true || ERROR_FLAG[IDSECZ] == true ||
       ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true ||
       ERROR_FLAG[IONTRB] == true ||
       ERROR_FLAG[FQSERR] == true || ERROR_FLAG[LOPOFS] == true)) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

  ID = 21;
  sprintf((proc_menu+ID)->name,
                 "Delay Spatial Structure Function (target)");
  proc_menu[ID].ID = ID;
  if (ANT_NUM >= 2 &&
      (ERROR_FLAG[APOSER] == true || ERROR_FLAG[EOPERR] == true ||
       ERROR_FLAG[TDSECZ] == true || ERROR_FLAG[IDSECZ] == true ||
       ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true ||
       ERROR_FLAG[IONTRB] == true ||
       ERROR_FLAG[FQSERR] == true || ERROR_FLAG[LOPOFS] == true)) {
    proc_menu[ID].flag = true;
  } else {
    proc_menu[ID].flag = false;
  }

/*
------------------------------------------
*/

  output_dev_num = 8;
  sprintf(pgdev[0], "/XS");
  sprintf(pgdev[1], "/PS");
  sprintf(pgdev[2], "/VPS");
  sprintf(pgdev[3], "/CPS");
  sprintf(pgdev[4], "/VCPS");
  sprintf(pgdev[5], "/GIF");
  sprintf(pgdev[6], "ASCII TEXT");
  sprintf(pgdev[7], "/NULL");
  if (TV_SWT == true) {
    PGDEV_ID = 0;
  } else {
    PGDEV_ID = 7;
  }

/*
------------------------------------------
*/

  if (TV_SWT == true) {
    cpgslct(cpgid1);
    cpgpap(1.5*pgpap_prm, 1.0);
    cpgsch(1.5*pgpap_prm/13.0);
    cpgsvp(0.0, 1.0, 0.0, 1.0);
    cpgswin(0.0, 1.0, 0.0, 1.0);

    cpgsci(0);
    cpgrect(0.0, 1.0, 0.0, 1.0);

    cmnt.ncol = 24;
    cmnt.xmin = 0.52;
    cmnt.xmax = 0.98;
    cmnt.ymin = 0.28;
    cmnt.ymax = 0.92;
    cmnt.pitch = 0.03;
    comment_init(&cmnt, comstr, true);
  }

/*
------------------------------------------
*/

  while (1) {
    for (i=0; i<sizeof(ascii_out); i++) {
      ascii_out[i] = 0;
    }

    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;

    if (TV_SWT == false) {

      printf("**** plot : OBSERVING RESULTS ****\n");
      for (i=0; i<menu_num; i++) {
        if (proc_menu[i].flag == true) {
          if (i < 9) {
            strtmp[0] = '1' + i; 
          } else {
            strtmp[0] = 'a' + i - 9; 
          }
          printf("%1s. %s\n", strtmp, (proc_menu + i)->name);
        }
      }
      printf("\n");
      printf("0. RETURN\n");
      printf("-. EXIT\n");
      printf("Select number: ");

      while (1) {
        if (fgets(string, sizeof(string), stdin) == NULL) {
          printf("ERROR: MENU_CONFIG: Invalid input.\n");
          return (__NG__);
        }
        if (string[0] == '-') {
          MENU_ID = -2;
          break;
        } else if (string[0] == '0') {
          MENU_ID = -1;
          break;
        } else if (string[0] - '0' >= 1 && string[0] - '0' <= 9) {
          MENU_ID = string[0] - '1';
          if (proc_menu[MENU_ID].flag == true) {
            break;
          }
        } else if (string[0] - 'a' >= 0 && string[0] - 'a' <= 11) {
          MENU_ID = string[0] - 'a' + 9;
          if (proc_menu[MENU_ID].flag == true) {
            break;
          }
          break;
        } else {
          printf("Wrong number. Select again : ");
        }
      }

      if (MENU_ID != -2 && MENU_ID != -1 && MENU_ID != 1 && MENU_ID != 5) {
        printf("PG Device? ");
        printf("(1./XS  2./PS  3./VPS  4./CPS  ");
        printf( "5./VCPS  /6.GIF  7.ASCII TEXT (CR->/NULL)) : ");
        while (1) {
          if (fgets(string, sizeof(string), stdin) == NULL) {
            printf("ERROR: MENU_CONFIG: Invalid input.\n");
            return (__NG__);
          }
          if (string[0] == '\n') {
            PGDEV_ID = 7;
            break;
          } else if (string[0] >= '1' && string[0] <= '7') {
            PGDEV_ID = string[0] - '1';
            break;
          } else {
            printf("? Input again : ");
          }
        }
      }

/*
------------
*/

    } else if (TV_SWT == true) {

      cpgbbuf();

      TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 0, 1);
      TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 0, 1);
      TV_menu_hatch(0.00, 1.00, 0.00, sel_pos+0.04, 0, 1);

      m_block = 2 + menu_num + output_dev_num;
      if ((bttn_box = (float **)malloc(m_block * sizeof(float *)))
                                                               == NULL) {
        printf("ERROR: menu_config: memory alloc for **bttn_box.\n");
        return (__NG__);
      }
      for (i=0; i<m_block; i++) {
        if ((bttn_box[i] = (float *)calloc(4, sizeof(float))) == NULL) {
          printf("ERROR: menu_config: memory alloc for **bttn_box.\n");
          free_memory_block_float(bttn_box, i);
          return (__NG__);
        }
      }

      I = 0;
      bttn_box[I][0] = 0.12;
      bttn_box[I][1] = 0.42;
      bttn_box[I][2] = 0.97;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      off_button(&idum, "RETURN\0", bttn_box[I]);

      I = 1;
      bttn_box[I][0] = 0.58;
      bttn_box[I][1] = 0.88;
      bttn_box[I][2] = 0.97;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      off_button(&idum, "EXIT\0", bttn_box[I]);

      MENU_SECTION  = 2;
      PGDEV_SECTION = MENU_SECTION + menu_num;

      for (i=0; i<menu_num; i++) {
        I = MENU_SECTION + i;
        if (proc_menu[i].flag == true) {
          if (i <= 7) {
            bttn_box[I][0] = 0.02 + 0.250 * (float)(i%2);
            bttn_box[I][1] = bttn_box[I][0] + 0.238;
            bttn_box[I][2] = 0.89 - 0.037 * (float)(i/2);
            bttn_box[I][3] = bttn_box[I][2] + pitch;
          } else {
            bttn_box[I][0] = 0.02;
            bttn_box[I][1] = bttn_box[I][0] + 0.480;
            bttn_box[I][2] = 0.89 - 0.037 * (float)(i - 4);
            bttn_box[I][3] = bttn_box[I][2] + pitch;
          }
          off_button(&idum, (proc_menu+i)->name, bttn_box[I]);
        } else {
          bttn_box[I][0] = -1.0;
          bttn_box[I][1] = -1.0;
          bttn_box[I][2] = -1.0;
          bttn_box[I][3] = -1.0;
        }
      }

      for (i=0; i<output_dev_num; i++) {
        I = PGDEV_SECTION +i;
        bttn_box[I][0] = 0.10 + 0.11 * (float)i;
        bttn_box[I][1] = 0.20 + 0.11 * (float)i;
        bttn_box[I][2] = 0.93;
        bttn_box[I][3] = bttn_box[I][2] + pitch;

        if (i == 0) {
          on_button(&idum, pgdev[i], bttn_box[I]);
        } else {
          off_button(&idum, pgdev[i], bttn_box[I]);
        }
        PGDEV_ID = 0;
      }

      while (1) {
        MENU_ID = -10;
        cpgcurs(cursor_pos, cursor_pos+1, string);

/*
-------------------
*/

        if (_button_chk(cursor_pos, bttn_box[0]) == true) {
          on_button(&idum, "RETURN\0", bttn_box[0]);
          MENU_ID = -1;
          break;
        }

        if (_button_chk(cursor_pos, bttn_box[1]) == true) {
          on_button(&idum, "EXIT\0", bttn_box[1]);
          MENU_ID = -2;
          break;
        }

/*
-------------------
*/

        for (i=0; i<menu_num; i++) {
          I = MENU_SECTION + i;
          if (_button_chk(cursor_pos, bttn_box[I]) == true) {
            if (i == 0 && ANT_NUM >= 1) {
              on_button(&idum, (proc_menu+i)->name, bttn_box[I]);
              MENU_ID = i;
              break;
            } else if (i >= 1 && i <= 2 && GRT_NUM >= 1) {
              on_button(&idum, (proc_menu+i)->name, bttn_box[I]);
              MENU_ID = i;
              break;
            } else if (i >= 3 && i <= 7 && ANT_NUM >= 2) {
              on_button(&idum, (proc_menu+i)->name, bttn_box[I]);
              MENU_ID = i;
              break;
            } else if (i >= 8 && i <= 15) {
/**** xxxxxxxxxxxxxxxxxxxxxxxxxxx
              if (SRT_NUM >= 1) {
****/
              if (SRT_NUM >= 0) {
                on_button(&idum, (proc_menu+i)->name, bttn_box[I]);
                MENU_ID = i;
                break;
              }
            } else if (i >= 16 && i <= 21 && ANT_NUM >= 2 &&
                 (ERROR_FLAG[APOSER] == true || ERROR_FLAG[EOPERR] == true ||
                  ERROR_FLAG[TDSECZ] == true || ERROR_FLAG[IDSECZ] == true ||
                  ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true ||
                  ERROR_FLAG[IONTRB] == true ||
                  ERROR_FLAG[FQSERR] == true || ERROR_FLAG[LOPOFS] == true)) {
              on_button(&idum, (proc_menu+i)->name, bttn_box[I]);
              MENU_ID = i;
              break;
            }
          }
        }
        if (MENU_ID != -10) {
          break;
        }

/*
-------------------
*/

        for (i=0; i<output_dev_num; i++) {
          I = PGDEV_SECTION + i;
          if (_button_chk(cursor_pos, bttn_box[I]) == true) {
            PGDEV_ID = i;
            for (j=0; j<output_dev_num; j++) {
              J = PGDEV_SECTION + j;
              if (j == i) {
                on_button(&idum, pgdev[j], bttn_box[J]);
              } else {
                off_button(&idum, pgdev[j], bttn_box[J]);
              }
            }
          }
        }
      }

      cpgebuf();
    }

/*
----------------------------------------------------
*/

    if (MENU_ID == -2) {
      sprintf(string, "See you!");
      if (TV_SWT == false) {
        printf("%s\n", string);
      } else if (TV_SWT == true) {
        comment_disp(&cmnt, comstr, string, true);
        free_memory_block_float(bttn_box, m_block);
      }
      return (_EXIT_);
    } else if (MENU_ID == -1) {
      if (TV_SWT == true) {
        free_memory_block_float(bttn_box, m_block);
      }
      return (__GO__);
    }

/*
----------------------------------------------------
*/

    if (PGDEV_ID == 6 &&
       (! (MENU_ID == 2 || MENU_ID == 3 || MENU_ID == 8 || MENU_ID == 20 || MENU_ID == 21))) {
      PGDEV_ID = 7;
      sprintf(string, "Sorry, the text output mode has not been");
      if (TV_SWT == true) {
        comment_disp(&cmnt, comstr, string, true);
      } else {
        printf("%s\n", string);
      }
      sprintf(string, "supported for the selected function.");
      if (TV_SWT == true) {
        comment_disp(&cmnt, comstr, string, true);
      } else {
        printf("%s\n", string);
      }
      sprintf(string, "PGPLOT device was changed to /NULL.");
      if (TV_SWT == true) {
        comment_disp(&cmnt, comstr, string, true);
      } else {
        printf("%s\n", string);
      }
    }

/*
----------------------------------------------------
*/

    if (TV_SWT == true) {
      if (MENU_ID !=  5
        && (! (PGDEV_ID ==  6 &&
               (MENU_ID ==  1 || MENU_ID ==  2 || MENU_ID ==  3 ||
                MENU_ID ==  8 || MENU_ID == 14 || MENU_ID == 15 ||
                MENU_ID == 20 || MENU_ID == 21) /*Not Implemented*/
                                                               ))) {
        cpgid2 = (int)cpgopen(pgdev[PGDEV_ID]);
        if (cpgid2 < 0) {
          cpgask(-1);
        }
      }
    }

/*
----------------------------------------------------
*/

    if (MENU_ID ==  0) {
      if (TV_SWT == true) {
        cpgslct(cpgid2);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
      }
      on_source_disp(ANT_NUM,  GRT_NUM,  nobs, TimUTC,   UT1_UTC,
                     nswt, ant_prm,  int_obs,  elevation_limit);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID ==  1) {

      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgslct(cpgid2);
        } else {
          cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        }
        cpgpap(pgpap_prm, 1.0);

        ant_pos_err_disp(GRT_NUM, ant_prm, ant_err, dz);
        if (TV_SWT == true) {
          cpgclos();
        } else {
          cpgend();
        }
      }

/*
--------
*/

      if (PGDEV_ID != 6) {
        PGDEV_ID = 0;
        if (TV_SWT == true) {
          cpgslct(cpgid1);
          TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 7, 4);
          TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 7, 4);
        }

        while (1) {
          I = END_ANT_I;
          J = END_ANT_J;
          if (END_ANT_I > GRT_NUM) {
            I = GRT_NUM;
          }
          if (END_ANT_J > GRT_NUM) {
            J = GRT_NUM;
          }
          proc_mode = station_select(ANT_NUM,
                         BGN_ANT_I, I, BGN_ANT_J, J, 
                         &NX, &NY, ant_prm, sel_pos, cursor_pos, TV_SWT, 1);
          if (proc_mode == __NG__) {
            if (TV_SWT == true) {
              free_memory_block_float(bttn_box, m_block);
            }
            return (__NG__);
          } else if (proc_mode == _EXIT_) {
            break;
          } else {
            if (ERROR_FLAG[TDSECZ] == true && ERROR_FLAG[TDSECZ] == true) {
              sprintf(string, "%10s : (X,Y,Z)[mm], Zt[ns], Zi[TECU]: ",
                      ant_prm[NX].IDC);
            } else if (ERROR_FLAG[TDSECZ] == true &&
                       ERROR_FLAG[TDSECZ] == false) {
              sprintf(string, "%10s : (X,Y,Z)[mm], Zt[ns]: ",
                      ant_prm[NX].IDC);
            } else if (ERROR_FLAG[TDSECZ] == false &&
                       ERROR_FLAG[TDSECZ] == true) {
              sprintf(string, "%10s : (X,Y,Z)[mm], Zi[TECU]: ",
                      ant_prm[NX].IDC);
            } else if (ERROR_FLAG[TDSECZ] == false &&
                       ERROR_FLAG[TDSECZ] == false) {
              sprintf(string, "%10s : (X,Y,Z)[mm]: ",
                      ant_prm[NX].IDC);
            }
            if (TV_SWT == false) {
              printf("%s", string);
            } else if (TV_SWT == true) {
              comment_disp(&cmnt, comstr, string, true);
            }

            if (ERROR_FLAG[TDSECZ] == true && ERROR_FLAG[TDSECZ] == true) {
              sprintf(string, "  (%6.2f,%6.2f,%6.2f), %6.2f, %6.2f",
                (float)((ant_prm[NX].ERR[0] - ant_prm[NX].XYZ[0])/ 1.0e-3),
                (float)((ant_prm[NX].ERR[1] - ant_prm[NX].XYZ[1])/ 1.0e-3),
                (float)((ant_prm[NX].ERR[2] - ant_prm[NX].XYZ[2])/ 1.0e-3),
                (float)(dz[NX].trp / 1.0e-9),
                (float)(dz[NX].tec / 1.0e16));
            } else if (ERROR_FLAG[TDSECZ] == true &&
                       ERROR_FLAG[TDSECZ] == false) {
              sprintf(string, "  (%6.2f,%6.2f,%6.2f), %6.2f",
                (float)((ant_prm[NX].ERR[0] - ant_prm[NX].XYZ[0])/ 1.0e-3),
                (float)((ant_prm[NX].ERR[1] - ant_prm[NX].XYZ[1])/ 1.0e-3),
                (float)((ant_prm[NX].ERR[2] - ant_prm[NX].XYZ[2])/ 1.0e-3),
                (float)(dz[NX].trp / 1.0e-9));
            } else if (ERROR_FLAG[TDSECZ] == false &&
                       ERROR_FLAG[TDSECZ] == true) {
              sprintf(string, "  (%6.2f,%6.2f,%6.2f), %6.2f",
                (float)((ant_prm[NX].ERR[0] - ant_prm[NX].XYZ[0])/ 1.0e-3),
                (float)((ant_prm[NX].ERR[1] - ant_prm[NX].XYZ[1])/ 1.0e-3),
                (float)((ant_prm[NX].ERR[2] - ant_prm[NX].XYZ[2])/ 1.0e-3),
                (float)(dz[NX].tec / 1.0e16));
            } else if (ERROR_FLAG[TDSECZ] == false &&
                       ERROR_FLAG[TDSECZ] == false) {
              sprintf(string, "  (%6.2f,%6.2f,%6.2f)",
                (float)((ant_prm[NX].ERR[0] - ant_prm[NX].XYZ[0])/ 1.0e-3),
                (float)((ant_prm[NX].ERR[1] - ant_prm[NX].XYZ[1])/ 1.0e-3),
                (float)((ant_prm[NX].ERR[2] - ant_prm[NX].XYZ[2])/ 1.0e-3));
            }
            if (TV_SWT == false) {
              printf("%s\n", string);
            } else if (TV_SWT == true) {
              comment_disp(&cmnt, comstr, string, true);
            }
          }
        }

      } else {
        ascii_out[0] = '!';
        sprintf(ascii_out+1, "aris_log/st.log");
        if ((log_fp = fopen(ascii_out+1, "w")) == NULL) {
          sprintf(string, "Warning: MENU_CONFIG: %s cannot open.\n",
                  ascii_out+1);
          if (TV_SWT == true) {
            cpgslct(cpgid1);
            comment_disp(&cmnt, comstr, string, true);
          } else {
            printf("%s", string);
          }
        } else {
          if (ERROR_FLAG[TDSECZ] == true && ERROR_FLAG[TDSECZ] == true) {
            fprintf(log_fp,
             "  ANTENNA   X[mm]  Y[mm]  Z[mm]    Zt[ns]  Zi[TECU]\n");
          } else if (ERROR_FLAG[TDSECZ] == true && ERROR_FLAG[TDSECZ] == false) {
            fprintf(log_fp,
             "  ANTENNA   X[mm]  Y[mm]  Z[mm]    Zt[ns]\n");
          } else if (ERROR_FLAG[TDSECZ] == false && ERROR_FLAG[TDSECZ] == true) {
            fprintf(log_fp,
             "  ANTENNA   X[mm]  Y[mm]  Z[mm]    Zi[TECU]\n");
          }
          fprintf(log_fp,
           "---------------------------------------------------\n");
          for (i=0; i<GRT_NUM; i++) {
            if (ERROR_FLAG[TDSECZ] == true && ERROR_FLAG[TDSECZ] == true) {
              sprintf(string, "%8s  (%6.2f,%6.2f,%6.2f), %6.2f, %6.2f",
                ant_prm[i].IDC,
                (float)((ant_prm[i].ERR[0] - ant_prm[i].XYZ[0])/ 1.0e-3),
                (float)((ant_prm[i].ERR[1] - ant_prm[i].XYZ[1])/ 1.0e-3),
                (float)((ant_prm[i].ERR[2] - ant_prm[i].XYZ[2])/ 1.0e-3),
                (float)(dz[i].trp / 1.0e-9),
                (float)(dz[i].tec / 1.0e16));
              fprintf(log_fp, "%s\n", string);
            } else if (ERROR_FLAG[TDSECZ] == true &&
                       ERROR_FLAG[TDSECZ] == false) {
              sprintf(string, "%8s  (%6.2f,%6.2f,%6.2f), %6.2f",
                ant_prm[i].IDC,
                (float)((ant_prm[i].ERR[0] - ant_prm[i].XYZ[0])/ 1.0e-3),
                (float)((ant_prm[i].ERR[1] - ant_prm[i].XYZ[1])/ 1.0e-3),
                (float)((ant_prm[i].ERR[2] - ant_prm[i].XYZ[2])/ 1.0e-3),
                (float)(dz[i].trp / 1.0e-9));
            } else if (ERROR_FLAG[TDSECZ] == false &&
                       ERROR_FLAG[TDSECZ] == true) {
              sprintf(string, "%8s  (%6.2f,%6.2f,%6.2f), %6.2f",
                ant_prm[i].IDC,
                (float)((ant_prm[i].ERR[0] - ant_prm[i].XYZ[0])/ 1.0e-3),
                (float)((ant_prm[i].ERR[1] - ant_prm[i].XYZ[1])/ 1.0e-3),
                (float)((ant_prm[i].ERR[2] - ant_prm[i].XYZ[2])/ 1.0e-3),
                (float)(dz[i].tec / 1.0e16));
            } else if (ERROR_FLAG[TDSECZ] == false &&
                       ERROR_FLAG[TDSECZ] == false) {
              sprintf(string, "%8s  (%6.2f,%6.2f,%6.2f)",
                ant_prm[i].IDC,
                (float)((ant_prm[i].ERR[0] - ant_prm[i].XYZ[0])/ 1.0e-3),
                (float)((ant_prm[i].ERR[1] - ant_prm[i].XYZ[1])/ 1.0e-3),
                (float)((ant_prm[i].ERR[2] - ant_prm[i].XYZ[2])/ 1.0e-3));
            }
          }
          fclose (log_fp);
          sprintf(string, "Output file name is %s.\n", ascii_out+1);
          if (TV_SWT == true) {
            cpgslct(cpgid1);
            comment_disp(&cmnt, comstr, string, true);
          } else {
            printf("%s", string);
          }
        }
      }

/*
----
*/

    } else if (MENU_ID ==  2) {
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgslct(cpgid2);
        } else {
          cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        }
        cpgpap(pgpap_prm, 1.0);
      } else {
        ascii_out[0] = '!';
        sprintf(ascii_out, "!aris_log/azel.log");
        sprintf(string, "Output file name is %s.\n", ascii_out+1);
        if (TV_SWT == true) {
          cpgslct(cpgid1);
          comment_disp(&cmnt, comstr, string, true);
        } else {
          printf("%s", string);
        }
      }
      azel_disp(GRT_NUM, nobs, ant_prm, int_obs, timUTC, elevation_limit,
                ascii_out);
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgclos();
        } else {
          cpgend();
        }
      }

/*
----
*/

    } else if (MENU_ID ==  3) {
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgslct(cpgid2);
        } else {
          cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        }
        cpgpap(pgpap_prm, 1.0);
        cpgsvp(0.15, 0.90, 0.15, 0.90);
      } else {
        ascii_out[0] = '!';
        sprintf(ascii_out, "!aris_log/uvw.log");
        sprintf(string, "Output file name is %s.\n", ascii_out+1);
        if (TV_SWT == true) {
          cpgslct(cpgid1);
          comment_disp(&cmnt, comstr, string, true);
        } else {
          printf("%s", string);
        }
      }
      uv_max = 0.0;
#ifdef __UV_SCALE_FIX__
      uv_max = UV_MAX;
#endif /* __UV_SCALE_FIX__ */
      if (uv_display(nobs, &uv_max, ANT_NUM,
                     BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J,
                     timUTC, ut1_utc,
                     ant_prm, src,
                     src_flag, int_obs, wave_length, nswt, true,
                     1.25, ascii_out, C__STRUCTURE) == __NG__) {
        if (TV_SWT == true) {
          free_memory_block_float(bttn_box, m_block);
        }
        return (__NG__);
      }
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgclos();
        } else {
          cpgend();
        }
      }

/********
      int    II;
      float  pga[100000];
      float  pgb[100000];
      float  pgc[100000];
      float  pgd[100000];
      for (i=10; i<nobs-10; i++) {
        pga[II] = (float)II;
        pgb[II] = (bluvw[0][i].w - bluvw[0][i-1].w) / wave_length[0];
        pgc[II] = int_obs[0][i].az / dpi * 180.0 + 90.0;
        pgd[II] = int_obs[0][i].el / dpi * 180.0;
        II++;
      }

      cpgpap(4.0, 2.0);

      cpgsvp(0.1, 0.9, 0.70, 0.85);
      cpgswin(0.0, (float)nobs, 0.0, 90.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      cpgpt(nobs, pga, pgd, 1);
      cpglab("Time", "Elevation", "");

      cpgsvp(0.1, 0.9, 0.50, 0.65);
      cpgswin(0.0, (float)nobs, -1.0, 1.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      cpglab("Time", "Fringe Frequency [Hz]", "");
      cpgpt(nobs, pga, pgb, 1);

      cpgsvp(0.1, 0.9, 0.30, 0.45);
      cpgswin(-180.0, 180.0, -1.0, 1.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      cpglab("Azimuth", "Fringe Frequency [Hz]", "");
      cpgpt(nobs, pgc, pgb, 1);

      cpgsvp(0.1, 0.9, 0.10, 0.25);
      cpgswin(0.0, 90.0, -1.0, 1.0);
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      cpgpt(nobs, pgd, pgb, 1);
      cpglab("Elevation", "Fringe Frequency [Hz]", "");

      cpgclos();
********/


/****
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgslct(cpgid2);
        } else {
          cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        }
        cpgpap(pgpap_prm, 1.0);
        cpgsvp(0.15, 0.90, 0.15, 0.90);
      } else {
        ascii_out[0] = '!';
        sprintf(ascii_out, "!aris_log/******.log");
        sprintf(string, "Output file name is %s.\n", ascii_out+1);
        if (TV_SWT == true) {
          cpgslct(cpgid1);
          comment_disp(&cmnt, comstr, string, true);
        } else {
          printf("%s", string);
        }
      }



      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgclos();
        } else {
          cpgend();
        }
      }
****/

/*
----
*/

    } else if (MENU_ID ==  4) {
/****
      int    ibase, K;
      double blen[3];
      float  pgx[100], pgy[100], pgz[100];
      I = 0;
      for (NX=BGN_ANT_I; NX<END_ANT_I-1; NX++) {
        for (NY=NX+1; NY<END_ANT_J; NY++) {

          ibase = baseline_number(ANT_NUM, NX, NY);
          K = ibase * nobs + 0;
          blen[0] = ant_prm[NX].XYZ[0] - ant_prm[NY].XYZ[0];
          blen[1] = ant_prm[NX].XYZ[1] - ant_prm[NY].XYZ[1];
          blen[2] = ant_prm[NX].XYZ[2] - ant_prm[NY].XYZ[2];

          proc_mode = 2;
          FRNG_NUM  = 2;
          if (fringe_disp(FRNG_NUM, ANT_NUM, NX, NY, ant_prm, nobs, nfrq,
                          timUTC, ut1_utc, frng, fringe_weight,
                          rms_phase, coherence) == -1) {
            printf("MENU_CONFIG: FRINGE_DISP: ERROR\n");
            exit (-1);
          }

          proc_mode = 3;
          FRNG_NUM = 1;
          if (fringe_disp(FRNG_NUM, ANT_NUM, NX, NY, ant_prm, nobs, nfrq,
                          timUTC, ut1_utc, &frng[2], &fringe_weight[2],
                          rms_phase, coherence) == -1) {
            printf("MENU_CONFIG: FRINGE_DISP: ERROR\n");
            exit (-1);
          }

          printf("%6.2f   %6.2f   %lf\n", rms_phase[0], rms_phase[1], vlen3(blen));
          pgx[I] = log10(vlen3(blen));
          pgy[I] = log10(rms_phase[0]);
          pgz[I] = log10(rms_phase[1]);
          I++;
        }
      }
      cpgslct(cpgid1);
      cpgpap(6.0, 1.0);
      cpgsvp(0.15, 0.9, 0.15, 0.9);
      cpgscf(2);
      cpgsch(1.5);
      cpgswin(log10(10.0), log10(10.0e3), log10(1.0), log10(200.0));
      cpgbox("BCNLTS", 0, 0, "BCNLTS", 0, 0);
      cpglab("Baseline length [m]", "RMS phase@BAND-3", "ARIS Simulation");
      cpgpt(I, pgx, pgy, 17);
      cpgsci(2);
      cpgpt(I, pgx, pgz, 23);
      cpgend();
      getchar();
****/




      if (TV_SWT == true) {
        cpgslct(cpgid1);
        TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 7, 4);
        TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 7, 4);
      }
      while (1) {
        if (TV_SWT == true) {
          cpgslct(cpgid1);
        }
        proc_mode = station_select(ANT_NUM,
                           BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J, 
                           &NX, &NY, ant_prm, sel_pos, cursor_pos, TV_SWT, 3);
        if (proc_mode == __NG__) {
          if (TV_SWT == true) {
            free_memory_block_float(bttn_box, m_block);
          }
          return (__NG__);
        } else if (proc_mode == _EXIT_) {
          break;
        } else {
          if (TV_SWT == true) {
            cpgslct(cpgid2);
          } else {
            cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
            cpgpap(pgpap_prm, 1.0);
          }
          if (proc_mode == 2) {
            FRNG_NUM = 2;
            if (fringe_disp(FRNG_NUM, ANT_NUM, NX, NY, ant_prm, nobs, nfrq,
                            timUTC, ut1_utc, frng, fringe_weight,
                            rms_phase, coherence) == -1) {
              printf("MENU_CONFIG: FRINGE_DISP: ERROR\n");
              exit (-1);
            }
          } else if (proc_mode == 3) {
            FRNG_NUM = 1;
            if (fringe_disp(FRNG_NUM, ANT_NUM, NX, NY, ant_prm, nobs, nfrq,
                            timUTC, ut1_utc, &frng[2], &fringe_weight[2],
                            rms_phase, coherence) == -1) {
              printf("MENU_CONFIG: FRINGE_DISP: ERROR\n");
              exit (-1);
            }
          }

          for (i=0; i<FRNG_NUM; i++) {
            if (i == 0) {
              sprintf(string, "RMS PHASE (TGT) : %6.2f [deg]", rms_phase[i]);
/****
              printf("COHERENCE : %5.3f\n", coherence[i]);
****/
            } else if (i == 1) {
              sprintf(string, "RMS PHASE (REF) : %6.2f [deg]", rms_phase[i]);
            }
            if (TV_SWT == true) {
              cpgslct(cpgid1);
              comment_disp(&cmnt, comstr, string, true);
            } else {
              printf("%s\n", string);
            }
          }
        }
      }
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID ==  5) {
#ifdef __FITS_SAVE_DEBUG__
      printf("#### __DEBUG__ : FITS_DATA_SAVE\n"); fflush(stdout);
#endif
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 7, 4);
        TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 7, 4);
      }
      while (1) {
        if (TV_SWT == true) {
          cpgslct(cpgid1);
        }
#ifdef __FITS_SAVE_DEBUG__
        printf("#### __DEBUG__ : fits_data_select in\n"); fflush(stdout);
#endif
        if (fits_data_select(fits_fname, fits_save_flag,
                             sel_pos, cursor_pos, TV_SWT) == _EXIT_) {
          break;
        } else {
#ifdef __FITS_SAVE_DEBUG__
          printf("#### __DEBUG__ : fits_data_select out\n"); fflush(stdout);
#endif
          for (i=0; i<6; i++) {
            timUTC[i] = TimUTC[i];
          }
          ut1_utc = UT1_UTC;

          obs_start_time_utc = (double)
                     (3600 * timUTC[0] + 60 * timUTC[1] + timUTC[2]) + ut1_utc;
          tim = (double *)calloc(nobs, sizeof(double));
          for (iobs=0; iobs<nobs; iobs++) {
            tim[iobs] = obs_start_time_utc + (double)iobs;
          }
          for (ns=0; ns<3; ns++) {
            if (fits_save_flag[ns] == true) {
#ifdef __FITS_SAVE_DEBUG__
              printf("#### __DEBUG__ : fitsidi_save in\n"); fflush(stdout);
#endif
              nseq = ns % 2;
              if (fitsidi_save(fits_fname[ns], nobs,
                    ANT_NUM, GRT_NUM, SRT_NUM,
                    BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J, ant_prm,
                    timUTC, ut1_utc, inttim, nfrq, band_width, nu[ns], EOP,
                    tim, bluvw[nseq], frng[ns], fringe_weight[ns],
                    ANT_NUM,
                    BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J, ant_prm,
                    src[nseq], srt) == -1) {
                printf("MENU_CONFIG: FITSIDI_SAVE: ERROR: STOP.\n");
                if (TV_SWT == true) {
                  free_memory_block_float(bttn_box, m_block);
                }
#ifdef __FITS_SAVE__
              printf("#### __DEBUG__ : fitsidi_save out\n");
#endif
                return (__NG__);
              }
            }
          }
          free (tim);
        }
      }

/*
----
*/

    } else if (MENU_ID ==  6) {
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgpap(pgpap_prm, 1.0);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        cpgpap(pgpap_prm, 1.0);
      }
      qlook_imager(ANT_NUM, BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J,
                   nswt, nobs, nfrq, bluvw, frng, fringe_weight, wave_length);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID ==  7) {

      sprintf(string, "Sorry, this function has not been implemented.");
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        comment_disp(&cmnt, comstr, string, true);
      } else {
        printf("%s\n", string);
      }

/********
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 7, 4);
        TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 7, 4);
      }
      while (1) {
        if (TV_SWT == true) {
          cpgslct(cpgid1);
        }
        proc_mode = station_select(ANT_NUM,
                           BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J, 
                           &NX, &NY, ant_prm, sel_pos, cursor_pos, TV_SWT, 3);
        if (proc_mode == __NG__) {
          if (TV_SWT == true) {
            free_memory_block_float(bttn_box, m_block);
          }
          return (__NG__);
        } else if (proc_mode == _EXIT_) {
          break;
        } else {
          if (TV_SWT == true) {
            cpgslct(cpgid2);
          } else {
            cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
            cpgpap(pgpap_prm, 1.0);
          }
          if (proc_mode == 2) {
            FRNG_NUM = 2;
            if (fringe_fitting(FRNG_NUM, ANT_NUM, NX, NY, ant_prm, nobs, nfrq,
                            timUTC, ut1_utc, frng, fringe_weight,
                            rms_phase) == -1) {
              printf("MENU_CONFIG: FRINGE_DISP: ERROR\n");
              exit (-1);
            }
          } else if (proc_mode == 3) {
            FRNG_NUM = 1;
            if (fringe_fitting(FRNG_NUM, ANT_NUM, NX, NY, ant_prm, nobs, nfrq,
                            timUTC, ut1_utc, &frng[2], &fringe_weight[2],
                            rms_phase) == -1) {
              printf("MENU_CONFIG: FRINGE_DISP: ERROR\n");
              exit (-1);
            }
          }

          for (i=0; i<FRNG_NUM; i++) {
            if (i == 0) {
              sprintf(string, "RMS PHASE (TGT) : %6.2f [deg]", rms_phase[i]);
              printf("COHERENCE : %5.3f\n", coherence[i]);
            } else if (i == 1) {
              sprintf(string, "RMS PHASE (REF) : %6.2f [deg]", rms_phase[i]);
            }
            if (TV_SWT == true) {
              cpgslct(cpgid1);
              comment_disp(&cmnt, comstr, string, true);
            } else {
              printf("%s\n", string);
            }
          }
        }
      }
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgclos();
      } else {
        cpgend();
      }
********/

/*
----
*/

    } else if (MENU_ID ==  8) {
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgslct(cpgid2);
        } else {
          cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        }
        cpgpap(pgpap_prm, 1.0);
      } else {
        ascii_out[0] = '!';
        sprintf(ascii_out, "!aris_log/orbit.log");
        sprintf(string, "Output file name is %s.\n", ascii_out+1);
        if (TV_SWT == true) {
          cpgslct(cpgid1);
          comment_disp(&cmnt, comstr, string, true);
        } else {
          printf("%s", string);
        }
      }
      orbit_disp(SRT_NUM, nobs, srt, timUTC, ut1_utc,
                 sep_angle_limit_from_earth_limb,
                 srt_link,
                 OBS_T, OBS_p,
                 TRK_NUM, GRT_NUM, int_obs,
                 trk_pos, src[0], sun, ascii_out);
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgclos();
        } else {
          cpgend();
        }
      }

/*
----
*/

    } else if (MENU_ID ==  9) {
      if (TV_SWT == true) {
        cpgslct(cpgid2);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        cpgpap(pgpap_prm, 1.0);
      }
      link_TRK_sim(SRT_NUM, nobs, srt, trk_pos,
                src[0], sun, srt_link,
                timUTC, ut1_utc, sep_angle_limit_from_earth_limb,
                OBS_T, OBS_p, TRK_NUM, TRACKING);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID == 10) {
      if (TV_SWT == true) {
        cpgslct(cpgid2);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        cpgpap(pgpap_prm, 1.0);
      }
      link_TRK_sim(SRT_NUM, nobs, srt, trk_pos,
                src[0], sun, srt_link,
                timUTC, ut1_utc, sep_angle_limit_from_earth_limb,
                OBS_T, OBS_p, TRK_NUM, ONBOARD);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID == 11) {
      if (TV_SWT == true) {
        cpgslct(cpgid2);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        cpgpap(pgpap_prm, 1.0);
      }
      data_link_schedule(SRT_NUM, nobs, srt, trk_pos,
                src[0], sun, srt_link,
                timUTC, ut1_utc, sep_angle_limit_from_earth_limb,
                OBS_T, OBS_p, TRK_NUM, ONBOARD);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID == 12) {
      if (TV_SWT == true) {
        cpgslct(cpgid2);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        cpgpap(pgpap_prm, 1.0);
      }
      link_SLR_sim(SRT_NUM, nobs, srt, trk_pos,
                src[0], sun, srt_link,
                timUTC, ut1_utc, sep_angle_limit_from_earth_limb,
                OBS_T, OBS_p, TRK_NUM, 0);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID == 13) {
      if (TV_SWT == true) {
        cpgslct(cpgid2);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        cpgpap(pgpap_prm, 1.0);
      }
      link_SLR_sim(SRT_NUM, nobs, srt, trk_pos,
                src[0], sun, srt_link,
                timUTC, ut1_utc, sep_angle_limit_from_earth_limb,
                OBS_T, OBS_p, TRK_NUM, 1);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID == 14) {
      sprintf(string, "Sorry, this function has not been implemented.");
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        comment_disp(&cmnt, comstr, string, true);
      } else {
        printf("%s\n", string);
      }

/****AAAA
      if (TV_SWT == true) {
        cpgslct(cpgid2);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        cpgpap(pgpap_prm, 1.0);
      }
      link_SLR_sim(SRT_NUM, nobs, srt, trk_pos,
                src[0], sun, srt_link,
                timUTC, ut1_utc, sep_angle_limit_from_earth_limb,
                OBS_T, OBS_p, TRK_NUM, 1);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }
AAAA****/

/*
----
*/

    } else if (MENU_ID == 15) {
      sprintf(string, "Sorry, this function has not been implemented.");
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        comment_disp(&cmnt, comstr, string, true);
      } else {
        printf("%s\n", string);
      }

/****AAAA
      if (TV_SWT == true) {
        cpgslct(cpgid2);
      } else {
        cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        cpgpap(pgpap_prm, 1.0);
      }
      link_SLR_sim(SRT_NUM, nobs, srt, trk_pos,
                src[0], sun, srt_link,
                timUTC, ut1_utc, sep_angle_limit_from_earth_limb,
                OBS_T, OBS_p, TRK_NUM, 1);
      if (TV_SWT == true) {
        cpgclos();
      } else {
        cpgend();
      }
AAAA****/

/*
----
*/

    } else if (MENU_ID == 16) {
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 7, 4);
        TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 7, 4);
      }
      while (1) {
        if (TV_SWT == true) {
          cpgslct(cpgid1);
        }
        proc_mode = station_select(ANT_NUM,
                           BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J, 
                           &NX, &NY, ant_prm, sel_pos, cursor_pos, TV_SWT, 2);
        if (proc_mode == __NG__) {
          if (TV_SWT == true) {
            free_memory_block_float(bttn_box, m_block);
          }
          return (__NG__);
        } else if (proc_mode == _EXIT_) {
          break;
        } else {
          if (TV_SWT == true) {
            cpgslct(cpgid2);
          } else {
            cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
            cpgpap(pgpap_prm, 1.0);
          }
          EPL_disp(NX, NY, nobs, ant_prm, nu, timUTC, int_obs);
          if (TV_SWT == false) {
            cpgend();
          }
        }
      }
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgclos();
      }

/*
----
*/

    } else if (MENU_ID == 17) {
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 7, 4);
        TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 7, 4);
      }
      while (1) {
        if (TV_SWT == true) {
          cpgslct(cpgid1);
        }
        proc_mode = station_select(ANT_NUM,
                           BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J, 
                           &NX, &NY, ant_prm, sel_pos, cursor_pos, TV_SWT, 2);
        if (proc_mode == __NG__) {
          if (TV_SWT == true) {
            free_memory_block_float(bttn_box, m_block);
          }
          return (__NG__);
        } else if (proc_mode == _EXIT_) {
          break;
        } else {
          if (TV_SWT == true) {
            cpgslct(cpgid2);
          } else {
            cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
            cpgpap(pgpap_prm, 1.0);
          }
          spectl_disp(NX, NY, nobs, ant_prm, nu, timUTC, int_obs);
          if (TV_SWT == false) {
            cpgend();
          }
        }
      }
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgclos();
      }

/*
----
*/

    } else if (MENU_ID == 18) {
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 7, 4);
        TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 7, 4);
      }
      while (1) {
        if (TV_SWT == true) {
          cpgslct(cpgid1);
        }
        proc_mode = station_select(ANT_NUM,
                           BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J, 
                           &NX, &NY, ant_prm, sel_pos, cursor_pos, TV_SWT, 2);
        if (proc_mode == __NG__) {
          if (TV_SWT == true) {
            free_memory_block_float(bttn_box, m_block);
          }
          return (__NG__);
        } else if (proc_mode == _EXIT_) {
          break;
        } else {
          if (TV_SWT == true) {
            cpgslct(cpgid2);
          } else {
            cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
            cpgpap(pgpap_prm, 1.0);
          }
          allanv_disp(NX, NY, nobs, ant_prm, nu, timUTC, int_obs);
          if (TV_SWT == false) {
            cpgend();
          }
        }
      }
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgclos();
      }

/*
----
*/

    } else if (MENU_ID == 19) {
      if (TV_SWT == true) {
        cpgslct(cpgid1);
        TV_menu_hatch(0.00, 1.00, 0.92, 1.00, 7, 4);
        TV_menu_hatch(0.00, 0.50, sel_pos+0.04, 0.92, 7, 4);
      }
      while (1) {
        if (TV_SWT == true) {
          cpgslct(cpgid1);
        }
        proc_mode = station_select(ANT_NUM,
                           BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J, 
                           &NX, &NY, ant_prm, sel_pos, cursor_pos, TV_SWT, 3);
        if (proc_mode == __NG__) {
          if (TV_SWT == true) {
            free_memory_block_float(bttn_box, m_block);
          }
          return (__NG__);
        } else if (proc_mode == _EXIT_) {
          break;
        } else {
          if (TV_SWT == true) {
            cpgslct(cpgid2);
          } else {
            cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
            cpgpap(pgpap_prm, 1.0);
          }
          if (proc_mode == 2) {
            FRNG_NUM = 2;
            if (coherence_disp(FRNG_NUM, ANT_NUM, NX, NY, ant_prm, nobs, nfrq,
                            timUTC, ut1_utc, frng, fringe_weight,
                            rms_phase) == -1) {
              printf("MENU_CONFIG: FRINGE_DISP: ERROR\n");
              exit (-1);
            }
          } else if (proc_mode == 3) {
            FRNG_NUM = 1;
            if (coherence_disp(FRNG_NUM, ANT_NUM, NX, NY, ant_prm, nobs, nfrq,
                            timUTC, ut1_utc, &frng[2], &fringe_weight[2],
                            rms_phase) == -1) {
              printf("MENU_CONFIG: FRINGE_DISP: ERROR\n");
              exit (-1);
            }
          }
        }
      }
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID == 20) {
      printf("AAAAAAAAAAAAAAAAA   %d   %s\n", PGDEV_ID, pgdev[PGDEV_ID]);
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgslct(cpgid2);
        } else {
          cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        }
        cpgpap(pgpap_prm, 1.0);
      } else {
        ascii_out[0] = '!';
        sprintf(ascii_out, "!aris_log/ssf.log");
        sprintf(string, "Output file name is %s.\n", ascii_out+1);
        if (TV_SWT == true) {
          cpgslct(cpgid1);
          comment_disp(&cmnt, comstr, string, true);
        } else {
          printf("%s", string);
        }
      }
      SSF_disp(ANT_NUM,   GRT_NUM,
               BGN_ANT_I, END_ANT_I,
               BGN_ANT_J, END_ANT_J,
               nobs,      nfrq, timUTC,  ut1_utc, nu,
               ant_prm,   frng, fringe_weight, int_obs, ascii_out, SSF_PHASE);
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgclos();
      } else {
        cpgend();
      }

/*
----
*/

    } else if (MENU_ID == 21) {
      if (PGDEV_ID != 6) {
        if (TV_SWT == true) {
          cpgslct(cpgid2);
        } else {
          cpgbeg(1, pgdev[PGDEV_ID], 1, 1);
        }
        cpgpap(pgpap_prm, 1.0);
      } else {
        ascii_out[0] = '!';
        sprintf(ascii_out, "!aris_log/ssf.log");
        sprintf(string, "Output file name is %s.\n", ascii_out+1);
        if (TV_SWT == true) {
          cpgslct(cpgid1);
          comment_disp(&cmnt, comstr, string, true);
        } else {
          printf("%s", string);
        }
      }
      SSF_disp(ANT_NUM,   GRT_NUM,
               BGN_ANT_I, END_ANT_I,
               BGN_ANT_J, END_ANT_J,
               nobs,      nfrq, timUTC,  ut1_utc, nu,
               ant_prm,   frng, fringe_weight, int_obs, ascii_out, SSF_DELAY);
      if (TV_SWT == true) {
        cpgslct(cpgid2);
        cpgclos();
      } else {
        cpgend();
      }
    }

/*
===========================================================
*/

    if (TV_SWT == true) {
      cpgslct(cpgid1);
      I = MENU_SECTION + MENU_ID;
      off_button(&idum, (proc_menu+MENU_ID)->name, bttn_box[I]);
    }
  }

  if (TV_SWT == true) {
    free_memory_block_float(bttn_box, m_block);
  }
  return (__GO__);
}
