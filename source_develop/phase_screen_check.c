#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <cpgplot.h>
#include <aris.h>

#define ATM_SECTION   10
#define CNTR_SECTION  30
#define NAME_SECTION  40

/****
#define __DEBUG__
****/

#define SINGLE_SCRN 0
#define MULTI__SCRN 1


int   ps_model_out(char *, char [][40], int , int , _Bool , int , char *);
int   phase_screen_cull(
             int       , int       , double   *,
             int       , int       , int       ,
             float    *, float    *, float    *);


int  phase_screen_check(int NMTRX, struct phase_screen_parameter wvc,
                        float  *cursor_pos, _Bool TV_SWT, int *pgid)
{
  int    i, j, I, NOD, idum;
  int    NX, NY;
  int    nx=0, ny=0;
  int    *IX=NULL, *IY=NULL;
  int    corner_position[2][2];
  int    COUNT_NOD;
  int    SCREEN_MODE=0;
  int    itmp;
  _Bool  CONT_SWT_S=false;
  _Bool  Bdum;
  double *DS;
  double *seed_dist;
  float  dmin=0.0, dmax=0.0;
  float  m_min=0.0, m_max=0.0;
  float  pmin, pmax, noise, err_x, err_y, delta_x=0.0, delta_y=0.0;
  float  *s_x_cnt=NULL, *s_y_cnt=NULL, *s_x_w=NULL, *s_y_w=NULL;
  float  pitch = 0.03;
  float  plot_size_y=0.0;
  double lftmp;

  FILE   *ifp, *ofp;
  float  bttn_box[60][4];
  char   string[100];
  char   param_char[20][40];
  char   fname[40];
  char   ps_model_prm[100];
  char   data_file[100];

  _Bool  PROC_SWT;
  _Bool  DISP_SWT;
  double ps_param[20];
  double diff_expon;

  float  tr[6];
  float  pgsvpxmin=0.0, pgsvpxmax=0.0, pgsvpymin=0.0, pgsvpymax=0.0;
  float  *m_dist;
  int    n_bandle=0;
  int    N128=128;
  int    NSCREEN=1, nscrn, NPLOT=0;
  float  box1=0.01, box2=0.99, box3=0.09, box4=0.62;
  float  plot_space=0.037;

/*
----------------------
*/

  sprintf(ps_model_prm, "aris_input/ps_model.prm");

/*
---------------------- Atmospheric Turbulence 2D
*/

  ps_param[0]  = wvc.H_d;
  ps_param[1]  = vlen2(wvc.v);
  ps_param[2]  = 180.0 / dpi * atan2(wvc.v[1], wvc.v[0]);
  ps_param[3]  = wvc.i_scale[0];
  ps_param[4]  = wvc.o_scale[0];
  ps_param[5]  = wvc.i_expon;
  ps_param[6]  = wvc.o_expon;
  ps_param[7]  = wvc.pixel;
  ps_param[8]  = 1024.0;
  ps_param[9]  = 1024.0;
  ps_param[11] = 100.0;
  ps_param[12] = 1.5;
  ps_param[10] = ps_param[12] * 0.02 / pow(ps_param[11], 0.5*ps_param[5]);

  for (i=0; i<13; i++) {
    sprintf(param_char[i], "%lf", ps_param[i]);
  }
  DISP_SWT = false;
  sprintf(fname, "screen.dat");

  if ((ifp = fopen(ps_model_prm, "r")) != NULL) {
    while (1) {
      if (fgets(string, sizeof(string), ifp) == NULL) {
        break;
      } else {
        string[strlen(string)-1] = '\0';
      }
      if (strncmp(string, "SCREEN HEIGHT       ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[0]);
        sscanf(string+20, "%lf", &ps_param[0]);
      } else if (strncmp(string, "WIND VELOCITY       ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[1]);
        sscanf(string+20, "%lf", &ps_param[1]);
      } else if (strncmp(string, "POSITION ANGLE      ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[2]);
        sscanf(string+20, "%lf", &ps_param[2]);
      } else if (strncmp(string, "INNER SCALE LENGTH  ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[3]);
        sscanf(string+20, "%lf", &ps_param[3]);
      } else if (strncmp(string, "OUTER SCALE LENGTH  ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[4]);
        sscanf(string+20, "%lf", &ps_param[4]);
      } else if (strncmp(string, "INNER EXPONENT      ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[5]);
        sscanf(string+20, "%lf", &ps_param[5]);
      } else if (strncmp(string, "OUTER EXPONENT      ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[6]);
        sscanf(string+20, "%lf", &ps_param[6]);
      } else if (strncmp(string, "PIXEL SIZE          ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[7]);
        sscanf(string+20, "%lf", &ps_param[7]);
      } else if (strncmp(string, "SCREEN WIDTH (X)    ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[8]);
        sscanf(string+20, "%lf", &ps_param[8]);
      } else if (strncmp(string, "SCREEN WIDTH (Y)    ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[9]);
        sscanf(string+20, "%lf", &ps_param[9]);
      } else if (strncmp(string, "RMS VALUE           ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[10]);
        sscanf(string+20, "%lf", &ps_param[10]);
      } else if (strncmp(string, "BASELINE OF THE RMS ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[11]);
        sscanf(string+20, "%lf", &ps_param[11]);
      } else if (strncmp(string, "BIAS                ", 20) == 0) {
        sscanf(string+20, "%s",  param_char[12]);
        sscanf(string+20, "%lf", &ps_param[12]);
      } else if (strncmp(string, "DISPLAY SWITCH      ", 20) == 0) {
        sscanf(string+20, "%d",  &idum);
        if (idum == 1) {
          DISP_SWT = false;
        } else if (idum == 1) {
          DISP_SWT = true;
        }
      } else if (strncmp(string, "SAVE FILE NAME      ", 20) == 0) {
        sscanf(string+20, "%s", fname);
      } else if (strncmp(string, "SCREEN GENERATION   ", 20) == 0) {
        sscanf(string+20, "%d", &SCREEN_MODE);
      }
    }
    fclose (ifp);
  }

/*
----------------------
*/

  if (TV_SWT == true) {
    cpgslct(pgid[0]);
    cpgpap(1.5*pgpap_prm, 1.0);
    cpgsvp(0.0, 1.0, 0.0, 1.0);
    cpgswin(0.0, 1.0, 0.0, 1.0);
    cpgsch(0.65);

    cpgsci(0);
    cpgrect(0.0, 1.0, 0.0, 1.0);
    cpgsci(1);

    for (i=0; i<13; i++) {
      I = ATM_SECTION + i;
      bttn_box[I][0] = 0.21 + 0.50 * (float)(i/7);
      bttn_box[I][1] = 0.44 + 0.50 * (float)(i/7);
      bttn_box[I][2] = 0.94 - 0.04*(float)(i%7);
      bttn_box[I][3] = bttn_box[I][2] + pitch;

      cpgsci(1);
      if (i == 0) {
        sprintf(string, "Screen Height [m]");
      } else if (i == 1) {
        sprintf(string, "Wind Velocity [m/s]");
      } else if (i == 2) {
        sprintf(string, "Position Angle [deg]");
      } else if (i == 3) {
        sprintf(string, "Inner Scale [m]");
      } else if (i == 4) {
        sprintf(string, "Outer Scale [m]");
      } else if (i == 5) {
        sprintf(string, "Inner Exponent");
      } else if (i == 6) {
        sprintf(string, "Outer Exponent");
      } else if (i == 7) {
        sprintf(string, "Pixel Size [m]");
      } else if (i == 8) {
        sprintf(string, "Width (x) [m]");
      } else if (i == 9) {
        sprintf(string, "Width (y) [m] (x \\(2244) y)");
      } else if (i == 10) {
        sprintf(string, "RMS Value");
      } else if (i == 11) {
        sprintf(string, "Baseline Length [m]");
      } else if (i == 12) {
        sprintf(string, "BIAS");
      }
      cpgtext(bttn_box[I][0]-0.19,
              0.6*bttn_box[I][2]+0.4*bttn_box[I][3], string);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015,
              0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
              0.0, 1.0, param_char[i]);
      cpgsci(1);
    }

/*
-------------------------------
*/

    TV_menu_hatch(0.00, 0.45, 0.855, 0.98, 7, 4);

/*
-------------------------------
*/

    I = NAME_SECTION;

    i = 13;
    bttn_box[I][0] = 0.21 + 0.50 * (float)(i/7);
    bttn_box[I][1] = 0.44 + 0.50 * (float)(i/7);
    bttn_box[I][2] = 0.94 - 0.04*(float)(i%7);
    bttn_box[I][3] = bttn_box[I][2] + pitch;

    cpgsci(1);
    cpgtext(bttn_box[I][0]-0.19,
            0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
            "Save File Name\0");
    cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
    cpgsci(0);
    cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
            0.0, 1.0, fname);
    cpgsci(1);

/*
-------------------------------
*/

    I = CNTR_SECTION;
    bttn_box[I][0] = 0.14;
    bttn_box[I][1] = 0.32;
    bttn_box[I][2] = 0.05;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    if (DISP_SWT == true) {
      _on_button(&Bdum, "DISPLAY SCREEN\0", bttn_box[I]);
      cpgsci(15);
      cpgrect(box1, box2, box3, box4);
    } else {
      _off_button(&Bdum, "DISPLAY SCREEN\0", bttn_box[I]);
      cpgsci(0);
      cpgrect(box1, box2, box3, box4);
    }
    I++;
    bttn_box[I][0] = 0.54;
    bttn_box[I][1] = 0.70;
    bttn_box[I][2] = 0.05;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    _off_button(&Bdum, "SAVE\0", bttn_box[I]);
    I++;
    bttn_box[I][0] = 0.75;
    bttn_box[I][1] = 0.91;
    bttn_box[I][2] = 0.05;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    _off_button(&Bdum, "Return\0", bttn_box[I]);


    I++;
    bttn_box[I][0] = 0.32;
    bttn_box[I][1] = 0.33;
    bttn_box[I][2] = 0.65;
    bttn_box[I][3] = bttn_box[I][2] + 0.01;
    if (SCREEN_MODE == SINGLE_SCRN) {
      on_button( &i, "", bttn_box[I]);
    } else {
      off_button(&i, "", bttn_box[I]);
    }
    cpgsci(1);
    cpgtext(bttn_box[I][0]+0.020, bttn_box[I][3]-0.010, "Single screen");

    cpgtext(0.020, bttn_box[I][3]-0.010, "Screen generation mode");
    cpgsfs(2);
    cpgrect(0.010, 0.650, 0.630, 0.679);
    cpgsfs(1);

    I++;
    bttn_box[I][0] = 0.49;
    bttn_box[I][1] = 0.50;
    bttn_box[I][2] = 0.65;
    bttn_box[I][3] = bttn_box[I][2] + 0.01;
    if (SCREEN_MODE == MULTI__SCRN) {
      on_button( &i, "", bttn_box[I]);
    } else {
      off_button(&i, "", bttn_box[I]);
    }
    cpgsci(1);
    cpgtext(bttn_box[I][0]+0.020, bttn_box[I][3]-0.010, "Multi screen");
  }

/*
-------------------------------
*/

  while (1) {
    PROC_SWT = false;
    if (TV_SWT == true) {
      cpgcurs(cursor_pos, cursor_pos+1, string);

      for (i=3; i<13; i++) {
        I = ATM_SECTION + i;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          char_copy(string, param_char[i]);
          tv_get_param("char", cursor_pos, bttn_box[I],
                       pitch, string, 0.0, 0.0);
          char_copy(param_char[i], string);
          str_init(string, sizeof(string));
          sscanf(param_char[i], "%lf", &ps_param[i]);
        }
      }

/*
-------------------------------
*/

      I = CNTR_SECTION;
      if (_button_chk(cursor_pos, bttn_box[I]) == true) {
        _toggle_button(&DISP_SWT, "DISPLAY SCREEN\0", bttn_box[I]);
        if (DISP_SWT == true) {
          cpgsci(15);
        } else if (DISP_SWT == false) {
          cpgsci(0);
        }
        cpgrect(box1, box2, box3, box4);
      }

      I = CNTR_SECTION + 1;
      if (_button_chk(cursor_pos, bttn_box[I]) == true) {
        _on_button(&Bdum, "Calculating...\0", bttn_box[I]);
        cpgebuf();
        PROC_SWT = true;
      }

      I = CNTR_SECTION + 2;
      if (_button_chk(cursor_pos, bttn_box[I]) == true) {
        _on_button(&Bdum, "EXIT\0", bttn_box[I]);
        break;
      }

/*
-------------------------------
*/

      I = CNTR_SECTION + 3;
      if (_button_chk(cursor_pos, bttn_box[I]) == true) {
        _on_button(&Bdum, "", bttn_box[I]);
        _off_button(&Bdum, "", bttn_box[I+1]);
        SCREEN_MODE = SINGLE_SCRN;
      }

      I = CNTR_SECTION + 4;
      if (_button_chk(cursor_pos, bttn_box[I]) == true) {
        _on_button(&Bdum, "", bttn_box[I]);
        _off_button(&Bdum, "", bttn_box[I-1]);
        SCREEN_MODE = MULTI__SCRN;
      }

/*
-------------------------------
*/

      I = NAME_SECTION;
      if (_button_chk(cursor_pos, bttn_box[I]) == true) {
        tv_get_param("char", cursor_pos, bttn_box[I],
                     pitch, fname, 0, 0);
      }

/*
-------------------------------
*/

    } else if (TV_SWT == false) {
      while (1) {
        printf("#### Parameter Input Menu ####\n");
        printf("1. Screen Height [m] (%f)\n",    (float)ps_param[0]);
        printf("2. Wind Velocity [m/s] (%f)\n",  (float)ps_param[1]);
        printf("3. Position Angle [deg] (%f)\n", (float)ps_param[2]);
        printf("4. Inner Scale [m] (%f)\n",      (float)ps_param[3]);
        printf("5. Outer Scale [m] (%f)\n",      (float)ps_param[4]);
        printf("6. Inner Exponent (%f)\n",       (float)ps_param[5]);
        printf("7. Outer Exponent (%f)\n",       (float)ps_param[6]);
        printf("8. Pixel Size [m] (%f)\n",       (float)ps_param[7]);
        printf("9. Width (x) [m] (%f)\n",        (float)ps_param[8]);
        printf("a. Width (y) [m] (%f) (x<=y)\n", (float)ps_param[9]);
        printf("b. RMS Value (%f)\n",            (float)ps_param[10]);
        printf("c. Baseline Length [m] (%f)\n",  (float)ps_param[11]);
        printf("d. BIAS (%f)\n",                 (float)ps_param[12]);
        printf("\n");
        printf("e. SAVE\n");
        printf("0. EXIT\n");
        printf("\n");
        printf("Input Number : ");
        if (fgets(string, sizeof(string), stdin) == NULL) {
          printf("ERROR: PHASE_SCREEN_CHECK: Invalid input.\n");
          ps_model_out(ps_model_prm, param_char, 0, 0,
                       DISP_SWT, SCREEN_MODE, fname);
          return (-1);
        }
        if (string[0] == '0') {
          break;
        } else if (string[0] == 'e') {
          if (DISP_SWT == true) {
            sprintf(string, "y");
          } else if (DISP_SWT == false) {
            sprintf(string, "n");
          }
          printf("DISPLAY MODE (y/n) [CR->%s] : ", string);
          if (fgets(string, sizeof(string), stdin) == NULL) {
            printf("ERROR: PHASE_SCREEN_CHECK: Invalid input.\n");
            ps_model_out(ps_model_prm, param_char, 0, 0,
                         DISP_SWT, SCREEN_MODE, fname);
            return (-1);
          }
          if (string[0] == 'y') {
            DISP_SWT = true;
          } else if (string[0] == 'n') {
            DISP_SWT = false;
          }
          printf("SAVE FILE NAME [CR->%s] : ", fname);
          if (fgets(string, sizeof(string), stdin) == NULL) {
            printf("ERROR: PHASE_SCREEN_CHECK: Invalid input.\n");
            ps_model_out(ps_model_prm, param_char, 0, 0,
                         DISP_SWT, SCREEN_MODE, fname);
            return (-1);
          }
          if (string[0] != '\n') {
            char_copy(fname, string);
          }

          PROC_SWT = true;
          break;
        } else if (string[0] >= '1' && string[0] <= '9' ||
                   string[0] >= 'a' && string[0] <= 'd') {
          if (string[0] >= '1' && string[0] <= '9') {
            i = string[0] - '1';
          } else if (string[0] >= 'a' && string[0] <= 'd') {
            i = string[0] - 'a' + 9;
          }
          printf("Input parameter : ");
          if (fgets(string, sizeof(string), stdin) == NULL) {
            printf("ERROR: PHASE_SCREEN_CHECK: Invalid input.\n");
            ps_model_out(ps_model_prm, param_char, 0, 0,
                         DISP_SWT, SCREEN_MODE, fname);
            return (-1);
          }
          sscanf(string, "%lf", &ps_param[i]);
          sprintf(param_char[i], "%lf", ps_param[i]);
        }
      }
      if (string[0] == '0') {
        break;
      }
    }

/*
-------------------------------
*/

    if (PROC_SWT == true) {

      if (ps_param[8] < ps_param[9]) {
        printf("CAUTION: PHASE_SCREEN_CHECK: ");
        printf("Width(y) is greater than Width(x).");
        printf("Please check Width (x) and Width (y).\n");
        lftmp = ps_param[9];
        ps_param[9] = ps_param[8];
        ps_param[8] = lftmp;
        sprintf(param_char[8], "%lf", ps_param[8]);
        sprintf(param_char[9], "%lf", ps_param[9]);
      }

      wvc.H_d        = ps_param[0];
      wvc.v[0]       = ps_param[1] * cos(ps_param[2] / 180.0 * dpi);
      wvc.v[1]       = ps_param[1] * sin(ps_param[2] / 180.0 * dpi);
      wvc.i_scale[0] = ps_param[3];
      wvc.o_scale[0] = ps_param[4];
      wvc.i_expon    = ps_param[5];
      wvc.o_expon    = ps_param[6];
      wvc.pixel      = ps_param[7];

      fit_power_two_number(wvc.i_scale[0], wvc.pixel, &(wvc.i_scale[0]), &itmp);
      fit_power_two_number(wvc.o_scale[0], wvc.pixel, &(wvc.o_scale[0]), &itmp);
      ps_param[3] = wvc.i_scale[0];
      ps_param[4] = wvc.o_scale[0];

/*
--------
*/

      NX = (int)lrint(ps_param[8] / ps_param[7]);
      NY = (int)lrint(ps_param[9] / ps_param[7]);

      nx = NX;
      while (nx%2 == 0)
        nx /= 2;
      if (nx != 1) {
        nx = 1;
        while (nx < NX)
          nx *= 2;
        NX = nx;
      }

      ny = NY;
      while (ny%2 == 0)
        ny /= 2;
      if (ny != 1) {
        ny = 1;
        while (ny < NY)
          ny *= 2;
        NY = ny;
      }

/*
--------
*/

      ps_param[8]    = ps_param[7] * (double)NX;
      ps_param[9]    = ps_param[7] * (double)NY;
      diff_expon     = wvc.i_expon - wvc.o_expon;

      if (ps_param[11] <= wvc.i_scale[0]) {
        wvc.i_coeffi   = ps_param[10] / pow(ps_param[11],   0.5*wvc.i_expon);
        wvc.o_coeffi   = wvc.i_coeffi * pow(wvc.i_scale[0], 0.5*diff_expon );
        wvc.c_coeffi   = wvc.o_coeffi * pow(wvc.o_scale[0], 0.5*wvc.o_expon);
      } else if (ps_param[11] >  wvc.i_scale[0] &&
                 ps_param[11] <= wvc.o_scale[0]) {
        wvc.o_coeffi   = ps_param[10] / pow(ps_param[11],   0.5*wvc.o_expon);
        wvc.i_coeffi   = wvc.o_coeffi / pow(wvc.i_scale[0], 0.5*diff_expon );
        wvc.c_coeffi   = wvc.o_coeffi * pow(wvc.o_scale[0], 0.5*wvc.o_expon);
      } else {
        wvc.c_coeffi   = ps_param[10];
        wvc.o_coeffi   = wvc.c_coeffi / pow(wvc.o_scale[0], 0.5*wvc.o_expon);
        wvc.i_coeffi   = wvc.o_coeffi / pow(wvc.i_scale[0], 0.5*diff_expon );
      }

/*
-------------------------------
*/

      fit_power_two_number(wvc.o_scale[0], wvc.pixel, &(wvc.o_scale[0]), &i);
      ps_param[4] = wvc.o_scale[0];

/*
-------------------------------
*/

      if (SCREEN_MODE == SINGLE_SCRN &&
          (double)NX * (double)NY > (double)NMTRX * (double)NMTRX) {
        SCREEN_MODE = MULTI__SCRN;
        off_button(&i, "", bttn_box[CNTR_SECTION+3]);
        on_button (&i, "", bttn_box[CNTR_SECTION+4]);
      }

      if (SCREEN_MODE == SINGLE_SCRN) {
        NSCREEN = 1;
      } else if (SCREEN_MODE == MULTI__SCRN) {
        NSCREEN = (int)rint((double)NX * wvc.pixel / wvc.o_scale[0]);
        if (NSCREEN == 0) {
          NSCREEN = 1;
          on_button( &i, "", bttn_box[CNTR_SECTION+3]);
          off_button(&i, "", bttn_box[CNTR_SECTION+4]);
        } else if (NSCREEN != 0) {
          NX /= NSCREEN;
        }
      }

#ifdef __DEBUG__
      printf("__DEBUG__: Phase_Screen_Check: NSCREEN: %d   NX:%d\n",
             NSCREEN, NX);
#endif /* __DEBUG__ */

/*
------------------
*/

      if (DISP_SWT == true) {
        if (NY > 128) {
          ny = N128;
          n_bandle = NY / N128;
          if (NY % N128 != 0) {
            n_bandle++;
          }
        } else {
          ny = NY;
          n_bandle = 1;
        }

        nx = NX / n_bandle;
        if (NX % n_bandle != 0) {
          nx++;
        }

        if (nx == ny) {
          pgsvpxmin = 0.30;
          pgsvpxmax = 0.70;
          plot_size_y = pgsvpxmax - pgsvpxmin;
          pgsvpymin = 0.38 - 0.50 * plot_size_y;
          pgsvpymax = 0.38 + 0.50 * plot_size_y;
        } else {
          pgsvpxmin = 0.12;
          pgsvpxmax = 0.92;
          plot_size_y = (pgsvpxmax - pgsvpxmin) * ((float)ny / (float)nx);
          if (plot_size_y < 0.005) {
            plot_size_y = 0.005;
          }
          pgsvpymin = 0.38 - 0.50 * plot_size_y;
          pgsvpymax = 0.38 + 0.50 * plot_size_y;
        }

        NPLOT = (int)((box4-box3) / (plot_size_y + plot_space));

        delta_x = wvc.pixel * (float)n_bandle;
        delta_y = wvc.pixel * (float)n_bandle;

        tr[0] = -0.5 * delta_x;
        tr[1] = 0.0;
        tr[2] = delta_y;
        tr[3] = -0.5 * delta_y;
        tr[4] = delta_x;
        tr[5] = 0.0;
      }

/*
------------------
*/

      NOD         = NX * NY;
      seed_dist   = (double *)calloc(NX+1, sizeof(double));
      DS          = (double *)calloc(NOD,  sizeof(double));
      CONT_SWT_S  = false;

/*
------------------
*/

      for (nscrn=0; nscrn<NSCREEN; nscrn++) {

#ifdef __DEBUG__
        printf("__DEBUG__: Phase_Screen_Check: nscrn/NSCREEN : (%3d/%3d)\n",
                nscrn+1, NSCREEN);
#endif /* __DEBUG__ */

        corner_position[0][0] = 0;
        corner_position[0][1] = 0;
        corner_position[1][0] = NX - 1;
        corner_position[1][1] = NY - 1;

#ifdef __DEBUG__
        printf("__DEBUG__: Phase_Screen_Check: (%d, %d, %d, %d)\n",
                          corner_position[0][0], corner_position[0][1],
                          corner_position[1][0], corner_position[1][1]);
               
#endif /* __DEBUG__ */

        if ((COUNT_NOD = turbulent_phase_screen
               (NX, 0, seed_dist, &wvc,
                CONT_SWT_S, 0, IX, IY, DS,
                corner_position, 0)) == -1) {
          free (seed_dist);
          free (DS);
          printf("ERROR: PHASE_SCREEN_CHECK: ");
          printf("generating phase screen failed ");
          printf("in TURBULENT_PHASE_SCREEN.\n");
          ps_model_out(ps_model_prm, param_char, 0, 0,
                       DISP_SWT, SCREEN_MODE, fname);
          return (-1);
        } else {
          NOD = COUNT_NOD;
          ps_param[3] = wvc.i_scale[1];
          ps_param[4] = wvc.o_scale[1];
          for (i=3; i<5; i++) {
            I = ATM_SECTION + i;
            sprintf(param_char[i], "%lf", ps_param[i]);
            if (TV_SWT == true) {
              cpgsvp(0.0, 1.0, 0.0, 1.0);
              cpgswin(0.0, 1.0, 0.0, 1.0);
              cpgsci(1);
              cpgrect(bttn_box[I][0], bttn_box[I][1],
                      bttn_box[I][2], bttn_box[I][3]);
              cpgsci(0);
              cpgptxt(bttn_box[I][1]-0.015,
                      0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
                      0.0, 1.0, param_char[i]);
              cpgsci(1);
            }
          }

          for (i=8; i<10; i++) {
            I = ATM_SECTION + i;
            sprintf(param_char[i], "%lf", ps_param[i]);
            if (TV_SWT == true) {
              cpgsvp(0.0, 1.0, 0.0, 1.0);
              cpgswin(0.0, 1.0, 0.0, 1.0);
              cpgsci(1);
              cpgrect(bttn_box[I][0], bttn_box[I][1],
                      bttn_box[I][2], bttn_box[I][3]);
              cpgsci(0);
              cpgptxt(bttn_box[I][1]-0.015,
                      0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
                      0.0, 1.0, param_char[i]);
              cpgsci(1);
            }
          }
        }

/*
--------
*/

        if (NSCREEN == 1) {
          sprintf(data_file, "%s", fname);
        } else if (NSCREEN >= 2) {
          sprintf(data_file, "%s-%3d", fname, nscrn);
        }
        for (i=0; i<strlen(data_file); i++) {
          if (data_file[i] == '\0') {
            break;
          }
          if (data_file[i] == ' ') {
            data_file[i] = '0';
          }
        }

        if (ps_model_out(data_file, param_char, NX, NY,
                         DISP_SWT, SCREEN_MODE, fname) == -1) {
          printf("ERROR: phase_screen_check: cannot open file: %s\n",
                 data_file);
          ps_model_out(ps_model_prm, param_char, 0, 0,
                       DISP_SWT, SCREEN_MODE, fname);
          return (-1);
        } else {
          ofp = fopen(data_file, "a");
          fprintf(ofp, "# END OF HEADER\n");
          for (i=0; i<NX; i++) {
            for (j=0; j<NY; j++) {
              fprintf(ofp, "%d,%d,%lf\n",
                      NX*nscrn+i, j, DS[NY*i+j]+ps_param[12]);
            }
          }
          fclose (ofp);
        }

/*
------------------------------------------------
*/

        if (DISP_SWT == true) {
          if ((m_dist = (float *)calloc(nx*ny, sizeof (float))) == NULL) {
            printf("WARNING: ");
            printf("Memory alloc error for showing phase screen map.\n");
            return 1;
          }
          phase_screen_cull(NX, NY, DS, nx, ny, n_bandle,
                            &m_min, &m_max, m_dist);

          if (NPLOT == 0 || nscrn % NPLOT == 0) {
            cpgsci(15);
            cpgrect(box1, box2, box3, box4);
            cpgsci(1);
          }

          if (NSCREEN == 1 || NPLOT == 0 || NPLOT == 1) {
            pgsvpymin = 0.38 - 0.50 * plot_size_y;
            pgsvpymax = 0.38 + 0.50 * plot_size_y;
          } else {
            pgsvpymax = 0.60
                      - (float)(nscrn%NPLOT) * (plot_size_y + plot_space);
            pgsvpymin = pgsvpymax - plot_size_y;
            sprintf(string, "%3d/%d", nscrn+1, NSCREEN);
            cpgtext(0.02, pgsvpymin, string);
          }

/*
-- The first screen color gradation is applied to the following screen.
*/

          if (nscrn % NPLOT == 0) {
            dmin = m_min;
            dmax = m_max - m_min;
          }

/****
          printf("(dmin, dmax) = (%f,   %f)\n", dmin, dmax);
          printf("dmin: ");
          fgets(string, sizeof(string), stdin);
          sscanf(string, "%f", &dmin);
          printf("dmax: ");
          fgets(string, sizeof(string), stdin);
          sscanf(string, "%f", &dmax);
****/


          for (i=0; i<nx*ny; i++) {
            m_dist[i] -= dmin;
            m_dist[i] /= dmax;
          }

/*
--
*/

          cpgsvp(pgsvpxmin, pgsvpxmax, pgsvpymin, pgsvpymax);
          cpgswin(0.0, delta_x*nx, 0.0, delta_y*ny);

          palett(2, 1.0, 0.5);
          cpgimag(m_dist, ny, nx, 1, ny, 1, nx, 0.0, 1.0, tr);
          if (NPLOT <= 4) {
            cpgbox("BN", 0, 0, "BN", 0, 0);
          } else {
            cpgbox("BN", 0, 0, "B", 0, 0);
          }

          free (m_dist);
        }

        CONT_SWT_S = true;
      }

/*
------------------------------------------------
*/

      free (seed_dist);
      free (DS);

/*
------------------------------------------------
*/

      ps_model_out(ps_model_prm, param_char, 0, 0,
                   DISP_SWT, SCREEN_MODE, fname);

/*
------------------------------------------------
*/

      if (TV_SWT == true) {
        cpgsvp(0.0, 1.0, 0.0, 1.0);
        cpgswin(0.0, 1.0, 0.0, 1.0);
        I = CNTR_SECTION + 1;
        _off_button(&Bdum, "SAVE\0", bttn_box[I]);
      }
    }
  }

/*
------------------------------------------------
*/

  return ( 1);
}


int  ps_model_out(char *prm_fname, char param_char[][40],
                  int NX, int NY,
                  _Bool DISP_SWT, int SCREEN_MODE,
                  char *save_fname)
{
  FILE   *ifp;

  if ((ifp = fopen(prm_fname, "w")) != NULL) {
    fprintf(ifp, "SCREEN HEIGHT       %s\n", param_char[0]);
    fprintf(ifp, "WIND VELOCITY       %s\n", param_char[1]);
    fprintf(ifp, "POSITION ANGLE      %s\n", param_char[2]);
    fprintf(ifp, "INNER SCALE LENGTH  %s\n", param_char[3]);
    fprintf(ifp, "OUTER SCALE LENGTH  %s\n", param_char[4]);
    fprintf(ifp, "INNER EXPONENT      %s\n", param_char[5]);
    fprintf(ifp, "OUTER EXPONENT      %s\n", param_char[6]);
    fprintf(ifp, "PIXEL SIZE          %s\n", param_char[7]);
    fprintf(ifp, "SCREEN WIDTH (X)    %s\n", param_char[8]);
    fprintf(ifp, "SCREEN WIDTH (Y)    %s\n", param_char[9]);
    if (NX != 0) {
      fprintf(ifp, "SCREEN WIDTH (NX)   %d\n", NX);
    }
    if (NY != 0) {
      fprintf(ifp, "SCREEN WIDTH (NY)   %d\n", NY);
    }
    fprintf(ifp, "RMS VALUE           %s\n", param_char[10]);
    fprintf(ifp, "BASELINE OF THE RMS %s\n", param_char[11]);
    fprintf(ifp, "BIAS                %s\n", param_char[12]);
    fprintf(ifp, "DISPLAY SWITCH      %1d\n", DISP_SWT);
    fprintf(ifp, "SCREEN GENERATION   %1d\n", SCREEN_MODE);
    fprintf(ifp, "SAVE FILE NAME      %s\n", save_fname);
    fclose (ifp);
    return ( 1);
  } else {
    printf("CAUTION: PHASE_SCREEN_CHECK: ");
    printf("Input parameters cannot be saved in %s.", prm_fname);
    return (-1);
  }
}


int   phase_screen_cull(
             int      NX, int      NY, double  *A,
             int      MX, int      MY, int    n_bandle,
             float *Bmin, float *Bmax,
             float   *B)
{

  int    i, j, k, l, nxy;
  int    imax=0, itmp, nmax;
  float  fmax, ftmp, f1, f2;
  float  *b;

  nxy = n_bandle * n_bandle;

  *Bmin = (float)A[0];
  *Bmax = (float)A[0];

  for (i=0; i<MX; i++) {
    for (j=0; j<MY; j++) {
      fmax = 0.0;
      for (k=0; k<n_bandle; k++) {
        for (l=0; l<n_bandle; l++) {
          itmp = NY * (n_bandle * i + k) + n_bandle * j + l;
          if (itmp < NX*NY) {
            ftmp = fabsf(*(A + itmp));
            if (ftmp > fmax) {
              fmax = ftmp;
              imax = itmp;
            }
          }
        }
      }
      b = B + MY * i + j;
      *b = *(A + imax);
      if (*b > *Bmax) {
        *Bmax = *b;
      }
      if (*b < *Bmin) {
        *Bmin = *b;
      }
    }
  }

  return (1);
}
