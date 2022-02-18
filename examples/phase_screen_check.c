#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <phase_screen.h>



int  main()
{
  static int    i, j, NOD;
  static int    NSHEET, NS;
  static int    NMTRX=8192;
  static int    NX, NY, NX_TMP;
  static int    *IX, *IY;
  static int    corner_position[2][2];
  static int    COUNT_NOD;
  static double *DS, *wvdist, *seed_dist, seed_end[2];
  static struct phase_screen_parameter atm;

  static FILE   *ifp, *ofp;
  static char   string[100];
  static char   param_char[20][40];
  static char   fname[40], fname_tmp[40];

  static int    PROC_SWT;
  static int    CONTINUE_STR_SWT, CONTINUE_END_SWT;
  static double ps_param[20];

/*
---------------------- Atmospheric Turbulence 2D
*/

#ifdef __RANDOM_SEED__
  seed_random(ON);
#endif /*__RANDOM_SEED__*/

/*
----------------------
*/

  atm.H_d        = 1000.0;
  atm.v[0]       = 10.0;
  atm.v[1]       =  0.0;
  atm.i_scale[0] = 1024.0;
  atm.o_scale[0] = 8192.0;
  atm.i_expon    = 5.0 / 3.0;
  atm.o_expon    = 2.0 / 3.0;

  ps_param[0]  = atm.H_d;
  ps_param[1]  = sqrt(atm.v[0]*atm.v[0] + atm.v[1]*atm.v[1]);
  ps_param[2]  = 180.0 / dpi * atan2(atm.v[1], atm.v[0]);
  ps_param[3]  = atm.i_scale[0];
  ps_param[4]  = atm.o_scale[0];
  ps_param[5]  = atm.i_expon;
  ps_param[6]  = atm.o_expon;
  ps_param[7]  = atm.pixel;
  ps_param[8]  = 1024.0;
  ps_param[9]  = 1024.0;
  ps_param[11] = 128.0;
  ps_param[12] = 1.5;
  ps_param[10] = ps_param[12] * 0.02 / pow(ps_param[11], 0.5*ps_param[5]);

  for (i=0; i<13; i++) {
    sprintf(param_char[i], "%lf", ps_param[i]);
  }
  sprintf(fname, "screen.dat\0");

  if ((ifp = fopen("ps_model.prm", "r")) != NULL) {
    while (1) {
      if (fgets(string, sizeof(string), ifp) == NULL) {
        break;
      } else {
        string[strlen(string)-1] = '\0';
      }
      if (strncmp(string, "SCREEN HEIGHT       ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[0]);
        sscanf(string+20, "%lf", &ps_param[0]);
      } else if (strncmp(string, "WIND VELOCITY       ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[1]);
        sscanf(string+20, "%lf", &ps_param[1]);
      } else if (strncmp(string, "POSITION ANGLE      ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[2]);
        sscanf(string+20, "%lf", &ps_param[2]);
      } else if (strncmp(string, "INNER SCALE LENGTH  ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[3]);
        sscanf(string+20, "%lf", &ps_param[3]);
      } else if (strncmp(string, "OUTER SCALE LENGTH  ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[4]);
        sscanf(string+20, "%lf", &ps_param[4]);
      } else if (strncmp(string, "INNER EXPONENT      ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[5]);
        sscanf(string+20, "%lf", &ps_param[5]);
      } else if (strncmp(string, "OUTER EXPONENT      ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[6]);
        sscanf(string+20, "%lf", &ps_param[6]);
      } else if (strncmp(string, "GRID                ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[7]);
        sscanf(string+20, "%lf", &ps_param[7]);
      } else if (strncmp(string, "SCREEN WIDTH (X)    ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[8]);
        sscanf(string+20, "%lf", &ps_param[8]);
      } else if (strncmp(string, "SCREEN WIDTH (Y)    ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[9]);
        sscanf(string+20, "%lf", &ps_param[9]);
      } else if (strncmp(string, "RMS VALUE           ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[10]);
        sscanf(string+20, "%lf", &ps_param[10]);
      } else if (strncmp(string, "BASELINE OF THE RMS ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[11]);
        sscanf(string+20, "%lf", &ps_param[11]);
      } else if (strncmp(string, "BIAS                ", 20) == 0) {
        sscanf(string+20, "%s",  &param_char[12]);
        sscanf(string+20, "%lf", &ps_param[12]);
      } else if (strncmp(string, "SAVE FILE NAME      ", 20) == 0) {
        sscanf(string+20, "%s", fname);
      }
    }
    fclose (ifp);
  }

/*
-------------------------------
*/

  while (1) {

    PROC_SWT = OFF;

    while (1) {
      printf("#### Parameter Input Menu ####\n");
      printf("1. Screen Height [m] (%f)\n",    (float)ps_param[0]);
      printf("2. Wind Velocity [m/s] (%f)\n",  (float)ps_param[1]);
      printf("3. Position Angle [deg] (%f)\n", (float)ps_param[2]);
      printf("4. Inner Scale [m] (%f)\n",      (float)ps_param[3]);
      printf("5. Outer Scale [m] (%f)\n",      (float)ps_param[4]);
      printf("6. Inner Exponent (%f)\n",       (float)ps_param[5]);
      printf("7. Outer Exponent (%f)\n",       (float)ps_param[6]);
      printf("8. Grid [m] (%f)\n",             (float)ps_param[7]);
      printf("9. Width (x) [m] (%f)\n",        (float)ps_param[8]);
      printf("a. Width (y) [m] (%f)\n",        (float)ps_param[9]);
      printf("b. RMS Value (%f)\n",            (float)ps_param[10]);
      printf("c. Baseline Length [m] (%f)\n",  (float)ps_param[11]);
      printf("d. BIAS (%f)\n",                 (float)ps_param[12]);
      printf("\n");
      printf("e. SAVE THE SCREEN DATA\n");
      printf("0. EXIT\n");
      printf("\n");
      printf("Input Character : ");
      fgets(string, sizeof(string), stdin);
      if (string[0] == '0') {
        break;
      } else if (string[0] == 'e') {
        printf("SAVE FILE NAME [CR->%s] : ", fname);
        fgets(string, sizeof(string), stdin);
        if (string[0] != '\n') {
          char_copy(fname, string);
        }

        PROC_SWT = ON;
        break;
      } else if (string[0] >= '1' && string[0] <= '9' ||
                 string[0] >= 'a' && string[0] <= 'd') {
        if (string[0] >= '1' && string[0] <= '9') {
          i = string[0] - '1';
        } else if (string[0] >= 'a' && string[0] <= 'd') {
          i = string[0] - 'a' + 9;
        }
        printf("Input parameter : ");
        fgets(string, sizeof(string), stdin);
        sscanf(string, "%lf", &ps_param[i]);
        sprintf(param_char[i], "%lf\0", ps_param[i]);
      }
    }
    if (string[0] == '0') {
      break;
    }

/*
-------------------------------
*/

    if (PROC_SWT == ON) {

      atm.H_d        = ps_param[0];
      atm.v[0]       = ps_param[1] * cos(ps_param[2] / 180.0 * dpi);
      atm.v[1]       = ps_param[1] * sin(ps_param[2] / 180.0 * dpi);
      atm.i_scale[0] = ps_param[3];
      atm.o_scale[0] = ps_param[4];
      atm.i_expon    = ps_param[5];
      atm.o_expon    = ps_param[6];
      atm.pixel      = ps_param[7];

      NX             = (int)rint(ps_param[8] / ps_param[7]);
      NY             = (int)rint(ps_param[9] / ps_param[7]);
      ps_param[8]    = ps_param[7] * (double)NX;
      ps_param[9]    = ps_param[7] * (double)NY;

      atm.i_coeffi   = ps_param[10] / pow(ps_param[11], 0.5*ps_param[5]);
      atm.o_coeffi   = 0.0;
      atm.c_coeffi   = 0.0;

/*
-------------------------------
*/

      NSHEET = NX / NMTRX + 1;
      if (NSHEET > 999) {
        printf("phase_screen_check: cannot process with %f.\n",
               ps_param[9]);
        return (-1);
      }

/*
-------------------------------
*/

      NOD = NY * NMTRX;
      seed_dist   = (double *)calloc(NMTRX+1, sizeof(double));
      DS          = (double *)calloc(NOD,     sizeof(double));
      IX          = (int    *)calloc(NOD,     sizeof(int));
      IY          = (int    *)calloc(NOD,     sizeof(int));

      for (NS=0; NS<NSHEET; NS++) {

        corner_position[0][0] = 0;
        corner_position[0][1] = NMTRX/2 - NY    /2;
        corner_position[1][1] = NMTRX/2 + NY    /2 - 1;
        if (NS != NSHEET - 1) {
          NX_TMP = NMTRX;
          corner_position[1][0] = NMTRX - 1;
        } else {
          NX_TMP = NX % NMTRX;
          corner_position[1][0] = NX_TMP - 1;
        }

        NOD      = NX_TMP * NY;

        if (NS == 0) {
          CONTINUE_STR_SWT = OFF;
          CONTINUE_END_SWT = ON;
        } else {
          CONTINUE_STR_SWT = ON;
          CONTINUE_END_SWT = ON;
        }

        seed_end[0] = gauss_dev();
        seed_end[1] = gauss_dev();

        if ((COUNT_NOD = turbulent_phase_screen
               (NMTRX, 0, wvdist, seed_dist, seed_end, &atm,
                CONTINUE_STR_SWT, CONTINUE_END_SWT, 0, IX, IY, DS,
                corner_position, 0)) == -1) {
          free (seed_dist);
          free (IX);
          free (IY);
          free (DS);
          return (-1);
        } else {
          NOD = COUNT_NOD;
          ps_param[3] = atm.i_scale[1];
          ps_param[4] = atm.o_scale[1];
          for (i=3; i<5; i++) {
            sprintf(param_char[i], "%lf\0", ps_param[i]);
          }
          sprintf(param_char[8], "%lf\0", ps_param[7] * (double)NX_TMP);
          sprintf(param_char[9], "%lf\0", ps_param[7] * (double)NY);
        }

/*
--------
*/

        if (NSHEET != 0) {
          sprintf(fname_tmp, "%s.%3d\0", fname, NS);
          for (i=0; i<strlen(fname_tmp); i++) {
            if (fname_tmp[i] == ' ') {
              fname_tmp[i] = '0';
            }
          }
        } else {
          sprintf(fname_tmp, "%s\0", fname);
        }
        printf("FILE SAVE (%3d/%3d) : %s [%d X %d]\n",
               NS+1, NSHEET, fname_tmp, NX_TMP, NY);

        if ((ofp = fopen(fname_tmp, "w")) == NULL) {
          printf("ERROR: phase_screen_check: cannot open file\n");
          return (-1);
        }

        fprintf(ofp, "# SCREEN HEIGHT       %s\n", param_char[0]);
        fprintf(ofp, "# WIND VELOCITY       %s\n", param_char[1]);
        fprintf(ofp, "# POSITION ANGLE      %s\n", param_char[2]);
        fprintf(ofp, "# INNER SCALE LENGTH  %s\n", param_char[3]);
        fprintf(ofp, "# OUTER SCALE LENGTH  %s\n", param_char[4]);
        fprintf(ofp, "# INNER EXPONENT      %s\n", param_char[5]);
        fprintf(ofp, "# OUTER EXPONENT      %s\n", param_char[6]);
        fprintf(ofp, "# GRID                %s\n", param_char[7]);
        fprintf(ofp, "# SCREEN WIDTH (X)    %s\n", param_char[8]);
        fprintf(ofp, "# SCREEN WIDTH (Y)    %s\n", param_char[9]);
        fprintf(ofp, "# RMS VALUE           %s\n", param_char[10]);
        fprintf(ofp, "# BASELINE OF THE RMS %s\n", param_char[11]);
        fprintf(ofp, "# BIAS                %s\n", param_char[12]);
        fprintf(ofp, "# SAVE FILE NAME      %s\n", fname_tmp);
        fprintf(ofp, "# END OF HEADER       \n");
        for (i=0; i<NX_TMP; i++) {
          for (j=0; j<NY; j++) {
/****
            fprintf(ofp, "%4d,%4d,%lf\n", i, j, DS[NY*i+j]+ps_param[12]);
****/
            fprintf(ofp, "%4d %4d %lf\n", i+NS*NMTRX, j, DS[NY*i+j]+ps_param[12]);
          }
          fprintf(ofp, "\n");
        }
        fclose (ofp);
      }

/*
------------------------------------------------
*/

      free (seed_dist);
      free (IX);
      free (IY);
      free (DS);
    }
  }

/*
------------------------------------------------
*/

  sprintf(param_char[8], "%lf\0", ps_param[7] * (double)NX);
  sprintf(param_char[9], "%lf\0", ps_param[7] * (double)NY);

  if ((ifp = fopen("ps_model.prm", "w")) != NULL) {
    fprintf(ifp, "SCREEN HEIGHT       %s\n", param_char[0]);
    fprintf(ifp, "WIND VELOCITY       %s\n", param_char[1]);
    fprintf(ifp, "POSITION ANGLE      %s\n", param_char[2]);
    fprintf(ifp, "INNER SCALE LENGTH  %s\n", param_char[3]);
    fprintf(ifp, "OUTER SCALE LENGTH  %s\n", param_char[4]);
    fprintf(ifp, "INNER EXPONENT      %s\n", param_char[5]);
    fprintf(ifp, "OUTER EXPONENT      %s\n", param_char[6]);
    fprintf(ifp, "GRID                %s\n", param_char[7]);
    fprintf(ifp, "SCREEN WIDTH (X)    %s\n", param_char[8]);
    fprintf(ifp, "SCREEN WIDTH (Y)    %s\n", param_char[9]);
    fprintf(ifp, "RMS VALUE           %s\n", param_char[10]);
    fprintf(ifp, "BASELINE OF THE RMS %s\n", param_char[11]);
    fprintf(ifp, "BIAS                %s\n", param_char[12]);
    fprintf(ifp, "SAVE FILE NAME      %s\n", fname);
    fclose (ifp);
  } else {
    printf("CAUTION: Input parameters cannot be saved. ");
  }

  return ( 1);

}
