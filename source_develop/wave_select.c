#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <aris.h>

int   wave_select(int  wave_id, double *wave_length, double *nu)
{
  int    i, WID;
  FILE   *fp;
  char   string[200];
  _Bool  PROC_SWT;

/*
--------
*/

  *wave_length = 0.0;
  *nu          = 0.0;

  if (wave_id == L_BAND) {
    WID =  0;
  } else if (wave_id == S_BAND) {
    WID =  1;
  } else if (wave_id == C_BAND) {
    WID =  2;
  } else if (wave_id == X_BAND) {
    WID =  3;
  } else if (wave_id == KU_BAND) {
    WID =  4;
  } else if (wave_id == K_BAND) {
    WID =  5;
  } else if (wave_id == Q_BAND) {
    WID =  6;
  } else if (wave_id == W_BAND) {
    WID =  7;
  } else if (wave_id == BAND03) {
    WID =  8;
  } else if (wave_id == BAND04) {
    WID =  9;
  } else if (wave_id == BAND05) {
    WID = 10;
  } else if (wave_id == BAND06) {
    WID = 11;
  } else if (wave_id == BAND07) {
    WID = 12;
  } else if (wave_id == BAND08) {
    WID = 13;
  } else if (wave_id == BAND09) {
    WID = 14;
  } else if (wave_id == BAND10) {
    WID = 15;
  } else {
    WID = -1;
  }

/*
--------
*/

  if ((fp=fopen("aris_input/fixed_parameter.prm", "r")) == NULL) {
    printf("ERROR: WAVE_SELCTION: ./aris_input/fixed_parameter.prm.\n");
    return -1;
  }
  PROC_SWT = false;
  while (1) {
    if (fgets(string, sizeof(string), fp) == NULL) {
      break;
    }
    if (strncmp(string, "WAVE LENGTH L_BAND", 18) == 0 && wave_id == L_BAND) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH S_BAND", 18) == 0 &&
               wave_id == S_BAND) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH C_BAND", 18) == 0 &&
               wave_id == C_BAND) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH X_BAND", 18) == 0 &&
               wave_id == X_BAND) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH KU_BAND", 19) == 0 &&
               wave_id == KU_BAND) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH K_BAND", 18) == 0 &&
               wave_id == K_BAND) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH Q_BAND", 18) == 0 &&
               wave_id == Q_BAND) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH W_BAND", 18) == 0 &&
               wave_id == W_BAND) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND01", 18) == 0 &&
               wave_id == BAND01) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND02", 18) == 0 &&
               wave_id == BAND02) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND03", 18) == 0 &&
               wave_id == BAND03) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND04", 18) == 0 &&
               wave_id == BAND04) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND05", 18) == 0 &&
               wave_id == BAND05) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND06", 18) == 0 &&
               wave_id == BAND06) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND07", 18) == 0 &&
               wave_id == BAND07) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND08", 18) == 0 &&
               wave_id == BAND08) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND09", 18) == 0 &&
               wave_id == BAND09) {
      PROC_SWT = true;
      break;
    } else if (strncmp(string, "WAVE LENGTH BAND10", 18) == 0 &&
               wave_id == BAND10) {
      PROC_SWT = true;
      break;
    }
  }
  fclose (fp);

  if (PROC_SWT == true) {
    i = 0;
    while (1) {
      if (string[i++] == ':') {
        break;
      }
    }
    sscanf(string+i, "%lf", wave_length);
    *nu = speed_of_light / *wave_length;
  }

/*
--------
*/

  return WID;
}
