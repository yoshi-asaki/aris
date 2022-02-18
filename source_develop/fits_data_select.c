#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


#define  SLCT_SECTION   10
#define  NAME_SECTION   20


int  fits_data_select(char fits_fname[][40], _Bool *fits_save_flag,
                      float sel_pos, float *cursor_pos, _Bool TV_SWT)
{
  char   string[100];
  int    i, j, I, J, ns;
  int    data_num;
  char   data_sel[3][40];
  float  bttn_box[30][4];
  float  pitch = 0.03;
  _Bool  Bdum;

  data_num = 3;
  for (i=0; i<data_num; i++) {
    fits_save_flag[i] = false;
    sprintf(fits_fname[i]+1, "ARIS-0%1d.FITS", i+1);
    fits_fname[i][0] = '!';
  }

  if (TV_SWT == false) {

    printf("Which data is selected : \n");
    printf("[1.target without P-R  2.Cal source  ");
    printf("3.target with P-R  0.RETURN] : ");
    if (fgets(string, sizeof(string), stdin) == NULL) {
      printf("ERROR: FITS_DATA_SELECT: Invalid input.\n");
      return (__NG__);
    }
    sscanf(string, "%d", &ns);
    if (ns == 0) {
      return (_EXIT_);
    } else {
      ns--;
    }
    fits_save_flag[ns] = true;

    printf("File name (CR->%s) : ", fits_fname[ns]+1);
    if (fgets(string, sizeof(string), stdin) == NULL) {
      printf("ERROR: FITS_DATA_SELECT: Invalid input.\n");
      return (__NG__);
    }
    if (string[0] != '\n') {
      char_copy(fits_fname[ns]+1, string);
    }
    printf("File name : %s\n", fits_fname[ns]+1);

  } else if (TV_SWT == true) {

    bttn_box[0][0] = 0.80;
    bttn_box[0][1] = 0.94;
    bttn_box[0][2] = sel_pos;
    bttn_box[0][3] = bttn_box[0][2] + pitch;
    _off_button(&Bdum, "SAVE\0", bttn_box[0]);
    bttn_box[1][0] = 0.80;
    bttn_box[1][1] = 0.94;
    bttn_box[1][2] = sel_pos - 0.035;
    bttn_box[1][3] = bttn_box[1][2] + pitch;
    _off_button(&Bdum, "QUIT\0", bttn_box[1]);

/*
-------------------------------
*/

    sprintf(data_sel[0], "Target without P-R");
    sprintf(data_sel[1], "Cal source");
    sprintf(data_sel[2], "Target with P-R");
    for (i=0; i<data_num; i++) {
      I = SLCT_SECTION + i;
      bttn_box[I][0] = 0.07;
      bttn_box[I][1] = 0.32;
      bttn_box[I][2] = sel_pos - (float)i * 0.035;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      _off_button(&fits_save_flag[i], data_sel[i], bttn_box[I]);
    }

/*
-------------------------------
*/

    for (i=0; i<data_num; i++) {
      I = NAME_SECTION + i;
      bttn_box[I][0] = 0.48;
      bttn_box[I][1] = 0.75;
      bttn_box[I][2] = bttn_box[SLCT_SECTION + i][2];
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      cpgsci(1);
      cpgtext(0.38, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3], "File Name\0");
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
              0.0, 1.0, fits_fname[i]+1);
      cpgsci(1);
    }

/*
-------------------------------
*/

    while (1) {
      cpgcurs(cursor_pos, cursor_pos+1, string);

      if (_button_chk(cursor_pos, bttn_box[0]) == true) {
        _on_button(&Bdum, "SAVE\0", bttn_box[0]);
        return (__GO__);
      }

      if (_button_chk(cursor_pos, bttn_box[1]) == true) {
        _on_button(&Bdum, "QUIT\0", bttn_box[1]);
        return (_EXIT_);
      }

      for (i=0; i<data_num; i++) {
        I = SLCT_SECTION + i;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          _toggle_button(&fits_save_flag[i], data_sel[i], bttn_box[I]);
        }
      }

      for (i=0; i<data_num; i++) {
        I = NAME_SECTION + i;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          tv_get_param("char", cursor_pos, bttn_box[I],
                       pitch, fits_fname[i]+1, 0, 0);
        }
      }
    }
  }
}
