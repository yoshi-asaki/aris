#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


void   srt_info_disp(int SRT_NUM,  float pitch,
                     float bttn_box[][4],    float y_offset,
                     char *error_source,     _Bool  *ERROR_FLAG,
                     struct char_srt_info    *ch_srt)
{
  int   i, I, idum;
  char  string[100];
  float x_pos = 0.305;

/*
------------
*/

  I = 0;

  bttn_box[I][0] = 0.02;
  bttn_box[I][1] = 0.20;
  bttn_box[I][2] = y_offset - 0.035;
  bttn_box[I][3] = bttn_box[I][2] + pitch;
  if (SRT_NUM >= 1) {
    if (*ERROR_FLAG == true) {
      on_button(&idum, error_source, bttn_box[I]);
    } else if (*ERROR_FLAG == false) {
      off_button(&idum, error_source, bttn_box[I]);
    }
  }
  I++;

  x_pos = 0.280;

  for (i=0; i<SRT_NUM; i++) {
    bttn_box[I][0] = x_pos;
    bttn_box[I][1] = x_pos + 0.095;
    bttn_box[I][2] = y_offset - 0.035 * (float)(i + 1);
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    tv_button_disp(bttn_box[I], ch_srt[i].apo);
    sprintf(string, "SRT%1d", i+1);
    cpgtext(x_pos-0.042, text_bottom(bttn_box[I][2], bttn_box[I][3]), string);
    if (i == 0) {
      cpgptxt(0.5 * (bttn_box[I][0] + bttn_box[I][1]),
              text_bottom(bttn_box[I][2], bttn_box[I][3]) + 0.03,
              0.0, 0.5, "Apo. [km]\0");
    }
    I++;

    bttn_box[I][0] = bttn_box[I - 1][1] + 0.005;
    bttn_box[I][1] = bttn_box[I    ][0] + 0.095;
    bttn_box[I][2] = bttn_box[I - 1][2];
    bttn_box[I][3] = bttn_box[I - 1][3];
    tv_button_disp(bttn_box[I], ch_srt[i].per);
    if (i == 0) {
      cpgptxt(0.5 * (bttn_box[I][0] + bttn_box[I][1]),
              text_bottom(bttn_box[I][2], bttn_box[I][3]) + 0.03,
              0.0, 0.5, "Peri. [km]\0");
    }
    I++;

    bttn_box[I][0] = bttn_box[I - 1][1] + 0.005;
    bttn_box[I][1] = bttn_box[I    ][0] + 0.095;
    bttn_box[I][2] = bttn_box[I - 1][2];
    bttn_box[I][3] = bttn_box[I - 1][3];
    tv_button_disp(bttn_box[I], ch_srt[i].inc);
    if (i == 0) {
      cpgptxt(0.5 * (bttn_box[I][0] + bttn_box[I][1]),
              text_bottom(bttn_box[I][2], bttn_box[I][3]) + 0.03,
              0.0, 0.5, "Incl. [deg]\0");
    }
    I++;

    bttn_box[I][0] = bttn_box[I - 1][1] + 0.005;
    bttn_box[I][1] = bttn_box[I    ][0] + 0.095;
    bttn_box[I][2] = bttn_box[I - 1][2];
    bttn_box[I][3] = bttn_box[I - 1][3];
    tv_button_disp(bttn_box[I], ch_srt[i].OMG);
    if (i == 0) {
      cpgptxt(0.5 * (bttn_box[I][0] + bttn_box[I][1]),
              text_bottom(bttn_box[I][2], bttn_box[I][3]) + 0.03,
              0.0, 0.5, "\\gW\\fn [deg]\0");
    }
    I++;

    bttn_box[I][0] = bttn_box[I - 1][1] + 0.005;
    bttn_box[I][1] = bttn_box[I    ][0] + 0.095;
    bttn_box[I][2] = bttn_box[I - 1][2];
    bttn_box[I][3] = bttn_box[I - 1][3];
    tv_button_disp(bttn_box[I], ch_srt[i].omg);
    if (i == 0) {
      cpgptxt(0.5 * (bttn_box[I][0] + bttn_box[I][1]),
              text_bottom(bttn_box[I][2], bttn_box[I][3]) + 0.03,
              0.0, 0.5, "\\gw\\fn [deg]\0");
    }
    I++;

    bttn_box[I][0] = bttn_box[I - 1][1] + 0.005;
    bttn_box[I][1] = bttn_box[I    ][0] + 0.190;
    bttn_box[I][2] = bttn_box[I - 1][2];
    bttn_box[I][3] = bttn_box[I - 1][3];
    tv_button_disp(bttn_box[I], ch_srt[i].t_0);
    if (i == 0) {
      cpgptxt(0.5 * (bttn_box[I][0] + bttn_box[I][1]),
              text_bottom(bttn_box[I][2], bttn_box[I][3]) + 0.03,
              0.0, 0.5, "t0(YYYYMMDDhhmmss)\0");
    }
    I++;

    bttn_box[I][0] = bttn_box[I - 1][1] + 0.005;
    bttn_box[I][1] = bttn_box[I    ][0] + 0.160;
    bttn_box[I][2] = bttn_box[I - 1][2];
    bttn_box[I][3] = bttn_box[I - 1][3];
    tv_button_disp(bttn_box[I], ch_srt[i].d_OMG);
    if (i == 0) {
      cpgptxt(0.5 * (bttn_box[I][0] + bttn_box[I][1]),
              text_bottom(bttn_box[I][2], bttn_box[I][3]) + 0.03,
              0.0, 0.5, "d\\gW\\fn/dt [deg/yr]\0");
    }
    I++;

    bttn_box[I][0] = bttn_box[I - 1][1] + 0.005;
    bttn_box[I][1] = bttn_box[I    ][0] + 0.160;
    bttn_box[I][2] = bttn_box[I - 1][2];
    bttn_box[I][3] = bttn_box[I - 1][3];
    tv_button_disp(bttn_box[I], ch_srt[i].d_omg);
    if (i == 0) {
      cpgptxt(0.5 * (bttn_box[I][0] + bttn_box[I][1]),
              text_bottom(bttn_box[I][2], bttn_box[I][3]) + 0.03,
              0.0, 0.5, "d\\gw\\fn/dt [deg/yr]\0");
    }
    I++;
  }

  return;
}
