#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


void  srtatt_button_disp(int SRT_NUM,
                         int srtatt_num, int *srtatt_code,
                         int trk_num, float y_pos, float pitch,
                         char srtatt_name[][10], float bttn_box[][4])
{
  int    i, I, iant, idum;
  float  pgxmin[SRTMAX], pgxmax[SRTMAX], pgymin[SRTMAX], pgymax[SRTMAX];
  char   string[20];

/*
--------------------
*/

  I = 0;
  for (iant=0; iant<SRT_NUM; iant++) {
    for (i=0; i<srtatt_num; i++) {
      bttn_box[I][0] = 0.660 + 0.105 * (float)i;
      bttn_box[I][1] = bttn_box[I][0] + 0.100;
      bttn_box[I][2] = y_pos + 0.016 - (float)iant * 0.031;
      bttn_box[I][3] = bttn_box[I][2] + 0.750 * pitch;
      I++;
    }

    pgxmin[iant] = bttn_box[srtatt_num*iant][0]              - 0.150;
    pgxmax[iant] = bttn_box[srtatt_num*iant+srtatt_num-1][1] + 0.000;
    pgymin[iant] = bttn_box[srtatt_num*iant][2]              - 0.002;
    pgymax[iant] = bttn_box[srtatt_num*iant][3]              + 0.002;

    TV_menu_hatch(pgxmin[iant], pgxmax[iant], pgymin[iant], pgymax[iant], 0, 1);
    cpgsci(1);
    if (SRT_NUM == 1) {
      sprintf(string, "SRT Attitude");
      cpgtext(pgxmin[iant]+0.021,
              y_pos + 0.750 * pitch - (float)iant * 0.031, string);
    } else {
      sprintf(string, "SRT [%d] Attitude", iant+1);
      cpgtext(pgxmin[iant]+0.001,
              y_pos + 0.750 * pitch - (float)iant * 0.031, string);
    }
    for (i=0; i<srtatt_num; i++) {
      if (i == srtatt_code[iant]) {
        on_button(&idum, srtatt_name[i], bttn_box[srtatt_num*iant+i]);
      } else {
        off_button(&idum, srtatt_name[i], bttn_box[srtatt_num*iant+i]);
      }
    }
  }

  if (trk_num == 0) {
    TV_menu_hatch(pgxmin[0], pgxmax[0], pgymin[SRT_NUM-1], pgymax[0], 7, 4);
  }

  return;
}
