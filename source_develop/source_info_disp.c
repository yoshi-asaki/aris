#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


void  source_info_disp(
        int   SEP_MODE,       int    POS_MODE,
        float  pitch,
        float bttn_box[][4],  float  y_pos,
        struct char_src_info  *ch_src)
{
  int    I;
  float  y_pos_d;
  char   string[100];

  y_pos_d = y_pos - 0.04;

  bttn_box[0][0] = 0.31;
  bttn_box[0][1] = 0.49;
  bttn_box[0][2] = y_pos;
  bttn_box[0][3] = bttn_box[0][2] + pitch;

  bttn_box[1][0] = 0.31;
  bttn_box[1][1] = 0.49;
  bttn_box[1][2] = y_pos_d;
  bttn_box[1][3] = bttn_box[1][2] + pitch;

  bttn_box[2][0] = 0.80;
  bttn_box[2][1] = 0.98;
  bttn_box[2][2] = y_pos;
  bttn_box[2][3] = bttn_box[2][2] + pitch;

  bttn_box[3][0] = 0.80;
  bttn_box[3][1] = 0.98;
  bttn_box[3][2] = y_pos_d;
  bttn_box[3][3] = bttn_box[3][2] + pitch;

  bttn_box[4][0] = 0.60;
  bttn_box[4][1] = 0.98;
  bttn_box[4][2] = y_pos + 0.04;
  bttn_box[4][3] = bttn_box[4][2] + pitch;

  cpgsci(1);
  if (SEP_MODE == SRC__RA__DEC) {
    I = 0;
    cpgtext(0.02, text_bottom(bttn_box[I][2], bttn_box[I][3]),
            "Source-1 RA (HH MM SS.)");
    tv_button_disp(bttn_box[I++], ch_src->tgt_ra);

    cpgtext(0.02, text_bottom(bttn_box[I][2], bttn_box[I][3]),
            "Source-1 DEC ([+/-]dd mm ss.)");
    tv_button_disp(bttn_box[I++], ch_src->tgt_dc);

    cpgtext(0.51, text_bottom(bttn_box[I][2], bttn_box[I][3]),
            "Source-2 RA (HH MM SS.)");
    tv_button_disp(bttn_box[I++], ch_src->ref_ra);

    cpgtext(0.51, text_bottom(bttn_box[I][2], bttn_box[I][3]),
            "Source-2 DEC ([+/-]dd mm ss.)");
    tv_button_disp(bttn_box[I  ], ch_src->ref_dc);

  } else {
    I = 0;
    if (       POS_MODE == SRC_POS1) {
      cpgtext(0.02, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              "Source-1 RA (HH MM SS.)");
      tv_button_disp(bttn_box[I++], ch_src->tgt_ra);

      cpgtext(0.02, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              "Source-1 DEC ([+/-]dd mm ss.)");
      tv_button_disp(bttn_box[I  ], ch_src->tgt_dc);
    } else if (POS_MODE == SRC_POS2) {
      cpgtext(0.02, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              "MID POINT RA (HH MM SS.)");
      tv_button_disp(bttn_box[I++], ch_src->mid_ra);

      cpgtext(0.02, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              "MID POINT DEC ([+/-]dd mm ss.)");
      tv_button_disp(bttn_box[I  ], ch_src->mid_dc);
    }

    I = 2;
    if (       SEP_MODE == SRC_dRA_dDEC) {
      cpgtext(0.65, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              "Delta RA  [deg]");
      tv_button_disp(bttn_box[I++], ch_src->dlt_ra);

      cpgtext(0.65, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              "Delta DEC [deg]");
      tv_button_disp(bttn_box[I], ch_src->dlt_dc);
    } else if (SEP_MODE == SRC_SEP_POSA) {
      cpgtext(0.60, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              "Separation Angle [deg]");
      tv_button_disp(bttn_box[I++], ch_src->sepang);

      cpgtext(0.60, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              "Position Angle [deg]");
      tv_button_disp(bttn_box[I], ch_src->posang);
    }
  }

  sprintf(string, "Separation Angle [deg]: %s", ch_src->sepang);
  cpgsci(4);
  cpgrect(bttn_box[4][0], bttn_box[4][1], bttn_box[4][2], bttn_box[4][3]); 
  cpgsci(7);
  cpgptxt(bttn_box[4][0]+0.015,
          text_bottom(bttn_box[4][2], bttn_box[4][3]),
          0.0, 0.0, string);
  cpgsci(1);

  return;
}
