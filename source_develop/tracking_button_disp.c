#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cpgplot.h>
#include <aris.h>

int   tracking_button_disp(
                int    TRK_SECTION,
                float  y_pos, float pitch,
                float  bttn_box[][4],
                int    TRK_NUM,
                int    *trk_priority,
                char   trk_name[][10],
                struct antenna_parameter *trk_pos)
{
  int    i, I, itrk, trk_num;
  char   string[100];

/*
---------------------------------------------------------
*/

  cpgsch(1.5*pgpap_prm/13.0);

  cpgsfs(2);
  cpgsci(1);
  cpgrect(0.010, 1.400, y_pos-0.065, y_pos+0.040);
  cpgsfs(1);

  cpgsci(1);
  cpgtext(0.035, y_pos + 0.55 * pitch, "Tracking Network\0");

  I = TRK_SECTION;
  bttn_box[I][0] = 0.190;
  bttn_box[I][1] = bttn_box[I][0] + 0.095;
  bttn_box[I][2] = y_pos + 0.005;
  bttn_box[I][3] = bttn_box[I][2] + pitch;
  I++;
  bttn_box[I][0] = 0.190 + 0.105;
  bttn_box[I][1] = bttn_box[I][0] + 0.095;
  bttn_box[I][2] = y_pos + 0.005;
  bttn_box[I][3] = bttn_box[I][2] + pitch;
  I++;

  y_pos -= 0.059;
  for (itrk=0; itrk<TRK_NUM; itrk++) {
    bttn_box[I][0] = 0.035 + (float)itrk * 0.097;
    bttn_box[I][1] = bttn_box[I][0] + 0.093;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    I++;
  }
  I = TRK_SECTION;
  off_button(&i, "ALL ON\0", bttn_box[I]);
  I++;
  off_button(&i, "RESET\0",  bttn_box[I]);
  I++;

  cpgsci(8);
  cpgtext(0.04, y_pos+0.037,
    "[Tracking station list: ./aris_input/tracking_network.prm]\0");
  cpgsci(1);
  trk_num = 0;
  for (itrk=0; itrk<TRK_NUM; itrk++) {
    if (trk_priority[itrk] >= 1) {
      sprintf(string, "%s(%d)", trk_name[itrk], trk_priority[itrk]);
      on_button(&i, string, bttn_box[I]);
      trk_num++;
    } else {
      off_button(&i, trk_name[itrk], bttn_box[I]);
    }
    I++;
  }

/*
--------------------
*/

  return trk_num;
}
