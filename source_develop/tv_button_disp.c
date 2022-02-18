#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <aris.h>


void tv_button_disp(float *bttn_box, char *txt_char)
{
  cpgsci(1);
  cpgrect(bttn_box[0], bttn_box[1], bttn_box[2], bttn_box[3]);
  cpgsci(0);
  cpgptxt(bttn_box[1]-0.015, text_bottom(bttn_box[2], bttn_box[3]),
          0.0, 1.0, txt_char);
  cpgsci(1);

  return;
}
