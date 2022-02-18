#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <cpgplot.h>
#include <aris.h>


void _off_button(_Bool *SWT, char *string, float *button_box)
{
  float  hx[2], hy[2];
  float  vx[2], vy[2];
  float  bx[2], by[2];
  float  step;

/*
--------
*/

  cpgsfs(1);
  step = 0.125 * (button_box[3] - button_box[2]);

  hx[0] = button_box[0];
  hx[1] = button_box[1];
  hy[0] = button_box[3] - step;
  hy[1] = button_box[3];

  vx[0] = button_box[0];
  vx[1] = button_box[0] + step;
  vy[0] = button_box[2];
  vy[1] = button_box[3];

  bx[0] = button_box[0] + step;
  bx[1] = button_box[1];
  by[0] = button_box[2];
  by[1] = button_box[3] - step;

  cpgscr(20, 0.40, 0.40, 0.40);
  cpgsci(20);
  cpgrect(bx[0], bx[1], by[0], by[1]);
  cpgscr(20, 0.20, 0.20, 0.20);
  cpgsci(20);
  cpgrect(hx[0], hx[1], hy[0], hy[1]);
  cpgrect(vx[0], vx[1], vy[0], vy[1]);

  cpgscr(20, 1.00, 1.00, 1.00);
  cpgsci(20);
  cpgptxt(0.50 * (button_box[0] + button_box[1]),
          0.70*button_box[2] + 0.30*button_box[3], 0.0, 0.5, string);
  *SWT = false;
}
