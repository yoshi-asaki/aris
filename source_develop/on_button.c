#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <cpgplot.h>
#include <aris.h>


void on_button(int *SWT, char *string, float *button_box)
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
  hy[0] = button_box[2];
  hy[1] = button_box[2] + step;

  vx[0] = button_box[1] - step;
  vx[1] = button_box[1];
  vy[0] = button_box[2];
  vy[1] = button_box[3];

  bx[0] = button_box[0];
  bx[1] = button_box[1] - step;
  by[0] = button_box[2] + step;
  by[1] = button_box[3];

  cpgscr(20, 0.95, 0.95, 0.95);
  cpgsci(20);
  cpgrect(bx[0], bx[1], by[0], by[1]);
  cpgscr(20, 0.75, 0.75, 0.75);
  cpgsci(20);
  cpgrect(hx[0], hx[1], hy[0], hy[1]);
  cpgrect(vx[0], vx[1], vy[0], vy[1]);

  cpgscr(20, 0.0, 0.0, 0.0);
  cpgsci(20);
  cpgptxt(0.50 * (button_box[0] + button_box[1]),
          0.60*button_box[2] + 0.40*button_box[3], 0.0, 0.5, string);
  cpgsci(1);
  *SWT = 1;
}
