#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <aris.h>

#define  NM  400

void pg_color_bar(float bias, float width,
                  float x1, float x2, float y1, float y2,
                  _Bool NEGATIVE_COMP_SWT, char *mode)
{
  int    i, N1;
  float  tx[2], ty[2];
  float  CR, CB, CG;
  float  E, BIAS, WIDTH;

/*
----------
*/

  ty[0] = 0.0;
  ty[1] = 1.0;
  WIDTH = width;

  cpgsvp(x1, x2, y1, y2);
  cpgswin(bias, bias+width, 0.0, 1.0);
  cpgsci(1);
  cpgbox("BCN", 0.0, 0, "BC", 0.0, 0);

  if (NEGATIVE_COMP_SWT == false) {
    E = (float)NM;
    for (i=0; i<NM; i++) {
      tx[0] = (float)i     / E;
      tx[1] = (float)(i+1) / E;
      set_color(tx[0], &CR, &CG, &CB, NEGATIVE_COMP_SWT, mode);
      cpgscr(10, CR, CG, CB);
      cpgsci(10);
      tx[0] = bias + tx[0] * width;
      tx[1] = bias + tx[1] * width;
      cpgrect(tx[0], tx[1], ty[0], ty[1]);
    }
  } else if (NEGATIVE_COMP_SWT == true) {
    if (bias >= 0.0) {
      N1 = 0;
      BIAS  = bias;
      WIDTH = width;
    } else {
      if (bias + width < 0.0) {
        N1 = NM + 1;
        BIAS  = fabs(bias);
      } else {
        N1 = (int)lrint((float)NM * (fabs(bias) / width));
        BIAS  = fabs(bias);
        WIDTH = width - BIAS;
      }
    }

    E = (float)N1;
    for (i=0; i<N1; i++) {
      tx[0] = (float)i     / E - 1.0;
      tx[1] = (float)(i+1) / E - 1.0;
      set_color(tx[0], &CR, &CG, &CB, NEGATIVE_COMP_SWT, mode);
      cpgscr(10, CR, CG, CB);
      cpgsci(10);
      tx[0] = tx[0] * BIAS;
      tx[1] = tx[1] * BIAS;
      cpgrect(tx[0], tx[1], ty[0], ty[1]);
    }

    if (bias < 0.0) {
      BIAS = 0.0;
    }
    E = (float)(NM - N1);
    for (i=N1; i<=NM; i++) {
      tx[0] = (float)(i - N1)     / E;
      tx[1] = (float)(i + 1 - N1) / E;
      set_color(tx[0], &CR, &CG, &CB, NEGATIVE_COMP_SWT, mode);
      cpgscr(10, CR, CG, CB);
      cpgsci(10);
      tx[0] = BIAS + tx[0] * WIDTH;
      tx[1] = BIAS + tx[1] * WIDTH;
      cpgrect(tx[0], tx[1], ty[0], ty[1]);
    }
  }

  cpgsci(1);
  return;
}
