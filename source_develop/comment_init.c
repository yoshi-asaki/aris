#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


void comment_init(struct comment_param  *cmnt, char comment[][NCOMLEN],
                  _Bool TV_SWT)
{
  int    i, j;

  for (i=0; i<cmnt->ncol; i++) {
    for (j=0; j<NCOMLEN; j++) {
      comment[i][j] = 0;
    }
  }
  if (TV_SWT == true) {
    cpgsci(4);
    cpgrect(cmnt->xmin, cmnt->xmax, cmnt->ymin, cmnt->ymax);
    cpgsci(1);
  }
}
