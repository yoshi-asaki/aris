#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


int      trk_priority_check(int TRK_NUM, int *trk_priority, 
                            char trk_name[][10])
{
  int   i, j, ncnt;
  int   order_tmp;

  order_tmp = 1;
  for (i=1; i<=TRK_NUM; i++) {
    ncnt = 0;
    for (j=0; j<TRK_NUM; j++) {
      if (trk_priority[j] == i) {
        trk_priority[j] = order_tmp;
        ncnt++;
      }
    }

    if (ncnt != 0) {
      order_tmp++;
    }
  }

  return 1;
}
