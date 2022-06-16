#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <aris.h>


int  antenna_selection(int *ANT_NUM, int *GRT_NUM, int *SRT_NUM,
                       int *wave_id,
                       double grt_elevation_limit,
                       struct antenna_parameter *ant_prm,
                       char   antenna_code[][10],
                       char   *antenna_file,
                       _Bool ERR_RESET_SWT)
{
  int    i, iant, itmp, I;

/*
----------
*/

  I = 0;
  for (iant=0; iant<*GRT_NUM; iant++) {
    array_config(NO_DEF_ARRAY,  ALL_ANT,    wave_id,         0,  &itmp,
                 ant_prm+I, antenna_code[iant], antenna_file,
                 true,  ERR_RESET_SWT);
    if (grt_elevation_limit > ant_prm[I].ELLIM) {
      ant_prm[I].ELLIM = grt_elevation_limit;
    }
    strncpy(antenna_code[I], (ant_prm+iant)->IDC,
            strlen((ant_prm+iant)->IDC));
    if ((i=strlen((ant_prm+iant)->IDC)) < 10) {
      antenna_code[I][i] = 0;
    }
    I++;
  }
  *GRT_NUM = I;
  *ANT_NUM = I;
  if (*SRT_NUM >= 1) {
    *ANT_NUM += array_config(NO_DEF_ARRAY,  ORBITING, wave_id,
                             *SRT_NUM,  &itmp,
                             ant_prm + *GRT_NUM, "", antenna_file,
                             false,   ERR_RESET_SWT);
  }

  for (iant=0; iant<*ANT_NUM; iant++) {
    ant_prm[iant].NOSTA = iant + 1;
  }

  return 1;
}
