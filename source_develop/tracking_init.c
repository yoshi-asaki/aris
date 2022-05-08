#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <aris.h>


int  tracking_init(int   trk_num,
                   int   *priority,
                   int   *trk_priority,
                   char  trk_name[][10],
                   struct antenna_parameter *trk_pos)
{
  int    i, itrk, TRK_NUM, idum;
  int    wave_id[SRC_NUM_P1];

/*
--------
*/

  TRK_NUM = array_config(0,        ALL_ANT, wave_id, 0, &idum,
                         trk_pos, "", "aris_input/tracking_network.prm",
                         false, true);
  for (i=0; i<trk_num; i++) {
    for (itrk=0; itrk<TRK_NUM; itrk++) {
      if (strncmp(trk_name[i], trk_pos[itrk].IDC, strlen(trk_name[i])) == 0) {
        if (priority[i] == 0) {
          trk_priority[itrk] = trk_num;
        } else {
          trk_priority[itrk] = priority[i];
        }
        break;
      }
    }
  }
  for (itrk=0; itrk<TRK_NUM; itrk++) {
    sprintf(trk_name[itrk], "%s", trk_pos[itrk].IDC);
  }
  trk_priority_check(TRK_NUM, trk_priority, trk_name);

  return TRK_NUM;
}
