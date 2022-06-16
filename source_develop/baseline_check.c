#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <aris.h>

_Bool  baseline_check(int I,         int J,     int ns,
                      int BGN_ANT_I, int END_ANT_I,
                      int BGN_ANT_J, int END_ANT_J,
                      struct antenna_parameter *ant_prm)
{
  if (I >= BGN_ANT_I && I <  END_ANT_I &&
      J >= BGN_ANT_J && J <  END_ANT_J) {
    if (ns == -1) {
      return true;
    } else {
      if (ant_prm[I].WID[ns] >= 0 && ant_prm[J].WID[ns] >= 0) {
        return true;
      } else {
        printf("000000000000000000000000    %d    %d\n", ant_prm[I].WID[ns], ant_prm[J].WID[ns]);
        return false;
      }
    }
  } else {
    return false;
  }
}
