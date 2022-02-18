#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __DEBUG__
****/

_Bool  weight_flag_chk(float src_flag, float *weight, float flag)
{
  int    i, N, nflag, npart, nchk;
  float  wt;

  nflag = (int)fabs(flag);
  N = 32;
  for (i=5; i>=0 ; i--) {
    npart = nflag / N;
    if (npart == 1) {
      break;
    } else {
      nflag %= N;
      N /= 2;
    }
  }
  nchk = i;
#ifdef __DEBUG__
  printf("__DEBUG__   %f  %f  %f  %d\n", src_flag, *weight, flag, nchk);
#endif /* __DEBUG__ */

  wt = *weight - src_flag;
  if (wt < 0.0) {
    nflag = (int)fabs(wt);
    N = 32;
    for (i=5; i>=nchk; i--) {
      npart = nflag / N;
      nflag %= N;
      N /= 2;
    }
    if (npart == 1) {
      return (true);
    }
  }

  return (false);
}
