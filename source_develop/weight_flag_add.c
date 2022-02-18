#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <aris.h>


void  weight_flag_add(float src_flag, float *weight, float flag)
{
  int    i, N, nflag, npart;
  float  wt, add_wt;
  _Bool  obs_flag[6];

  for (i=0; i<6; i++) {
    obs_flag[i] = false;
  }

  nflag = (int)fabs(flag);
  N = 32;
  for (i=5; i>=0; i--) {
    npart = nflag / N;
    if (npart == 1) {
      obs_flag[i] = true;
    }
    nflag %= N;
    N /= 2;
  }

  wt = *weight - src_flag;
  if (wt < 0.0) {
    nflag = (int)fabs(wt);
    N = 32;
    for (i=5; i>=0; i--) {
      npart = nflag / N;
      if (npart == 1) {
        obs_flag[i] = true;
      }
      nflag %= N;
      N /= 2;
    }
  }

  add_wt = 0.0;
  for (i=0; i<6; i++) {
    if (obs_flag[i] == true) {
      add_wt += pow(2.0, (float)i);
    }
  }

  if (add_wt > 0.0) {
    *weight = src_flag - add_wt;
  } else {
    *weight = src_flag + wt;
  }

  return;
}
