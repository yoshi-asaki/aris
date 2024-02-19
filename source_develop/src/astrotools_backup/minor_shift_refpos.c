#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"


#define __POS_MID__
#ifndef __POS_MID__
  #define __POS_REF__
#endif

int  minor_shift_refpos(double *cal0, double *calE,
                        double *tgt0, double *tgtE,
                        double *q,    double *A)
{
  double ds[3], e[3];
  double s0[3], t;
  double DS;
  double cos_half_fai, sin_half_fai, sin_2_half_fai, in_pro;

/*
-------------------------------
*/

#ifdef __POS_MID__
  s0[0] = tgt0[0] + cal0[0];
  s0[1] = tgt0[1] + cal0[1];
  s0[2] = tgt0[2] + cal0[2];
  t = vlen3(s0);
  s0[0] /= t;
  s0[1] /= t;
  s0[2] /= t;
#else
  s0[0] = cal0[0];
  s0[1] = cal0[1];
  s0[2] = cal0[2];
#endif

  ds[0] = calE[0] - cal0[0];
  ds[1] = calE[1] - cal0[1];
  ds[2] = calE[2] - cal0[2];

  e[0] = s0[1] * ds[2] - s0[2] * ds[1];
  e[1] = s0[2] * ds[0] - s0[0] * ds[2];
  e[2] = s0[0] * ds[1] - s0[1] * ds[0];
  t = vlen3(e);

  if (t == 0) {
    tgtE[0] = tgt0[0];
    tgtE[1] = tgt0[1];
    tgtE[2] = tgt0[2];
  } else {
    e[0] /= t;
    e[1] /= t;
    e[2] /= t;

    DS  = ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2];
    in_pro = cal0[0] * e[0] + cal0[1] * e[1] + cal0[2] * e[2];

    sin_2_half_fai = DS / 4.0 / (1.0 - in_pro*in_pro);
    cos_half_fai   = sqrt(1.0 - sin_2_half_fai);
    sin_half_fai   = sqrt(sin_2_half_fai);

    q[0] = e[0] * sin_half_fai;
    q[1] = e[1] * sin_half_fai;
    q[2] = e[2] * sin_half_fai;
    q[3] = cos_half_fai;
    q2xyz(q, A, tgt0, tgtE);
  }

  return 1;
}
