#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

void attitude_Q(int    BODY_X_SUN,
                struct source_parameter src,
                struct source_parameter sun,
                double *q)
{
  double srt_x[3], srt_y[3], srt_z[3];
  double r, A[9];

/*
--------
*/

  if (BODY_X_SUN == 0) {
    return;
  } else {
    srt_z[0] = src.s[0];
    srt_z[1] = src.s[1];
    srt_z[2] = src.s[2];

    if (BODY_X_SUN > 0) {
      srt_y[0] = srt_z[1] * sun.s[2] - srt_z[2] * sun.s[1];
      srt_y[1] = srt_z[2] * sun.s[0] - srt_z[0] * sun.s[2];
      srt_y[2] = srt_z[0] * sun.s[1] - srt_z[1] * sun.s[0];
    } else if (BODY_X_SUN < 0) {
      srt_y[0] = sun.s[1] * srt_z[2] - sun.s[2] * srt_z[1];
      srt_y[1] = sun.s[2] * srt_z[0] - sun.s[0] * srt_z[2];
      srt_y[2] = sun.s[0] * srt_z[1] - sun.s[1] * srt_z[0];
    }

    r = vlen3(srt_y);
    srt_y[0] /= r;
    srt_y[1] /= r;
    srt_y[2] /= r;
    srt_x[0] = srt_y[1] * srt_z[2] - srt_y[2] * srt_z[1];
    srt_x[1] = srt_y[2] * srt_z[0] - srt_y[0] * srt_z[2];
    srt_x[2] = srt_y[0] * srt_z[1] - srt_y[1] * srt_z[0];

    xyz2q(srt_x, srt_y, srt_z, A, q);
 
  }

  return;
}
