#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


void   position2uvw(double *u, double *v, double *w,
                    double *geometrical_delay,
                    double *geometrical_delay_rate,
                    double *X_xyz, double *Y_xyz,
                    double RA, double DC, double *s)
{
  double baseline[3];

  baseline[0] = X_xyz[0] - Y_xyz[0];
  baseline[1] = X_xyz[1] - Y_xyz[1];
  baseline[2] = X_xyz[2] - Y_xyz[2];
  *geometrical_delay_rate = OMEGA / speed_of_light *
                (-baseline[1] * s[0] + baseline[0] * s[1]);

  drotate(baseline, -0.5*dpi - RA, "z");
  drotate(baseline, -0.5*dpi + DC, "x");

  /*
     $B4QB,<T$O!"E75e>e$NE7BN%$%a!<%8$rN"B&(B($BFbB&(B)$B$+$i8+$F$$$k$3$H$K$J$k$,!"(B
     $BE7BN0LCV$KCV$/@\J?LL>e$G$O?eJ?J}8~$NA}2CJ}8~$,4QB,<T$N:8<j$NJ}8~$K(B
     $B$J$C$F$$$k!#(B
  */

  *u = baseline[0];
  *v = baseline[1];
  *w = baseline[2];
  *geometrical_delay   =  *w / speed_of_light;

}
