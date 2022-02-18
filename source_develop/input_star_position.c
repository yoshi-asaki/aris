#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __DEBUG__
****/


void input_star_position(char *RA,  char *DC,  double *ra, double *dc)
{
  int    i;
  int    rah=0, ram=0, dcd=0, dcm=0;
  double ras=0, dcs=0;

/*
----------------------------
*/

  sscanf(RA, "%d %d %lf", &rah, &ram, &ras);
#ifdef __DEBUG__
  printf("__DEBUG__ : RA  %d %d %lf\n", rah, ram, ras);
#endif
  *ra = ((double)rah + (double)ram/60.0 + (double)ras/3600.0) * 15.0;
#ifdef __DEBUG__
  printf("__DEBUG__ : RA[deg] = %16.11lf\n", *ra);
#endif
  *ra *= (dpi / 180.0);

/*
----
*/

  sscanf(DC, "%d %d %lf", &dcd, &dcm, &dcs);
#ifdef __DEBUG__
  printf("__DEBUG__ : DC  %d %d %lf\n", dcd, dcm, dcs);
#endif
  i = 0;
  while (1) {
    if (DC[i] != ' ') {
      break;
    }
    i++;
  }
  if (DC[i] == '-') {
    *dc = -1.0*(double)dcd + (double)dcm/60.0 + (double)dcs/3600.0;
    *dc *= -1.0;
  } else {
    *dc =      (double)dcd + (double)dcm/60.0 + (double)dcs/3600.0;
  }
#ifdef __DEBUG__
  printf("__DEBUG__ : DC[deg] = %16.11lf\n", *dc);
#endif
  *dc *= (dpi / 180.0);

/*
----
*/

  return;
}
