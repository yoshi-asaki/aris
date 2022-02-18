#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


_Bool  output_star_position(char *RA,  char *DC,  double *s)
{
  int    i;
  _Bool  disp_flag;
  int    rahh, ramm, dcdd, dcmm;
  double rass, dcss;

/*
----------------------------
*/

  xyz2radec(s, &rahh, &ramm, &rass, &dcdd, &dcmm, &dcss);

  while (1) {
    sprintf(RA, "%2d %2d %9.6lf", rahh, ramm, rass);
    if (RA[0] == ' ') {
      RA[0] = '0';
    }
    if (RA[3] == ' ') {
      RA[3] = '0';
    }
    if (RA[6] == ' ') {
      RA[6] = '0';
    }
    while ((i=strlen(RA)-1) >= 10) {
      if (RA[i] == '0') {
        RA[i] = '\0';
      } else {
        break;
      }
    }

    disp_flag = false;
    sscanf(RA, "%2d %2d %lf", &rahh, &ramm, &rass);
    if (rass >= 60.0) {
      rass -= 60.0;
      ramm += 1;
      disp_flag = true;
    }
    if (ramm >= 60) {
      ramm -= 60;
      rahh += 1;
      disp_flag = true;
    }
    if (rahh >= 24) {
      rahh -= 24;
      disp_flag = true;
    }
    if (disp_flag == false) {
      break;
    }
  }

/*
----------------------------
*/

  while (1) {
    if (dcdd >= 0) {
      DC[0] = '+';
    } else if (dcdd < 0) {
      DC[0] = '-';
    }

    if (abs(dcdd) < 10) {
      sprintf(DC+1, "0%1d %2d %9.6lf", abs(dcdd), abs(dcmm), fabs(dcss));
    } else if (abs(dcdd) >= 10) {
      sprintf(DC+1,  "%2d %2d %9.6lf", abs(dcdd), abs(dcmm), fabs(dcss));
    }
    if (DC[4] == ' ') {
      DC[4] = '0';
    }
    if (DC[7] == ' ') {
      DC[7] = '0';
    }
    while ((i=strlen(DC)-1) >= 11) {
      if (DC[i] == '0') {
        DC[i] = '\0';
      } else {
        break;
      }
    }

    disp_flag = false;
    sscanf(DC, "%3d %2d %lf", &dcdd, &dcmm, &dcss);
    if (dcss >= 60.0) {
      dcss -= 60.0;
      dcmm += 1;
      disp_flag = true;
    }
    if (dcmm >= 60) {
      dcmm -= 60;
      dcdd += 1;
      disp_flag = true;
    }
    if (dcdd >= 90) {
      printf("ERROR: ");
      printf("output_star_position: ");
      printf("something is wrong with Declination: %s\n", DC);
      return false;
    }
    if (disp_flag == false) {
      break;
    }
  }

/*
----------------------------
*/

  return true;
}
