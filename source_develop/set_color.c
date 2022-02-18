#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <aris.h>


void set_color(float strength, float *CR, float *CG, float *CB,
               _Bool NEGATIVE_COMP_SWT, char *mode)
{
  if (strncmp(mode, "clr", 3) == 0) {
    if (NEGATIVE_COMP_SWT == false) {
      if (strength == 0.0) {
        *CR = 0.0;
        *CG = 0.0;
        *CB = 0.0;
      } else {
        *CR = strength;
        *CG = -2.0*pow((strength-0.5), 2.0) + 0.6;
        *CB = -0.5/0.64*pow((strength-0.3), 2.0) + 0.5;
      }
    } else if (NEGATIVE_COMP_SWT == true) {
      if (strength >= 0.0) {
        strength += 1.0;
        strength /= 2.0;
        *CR = strength;
        *CG = -2.0*pow((strength-0.5), 2.0) + 0.6;
        *CB = -0.5/0.64*pow((strength-0.3), 2.0) + 0.5;
      } else if (strength < 0.0) {
        *CR = 0.5 * (1.0 + strength);
        *CG = 0.5 * (1.0 + strength);
        *CB = 1.0 + strength;
      }
    }
  } else if (strncmp(mode, "w2b", 3) == 0) {
    if (NEGATIVE_COMP_SWT == false) {
      if (strength == 0.0) {
        *CR = 1.0;
        *CG = 1.0;
        *CB = 1.0;
      } else {
        *CR = pow(strength-1.0, 2.0);
        *CG = pow(strength-1.0, 2.0);
        *CB = pow(strength-1.0, 2.0);
      }
    } else {
      strength += 1.0;
      strength /= 2.0;
      *CR = pow(strength-1.0, 2.0);
      *CG = pow(strength-1.0, 2.0);
      *CB = pow(strength-1.0, 2.0);
    }
  } else if (strncmp(mode, "b2w", 3) == 0) {
    if (NEGATIVE_COMP_SWT == false) {
      if (strength == 0.0) {
        *CR = 0.0;
        *CG = 0.0;
        *CB = 0.0;
      } else {
/****
        *CR = 1.0 - pow(strength-1.0, 2.0);
        *CG = 1.0 - pow(strength-1.0, 2.0);
        *CB = 1.0 - pow(strength-1.0, 2.0);
****/
        *CR = strength;
        *CG = strength;
        *CB = strength;
      }
    } else {
      strength += 1.0;
      strength /= 2.0;
/****
      *CR = 1.0 - pow(strength-1.0, 2.0);
      *CG = 1.0 - pow(strength-1.0, 2.0);
      *CB = 1.0 - pow(strength-1.0, 2.0);
****/
      *CR = strength;
      *CG = strength;
      *CB = strength;
    }
  }

/********
  if (*CR < 0.0) {
    *CR = 0.0;
  }
  if (*CG < 0.0) {
    *CG = 0.0;
  }
  if (*CB < 0.0) {
    *CB = 0.0;
  }
  if (*CR > 1.0) {
    *CR = 1.0;
  }
  if (*CG > 1.0) {
    *CG = 1.0;
  }
  if (*CB > 1.0) {
    *CB = 1.0;
  }
********/
  return;
}
