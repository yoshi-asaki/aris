#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <phase_screen.h>

#define IM1   1073741823  /* 2^{30} - 1 */

/****
  random_val1 : -1.0  -  1.0
  rand(): 0  ---- 2^{31}-1 (2147483647)
****/

double  random_val1()
{
  return (
      (double)(rand() - IM1) / (double)IM1
  );
}
