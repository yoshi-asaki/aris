#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define IM1   1073741823  /* 2^{30} - 1 */

/****
  random_val1 : -1.0  -  1.0
  rand(): 0  ---- 2^{31}-1 (2147483647)
****/

double  random_val1()
{
  unsigned int  im = IM1;
/**
  return (
      (double)((double)rand() - (double)IM1) / (double)IM1
  );
**/
/**
  srand((unsigned)time(NULL));
**/
  return (
      (double)((double)rand() - (double)im) / (double)im
  );
}
