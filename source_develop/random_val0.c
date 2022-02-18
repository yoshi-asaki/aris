#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define IM0   2147483647  /* 2^{31} - 1 */

/****
  random_val0 :  0.0  -  1.0
  rand(): 0  ---- 2^{31}-1 (2147483647)
****/

double  random_val0()
{
  unsigned int  im = IM0;
/**
  return (
      (double)rand() / (double)IM0
  );
**/
/**
  srand((unsigned)time(NULL));
**/
  return (
      (double)rand() / (double)im
  );
}
