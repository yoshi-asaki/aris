#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <aris.h>

/****
#define __EL_INDEPENDENCY__
****/
#define __EL_INDEPENDENCY__

double airmass(double el)
{
#ifdef __EL_INDEPENDENCY__
  return (1.0);
#else
  return (1.0/sin(el));
#endif
}
