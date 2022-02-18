#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <aris.h>

void  in__src_proc(int INPUT_MODE, int *SEP_MODE, int *POS_MODE)
{
  *SEP_MODE = INPUT_MODE / 2;
  *POS_MODE = INPUT_MODE % 2;

  return;
}
