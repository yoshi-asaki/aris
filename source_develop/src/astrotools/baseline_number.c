#include <stdio.h>
#include <stdlib.h>

int  baseline_number(int ANT_NUM, int NX, int NY)
{
  int    i, ibase;

  ibase = 0;
  for (i=0; i<NX; i++) {
    ibase += ANT_NUM - i - 1;
  }
  ibase += NY - NX - 1;

  return (ibase);
}
