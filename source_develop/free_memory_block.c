#include <stdio.h>
#include <stdlib.h>
#include <aris.h>


void free_memory_block(void  *ptr[],  int M)
{
  int    i;

  for (i=0; i<M; i++) {
    free (ptr[i]);
  }
  return;
}
