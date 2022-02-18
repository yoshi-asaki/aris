#include <stdio.h>
#include <stdlib.h>
#include <aris.h>


void free_memory_block_double(double **dp, int M)
{
  int   i;

  for (i=0; i<M; i++) {
    free (dp[i]);
  }
  free (dp);
}
