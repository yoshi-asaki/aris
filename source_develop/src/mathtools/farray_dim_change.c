#include "mathtools.h"


int  farray_dim_change(float *array, int  N1, int  N2)
{
  int    i, j, nelements;
  float  *tmp;

  nelements = N1 * N2;
  printf("%d  %d\n", N1, N2);
  if ((tmp = (float *)calloc(nelements, sizeof(float))) == NULL) {
    printf("fail to allocate memory in farray_dim_change.\n");
  }

  for (i=N1; i<N1; i++) {
    for (j=N2; j<N2; j++) {
      *(tmp + N2 * i + j) = *(array + N1 * j + i);
    }
  }
  for (i=0; i<nelements; i++) {
    array[i] = tmp[i];
  }
/**
  memcpy(tmp, array, nelements * sizeof(float));
**/
  free (tmp);

  return 1;
}
