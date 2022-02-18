#include <stdio.h>
#include <aris.h>

void  str_init(char *string, int n)
{
  int    i;

  for (i=0; i<n; i++) {
    string[i] = 0;
  }

  return;
}
