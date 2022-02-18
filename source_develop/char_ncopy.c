#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <aris.h>


void char_ncopy(char *new_string, char *original_string, int  n)
{
  int    i;

  for (i=0; i<n; i++) {
    if (original_string[i] == '\n' || original_string[i] == '\0') {
      new_string[i] = '\0';
      break;
    }
    new_string[i] = original_string[i];
  }
  new_string[n] = '\0';

  return;
}
