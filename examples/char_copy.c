#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <phase_screen.h>


void char_copy(char *new_string, char *original_string)
{
  static int    i;

  i = 0;
  while (original_string[i] != '\n' && original_string[i] != '\0') {
    new_string[i] = original_string[i];
    i++;
  }
  new_string[i] = '\0';

  return;
}
