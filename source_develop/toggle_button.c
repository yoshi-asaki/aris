#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


void toggle_button(int *ERROR_FLAG, char *error_source, float *button_box)
{
  if (*ERROR_FLAG == 1) {
    off_button(ERROR_FLAG, error_source, button_box);
  } else if (*ERROR_FLAG == 0) {
    on_button(ERROR_FLAG, error_source, button_box);
  }
}
