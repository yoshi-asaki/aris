#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


void _toggle_button(_Bool *SWT, char *error_source, float *button_box)
{
  if (*SWT == true) {
    _off_button(SWT, error_source, button_box);
  } else if (*SWT == false) {
    _on_button (SWT, error_source, button_box);
  }
}
