#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"


double EPSIRON(double t)
{
  double epsiron;
  double dpi = 3.141592653589793238462643;
  double DPI;

  DPI = dpi / 180.0;

  epsiron = 23.44253 - 0.00013*t + 0.00256 * cos((249.0 - 19.3*t)*DPI)
             + 0.00015 * cos((198.0 + 720.0*t)*DPI);

  return epsiron * DPI;
}
