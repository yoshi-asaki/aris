#include <math.h>
#include <phase_screen.h>

double gauss_dev(void )
{
  static double fac, rsq, v1, v2;

  do {
    v1 = random_val1();
    v2 = random_val1();
    rsq = v1*v1 + v2*v2;
  } while (rsq >= 1.0 || rsq == 0.0);
  fac = sqrt(-2.0 * log(rsq)/rsq);

  return v2*fac;
}
