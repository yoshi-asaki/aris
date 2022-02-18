#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <aris.h>

double slew_time(double az1, double az2, double az_lim,
                 double azsv, double azsa,
                 double el1, double el2, double el_lim,
                 double elsv, double elsa)
{
  double d_az, t_az, d_el, t_el;

  d_az = fabs(az1 - az2);
  if (d_az <= az_lim) {
    t_az = sqrt(d_az / azsa);
  } else {
    t_az = (d_az + az_lim) / azsv;
  }

  d_el = fabs(el1 - el2);
  if (d_el <= el_lim) {
    t_el = sqrt(d_el / elsa);
  } else {
    t_el = (d_el + el_lim) / elsv;
  }

  if (t_az > t_el) {
    return t_az;
  } else {
    return t_el;
  }
}
