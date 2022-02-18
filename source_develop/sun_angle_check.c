#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

double  sun_angle_check(int  *TimUTC, double obs_duration,
                        struct source_parameter *src,
                        struct source_parameter *sun)
{
  double inner_sun_angle;

/*
--------
*/

  sun_position(TimUTC, 0.5 * obs_duration, &(sun->RA), &(sun->DC));
  *(sun->s    ) = cos(sun->DC) * cos(sun->RA);
  *(sun->s + 1) = cos(sun->DC) * sin(sun->RA);
  *(sun->s + 2) = sin(sun->DC);
  inner_sun_angle = acos(inner_product3(sun->s, src->s2k));

  return inner_sun_angle;
}
