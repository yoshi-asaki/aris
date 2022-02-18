#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <aris.h>


void  position_on_screen(
             double  past_time,
             double *ant_off,
             struct phase_screen_parameter atm,
             double AZ,           double EL,
             double *soffx,       double  *soffy)
{
  double ar, ai, fr, fi, rho, theta_el, dz, fai;
  double z, Z;
  double XOFF, YOFF;

  XOFF = atm.v[0] * past_time;
  YOFF = atm.v[1] * past_time;

  spherical_geometry(atm.H_d, EL, &z, &rho, &theta_el, &Z, &dz);

  fr = rho * cos(AZ) + ant_off[0] - XOFF;
  fi = rho * sin(AZ) + ant_off[1] - YOFF;

  fai = atan2(atm.v[1], atm.v[0]) - dpi;
  ar = cos(fai);
  ai = sin(fai);

  *soffx =  fr * ar + fi * ai;
  *soffy = -fr * ai + fi * ar;
}
