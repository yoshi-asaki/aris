#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


int    GRT_TEC_TID(
         int GRT_NUM,         int nobs,
         double *ion_ds[],
         struct antenna_parameter      *ant_prm,
         struct phase_screen_parameter *ion,
         struct st_observable  *int_obs[], struct TID TID)
{
  int    i, I, iobs, ns, iant;
  double rho, el, arc, R, r, dz, ehta, S[3];
  double ant_XYZ[3];
  double z, Z;
  double tdum, pdum;

/*
==============================================================
*/

  for (iant=0; iant<GRT_NUM; iant++) {
    R = ion[iant].H_d + earth_radius;

    ant_XYZ[0] = ant_prm[iant].XYZ[0];
    ant_XYZ[1] = ant_prm[iant].XYZ[1];
    ant_XYZ[2] = ant_prm[iant].XYZ[2];

    for (iobs=0; iobs<nobs; iobs++) {

      I = iant * nobs + iobs;

      for (ns=0; ns<SRC_NUM; ns++) {
        if (int_obs[ns][I].el >= ant_prm[iant].ELLIM) {
          spherical_geometry(ion[iant].H_d, int_obs[ns][I].el, &z,
                             &rho, &el, &Z, &dz);
          r = sin(dz) / sin(z) * R;
          S[0] = r * cos(int_obs[ns][I].az) * cos(int_obs[ns][I].el);
          S[1] = r * sin(int_obs[ns][I].az) * cos(int_obs[ns][I].el);
          S[2] = r * sin(int_obs[ns][I].el);

          drotate(S, ant_prm[iant].LLH[1]-0.5*dpi, "x");
          drotate(S, ant_prm[iant].LLH[0]-0.5*dpi, "z");

          S[0] += ant_XYZ[0];
          S[1] += ant_XYZ[1];
          S[2] += ant_XYZ[2];

          arc = R * (0.5 * dpi - atan2(S[2], vlen2(S)));
          ehta = 2.0 * dpi * (TID.v[1] * (double)iobs - arc) / TID.lambda;

          *(ion_ds[ns] + I) += TID.amp * sin(ehta);
        }
      }
    }
  }

  return 1;
}
