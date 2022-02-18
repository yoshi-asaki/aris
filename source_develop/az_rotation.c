#include <stdio.h>
#include <math.h>
#include <aris.h>

int  az_rotation(int GRT_NUM,
                 struct data_number data_num,
                 struct st_observable *int_obs,
                 struct antenna_parameter *ant_prm)
{
  int    iant, iobs, I;
  double az_diff, az_tmp=0.0;

  for (iant=0; iant<GRT_NUM; iant++) {
    for (iobs=data_num.sobs; iobs<data_num.eobs; iobs++) {
      I = iant * data_num.nobs + iobs;
      if (int_obs[I].el >= ant_prm[iant].ELLIM) {
        az_tmp = int_obs[I].az;
        break;
      }
    }

    for (++iobs; iobs<data_num.eobs; iobs++) {
      I = iant * data_num.nobs + iobs;
      az_diff = diff(az_tmp, int_obs[I].az);
      if (fabs(az_diff) > dpi) {
        int_obs[I].az += 2.0 * dpi * (az_diff / fabs(az_diff));
      }
      az_tmp = int_obs[I].az;
    }
  }

  return 1;
}
