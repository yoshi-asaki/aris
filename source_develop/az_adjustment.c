#include <stdio.h>
#include <math.h>
#include <aris.h>

int  az_adjustment(int GRT_NUM,
                   struct data_number data_num,
                   struct st_observable *int_obs[],
                   struct antenna_parameter *ant_prm)
{
  int    iant, iobs, I;
  double az_diff;

  for (iant=0; iant<GRT_NUM; iant++) {
    I = iant * data_num.nobs;
    for (iobs=data_num.sobs; iobs<data_num.eobs; iobs++) {
      if (int_obs[0][I].el >= ant_prm[iant].ELLIM &&
          int_obs[1][I].el >= ant_prm[iant].ELLIM) {
        break;
      }
      I++;
    }

    I = iant * data_num.nobs + iobs;
    for (++iobs; iobs<data_num.eobs; iobs++) {
      az_diff = diff(int_obs[0][I].az, int_obs[1][I].az);
      if (fabs(az_diff) > dpi) {
        int_obs[1][I].az += 2.0 * dpi * (az_diff / fabs(az_diff));
      }
      I++;
    }
  }

  return 1;
}
