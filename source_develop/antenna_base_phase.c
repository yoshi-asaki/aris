#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

int   antenna_base_phase(int  ANT_NUM,   int  nfrq,
                      int  BGN_ANT_I, int  END_ANT_I,
                      int  BGN_ANT_J, int  END_ANT_J,
                      struct antenna_parameter *ant_prm,
                      struct data_number data_num,
                      int nobs, int nswt,
                      double *nu,
                      struct fringe *frng[], float *fringe_weight[])
{
  int    i, I, J, ibase, iobs, IOBS, NOBS;
  int    iant, jant, ifrq;
  double ar, ai, amp, phs;
  double ar_tmp, ai_tmp, w_ave;

  int    refant;
  double *bl_base_fringe;
  double *bl_weight;
  double *ant_base_phase, *ant_base_amp;
  double improve_factor;
  double xtmp, ytmp, ztmp;
  double dist, DIST;

/*
--------------------------------
*/

  xtmp = 0.0;
  ytmp = 0.0;
  ztmp = 0.0;
  i = 0;
  for (iant=0; iant<ANT_NUM; iant++) {
    if (iant >= BGN_ANT_I && iant < END_ANT_I) {
      xtmp += ant_prm[iant].XYZ[0];
      ytmp += ant_prm[iant].XYZ[1];
      ztmp += ant_prm[iant].XYZ[2];
      i++;
    }
  }
  xtmp /= (double)i;
  ytmp /= (double)i;
  ztmp /= (double)i;
  refant = 0;
  DIST = 1.0e22;
  for (iant=0; iant<ANT_NUM; iant++) {
    if (iant >= BGN_ANT_I && iant < END_ANT_I) {
      dist = pow(xtmp - ant_prm[iant].XYZ[0], 2.0)
           + pow(ytmp - ant_prm[iant].XYZ[1], 2.0)
           + pow(ztmp - ant_prm[iant].XYZ[2], 2.0);
      if (dist < DIST) {
        dist = DIST;
        refant = iant;
      }
    }
  }

/*
--------
*/

  bl_base_fringe = (double *)calloc(ANT_NUM*(ANT_NUM-1), sizeof(double));
  bl_weight      = (double *)calloc(ANT_NUM*(ANT_NUM-1), sizeof(double));
  ant_base_amp   = (double *)calloc(ANT_NUM, sizeof(double));
  ant_base_phase = (double *)calloc(ANT_NUM, sizeof(double));

  for (iobs=0; iobs<nobs; iobs++) {
    for (ifrq=0; ifrq<nfrq; ifrq++) {
      for (iant=0; iant<ANT_NUM; iant++) {
        for (jant=iant+1; jant<ANT_NUM; jant++) {
          ibase = baseline_number(ANT_NUM, iant, jant);
          I = (ibase * nobs + iobs) * nfrq + ifrq;
          if (iant >= BGN_ANT_I && iant < END_ANT_I &&
              jant >= BGN_ANT_J && jant < END_ANT_J) {
            J = ibase * nobs + iobs;
            if (fringe_weight[1][J] > 0.0) {
              bl_base_fringe[2*ibase    ] = frng[1][I].rl;
              bl_base_fringe[2*ibase + 1] = frng[1][I].im;
              bl_weight[ibase]            =  1.0;
            }
          } else {
            bl_base_fringe[2*ibase    ] = frng[1][I].rl;
            bl_base_fringe[2*ibase + 1] = frng[1][I].im;
            bl_weight[ibase]            =  1.0;
          }
        }
      }
      improve_factor = baseline_base2antenna_base_solution
                             ((size_t)ANT_NUM, refant, bl_base_fringe,
                              bl_weight, ant_base_phase, ant_base_amp, 0, 0);
      for (iant=0; iant<ANT_NUM; iant++) {
        for (jant=iant+1; jant<ANT_NUM; jant++) {
          ibase = baseline_number(ANT_NUM, iant, jant);
          I = (ibase * nobs + iobs) * nfrq + ifrq;
          if (iant >= BGN_ANT_I && iant < END_ANT_I &&
              jant >= BGN_ANT_J && jant < END_ANT_J) {
            J = ibase * nobs + iobs;
            if (fringe_weight[1][J] > 0.0) {
              xtmp = ant_base_phase[iant] - ant_base_phase[jant];
              ytmp = sqrt(frng[1][I].rl*frng[1][I].rl
                        + frng[1][I].im*frng[1][I].im);
              frng[1][I].rl = ytmp * cos(xtmp);
              frng[1][I].im = ytmp * sin(xtmp);
            } else {
              frng[1][I].rl = 0.0;
              frng[1][I].im = 0.0;
            }
          }
        }
      }
    }
  }
  free (bl_base_fringe);
  free (bl_weight);
  free (ant_base_amp  );
  free (ant_base_phase);

  return 1;
}
