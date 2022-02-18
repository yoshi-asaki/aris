#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

int   phase_reference(int  ANT_NUM,   int  nfrq,
                      int  BGN_ANT_I, int  END_ANT_I,
                      int  BGN_ANT_J, int  END_ANT_J,
                      struct data_number data_num,
                      int nobs, int nswt,
                      double *nu,
                      struct fringe *frng[], float *fringe_weight[])
{
  int    i, I, J, ibase, iobs, IOBS, NOBS, iave;
  int    iant, jant, ifrq;
  int    *ref_obs;
  _Bool  AVE_SWT;
  double *ref_phs, dref, phs_jamp;
  double ar, ai, amp, phs, alpha;
  double ar_tmp, ai_tmp, w_ave;

/*
--------------------------------
*/

  alpha = nu[0] / nu[1];

/*
--------------------------------
*/

  if (nswt == 0) {
    for (iant=0; iant<ANT_NUM; iant++) {
      for (jant=iant+1; jant<ANT_NUM; jant++) {
        if (iant >= BGN_ANT_I && iant < END_ANT_I &&
            jant >= BGN_ANT_J && jant < END_ANT_J) {
          ibase = baseline_number(ANT_NUM, iant, jant);

          for (iobs=0; iobs<nobs; iobs++) {
            J = ibase * nobs + iobs;
            if (fringe_weight[0][J] > 0.0 && fringe_weight[1][J] > 0.0) {
              for (ifrq=0; ifrq<nfrq; ifrq++) {
                I = (ibase * nobs + iobs) * nfrq + ifrq;
                amp = sqrt(frng[1][I].rl*frng[1][I].rl
                         + frng[1][I].im*frng[1][I].im);
                if (nu[0] != nu[1]) {
                  phs = alpha * atan2(frng[1][I].im, frng[1][I].rl);
                  ar = amp * cos(phs);
                  ai = amp * sin(phs);
                  frng[2][I].rl =  frng[0][I].rl * ar
                                 + frng[0][I].im * ai;
                  frng[2][I].im = -frng[0][I].rl * ai
                                 + frng[0][I].im * ar;
                } else {
                  frng[2][I].rl =  frng[0][I].rl * frng[1][I].rl
                                 + frng[0][I].im * frng[1][I].im;
                  frng[2][I].im = -frng[0][I].rl * frng[1][I].im
                                 + frng[0][I].im * frng[1][I].rl;
                }
                frng[2][I].rl /= amp;
                frng[2][I].im /= amp;
              }
              fringe_weight[2][J] = 101.0;
            } else {
              fringe_weight[2][J] = 0.0;
            }
          }
        }
      }
    }

/*
--------------------------------
*/

  } else {
    NOBS = nobs / nswt + 1;
    if ((ref_phs = (double *)calloc(NOBS, sizeof(double))) == NULL) {
      printf("ERROR: PHASE_REFERENCE.\n");
      printf("ERROR: fail in allocating memories for ref_phs.\n");
      return -1;
    }
    if ((ref_obs = (int    *)calloc(NOBS, sizeof(int))) == NULL) {
      printf("ERROR: PHASE_REFERENCE.\n");
      printf("ERROR: fail in allocating memories for ref_obs.\n");
      free (ref_phs);
      return -1;
    }

    for (iant=0; iant<ANT_NUM; iant++) {
      for (jant=iant+1; jant<ANT_NUM; jant++) {
        if (iant >= BGN_ANT_I && iant < END_ANT_I &&
            jant >= BGN_ANT_J && jant < END_ANT_J) {
          ibase = baseline_number(ANT_NUM, iant, jant);

          iobs = 0;
          IOBS = 0;
          while (iobs < nobs) {
            AVE_SWT = true;
            ar = (double)0.0;
            ai = (double)0.0;
            iave = 0;

            while (1) {
              ar_tmp = 0.0;
              ai_tmp = 0.0;
              J = ibase * nobs + iobs;
              if (fringe_weight[1][J] > 0.0) {
                for (ifrq=0; ifrq<nfrq; ifrq++) {
                  I = (ibase * nobs + iobs) * nfrq + ifrq;
                  ar_tmp += frng[1][I].rl;
                  ai_tmp += frng[1][I].im;
                }
              } else {
                break;
              }

              ar += ar_tmp / (double)nfrq;
              ai += ai_tmp / (double)nfrq;
              iave++;
              iobs++;
              if (iobs >= nobs - 1) {
                AVE_SWT = false;
                break;
              }
            }

            if (AVE_SWT == false) {
              break;
            } else if (AVE_SWT == true) {
              if (iave != 0) {
                ar /= (double)iave;
                ai /= (double)iave;
                ref_phs[IOBS] = atan2(ai, ar);
                ref_obs[IOBS] = iobs - iave / 2;
                IOBS++;
              }

              while (1) {
                J = ibase * nobs + iobs;
                if (fringe_weight[1][J] <= 0.0) {
                  iobs++;
                  if (iobs >= nobs - 1) {
                    break;
                  }
                } else {
                  break;
                }
              }
            }
          }
          NOBS = IOBS;

/*
-----------------------------------------------------------
*/

          while (1) {
            AVE_SWT = false;
            for (iobs=ref_obs[0]+nswt/2; iobs<nobs; iobs++) {
              I = ibase * nobs + iobs;
              if (fringe_weight[0][I] > 0.0) {
                break;
              }
            }

            if (! (iobs < ref_obs[1])) {
              for (IOBS=0; IOBS<NOBS; IOBS++) {
                if (ref_obs[IOBS] > iobs) {
                  break;
                }
              }
              J = 0;
              for (I=IOBS; I<NOBS; I++) {
                ref_phs[J] = ref_phs[I];
                ref_obs[J] = ref_obs[I];
                J++;
              }
              NOBS -= IOBS;
              break;
            } else {
              break;
            }
          }

/*
-----------------------------------------------------------
*/

          for (IOBS=0; IOBS<NOBS-1; IOBS++) {
            if (ref_obs[IOBS+1] - ref_obs[IOBS] <= 2 * nswt) {
              if (ref_phs[IOBS+1] - ref_phs[IOBS] >= dpi) {
                phs_jamp = -2.0 * dpi;
              } else if (ref_phs[IOBS+1] - ref_phs[IOBS] <= -dpi) {
                phs_jamp =  2.0 * dpi;
              } else {
                phs_jamp =  0.0;
              }
              for (i=IOBS+1; i<NOBS; i++) {
                ref_phs[i] += phs_jamp;
              }
            }
          }
          if (nu[0] != nu[1]) {
            for (IOBS=0; IOBS<NOBS; IOBS++) {
              ref_phs[IOBS] *= alpha;
            }
          }

          for (IOBS=0; IOBS<NOBS-1; IOBS++) {
            if (ref_obs[IOBS+1] - ref_obs[IOBS] <= 2 * nswt) {
              if (ref_phs[IOBS+1] - ref_phs[IOBS] >= dpi) {
                phs_jamp = -2.0*dpi;
              } else if (ref_phs[IOBS+1] - ref_phs[IOBS] <= -dpi) {
                phs_jamp =  2.0*dpi;
              } else {
                phs_jamp =  0.0;
              }
              dref = (ref_phs[IOBS+1] - ref_phs[IOBS] + phs_jamp)
                                  / (double)nswt;

              for (iobs=ref_obs[IOBS]; iobs<ref_obs[IOBS+1]; iobs++) {
                J = ibase * nobs + iobs;
                if (fringe_weight[0][J] > 0.0) {
                  for (ifrq=0; ifrq<nfrq; ifrq++) {
                    I = (ibase * nobs + iobs) * nfrq + ifrq;
                    phs_jamp
                     = ref_phs[IOBS] + dref * (double)(iobs - ref_obs[IOBS]);
                    ar = cos(phs_jamp);
                    ai = sin(phs_jamp);
                    frng[2][I].rl =  frng[0][I].rl*ar + frng[0][I].im*ai;
                    frng[2][I].im = -frng[0][I].rl*ai + frng[0][I].im*ar;
                  }
                  fringe_weight[2][J] = 101.0;
                }
              }
            }
          }
        }
      }
    }
    free (ref_phs);
    free (ref_obs);
  }

  return 1;
}
