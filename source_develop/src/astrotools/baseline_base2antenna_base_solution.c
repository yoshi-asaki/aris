#include "astrotools.h"


double   baseline_base2antenna_base_solution(
            int    ant_num,
            int    refant,
            double *bl_base_fringe,
            double *bl_weight,
            double *ant_base_phase,
            double *ant_base_amp,
            int    SOL_MODE,
            int    SOL_TYPE)
{

/****
Assuming for the phase that

- the number of antennas is ant_num,
- the number of baselines shoud be (ant_num X (ant_num - 1)) / 2,
- baseline order should be (1,2), (1,3),..,(1,N), (2,3), (2,4),.., (2,N), 
  (3,4),.., (N,N-1),
- if baseline_base data is flagged out, the real and imaginary parts 
  are all zero.

SOL_MODE
  0 : Only phase solution
  1 : Phase and amplitude (closure amplitude) solution
  2 : Only phase solution (Amplitude is baseline-base value.)

SOL_TYPE
  0<= : Antenna base solution
  0>  : Baseline base solution with respect to the reference antenna
****/

  int     tgtant;
  int     ibl, nbl, icnt, jcnt, iant, jant;
  int     ibl1, ibl2, ibl3, ibl4;
  int     ant1, ant2;
  int     N1, N2;
  int     ncnt, mcnt;
  int     NCNT;
  double  *wr, *wi;
  double  AR, AI;
  int     *bl_cnt, *dt_cnt[2];
  double  theta;
  double  closure_amp;
  double  *amp;
  double  phase_SNR_improvement_factor = 0.0;

/*
----
*/

  NCNT = 0;

  for (iant=0; iant<ant_num; iant++) {
    ant_base_phase[iant] = 0.0;
  }

  nbl = (ant_num * (ant_num - 1)) / 2;

  wr = (double *)calloc(nbl, sizeof(double));
  wi = (double *)calloc(nbl, sizeof(double));
  amp = (double *)calloc(nbl, sizeof(double));
  bl_cnt    = (int    *)calloc(ant_num, sizeof(int));
  dt_cnt[0] = (int    *)calloc(ant_num, sizeof(int));
  dt_cnt[1] = (int    *)calloc(ant_num, sizeof(int));
  for (ibl=0; ibl<nbl; ibl++) {
    amp[ibl] = sqrt(bl_base_fringe[2*ibl]   * bl_base_fringe[2*ibl]
                  + bl_base_fringe[2*ibl+1] * bl_base_fringe[2*ibl+1]);
  }

  for (tgtant=0; tgtant<ant_num; tgtant++) {
    if (tgtant != refant) {

      for (iant=0; iant<ant_num; iant++) {
        bl_cnt[iant]    = 0;
        dt_cnt[0][iant] = 0;
        dt_cnt[1][iant] = 0;
      }

      icnt = 0;
      for (ibl=0; ibl<nbl; ibl++) {
        if (bl_weight[ibl] > 0.0) {
          N1 = 2 * ibl;
          N2 = N1 + 1;
          baseline2antenna_number(ibl, ant_num, &ant1, &ant2);
          if (      ant1 == tgtant && ant2 == refant) {
            wr[icnt] =  bl_base_fringe[N1];
            wi[icnt] =  bl_base_fringe[N2];
            dt_cnt[bl_cnt[tgtant]][tgtant] = icnt;
            bl_cnt[tgtant] += 10;
            icnt++;

          } else if (ant1 == refant && ant2 == tgtant) {
            wr[icnt] =  bl_base_fringe[N1];
            wi[icnt] = -bl_base_fringe[N2];
            dt_cnt[bl_cnt[tgtant]][tgtant] = icnt;
            bl_cnt[tgtant] += 10;
            icnt++;

          } else if (SOL_TYPE >= 0) {
            if (ant1 == refant && ant2 != tgtant) {
              wr[icnt] =  bl_base_fringe[N1];
              wi[icnt] = -bl_base_fringe[N2];
              dt_cnt[bl_cnt[ant2  ]][ant2  ] = icnt;
              bl_cnt[ant2  ] += 1;
              icnt++;

            } else if (ant1 != tgtant && ant2 == refant) {
              wr[icnt] =  bl_base_fringe[N1];
              wi[icnt] =  bl_base_fringe[N2];
              dt_cnt[bl_cnt[ant1  ]][ant1  ] = icnt;
              bl_cnt[ant1  ] += 1;
              icnt++;

            } else if (ant1 == tgtant && ant2 != refant) {
              wr[icnt] =  bl_base_fringe[N1];
              wi[icnt] =  bl_base_fringe[N2];
              dt_cnt[bl_cnt[ant2  ]][ant2  ] = icnt;
              bl_cnt[ant2  ] += 1;
              icnt++;

            } else if (ant1 != refant && ant2 == tgtant) {
              wr[icnt] =  bl_base_fringe[N1];
              wi[icnt] = -bl_base_fringe[N2];
              dt_cnt[bl_cnt[ant1  ]][ant1  ] = icnt;
              bl_cnt[ant1  ] += 1;
              icnt++;
            }
          }
        }
      }

/*
----
*/

      AR   = 0.0;
      AI   = 0.0;
      mcnt = 0;
      ncnt = 0;

      for (iant=0; iant<ant_num; iant++) {
        if (bl_cnt[iant] == 10) {
          icnt = dt_cnt[0][iant];
          theta = atan2(wi[icnt], wr[icnt]);
          AR += cos(theta);
          AI += sin(theta);
          mcnt++;

        } else if (bl_cnt[iant] == 2) {
          icnt = dt_cnt[0][iant];
          jcnt = dt_cnt[1][iant];
          theta = atan2(wr[icnt] * wi[jcnt] + wi[icnt] * wr[jcnt], 
                        wr[icnt] * wr[jcnt] - wi[icnt] * wi[jcnt]);
          AR += cos(theta);
          AI += sin(theta);
          ncnt++;

        } else if (bl_cnt[iant] != 0 && bl_cnt[iant] != 1) {
          printf("#### There is something wrong in bl_cnt (%d, %d). Stop!\n",
                  bl_cnt[iant], iant);
          exit (-1);
        }
      }
      if (AR != 0.0 && AI != 0.0) {
        ant_base_phase[tgtant] = atan2(AI, AR);
        ant_base_amp[tgtant]   = 1.0;
        phase_SNR_improvement_factor += (double)(mcnt + 2 * ncnt);
        NCNT++;
      } else {
        ant_base_phase[tgtant] = 0.0;
        ant_base_amp[tgtant]   = 0.0;
      }

/*
----
*/

      if (SOL_MODE != 0) {
        if (refant < tgtant) {
          ibl1 = baseline_number(ant_num, refant, tgtant);
        } else {
          ibl1 = baseline_number(ant_num, tgtant, refant);
        }

        if (SOL_MODE == 1) {
          icnt = 0;
          closure_amp = 0.0;
          for (iant=0; iant<ant_num; iant++) {
            for (jant=0; jant<ant_num; jant++) {
              if (iant != jant) {
                if (iant != tgtant && iant != refant &&
                    jant != tgtant && jant != refant) {

                  if (iant < refant) {
                    ibl2 = baseline_number(ant_num, iant, refant);
                  } else {
                    ibl2 = baseline_number(ant_num, refant, iant);
                  }
                  if (jant < tgtant) {
                    ibl3 = baseline_number(ant_num, jant, tgtant);
                  } else {
                    ibl3 = baseline_number(ant_num, tgtant, jant);
                  }
                  if (iant < jant) {
                    ibl4 = baseline_number(ant_num, iant, jant);
                  } else {
                    ibl4 = baseline_number(ant_num, jant, iant);
                  }

                  if (bl_weight[ibl1] > 0.0 &&
                      bl_weight[ibl2] > 0.0 &&
                      bl_weight[ibl3] > 0.0 &&
                      bl_weight[ibl4] > 0.0 &&
                      amp[ibl1] > 0.0 &&
                      amp[ibl2] > 0.0 &&
                      amp[ibl3] > 0.0 &&
                      amp[ibl4] > 0.0) {
                    closure_amp +=
                       amp[ibl2] * amp[ibl3] / amp[ibl4] / amp[ibl1];
                    icnt++;
                  }
                }
              }
            }
          }
          ant_base_amp[tgtant] = closure_amp / (double)icnt;

/*
----
*/

        } else if (SOL_MODE == 2) {
          ant_base_amp[tgtant] = amp[ibl1];
        }
      }

/*
----
*/

    } else {
      ant_base_phase[tgtant] = 0.0;
      ant_base_amp[tgtant]   = 1.0;
    }
  }

/*
----
*/

  free (wr);
  free (wi);
  free (amp);
  free (bl_cnt);
  free (dt_cnt[0]);
  free (dt_cnt[1]);

  return (sqrt(phase_SNR_improvement_factor / (double)NCNT));
}
