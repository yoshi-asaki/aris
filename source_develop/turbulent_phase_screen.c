#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <stdbool.h>
#include <aris.h>

#define __INDEPENDENT__   0
#define ___DEPENDENT___   1

/**
#define __TPS_DEBUG0__
#define __TPS_DEBUG1__
#define __TPS_DEBUG2__
#define __TPS_DEBUG3__
#define __TPS_DEBUG4__
#define __TPS_DEBUG5__
#define __TPS_DEBUG6__
#define __TPS_DEBUG7__
**/

/****
    f(x)  ~  x^2 + (N/x)^2   --> minimum
                 |
                 V
    df/dx = 2 (x^3 - N^2/2) / x^2

     x = (N^2 / 2)^{1/3}  for minimum f(x)
****/


void   segment_calc(double *, double , int , int , int , int , int , int );
int    area_sort(int , int , int *, int *, int *, int *, int );
double excess_value(double );


int    turbulent_phase_screen
                   (int NDAT, int NDAT_SEGM,
                    double *seed_dist,
                    struct phase_screen_parameter *atm,
                    _Bool  CONTINUE_STR_SWT,
                    int    Num_of_Samp, int *NX, int  *NY, double *DS,
                    int    corner_position[][2], int Y_SHIFT)
{
  int    i, j, k, I, ndiv, iarea, N1, N2;
  int    NLOG, nlog, ilog;
  int    NC, M, NS=0, Num_of_Area;
  int    n_inner_log,  n_outer_log;
  int    nh_inner_log, nh_outer_log;
  int    n_outer;
  int    nsegm_P1, narea_P1;
  int    *nx, *ny;
  int    nsx, nsy, NDAT_AREA;
  int    RAND_MODE;
  _Bool  *AREA_SWT;

  double *wx;
  double *amp, amp0;
  double amp_factor=1.0;
  double D, d_expon;
  double wd[2][2] = {0.0, 0.0, 0.0, 0.0};
  double *init_dist;
  double lftmp;

  double *seg_dist;
  double *x_lods, *y_lods;
  double *new_seed;

  int    nxmin, nxmax, nymin, nymax, COUNT_Num_of_Samp=0;
  int    ix, iy, y_len=0.0;

/*
-------------------------------------
*/

  M = NDAT + 1;
  d_expon = atm->i_expon - atm->o_expon;

/*
-------------------------------------
*/

  if (Num_of_Samp != 0) {
    for (i=0; i<Num_of_Samp; i++) {
      if (NX[i] < 0 || NX[i] >= NDAT || NY[i] < 0 || NY[i] >= NDAT) {
        printf("ERROR: Turbulent_Phase_Screen: NX=%d, NY=%d (%d)\n",
                NX[i], NY[i], NDAT);
        exit (-1);
      }
    }
  }

/*
-------------------------------------
*/

  fit_power_two_number(atm->i_scale[0], atm->pixel, &lftmp, &n_inner_log);
  fit_power_two_number(atm->o_scale[0], atm->pixel, &lftmp, &n_outer_log);

  nh_inner_log = n_inner_log - 1;
  nh_outer_log = n_outer_log - 1;

  atm->i_scale[1] = atm->pixel * pow(2.0, (double)n_inner_log);
  atm->o_scale[1] = atm->pixel * pow(2.0, (double)n_outer_log);

  n_outer = 1;
  for (ilog=0; ilog<n_outer_log; ilog++) {
    n_outer *= 2;
  }

  NLOG = 0;
  i = 1;
  while (i <= NDAT) {
    i *= 2;
    NLOG++;
  }
  if ((amp = (double *)calloc(NLOG, sizeof(double))) == NULL) {
    printf("fail to allocate memories for amp.\n");
    exit (-1);
  }

/*
-------------------------------------------
*/

  amp0 = 9.0 * pow(atm->i_scale[1], d_expon)
             * pow(atm->o_scale[1], atm->o_expon);

  for (ilog=0; ilog<NLOG; ilog++) {
    I = NLOG -1 - ilog;
    if (I >= n_outer_log) {
      amp[ilog] = amp0;
#ifdef __TPS_DEBUG0__
      printf("__TPS_DEBUG0[0]__ %3d   %lf\n", ilog, amp[ilog]);
#endif /*__TPS_DEBUG0__*/
    } else if (I < n_outer_log  && I > nh_outer_log) {
      D = atm->pixel * pow(2.0, (double)I);
      amp[ilog] = 3.0 * pow(atm->i_scale[1], d_expon)
              * (4.0*pow(D, atm->o_expon)
                      - pow(atm->o_scale[1], atm->o_expon));
#ifdef __TPS_DEBUG0__
      printf("__TPS_DEBUG0[1]__ %3d   %lf\n", ilog, amp[ilog]);
#endif /*__TPS_DEBUG0__*/
    } else if (I <= nh_outer_log && I > n_inner_log) {
      D = atm->pixel * pow(2.0, (double)I);
      amp[ilog] = 3.0 * pow(atm->i_scale[1], d_expon)
                   * pow(D, atm->o_expon) * (4.0 - pow(2.0, atm->o_expon));
#ifdef __TPS_DEBUG0__
      printf("__TPS_DEBUG0[2]__ %3d   %lf\n", ilog, amp[ilog]);
#endif /*__TPS_DEBUG0__*/
    } else if (I <= n_inner_log  && I > nh_inner_log) {
      D = atm->pixel * pow(2.0, (double)I);
      amp[ilog] = 3.0
        * (4.0 * pow(D, atm->i_expon)
               - pow(atm->i_scale[1], d_expon) * pow(2.0*D, atm->o_expon));
#ifdef __TPS_DEBUG0__
      printf("__TPS_DEBUG0[3]__ %3d   %lf\n", ilog, amp[ilog]);
#endif /*__TPS_DEBUG0__*/
    } else if (I <= nh_inner_log) {
      D = atm->pixel * pow(2.0, (double)I);
      amp[ilog] = 3.0 * pow(D, atm->i_expon) * (4.0 - pow(2.0, atm->i_expon));
#ifdef __TPS_DEBUG0__
      printf("__TPS_DEBUG0[4]__ %3d   %lf\n", ilog, amp[ilog]);
#endif /*__TPS_DEBUG0__*/
    }
  }

  if (atm->i_coeffi != 0.0) {
    amp_factor = atm->i_coeffi;
  } else if (atm->o_coeffi != 0.0) {
    amp_factor = atm->o_coeffi / pow(atm->i_scale[1], 0.5*atm->i_expon);
  } else if (atm->c_coeffi != 0.0) {
    amp_factor = atm->c_coeffi / pow(atm->i_scale[1], 0.5*d_expon)
                               / pow(atm->o_scale[1], 0.5*atm->o_expon);
  }

  for (ilog=0; ilog<NLOG; ilog++) {
    amp[ilog] = amp_factor / 2.0 * sqrt(amp[ilog]);
#ifdef __TPS_DEBUG0__
    printf("__TPS_DEBUG0[-]__ %3d   %lf\n", ilog, amp[ilog]);
#endif /*__TPS_DEBUG0__*/
  }
/****
  getchar();
****/

/*
---------------------------------------------
*/

  if (CONTINUE_STR_SWT == true) {
    NS = 1;
  } else if (CONTINUE_STR_SWT == false) {
    NS = 0;
  }

  wd[0][0] = excess_value(amp[0]);
  wd[0][1] = excess_value(amp[0]);
  wd[1][0] = excess_value(amp[0]);
  wd[1][1] = excess_value(amp[0]);

#ifdef __TPS_DEBUG1__
  printf("__TPS_DEBUG1__  %lf   %lf  %lf  %lf\n",
         wd[0][0], wd[0][1], wd[1][0], wd[1][1]);
#endif /* __TPS_DEBUG1__ */

/*
-----------------------------------------
*/

  if (((corner_position[0][0] < 0 || corner_position[0][0] >= NDAT) ||
       (corner_position[0][1] < 0 || corner_position[0][1] >= NDAT) ||
       (corner_position[1][0] < 0 || corner_position[1][0] >= NDAT) ||
       (corner_position[1][1] < 0 || corner_position[1][1] >= NDAT)) &&
        Num_of_Samp == 0) {
    printf("ERROR: Turbulent_Phase_Screen: ");
    printf("range of corner_position is incorrect.\n");
    printf("ERROR:        %d : (%d, %d), (%d, %d)\n", NDAT,
               corner_position[0][0], corner_position[0][1],
               corner_position[1][0], corner_position[1][1]);
    free (amp);
    return (-1);

/*
-----------------------------------------
*/

  } else if (corner_position[0][0] == 0      &&
             corner_position[0][1] == 0      &&
             corner_position[1][0] == NDAT-1 &&
             corner_position[1][1] == NDAT-1) {
    if ((wx = (double *)calloc(M*M, sizeof(double))) == NULL) {
      printf("ERROR: Turbulent_Phase_Screen: ");
      printf("memory allocation for wx\n");
      exit (-1);
    }
    *(wx               ) = wd[0][0];
    *(wx  +      M  - 1) = wd[0][1];
    *(wx  + NDAT*M     ) = wd[1][0];
    *(wx  + M*M     - 1) = wd[1][1];
    if (CONTINUE_STR_SWT == true) {
      memcpy(wx, seed_dist, M * sizeof(double));
    }

    ndiv = 1;
    for (ilog=1; ilog<NLOG; ilog++) {
      NC = NDAT / ndiv;
      if (NLOG - ilog - 1> n_outer_log) {
        RAND_MODE = __INDEPENDENT__;
      } else {
        RAND_MODE = ___DEPENDENT___;
      }
      segment_calc(wx, amp[ilog], ndiv, NC, M, NS, 0, RAND_MODE);
      ndiv *= 2;
    }
    memcpy(seed_dist, wx + NDAT*M, M * sizeof(double));

    if (Num_of_Samp != 0) {
      for (i=0; i<Num_of_Samp; i++) {
        DS[i] = *(wx + NX[i]*M + NY[i]);
      }
    } else if (Num_of_Samp == 0) {
      for (i=0; i<NDAT; i++) {
        memcpy(DS + i*NDAT, wx + i*M, NDAT * sizeof(double));
      }
    }
    free (wx);
    COUNT_Num_of_Samp = NDAT * NDAT;

/*
----------------------------------------------------------
*/

  } else {
    if (NDAT_SEGM != 0) {
      i = 1;
      while (2 * i <= NDAT_SEGM) {
        i *= 2;
      }
      NDAT_SEGM = i;
      NDAT_AREA = NDAT / NDAT_SEGM;

    } else {

/****
     In the below, the number of (NDAT_AREA+1) x M is needed for a matrix,
     but it happens that the number is greater than 2^31. The following
     condition is set to avoid such a situation.
****/

      NDAT_AREA = NDAT;
      while ((double)(NDAT_AREA + 1) * (double)M >= pow(2.0, 30.0)) {
        NDAT_AREA /= 2;
      }

/****
     In the below, it is intended that segment data size is eqaual to the
     Y-axis data size. (N1 is the number of the Y-axis data size.
****/

      N1 = corner_position[1][1] - corner_position[0][1] + 1;
      while (NDAT_AREA > N1) {
        NDAT_AREA /= 2;
      }
      NDAT_SEGM = NDAT / NDAT_AREA;




/****
      NDAT_AREA = NDAT;
      NDAT_SEGM = NDAT / NDAT_AREA;
****/

      double a1, a2;
      a1 = log10(pow(2, 30) / (double)(NDAT_AREA + 1)) / log10(2);
      a2 = log10(pow(2, 30) / pow((double)(NDAT_SEGM + 1), 2))
                            / log10(2);
#ifdef __TPS_DEBUG2__
      printf("__TPS_DEBUG2__  %d  %d   %lf  %lf\n",
             NDAT_AREA, NDAT_SEGM, a1, a2);
#endif /* __TPS_DEBUG2__ */
      if (a1 < 0.0 && a2 > 0.0) {
        while (a1 < 0.0 || a2 < 0.0) {
          NDAT_AREA /= 2;
          a1 += 1.0;
          NDAT_SEGM *= 2;
          a2 -= 1.0;
        }
      } else if (a1 > 0.0 && a2 < 0.0) {
        while (a1 < 0.0 || a2 < 0.0) {
          NDAT_AREA *= 2;
          a1 -= 1.0;
          NDAT_SEGM /= 2;
          a2 += 1.0;
        }
      }
/**
      printf("__DEBUG__    %d  %d      %lf  %lf\n", NDAT_AREA, NDAT_SEGM, a1, a2);
**/



    }
/**
    getchar();
**/

/*
-----------------------------------------
*/

    if (Num_of_Samp != 0) {
      if ((nx = (int *)calloc(Num_of_Samp, sizeof(int))) == NULL) {
        printf("ERROR: Turbulent_Phase_Screen: calloc failed for nx.\n");
        exit (-1);
      }
      if ((ny = (int *)calloc(Num_of_Samp, sizeof(int))) == NULL) {
        printf("ERROR: Turbulent_Phase_Screen: calloc failed for ny.\n");
        exit (-1);
      }
      if ((Num_of_Area
          = area_sort(NDAT, NDAT_AREA, nx, ny, NX, NY, Num_of_Samp)) == -1) {
        printf("ERROR: Turbulent_Phase_Screen: Not enough width N: %d\n", NDAT);
        exit (-1);
      }
    } else if (Num_of_Samp == 0) {
      y_len =  corner_position[1][1] - corner_position[0][1] + 1;
      nxmin =  corner_position[0][0]      / NDAT_SEGM;
      nxmax = (corner_position[1][0] + 1) / NDAT_SEGM;
      if ((corner_position[1][0] + 1) % NDAT_SEGM != 0) {
        nxmax++;
      }
      nymin =  corner_position[0][1]      / NDAT_SEGM;
      nymax = (corner_position[1][1] + 1) / NDAT_SEGM;
      if ((corner_position[1][1] + 1) % NDAT_SEGM != 0) {
        nymax++;
      }

      Num_of_Area = (nxmax - nxmin) * (nymax - nymin);
#ifdef __TPS_DEBUG3__
      printf("__TPS_DEBUG3__  %d  [%7d  %7d]    [%7d  %7d]    (%d)    [%7d  %7d]\n",
             Num_of_Area,
             nxmin,     nxmax,
             nymin,     nymax,      Y_SHIFT,
             NDAT_AREA, NDAT_SEGM);
      fflush(stdout);
#endif /* __TPS_DEBUG3__ */

      if ((nx = (int *)calloc(Num_of_Area, sizeof(int))) == NULL) {
        printf("ERROR: Turbulent_Phase_Screen: calloc failed for nx.\n");
        exit (-1);
      }
      if ((ny = (int *)calloc(Num_of_Area, sizeof(int))) == NULL) {
        printf("ERROR: Turbulent_Phase_Screen: calloc failed for ny.\n");
        exit (-1);
      }
      I = 0;
      for (i=nxmin; i<nxmax; i++) {
        for (j=nymin; j<nymax; j++) {
          nx[I] = i;
          ny[I] = j;
          I++;
        }
      }
      COUNT_Num_of_Samp = 0;
    }
    narea_P1 = NDAT_AREA + 1;
    nsegm_P1 = NDAT_SEGM + 1;

#ifdef __TPS_DEBUG4__
    printf("__TPS_DEBUG4__   %d  %d  %d\n", narea_P1, M, narea_P1 * M);
#endif /* __TPS_DEBUG4__ */
    if ((x_lods = (double *)calloc(M * narea_P1, sizeof(double)))
                == NULL) {
      printf("fail to allocate memories for x_lods.\n");
      exit (-1);
    }
    if ((y_lods = (double *)calloc(M * narea_P1, sizeof(double)))
                == NULL) {
      printf("fail to allocate memories for y_lods.\n");
      exit (-1);
    }
    if ((init_dist = (double *)calloc(narea_P1*narea_P1, sizeof(double)))
                == NULL) {
      printf("fail in allocating memories for init_dist.\n");
      exit (-1);
    }

    if (CONTINUE_STR_SWT == true) {
      *(x_lods + NDAT_AREA*M                  ) = wd[1][0];
      *(x_lods + NDAT_AREA*M          + M - 1 ) = wd[1][1];
      *(y_lods + NDAT_AREA*M                  ) = wd[1][0];
      *(y_lods + NDAT_AREA*M          + M - 1 ) = wd[1][1];

      *(init_dist   + narea_P1*NDAT_AREA      ) = wd[1][0];
      *(init_dist   + narea_P1*narea_P1    - 1) = wd[1][1];

      memcpy(y_lods, seed_dist, M * sizeof(double));
      for (i=0; i<narea_P1; i++) {
        *(x_lods       + i*M) = *(seed_dist + i*NDAT_SEGM);
        *(init_dist      + i) = *(seed_dist + i*NDAT_SEGM);
      }
    } else if (CONTINUE_STR_SWT == false) {
      *(x_lods                                ) = wd[0][0];
      *(x_lods                        + M - 1 ) = wd[1][0];
      *(x_lods + NDAT_AREA*M                  ) = wd[0][1];
      *(x_lods + NDAT_AREA*M          + M - 1 ) = wd[1][1];

      *(y_lods                                ) = wd[0][0];
      *(y_lods                        + M - 1 ) = wd[0][1];
      *(y_lods + NDAT_AREA*M                  ) = wd[1][0];
      *(y_lods + NDAT_AREA*M          + M - 1 ) = wd[1][1];

      *(init_dist                             ) = wd[0][0];
      *(init_dist               + narea_P1 - 1) = wd[0][1];
      *(init_dist + narea_P1*NDAT_AREA        ) = wd[1][0];
      *(init_dist + narea_P1*narea_P1      - 1) = wd[1][1];
    }
#ifdef __TPS_DEBUG5__
    printf("__TPS_DEBUG5__   %lf   %lf   %lf   %lf\n",
      *(init_dist                             ),
      *(init_dist   + narea_P1             - 1),
      *(init_dist   + narea_P1*NDAT_AREA      ),
      *(init_dist   + narea_P1*narea_P1    - 1));
#endif /*__TPS_DEBUG5__*/

    nlog = 0;
    i = NDAT_AREA;
    while (i > 1) {
      i /= 2;
      nlog++;
    }

    ndiv = 1;
    for (ilog=1; ilog<=nlog; ilog++) {
      NC = NDAT_AREA / ndiv;
      if (NLOG - ilog > n_outer_log) {
        RAND_MODE = __INDEPENDENT__;
      } else {
        RAND_MODE = ___DEPENDENT___;
      }
#ifdef __TPS_DEBUG6__
      printf(
        "__TPS_DEBUG6__:__AREA__  %3d %3d   [%7d  %3d]   %7d    %1d    [%3d   %lf]\n",
             NLOG, nlog, n_outer, n_outer_log,
             NDAT/ndiv, RAND_MODE, ilog, amp[ilog]);
#endif /* __TPS_DEBUG6__ */
      segment_calc(init_dist, amp[ilog], ndiv, NC,
                   narea_P1, NS, 0, RAND_MODE);
      ndiv *= 2;
    }

    for (i=0; i<narea_P1; i++) {
      for (j=0; j<narea_P1; j++) {
        *(x_lods + j*M + NDAT_SEGM*i) = *(init_dist + i*narea_P1 + j);
        *(y_lods + i*M + NDAT_SEGM*j) = *(init_dist + i*narea_P1 + j);
      }
    }
/****
    getchar();
****/

    free (init_dist);

/*
--------------------------------------------------------------
*/

    for (iarea=0; iarea<Num_of_Area; iarea++) {

      if ((seg_dist = (double *)calloc(nsegm_P1*nsegm_P1,
                                       sizeof(double)))
                             == NULL) {
        printf("fail in allocating memories for seg_dist.\n");
        fflush(stdout);
        exit (-1);
      }

      *(seg_dist                      )
              = *(y_lods +  nx[iarea]   *M +  ny[iarea]   *NDAT_SEGM);
      *(seg_dist + nsegm_P1 - 1        )
              = *(y_lods +  nx[iarea]   *M + (ny[iarea]+1)*NDAT_SEGM);
      *(seg_dist + NDAT_SEGM*nsegm_P1       )
              = *(y_lods + (nx[iarea]+1)*M +  ny[iarea]   *NDAT_SEGM);
      *(seg_dist + nsegm_P1*nsegm_P1 - 1)
              = *(y_lods + (nx[iarea]+1)*M + (ny[iarea]+1)*NDAT_SEGM);

/*
--------
*/

      nsx = 0;
      if (CONTINUE_STR_SWT == true && nx[iarea] == 0) {
        nsx = 1;
      } else {
        for (j=0; j<iarea; j++) {
          if (nx[iarea] == nx[j] + 1 && ny[iarea] == ny[j]) {
            nsx = 1;
            break;
          }
        }
      }

      nsy = 0;
      for (j=0; j<iarea; j++) {
        if (ny[iarea] == ny[j] + 1 && nx[iarea] == nx[j]) {
          nsy = 1;
          break;
        }
      }

      if (nsx == 1) {
        for (i=1; i<nsegm_P1-1; i++) {
          *(seg_dist + i)
               = *(y_lods + M*nx[iarea] + ny[iarea]*NDAT_SEGM + i);
        }
      }
      if (nsy == 1) {
        for (i=1; i<nsegm_P1-1; i++) {
          *(seg_dist + i*nsegm_P1)
               = *(x_lods + M*ny[iarea] + nx[iarea]*NDAT_SEGM + i);
        }
      }

/*
--------
*/

      ndiv = 1;
      for (ilog=nlog+1; ilog<NLOG; ilog++) {
        NC = NDAT_SEGM / ndiv;
        if (NLOG - ilog > n_outer_log) {
          RAND_MODE = __INDEPENDENT__;
        } else {
          RAND_MODE = ___DEPENDENT___;
        }
#ifdef __TPS_DEBUG7__
        printf("__TPS_DEBUG7__:__SEGM__  %7d  %7d    %3d %3d   [%7d  %3d]   %7d    %1d    [%3d   %lf]\n",
               Num_of_Area, iarea,
               NLOG,        nlog,
               n_outer,     n_outer_log,
               NC,          RAND_MODE,
               ilog,        amp[ilog]);
#endif /* __TPS_DEBUG7__ */
        segment_calc(seg_dist, amp[ilog], ndiv, NC, nsegm_P1, nsx, nsy, RAND_MODE);
        ndiv *= 2;
      }

      if (Num_of_Samp != 0) {
        for (i=0; i<Num_of_Samp; i++) {
          if (NX[i] >= NDAT_SEGM*nx[iarea] && NX[i] <  NDAT_SEGM*(nx[iarea]+1) &&
              NY[i] >= NDAT_SEGM*ny[iarea] && NY[i] <  NDAT_SEGM*(ny[iarea]+1)) {
            DS[i] = *(seg_dist + (NX[i]%NDAT_SEGM)*nsegm_P1 + NY[i]%NDAT_SEGM);
            COUNT_Num_of_Samp++;
          }
        }
      } else if (Num_of_Samp == 0) {
        for (i=0; i<NDAT_SEGM; i++) {
          ix = NDAT_SEGM * nx[iarea] + i;
          if (ix  >= corner_position[0][0] &&
              ix  <= corner_position[1][0]) {
            for (j=0; j<NDAT_SEGM; j++) {
              iy = NDAT_SEGM * ny[iarea] + j;
              if (iy >= corner_position[0][1] &&
                  iy <= corner_position[1][1]) {
                k = (ix - corner_position[0][0]) * y_len
                         + (iy - corner_position[0][1]);
                DS[k] = *(seg_dist + i*nsegm_P1 + j);
                COUNT_Num_of_Samp++;
              }
            }
          }
        }
      }
/****
      getchar();
****/

/*
--------
*/

      memcpy(y_lods + M*(nx[iarea]+1) + NDAT_SEGM*ny[iarea] + 1,
             seg_dist +  NDAT_SEGM*nsegm_P1 + 1,
             sizeof(double) * (nsegm_P1 - 2));
      for (i=1; i<nsegm_P1-1; i++) {
        *(x_lods + M*(ny[iarea]+1) + NDAT_SEGM*nx[iarea] + i) =
                    *(seg_dist + (i+1)*nsegm_P1 - 1);
/*****
        *(y_lods + M*(nx[iarea]+1) + NDAT_SEGM*ny[iarea] + i) =
                    *(seg_dist +  NDAT_SEGM*nsegm_P1 + i);
*****/
      }
      free (seg_dist);
    }
    free (x_lods);

/*
-------------- Seed of the next screen -----------------
*/

    memcpy(seed_dist, y_lods + NDAT_AREA*M, M * sizeof(double));
    free (y_lods);

    if ((AREA_SWT = (_Bool *)calloc(NDAT_AREA, sizeof(_Bool))) == NULL) {
      printf("ERROR: Turbulent_Phase_Screen: \n");
      printf("failed in allocating memories for AREA_SWT\n");
      exit (-1);
    }

    for (i=0; i<NDAT_AREA; i++) {
      AREA_SWT[i] = true;
    }
    for (i=0; i<Num_of_Area; i++) {
      if (nx[i] == NDAT_AREA - 1) {
        AREA_SWT[ny[i]] = false;
      }
    }

    for (I=0; I<NDAT_AREA; I++) {
      if (AREA_SWT[I] == true) {
        ndiv = 1;
        for (ilog=nlog+1; ilog<NLOG; ilog++) {
          NC = NDAT_SEGM / ndiv;
          for (j=0; j<ndiv; j++) {
            *(seed_dist + NDAT_SEGM * I + j*NC + NC/2) =
              0.5 * (*(seed_dist + NDAT_SEGM * I + j*NC) +
                     *(seed_dist + NDAT_SEGM * I + (j+1)*NC))
                                      + excess_value(amp[ilog]);
          }
          ndiv *= 2;
        }
      }
    }

    free (AREA_SWT);
    free (nx);
    free (ny);

/*
--------------------------------------------------------
*/

  }

/*
------ Y Direction Shift of the Seed for the next ------
*/

  if (Y_SHIFT != 0) {

    if ((new_seed = (double *)calloc(M, sizeof(double))) == NULL) {
      printf("fail to allocate memories for new_seed.\n");
      exit (-1);
    }

    if (Y_SHIFT > 0) {
      new_seed[0]   = seed_dist[M-1];
      new_seed[M-1] = excess_value(amp[0]);
    } else if (Y_SHIFT < 0) {
      new_seed[0]   = excess_value(amp[0]);
      new_seed[M-1] = seed_dist[M-1];
    }

    ndiv = 1;
    for (ilog=0; ilog<NLOG; ilog++) {
      NC = NDAT / ndiv;
      for (j=0; j<ndiv; j++) {
        *(new_seed + j*NC + NC/2) =
          0.5 * (*(new_seed + j*NC) + *(new_seed + (j+1)*NC))
                                  + excess_value(amp[ilog]);
      }
      ndiv *= 2;
    }

    I = abs(Y_SHIFT);
    if (Y_SHIFT > 0) {
      memcpy(seed_dist,            seed_dist + I, (M - I) * sizeof(double));
      memcpy(seed_dist + (M - I),  new_seed + 1,        I * sizeof(double));
    } else if (Y_SHIFT < 0) {
      memcpy(seed_dist + I, seed_dist,            (M - I) * sizeof(double));
      memcpy(seed_dist,     new_seed + (M - I) - 1,     I * sizeof(double));
    }

    free (new_seed);
  }

/*
-----------------------------------
*/

  free (amp);

  return COUNT_Num_of_Samp;
}




void  segment_calc(double *wx, double amp, int ndiv, int NC, int M,
                   int  nsx,   int  nsy,   int  RAND_MODE)
{
  int    i, j;
  double *WX1, *WX2;

/*
-----------------
*/

  if (RAND_MODE == ___DEPENDENT___) {
    for (i=nsx; i<=ndiv; i++) {
      WX1 = wx + i*NC*M;
      for (j=0; j<ndiv; j++) {
        *(WX1 + j*NC+NC/2) = 0.5 * (*(WX1 + j*NC) + *(WX1 + (j+1)*NC))
                           + excess_value(amp);
      }
    }

    for (i=nsy; i<=ndiv; i++) {
      WX1 = wx + i*NC;
      for (j=0; j<ndiv; j++) {
        *(WX1 + (j*NC +NC/2)*M) = 0.5 * (*(WX1 + (j+1)*NC*M) + *(WX1 + j*NC*M))
                           + excess_value(amp);
      }
    }

    for (i=0; i<ndiv; i++) {
      WX1 = wx + i*NC*M;
      for (j=0; j<ndiv; j++) {
        WX2 = WX1 + j*NC;
        *(WX2      + NC/2 *M        +NC/2) = 0.25 *
            (*(WX2                  +NC/2)
           + *(WX2 + NC   *M        +NC/2)
           + *(WX2 + NC/2 *M             )
           + *(WX2 + NC/2 *M        +NC  ))
           + excess_value(amp);
      }
    }

/*
-----------------
*/

  } else if (RAND_MODE == __INDEPENDENT__) {
    for (i=nsx; i<=ndiv; i++) {
      WX1 = wx + i*NC*M;
      for (j=0; j<ndiv; j++) {
        *(WX1 + j*NC+NC/2)
             = excess_value(amp);
      }
    }

    for (i=nsy; i<=ndiv; i++) {
      WX1 = wx + i*NC;
      for (j=0; j<ndiv; j++) {
        *(WX1 + (j*NC +NC/2)*M)
             = excess_value(amp);
      }
    }

    for (i=0; i<ndiv; i++) {
      WX1 = wx + i*NC*M;
      for (j=0; j<ndiv; j++) {
        WX2 = WX1 + j*NC;
        *(WX2      + NC/2 *M        +NC/2)
             = excess_value(amp);
      }
    }
  }
}


int    area_sort(int N, int n_area,
                 int *nx, int *ny, int *NX, int *NY, int Num_of_Samp)
{
  int i, j, ntmp, n1, n2, nw;

  for (i=0; i<Num_of_Samp; i++) {
    nx[i] = NX[i] / (N / n_area);
    ny[i] = NY[i] / (N / n_area);
    if (nx[i] < 0 || ny[i] < 0) {
      printf("ERROR: area_sort: (nx, ny)=(%d, %d): ", nx[i], ny[i]);
      printf("screen width is not enough.\n");
      return -1;
    }
  }

  for (i=0; i<Num_of_Samp; i++) {
    n1 = nx[i] * n_area + ny[i];
    ntmp = i;
    for (j=i+1; j<Num_of_Samp; j++) {
      n2 = nx[j] * n_area + ny[j];
      if (n1 > n2) {
        n1 = n2;
        ntmp = j;
      }
    }
    nw = nx[i];
    nx[i] = nx[ntmp];
    nx[ntmp] = nw;
    nw = ny[i];
    ny[i] = ny[ntmp];
    ny[ntmp] = nw;
  }

  j = 1;
  for (i=1; i<Num_of_Samp; i++) {
    if (! (nx[i] == nx[j-1] && ny[i] == ny[j-1])) {
      nx[j] = nx[i];
      ny[j] = ny[i];
      j++;
    }
  }

  return (j);
}


double excess_value(double amp)
{
  return amp * random_val1();
/**
  return amp * gauss_dev();
  return amp * random_val1();
**/
}
