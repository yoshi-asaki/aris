#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <phase_screen.h>


/**
#define __DEBUG1__
**/


/****
    f(x)  ~  x^2 + (N/x)^2   --> minimum
                 |
                 V
    df/dx = 2 (x^3 - N^2/2) / x^2

     x = (N^2 / 2)^{1/3}  for minimum f(x)
****/


void   segment_calc(double *, double , int , int , int , int , int );
int    area_sort(int , int , int *, int *, int *, int *, int );


int    turbulent_phase_screen
                   (int N, int NSEG, double *two_D_dist,
                    double *seed_dist, double *seed_end,
                    struct phase_screen_parameter *atm,
                    int    CONTINUE_STR_SWT, int CONTINUE_END_SWT,
                    int    NOD, int *NX, int  *NY, double *DS,
                    int    corner_position[][2], int Y_SHIFT)
{
  static int    i, j, k, I, nlog, NLOG, ndiv;
  static int    NC, M, NS, NOA;
  static int    n_inner,  n_outer;
  static int    nh_inner, nh_outer;
  static int    nseg_P1,  narea_P1;
  static int    *nx, *ny;
  static int    nsx, nsy, N_AREA, *AREA_SWT;

  static double *wx;
  static double *amp;
  static double amp_factor;
  static double D, d_expon;
  static double wd[2][2];
  static double *init_dist;

  static double *seg_dist;
  static double *x_lods, *y_lods;
  static double *new_seed;

  static int    nxmin, nxmax, nymin, nymax, COUNT_NOD;
  static int    ix, iy, y_len;

/*
-------------------------------------
*/

  if (NOD != 0) {
    for (i=0; i<NOD; i++) {
      if (NX[i] < 0 || NX[i] >= N || NY[i] < 0 || NY[i] >= N) {
        printf("ERROR: turbulent_phase_screen: out of range NX or NY.\n");
        exit (-1);
      }
    }
  }

/*
-------------------------------------
*/

  d_expon = atm->i_expon - atm->o_expon;
  M = N + 1;

/*
-------------------------------------
*/

  i = (int)(log10(atm->i_scale[0]/atm->pixel) / log10(2.0));
  j = (int)(log10(atm->i_scale[0]/atm->pixel) / log10(2.0)) + 1;
  if (fabs(pow(2.0, (double)i) * atm->pixel - atm->i_scale[0]) < 
      fabs(pow(2.0, (double)j) * atm->pixel - atm->i_scale[0])) {
    n_inner = i;
  } else {
    n_inner = j;
  }

  i = (int)(log10(atm->o_scale[0]/atm->pixel) / log10(2.0));
  j = (int)(log10(atm->o_scale[0]/atm->pixel) / log10(2.0)) + 1;
  if (fabs(pow(2.0, (double)i) * atm->pixel - atm->o_scale[0]) < 
      fabs(pow(2.0, (double)j) * atm->pixel - atm->o_scale[0])) {
    n_outer = i;
  } else {
    n_outer = j;
  }

  nh_inner = n_inner - 1;
  nh_outer = n_outer - 1;

  atm->i_scale[1] = atm->pixel * pow(2.0, (double)n_inner);
  atm->o_scale[1] = atm->pixel * pow(2.0, (double)n_outer);

  nlog = 0;
  i = 1;
  while (i <= N) {
    i *= 2;
    nlog++;
  }
  if ((amp = (double *)calloc(nlog, sizeof(double))) == NULL) {
    printf("fail in allocating memories for amp.\n");
    exit (-1);
  }

/*
-------------------------------------------
*/

  for (i=0; i<nlog; i++) {
    I = nlog -1 - i;
    if (I > n_outer) {
      amp[i] = 9.0 * pow(atm->i_scale[1], d_expon)
                   * pow(atm->o_scale[1], atm->o_expon);
    } else if (I <= n_outer  && I > nh_outer) {
      D = atm->pixel * pow(2.0, (double)I);
      amp[i] = 3.0 * pow(atm->i_scale[1], d_expon)
              * (4.0*pow(D, atm->o_expon)
                      - pow(atm->o_scale[1], atm->o_expon));
    } else if (I <= nh_outer && I > n_inner) {
      D = atm->pixel * pow(2.0, (double)I);
      amp[i] = 3.0 * pow(atm->i_scale[1], d_expon)
                   * pow(D, atm->o_expon) * (4.0 - pow(2.0, atm->o_expon));
    } else if (I <= n_inner  && I > nh_inner) {
      D = atm->pixel * pow(2.0, (double)I);
      amp[i] = 3.0
        * (4.0 * pow(D, atm->i_expon)
               - pow(atm->i_scale[1], d_expon) * pow(2.0*D, atm->o_expon));
    } else if (I <= nh_inner) {
      D = atm->pixel * pow(2.0, (double)I);
      amp[i] = 3.0 * pow(D, atm->i_expon) * (4.0 - pow(2.0, atm->i_expon));
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

  for (i=0; i<nlog; i++) {
    amp[i] = amp_factor / 2.0 * sqrt(amp[i]);
  }

/*
---------------------------------------------
*/

  if (CONTINUE_STR_SWT == ON) {
    NS = 1;
  } else if (CONTINUE_STR_SWT == OFF) {
    wd[0][0] = amp[0] * random_val1();
    wd[0][1] = amp[0] * random_val1();
    NS = 0;
  }

  if (CONTINUE_END_SWT == OFF) {
    wd[1][0] = amp[0] * random_val1();
    wd[1][1] = amp[0] * random_val1();
  } else if (CONTINUE_END_SWT == ON) {
    wd[1][0] = amp[0] * seed_end[0];
    wd[1][1] = amp[0] * seed_end[1];
  }

/*
-----------------------------------------
*/

  if (((corner_position[0][0] < 0 || corner_position[0][0] >= N) ||
       (corner_position[0][1] < 0 || corner_position[0][1] >= N) ||
       (corner_position[1][0] < 0 || corner_position[1][0] >= N) ||
       (corner_position[1][1] < 0 || corner_position[1][1] >= N)) &&
        NOD == 0) {
    printf("ERROR: turbulent_phase_screen: ");
    printf("Range of corner_position is wrong.\n");
    free (amp);
    return (-1);

/*
-----------------------------------------
*/

  } else if (corner_position[0][0] == 0   && corner_position[0][1] == 0 &&
             corner_position[1][0] == N-1 && corner_position[1][1] == N-1) {
    if ((wx = (double *)calloc(M*M, sizeof(double))) == NULL) {
      printf("ERROR: turbulent_screen: memory allocation for wx\n");
      exit (-1);
    }
    ndiv = 1;
    *(wx            ) = wd[0][0];
    *(wx  +   M  - 1) = wd[0][1];
    *(wx  + N*M     ) = wd[1][0];
    *(wx  + M*M  - 1) = wd[1][1];
    if (CONTINUE_STR_SWT == ON) {
      memcpy(wx, seed_dist, M * sizeof(double));
    }

    for (i=1; i<nlog; i++) {
      NC = N / ndiv;
      segment_calc(wx, amp[i], ndiv, NC, M, NS, 0);
      ndiv *= 2;
    }

    for (i=0; i<N; i++) {
      memcpy(two_D_dist + i*N, wx + i*M, N * sizeof(double));
    }
    memcpy(seed_dist, wx + N*M, M * sizeof(double));

    if (NOD != 0) {
      for (i=0; i<NOD; i++) {
        DS[i] = *(wx + NX[i]*M + NY[i]);
      }
    } else if (NOD == 0) {
      for (i=0; i<N; i++) {
        memcpy(DS + i*N, wx + i*M, N * sizeof(double));
      }
    }
    free (wx);

/*
----------------------------------------------------------
*/

  } else {

    if (NSEG != 0) {
      i = 1;
      while (2 * i <= NSEG) {
        i *= 2;
      }
      NSEG = i;
      N_AREA = N / NSEG;
    } else {
      N_AREA = (int)rint(pow(0.5 * (double)(N*N), 1.0/3.0));
      i = (int)(log10((double)N_AREA) / log10(2.0));
      j = (int)(log10((double)N_AREA) / log10(2.0)) + 1;
      if (fabs(pow(2.0, (double)i)  - N_AREA) < 
          fabs(pow(2.0, (double)j)  - N_AREA)) {
        N_AREA = (int)pow(2.0, (double)i);
      } else {
        N_AREA = (int)pow(2.0, (double)j);
      }
      NSEG     = N / N_AREA;
    }
    if ((AREA_SWT = (int *)calloc(N_AREA, sizeof(int))) == NULL) {
      printf("fail in allocating memories for AREA_SWT\n");
      exit (-1);
    }

/*
-----------------------------------------
*/

    if (NOD != 0) {
      if ((nx = (int *)calloc(NOD, sizeof(int))) == NULL) {
        printf("ERROR: turbulent_screen: calloc failed.\n");
        exit (-1);
      }
      if ((ny = (int *)calloc(NOD, sizeof(int))) == NULL) {
        printf("ERROR: turbulent_screen: calloc failed.\n");
        exit (-1);
      }
      if ((NOA = area_sort(N, N_AREA, nx, ny, NX, NY, NOD)) == -1) {
        printf("ERROR: turbulent_screen: Not enough width N: %d\n", N);
        exit (-1);
      }
    } else if (NOD == 0) {
      y_len = corner_position[1][1] - corner_position[0][1] + 1;
      nxmin =  corner_position[0][0]      / NSEG;
      nxmax = (corner_position[1][0] + 1) / NSEG;
      if ((corner_position[1][0] + 1) % NSEG != 0) {
        nxmax++;
      }
      nymin =  corner_position[0][1]      / NSEG;
      nymax = (corner_position[1][1] + 1) / NSEG;
      if ((corner_position[1][1] + 1) % NSEG != 0) {
        nymax++;
      }

      NOA = (nxmax - nxmin) * (nymax - nymin);
#ifdef __DEBUG1__
      printf("__DEBUG1__  %d  %d  %d  %d  %d  %d\n",
            NOA, nxmin, nxmax, nymin, nymax, Y_SHIFT);
      fflush(stdout);
#endif /* __DEBUG1__ */

      if ((nx = (int *)calloc(NOA, sizeof(int))) == NULL) {
        printf("ERROR: turbulent_screen: calloc failed.\n");
        exit (-1);
      }
      if ((ny = (int *)calloc(NOA, sizeof(int))) == NULL) {
        printf("ERROR: turbulent_screen: calloc failed.\n");
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
      COUNT_NOD = 0;
    }
    narea_P1 = N_AREA + 1;
    nseg_P1  = NSEG   + 1;

    if ((x_lods = (double *)calloc(M * narea_P1, sizeof(double))) == NULL) {
      printf("fail in allocating memories for x_lods.\n");
      exit (-1);
    }
    if ((y_lods = (double *)calloc(M * narea_P1, sizeof(double))) == NULL) {
      printf("fail in allocating memories for y_lods.\n");
      exit (-1);
    }
    if ((init_dist = (double *)calloc(narea_P1*narea_P1, sizeof(double)))
                == NULL) {
      printf("fail in allocating memories for init_dist.\n");
      exit (-1);
    }

    *(x_lods                    ) = wd[0][0];
    *(x_lods            + M - 1 ) = wd[1][0];
    *(x_lods + N_AREA*M         ) = wd[0][1];
    *(x_lods + N_AREA*M + M - 1 ) = wd[1][1];

    *(y_lods                    ) = wd[0][0];
    *(y_lods            + M - 1 ) = wd[0][1];
    *(y_lods + N_AREA*M         ) = wd[1][0];
    *(y_lods + N_AREA*M + M - 1 ) = wd[1][1];

    if (CONTINUE_STR_SWT == ON) {
      memcpy(y_lods, seed_dist, M * sizeof(double));
      for (i=0; i<=N_AREA; i++) {
        *(x_lods + i*M       ) = *(seed_dist + i*NSEG);
      }
      for (i=0; i<narea_P1; i++) {
        *(init_dist + i) = *(seed_dist + i*NSEG);
      }
    } else {
      *(init_dist                          ) = wd[0][0];
      *(init_dist            + narea_P1 - 1) = wd[0][1];
      *(init_dist   + narea_P1*N_AREA      ) = wd[1][0];
      *(init_dist   + narea_P1*narea_P1 - 1) = wd[1][1];
    }

    NLOG = 0;
    i = N_AREA;
    while (i > 1) {
      i /= 2;
      NLOG++;
    }

    ndiv = 1;
    for (i=1; i<=NLOG; i++) {
      NC = N_AREA / ndiv;
      segment_calc(init_dist, amp[i], ndiv, NC, narea_P1, NS, 0);
      ndiv *= 2;
    }

    for (i=0; i<narea_P1; i++) {
      for (j=0; j<narea_P1; j++) {
        *(x_lods + j*M + NSEG*i) = *(init_dist + i*narea_P1 + j);
        *(y_lods + i*M + NSEG*j) = *(init_dist + i*narea_P1 + j);
      }
    }
    free (init_dist);

/*
--------------------------------------------------------------
*/

    if ((seg_dist = (double *)calloc(nseg_P1*nseg_P1, sizeof(double)))
                           == NULL) {
      printf("fail in allocating memories for seg_dist.\n");
      fflush(stdout);
    }

    for (I=0; I<NOA; I++) {

      for (i=0; i<nseg_P1; i++) {
        for (j=0; j<nseg_P1; j++) {
          *(seg_dist + i*nseg_P1 + j) = 0.0;
        }
      }

      *(seg_dist                      )
              = *(y_lods +  nx[I]   *M +  ny[I]   *NSEG);
      *(seg_dist + nseg_P1 - 1        )
              = *(y_lods +  nx[I]   *M + (ny[I]+1)*NSEG);
      *(seg_dist + NSEG*nseg_P1       )
              = *(y_lods + (nx[I]+1)*M +  ny[I]   *NSEG);
      *(seg_dist + nseg_P1*nseg_P1 - 1)
              = *(y_lods + (nx[I]+1)*M + (ny[I]+1)*NSEG);

/*
--------
*/

      nsx = 0;
      if (CONTINUE_STR_SWT == ON && nx[I] == 0) {
        nsx = 1;
      } else {
        for (j=0; j<I; j++) {
          if (nx[I] == nx[j] + 1 && ny[I] == ny[j]) {
            nsx = 1;
            break;
          }
        }
      }

      nsy = 0;
      for (j=0; j<I; j++) {
        if (ny[I] == ny[j] + 1 && nx[I] == nx[j]) {
          nsy = 1;
          break;
        }
      }

      if (nsx == 1) {
        for (i=1; i<nseg_P1-1; i++) {
          *(seg_dist + i)
               = *(y_lods + M*nx[I] + ny[I]*NSEG + i);
        }
      }
      if (nsy == 1) {
        for (i=1; i<nseg_P1-1; i++) {
          *(seg_dist + i*nseg_P1)
               = *(x_lods + M*ny[I] + nx[I]*NSEG + i);
        }
      }

/*
--------
*/

      ndiv = 1;
      for (i=NLOG+1; i<nlog; i++) {
        NC = NSEG / ndiv;
        segment_calc(seg_dist, amp[i], ndiv, NC, nseg_P1, nsx, nsy);
        ndiv *= 2;
      }

      if (NOD != 0) {
        for (i=0; i<NOD; i++) {
          if (NX[i] >= NSEG*nx[I] && NX[i] <  NSEG*(nx[I]+1) &&
              NY[i] >= NSEG*ny[I] && NY[i] <  NSEG*(ny[I]+1)) {
            DS[i] = *(seg_dist + (NX[i]%NSEG)*nseg_P1 + NY[i]%NSEG);
            COUNT_NOD++;
          }
        }
      } else if (NOD == 0) {
        for (i=0; i<NSEG; i++) {
          ix = NSEG * nx[I] + i;
          if (ix  >= corner_position[0][0] &&
              ix  <= corner_position[1][0]) {
            for (j=0; j<NSEG; j++) {
              iy = NSEG * ny[I] + j;
              if (iy >= corner_position[0][1] &&
                  iy <= corner_position[1][1]) {
                k = (ix - corner_position[0][0]) * y_len
                         + (iy - corner_position[0][1]);
                DS[k] = *(seg_dist + i*nseg_P1 + j);
                COUNT_NOD++;
              }
            }
          }
        }
      }

      for (i=1; i<nseg_P1-1; i++) {
        *(x_lods + M*(ny[I]+1) + NSEG*nx[I] + i) =
                    *(seg_dist + (i+1)*nseg_P1 - 1);
        *(y_lods + M*(nx[I]+1) + NSEG*ny[I] + i) =
                    *(seg_dist +  NSEG*nseg_P1 + i);
      }
    }
    free (seg_dist);
    free (x_lods);

/*
-------------- Seed of the next screen -----------------
*/

    memcpy(seed_dist, y_lods + N_AREA*M, M * sizeof(double));
    free (y_lods);
    for (i=0; i<N_AREA; i++) {
      AREA_SWT[i] = ON;
    }
    for (i=0; i<NOA; i++) {
      if (nx[i] == N_AREA - 1) {
        AREA_SWT[ny[i]] = OFF;
      }
    }

    for (I=0; I<N_AREA; I++) {
      if (AREA_SWT[I] == ON) {
        ndiv = 1;
        for (i=NLOG+1; i<nlog; i++) {
          NC = NSEG / ndiv;
          for (j=0; j<ndiv; j++) {
            *(seed_dist + NSEG * I + j*NC + NC/2) =
              0.5 * (*(seed_dist + NSEG * I + j*NC) +
                     *(seed_dist + NSEG * I + (j+1)*NC))
                                      + amp[i]*random_val1();
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
      printf("fail in allocating memories for new_seed.\n");
      exit (-1);
    }

    if (Y_SHIFT > 0) {
      new_seed[0]   = seed_dist[M-1];
      new_seed[M-1] = amp[0] * random_val1();
    } else if (Y_SHIFT < 0) {
      new_seed[0]   = amp[0] * random_val1();
      new_seed[M-1] = seed_dist[M-1];
    }

    ndiv = 1;
    for (i=0; i<nlog; i++) {
      NC = N / ndiv;
      for (j=0; j<ndiv; j++) {
        *(new_seed + j*NC + NC/2) =
          0.5 * (*(new_seed + j*NC) + *(new_seed + (j+1)*NC))
                                  + amp[i]*random_val1();
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

  return COUNT_NOD;
}




void  segment_calc(double *wx, double amp, int ndiv, int NC, int M,
                   int  nsx, int  nsy)
{
  static int    i, j;
  static double *WX1, *WX2;

  for (i=nsx; i<=ndiv; i++) {
    WX1 = wx + i*NC*M;
    for (j=0; j<ndiv; j++) {
      *(WX1 + j*NC+NC/2) = 0.5 * (*(WX1 + j*NC) + *(WX1 + (j+1)*NC))
                         + amp*random_val1();
    }
  }

  for (i=nsy; i<=ndiv; i++) {
    WX1 = wx + i*NC;
    for (j=0; j<ndiv; j++) {
      *(WX1 + (j*NC +NC/2)*M) = 0.5 * (*(WX1 + (j+1)*NC*M) + *(WX1 + j*NC*M))
                         + amp*random_val1();
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
         + amp*random_val1();
    }
  }
}


int    area_sort(int N, int n_area,
                 int *nx, int *ny, int *NX, int *NY,
                 int NOD)
{
  static int i, j, ntmp, n1, n2, nw;

  for (i=0; i<NOD; i++) {
    nx[i] = NX[i] / (N / n_area);
    ny[i] = NY[i] / (N / n_area);
    if (nx[i] < 0 || ny[i] < 0) {
      printf("ERROR: area_sort: (nx, ny)=(%d, %d): ", nx[i], ny[i]);
      printf("screen width is not enough.\n");
      return -1;
    }
  }

  for (i=0; i<NOD; i++) {
    n1 = nx[i] * n_area + ny[i];
    ntmp = i;
    for (j=i+1; j<NOD; j++) {
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
  for (i=1; i<NOD; i++) {
    if (! (nx[i] == nx[j-1] && ny[i] == ny[j-1])) {
      nx[j] = nx[i];
      ny[j] = ny[i];
      j++;
    }
  }

  return (j);
}
