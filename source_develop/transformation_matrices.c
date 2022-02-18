#include <stdio.h>
#include <math.h>
#include <aris.h>

/****
#define __DEBUG__
****/

int  transformation_matrices(
            int    *TimUTC,
            double Qt[][3],     double W[][3],
            double d_Delta_psi, double d_Delta_epsiron,
            double WX,          double WY)
{
  int    i, j, k;
  double tueta, theta, z;
  double epsiron, Delta_epsiron, Delta_psi;
  double R1[3][3], R2[3][3], R3[3][3], R4[3][3], R5[3][3];

/****************
    TRS -> CRS   : [Q] = [P][N]
      [P] = Rz(tueta)    Ry(-theta) Rz(z)
      [N] = Rx(-epsiron) Rz(D psi)  Rx(epsiron + D epsiron)

      [Q]^{-1} = [Q]^{T}
               = Rx(-epsiron - D epsiron)
                 Rz(-D psi)
                 Rx(epsiron)
                 Rz(z)
                 Ry(theta)
                 Rz(-tueta)
****************/

/*
--------
*/

  nutation_calc(TimUTC, 0.0, &epsiron, &Delta_psi, &Delta_epsiron);
#ifdef __DEBUG__
  printf("TRANSFORMATION_MATRICES: __DEBUG__   %lf  %lf  %lf\n",
          epsiron, Delta_psi, Delta_epsiron); fflush(stdout);
#endif /* __DEBUG__ */

  drotation_matrix_set(R1[0], -epsiron
                           - (Delta_epsiron + d_Delta_epsiron), "x");
  drotation_matrix_set(R2[0], -(Delta_psi + d_Delta_psi),       "z");
  drotation_matrix_set(R3[0],  epsiron,                         "x");

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      R4[i][j] = 0.0;
      for (k=0; k<3; k++) {
        R4[i][j] += R1[i][k] * R2[k][j];
      }
    }
  }
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      R5[i][j] = 0.0;
      for (k=0; k<3; k++) {
        R5[i][j] += R4[i][k] * R3[k][j];
      }
    }
  }

  precession_calc(TimUTC, 0.0, &tueta, &theta, &z);
  drotation_matrix_set(R1[0],       z,  "z");
  drotation_matrix_set(R2[0],   theta,  "y");
  drotation_matrix_set(R3[0],  -tueta,  "z");

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      R4[i][j] = 0.0;
      for (k=0; k<3; k++) {
        R4[i][j] += R5[i][k] * R1[k][j];
      }
    }
  }
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      R5[i][j] = 0.0;
      for (k=0; k<3; k++) {
        R5[i][j] += R4[i][k] * R2[k][j];
      }
    }
  }
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      Qt[i][j] = 0.0;
      for (k=0; k<3; k++) {
        Qt[i][j] += R5[i][k] * R3[k][j];
      }
    }
  }

/*
--------
*/

  drotation_matrix_set(R1[0], WY, "x");
  drotation_matrix_set(R2[0], WX, "y");

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      W[i][j] = 0.0;
      for (k=0; k<3; k++) {
        W[i][j] += R1[i][k] * R2[k][j];
      }
    }
  }

/*
--------
*/

  return 1;
}
