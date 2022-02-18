#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"

/****
#define __DEBUG__
****/

void Q_param(double *q, double *e,  double phi)
{
/* Q parameter */
  q[0] = e[0] * sin(0.5 * phi);
  q[1] = e[1] * sin(0.5 * phi);
  q[2] = e[2] * sin(0.5 * phi);
  q[3] = cos(0.5 * phi);

  return;
}


void Euler_parameter(double *s1, double *s2,
                     double *e,  double *phi)
{
  double t;

/* Euler vector */
  e[0] = s1[1]*s2[2] - s1[2]*s2[1];
  e[1] = s1[2]*s2[0] - s1[0]*s2[2];
  e[2] = s1[0]*s2[1] - s1[1]*s2[0];
  t = vlen3(e);
  e[0] /= t;
  e[1] /= t;
  e[2] /= t;

/* Euler angle */
  *phi = sepang(s1, s2);

  return;
}


void Q_parameter(double *s1, double *s2,
                 double *e,  double *phi,
                 double *q)
{

  Euler_parameter(s1, s2, e, phi);
  Q_param(q, e, *phi);

  return;
}


void coordinate_rotation_parameter
       (double *e, double phi, double *q, double *A)
{
  Q_param(q, e, phi);
  q2matrix(A, q);

  return;
}



void vector_rotation(double *s, double *e, double phi)
{
  double q[4], A[9], t[3];

  Q_param(q, e, phi);
  q2xyz(q, A, s, t);
  memcpy(s, t, 3*sizeof(double));

  return;
}


void coordinate_rotation(double *e, double phi, double *q, double *A,
                         double *s, double *t)
{
  Q_param(q, e, phi);

#ifdef __DEBUG__
  printf("Q-check: %lf\n", q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
#endif /* __DEBUG__ */

  q2xyz(q, A, s, t);

  return;
}


void matrix_rotation(double *A, double *s, double *t)
{
  t[0] = *(A+3*0+0) *s[0] + *(A+3*1+0) *s[1] + *(A+3*2+0) *s[2];
  t[1] = *(A+3*0+1) *s[0] + *(A+3*1+1) *s[1] + *(A+3*2+1) *s[2];
  t[2] = *(A+3*0+2) *s[0] + *(A+3*1+2) *s[1] + *(A+3*2+2) *s[2];

  return;
}


void q2xyz(double *q, double *A, double *s, double *t)
{
  q2matrix(A, q);
  matrix_rotation(A, s, t);

  return;
}


void  xyz2radec(double *s, int *rahh, int *ramm, double *rass,
                           int *dcdd, int *dcmm, double *dcss)
{
  double ra, dc;
  double dpi = 3.141592653589793238462643;

  xyz2radec_rad(s, &ra, &dc);
  ra *= ( 12.0 / dpi);
  dc *= (180.0 / dpi);
  if (ra < 0.0) {
    ra += 24.0;
  }

  *rahh = (int)ra;
  *ramm = (int)((ra - (double)*rahh) * 60.0);
  *rass = (ra - (double)*rahh - (double)*ramm/60.0) * 3600.0;

  *dcdd = (int)dc;
  *dcmm = (int)((dc - (double)*dcdd) * 60.0);
  *dcss = (dc - (double)*dcdd - (double)*dcmm/60.0) * 3600.0;

  return;
}


void  xyz2radec_rad(double *s, double *ra, double *dc)
{
  *ra = atan2(s[1], s[0]);
  *dc = atan2(s[2], vlen2(s));

  return;
}



void  xyz2q(double *x, double *y, double *z,
            double *A, double *q)
{

  *(A+3*0+0) = x[0]; *(A+3*0+1) = x[1]; *(A+3*0+2) = x[2];
  *(A+3*1+0) = y[0]; *(A+3*1+1) = y[1]; *(A+3*1+2) = y[2];
  *(A+3*2+0) = z[0]; *(A+3*2+1) = z[1]; *(A+3*2+2) = z[2];

  q[3] = 0.5 * sqrt(1.0 + *(A+3*0+0) + *(A+3*1+1) + *(A+3*2+2));
  q[0] = 1.0 / 4.0 / q[3] * (*(A+3*1+2) - *(A+3*2+1));
  q[1] = 1.0 / 4.0 / q[3] * (*(A+3*2+0) - *(A+3*0+2));
  q[2] = 1.0 / 4.0 / q[3] * (*(A+3*0+1) - *(A+3*1+0));

  return;
}


void  xyz2Euler_param(double *x, double *y, double *z,
                      double *e, double *phi)
{
  double q[4], A[9];

  xyz2q(x, y, z, A, q);
  *phi = acos(q[3]);
  e[0] = q[0] / sin(*phi);
  e[1] = q[1] / sin(*phi);
  e[2] = q[2] / sin(*phi);
  *phi *= 2.0;

  return;
}


void   radec_rad2xyz(double *s, double *ra, double *dc)
{
  s[0] = cos(*dc) * cos(*ra);
  s[1] = cos(*dc) * sin(*ra);
  s[2] = sin(*dc);

  return;
}


void   q2matrix(double *A, double *q)
{
  *(A+3*0+0) = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  *(A+3*0+1) = 2.0 * (q[0]*q[1] + q[2]*q[3]);
  *(A+3*0+2) = 2.0 * (q[0]*q[2] - q[1]*q[3]);
  *(A+3*1+0) = 2.0 * (q[0]*q[1] - q[2]*q[3]);
  *(A+3*1+1) = -q[0]*q[0] + q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  *(A+3*1+2) = 2.0 * (q[1]*q[2] + q[0]*q[3]);
  *(A+3*2+0) = 2.0 * (q[0]*q[2] + q[1]*q[3]);
  *(A+3*2+1) = 2.0 * (q[1]*q[2] - q[0]*q[3]);
  *(A+3*2+2) = -q[0]*q[0] - q[1]*q[1] + q[2]*q[2] + q[3]*q[3];

  return;
}
