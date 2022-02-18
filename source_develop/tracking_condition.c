#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

#define SRT_OFF_BORESIGHT_LIM  110
#define SRT_DAZDT_LIM          2.4
#define SRT_DELDT_LIM          2.4


/****
#define __DEBUG__
****/

int  tracking_condition(int  *TimUTC, double UT1_UTC,
                   struct antenna_parameter trk_pos,
                   double OBS_T,     double OBS_p,   _Bool  ATM_SWT,
                   double *srt_pos,  double *V_xyz,
                   double *AZ,       double *EL,
                   double *dELdt,    double *dAZdt,
                   struct srt_data_link *srt_link,
                   struct source_parameter src,
                   struct source_parameter sun,
                   double *srt_AZ,      double *srt_EL,
                   double *srt_dAZdt,   double *srt_dELdt,
                   double *ptrk,        double *pdis,
                   int    BODY_X_SUN)
{
  int    I, J, N;
  _Bool  srt_link_status;
  double r, R, omega, VR, dist;
  double theta, alpha, delta;
  double s[3], v[3], vtrk[3];
  double obj_rot_axis[3];
  double a[3], srt_x[3], srt_y[3], srt_z[3];

  theta = LST(TimUTC, 0.0, UT1_UTC, 0.0);
  ptrk[0] = trk_pos.XYZ[0];
  ptrk[1] = trk_pos.XYZ[1];
  ptrk[2] = trk_pos.XYZ[2];
  drotate(ptrk, theta, "z");

  pdis[0] = srt_pos[0] - ptrk[0];
  pdis[1] = srt_pos[1] - ptrk[1];
  pdis[2] = srt_pos[2] - ptrk[2];
  dist = vlen3(pdis);
  s[0] = pdis[0] / dist;
  s[1] = pdis[1] / dist;
  s[2] = pdis[2] / dist;
  xyz2radec_rad(s, &alpha, &delta);

  R = vlen2(ptrk);
  vtrk[0] = -R * OMEGA * sin(theta);
  vtrk[1] =  R * OMEGA * cos(theta);
  vtrk[2] = 0.0;
  v[0] = V_xyz[0] - vtrk[0];
  v[1] = V_xyz[1] - vtrk[1];
  v[2] = V_xyz[2] - vtrk[2];

  VR = v[0]*s[0] + v[1]*s[1] + v[2]*s[2];
  v[0] -= VR * s[0];
  v[1] -= VR * s[1];
  v[2] -= VR * s[2];

  obj_rot_axis[0] = s[1] * v[2] - s[2] * v[1];
  obj_rot_axis[1] = s[2] * v[0] - s[0] * v[2];
  obj_rot_axis[2] = s[0] * v[1] - s[1] * v[0];
  r = vlen3(obj_rot_axis);
  obj_rot_axis[0] /= r;
  obj_rot_axis[1] /= r;
  obj_rot_axis[2] /= r;
  omega = vlen3(v) / dist;

  azel_position(TimUTC, UT1_UTC,
                trk_pos.LLH[0], trk_pos.LLH[1], trk_pos.LLH[2],
                OBS_T, OBS_p, ATM_SWT, alpha, delta,
                AZ, EL, dAZdt, dELdt, omega, obj_rot_axis);

/*
--------
*/

  srt_link_status = true;
  if (BODY_X_SUN != 0) {
    srt_z[0] = src.s[0];
    srt_z[1] = src.s[1];
    srt_z[2] = src.s[2];

    if (BODY_X_SUN > 0) {
      srt_y[0] = srt_z[1] * sun.s[2] - srt_z[2] * sun.s[1];
      srt_y[1] = srt_z[2] * sun.s[0] - srt_z[0] * sun.s[2];
      srt_y[2] = srt_z[0] * sun.s[1] - srt_z[1] * sun.s[0];
    } else if (BODY_X_SUN < 0) {
      srt_y[0] = sun.s[1] * srt_z[2] - sun.s[2] * srt_z[1];
      srt_y[1] = sun.s[2] * srt_z[0] - sun.s[0] * srt_z[2];
      srt_y[2] = sun.s[0] * srt_z[1] - sun.s[1] * srt_z[0];
    }

    r = vlen3(srt_y);
    srt_y[0] /= r;
    srt_y[1] /= r;
    srt_y[2] /= r;
    srt_x[0] = srt_y[1] * srt_z[2] - srt_y[2] * srt_z[1];
    srt_x[1] = srt_y[2] * srt_z[0] - srt_y[0] * srt_z[2];
    srt_x[2] = srt_y[0] * srt_z[1] - srt_y[1] * srt_z[0];

    a[0] = -(s[0] * srt_x[0] + s[1] * srt_x[1] + s[2] * srt_x[2]);
    a[1] = -(s[0] * srt_y[0] + s[1] * srt_y[1] + s[2] * srt_y[2]);
    a[2] = -(s[0] * srt_z[0] + s[1] * srt_z[1] + s[2] * srt_z[2]);
    r = vlen3(a);
    a[0] /= r;
    a[1] /= r;
    a[2] /= r;
    vector_rotation(a, srt_link->FOV_rot_axis, (*srt_link).FOV_rot_angle);

    *srt_AZ = atan2(a[1], a[0]);
    *srt_EL = atan2(a[2], vlen2(a));

    J = (int)(90.0 - *srt_EL / dpi * 180.0);
    if (J > SRT_OFF_BORESIGHT_LIM) {
      srt_link_status = false;
    } else {
      I = (int)(*srt_AZ * 180.0 / dpi);
      if (I < 0) {
        I += 360;
      }
      if (*(srt_link->mask + srt_link->nzenith*I + J) == true) {
        srt_link_status = false;
      } else {
        obj_rot_axis[0] *= -1.0;
        obj_rot_axis[1] *= -1.0;
        obj_rot_axis[2] *= -1.0;
        vector_rotation(obj_rot_axis,
                        srt_link->FOV_rot_axis, (*srt_link).FOV_rot_angle);
        omega *= -1.0;
        azel_rot_speed(a, obj_rot_axis, omega,
                       *srt_AZ, *srt_EL, srt_dAZdt, srt_dELdt);
        if (*srt_dAZdt > SRT_DAZDT_LIM || *srt_dELdt > SRT_DELDT_LIM) {
          srt_link_status = false;
        }
      }
    }

#ifdef __DEBUG__
    if (*srt_dAZdt * 180.0 / dpi >= 1.0 ||
        *srt_dELdt * 180.0 / dpi >= 1.0) {
      printf("__DEBUG__   (%f, %f)    (%f, %f)\n",
              *srt_AZ    * 180.0 / dpi, *srt_EL    * 180.0 / dpi,
               *srt_dAZdt * 180.0 / dpi, *srt_dELdt * 180.0 / dpi);
    }
#endif /* __DEBUG__ */
 
  }

/*
--------
*/

  N = 0;
  if (*EL < trk_pos.ELLIM) {
    N += 4;
  }
  if (*dAZdt > trk_pos.AZSV || *dELdt > trk_pos.ELSV) {
    N += 2;
  }
  if (srt_link_status == false) {
    N += 1;
  }

  return N;
}
