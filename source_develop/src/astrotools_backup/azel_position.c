#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"

void azel_position(int  *TimUTC, double UT1_UTC,
                  double Lambda, double fai,      double h,
                  double T,      double p,        _Bool  ATM_SWT,
                  double RA,     double DC,
                  double *AZ,    double *EL,
                  double *dAZdt, double *dELdt,
                  double omega,  double *z)
{
  double x[3], y[3], e[3], q[4], A[3][3];
  double A1[3][3], A2[3][3];
  double U[3], dUdt[3], V;
  double rotz, rotx;
  double dpi = 3.141592653589793238462643;


  x[0] = cos(DC) * cos(RA);
  x[1] = cos(DC) * sin(RA);
  x[2] = sin(DC);
  rotz = -LST(TimUTC, 0.0, UT1_UTC, Lambda) - 0.5*dpi;
  rotx = fai - 0.5*dpi;

  drotate(x, rotz, "z");
  drotate(x, rotx, "x");
  *AZ = atan2(x[1], x[0]);
  *EL = atan2(x[2], vlen2(x));

  if (omega == 0.0) {
    if (ATM_SWT == true) {
      *EL += atmospheric_effect(*EL, T, p);
    }
    return;
  }

/*
-------------------------------
*/

  drotate(z, rotz, "z");
  drotate(z, rotx, "x");

  azel_rot_speed(x, z, omega, *AZ, *EL, dAZdt, dELdt);

/*
-------------------------------
*/

  if (ATM_SWT == true) {
    *EL += atmospheric_effect(*EL, T, p);
  }
  return;
}




double atmospheric_effect(double z, double T, double p)
{
  double R;
  double Ct, Cp;
  double dpi = 3.141592653589793238462643;

  z = 0.5 * dpi - z;
  Ct = 273.15 / (273.15 + T);
  Cp = p / 1013.25;
  R = 60.0615 * tan(z) - 0.0841* pow(tan(z), 3.0) * Ct * Cp;
  R *= dpi / 180.0 / 3600.0;
  return R;
}



void azel_rot_speed(double *x,     double *z, double omega,
                    double AZ,     double EL,
                    double *dAZdt, double *dELdt)
{
  double y[3], q[4], A[3][3];
  double U[3], dUdt[3], V;

  y[0] = z[1] * x[2] - z[2] * x[1];
  y[1] = z[2] * x[0] - z[0] * x[2];
  y[2] = z[0] * x[1] - z[1] * x[0];
  xyz2q(x, y, z, A[0], q);

  U[0]  = A[0][0];
  U[1]  = A[0][1];
  U[2]  = A[0][2];

  dUdt[0] = omega * A[1][0];
  dUdt[1] = omega * A[1][1];
  dUdt[2] = omega * A[1][2];

  V = vlen2(U);
  *dAZdt = pow(cos(AZ), 2.0)
             * (dUdt[0]*U[1] - U[0]*dUdt[1]) / pow(U[0], 2.0);
  *dELdt = pow(cos(EL), 2.0)
             * (dUdt[2]*V - U[2]*(U[0]*dUdt[0]+U[1]*dUdt[1])/V) / (V*V);

  return;
}
