#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>

double force_model(double *, int  *,  double );

#define __TRUE_ANOMALY__
#ifndef __TRUE_ANOMALY__
  #define __ECCENTRIC_ANOMALY__
#endif


void spacecraft_position(struct srt_orbit_parameter srt,
                         int    *TimUTC, double UT1_UTC, double   h,
                         double *P_xyz,  double *Pe_xyz, double *Oe_rac,
                         double *V_xyz,  double *init_l)
{
  double l, f, r, u, q;
  double E, E1, E2;
  double cos_u, sin_u;
  double omega, Omega;
  double d_day, delta;
  double fai;

/*
--------
*/

  d_day = MJD(TimUTC[0], TimUTC[1], TimUTC[2],
              TimUTC[3], TimUTC[4], TimUTC[5], UT1_UTC) - srt.t0;
  l = srt.n * 86400.0 * d_day;
  fai = srt.initial_phase + 2.0 * dpi * d_day;

/*
-------- SATELLITE POSITION --------------------------
*/

  delta = dpi / 180.0 / 3600.0 * 1.0e-3;
  cos_u = cos(*init_l);
  sin_u = sin(*init_l);

  while (1) {
    u = *init_l - (*init_l - srt.e * sin_u - l) / (1.0 - srt.e * cos_u);
    if (fabs(u - *init_l) < delta) {
      *init_l = u;
      break;
    } else {
      *init_l = u;
      cos_u = cos(*init_l);
      sin_u = sin(*init_l);
    }
  }

  r = srt.a * (1.0 - srt.e*cos_u);
  f = atan2(srt.ef*sin_u, cos_u-srt.e);
  q = sqrt(GM * srt.a) / r;

  P_xyz[0] = srt.a          * (cos_u - srt.e);
  P_xyz[1] = srt.a * srt.ef *  sin_u;
  P_xyz[2] = 0.0;

  V_xyz[0] = -q          * sin_u;
  V_xyz[1] =  q * srt.ef * cos_u;
  V_xyz[2] = 0.0;

/*
-------- ERROR VECTOR --------------------------------
*/

  E1 = srt.ODDA * sin(0.50 * l          );
  E2 = srt.ODDP * sin(0.50 * l + dpi/2.0);
  E = sqrt(E1*E1 + E2*E2);

  Oe_rac[0] = 0.0;
  Oe_rac[1] = E;
  Oe_rac[2] = 0.0;

#ifdef __TRUE_ANOMALY__
  drotate(Oe_rac,  f + fai,  "x");
#elif defined __ECCENTRIC_ANOMALY__
  drotate(Oe_rac,  l + fai,  "x");
#endif /* ANOMALY */

  drotate(Oe_rac, atan2(E2, fabs(E1)),     "y");
  Pe_xyz[0] = Oe_rac[0];
  Pe_xyz[1] = Oe_rac[1];
  Pe_xyz[2] = Oe_rac[2];
  drotate(Pe_xyz, f,                       "z");

/*
------------------------------------------------------
*/

  omega = srt.omega + srt.d_omega * d_day;
  Omega = srt.Omega + srt.d_Omega * d_day;

/*
------------------------------------------------------
*/

  drotate(P_xyz, omega,            "z");
  drotate(P_xyz, srt.inclination,  "x");
  drotate(P_xyz, Omega,            "z");

  drotate(V_xyz, omega,            "z");
  drotate(V_xyz, srt.inclination,  "x");
  drotate(V_xyz, Omega,            "z");

  drotate(Pe_xyz, omega,           "z");
  drotate(Pe_xyz, srt.inclination, "x");
  drotate(Pe_xyz, Omega,           "z");

/*
------------------------------------------------------
*/

  return;
}


void orbit_propagation(struct srt_orbit_parameter srt,
                       int    *TimUTC, double UT1_UTC, double  h,
                       double *P_xyz,  double *Pe_xyz, double  *Oe_rac,
                       double *V_xyz,  double *init_l, int     RK_stage,
                       double *va,     double *vb,     double  *vc)
{
  int    i, I, J, N;
  double t, r, F;
  double p[3], v[3];
  double ka[13][3], kv[13][3];
  double vb7[13], vb8[13];

/*
---- Initial
*/

  F = force_model(P_xyz, TimUTC, UT1_UTC);
  for (N=0; N<3; N++) {
    kv[0][N] = f1(t, P_xyz[N], V_xyz[N], F);
    ka[0][N] = f2(t, P_xyz[N], V_xyz[N], F);
  }

/*
---- Loop
*/

  for (J=1; J<RK_stage; J++) {
    for (N=0; N<3; N++) {
      p[N] = 0.0;
      v[N] = 0.0;
      for (i=0; i<J; i++) {
        p[N] += h * *(va + RK_stage*J + i) * *(kv[0] + 3*i + N);
        v[N] += h * *(va + RK_stage*J + i) * *(ka[0] + 3*i + N);
      }
      p[N] += P_xyz[N];
      v[N] += V_xyz[N];
    }
    F = force_model(p, TimUTC, UT1_UTC);
    for (N=0; N<3; N++) {
      *(kv[0] + 3*J + N) = f1(t+h*vc[J], p[N], v[N], F);
      *(ka[0] + 3*J + N) = f2(t+h*vc[J], p[N], v[N], F);
    }
  }

/*
---- Aquisition
*/

  for (N=0; N<3; N++) {
    p[N] = 0.0;
    v[N] = 0.0;
    for (J=0; J<RK_stage; J++) {
      p[N] += vb[J] * *(kv[0] + 3*J + N);
      v[N] += vb[J] * *(ka[0] + 3*J + N);
    }
    P_xyz[N] += h * p[N];
    V_xyz[N] += h * v[N];
  }

  return;
}


double f1(double t, double x, double v, double f)
{
  return v;
}

double f2(double t, double x, double v, double f)
{
  return f * x;
}


double force_model(double *p, int  *TimUTC,  double UT1_UTC)
{
  double r;

  r = vlen3(p);
  return -GM / r / r / r;
}


void  RK07_init(double *va, double *vb, double *vc)
{
  int    i, RK_stage=13;
  double vb7[13], vb8[13];

/*
---------------------
*/

  i = RK_stage *  0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  1;
  *(va+i++) = 1.0 / 18.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  2;
  *(va+i++) = 1.0 / 48.0;
  *(va+i++) = 1.0 / 16.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  3;
  *(va+i++) = 1.0 / 32.0;
  *(va+i++) = 0.0;
  *(va+i++) = 3.0 / 32.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  4;
  *(va+i++) = 5.0 / 16.0;
  *(va+i++) = 0.0;
  *(va+i++) = -75.0 / 64.0;
  *(va+i++) =  75.0 / 64.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  5;
  *(va+i++) = 3.0 / 80.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 3.0 / 16.0;
  *(va+i++) = 3.0 / 20.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  6;
  *(va+i++) = 29443841.0 / 614563906.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) =  77736538.0 / 692538347.0;
  *(va+i++) = -28693883.0 / 1125000000.0;
  *(va+i++) = 23124283.0 / 1800000000.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  7;
  *(va+i++) = 16016141.0 / 946692911.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 61564180.0 / 158732637.0;
  *(va+i++) = 22789713.0 / 633445777.0;
  *(va+i++) = 545815736.0 / 2771057229.0;
  *(va+i++) = -180193667.0 / 1043307555.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  8;
  *(va+i++) = 39632708.0 /573591083.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = -433636366.0 / 683701615.0;
  *(va+i++) = -421739975.0 / 2616292301.0;
  *(va+i++) = 100302831.0 / 723423059.0;
  *(va+i++) = 790204164.0 / 839813087.0;
  *(va+i++) = 800635310.0 / 3783071287.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage *  9;
  *(va+i++) = 246121993.0 / 1340847787.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = -37695042795.0 / 15268766246.0;
  *(va+i++) = -309121744.0 / 1061227803.0;
  *(va+i++) = -12992083.0 / 490766935.0;
  *(va+i++) = 6005943493.0 / 2108947869.0;
  *(va+i++) = 393006217.0 / 1396673457.0;
  *(va+i++) = 123872331.0 / 1001029789.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage * 10;
  *(va+i++) = -1028468189.0 / 846180014.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = 8478235783.0 / 508512852.0;
  *(va+i++) = 1311729495.0 / 1432422823.0;
  *(va+i++) = -10304129995.0 / 1701304382.0;
  *(va+i++) = -48777925059.0 / 3047939560.0;
  *(va+i++) = 15336726248.0 / 1032824649.0;
  *(va+i++) = -45442868181.0 / 3398467696.0;
  *(va+i++) = 3065993473.0 / 597172653.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage * 11;
  *(va+i++) = 185892177.0 / 718116043.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = -3185094517.0 / 667107341.0;
  *(va+i++) = -477755414.0 / 1098053517.0;
  *(va+i++) = -703635378.0 / 230739211.0;
  *(va+i++) = 5731566787.0 / 1027545527.0;
  *(va+i++) = 5232866602.0 / 850066563.0;
  *(va+i++) = -4093664535.0 / 808688257.0;
  *(va+i++) = 3962137247.0 / 1805957418.0;
  *(va+i++) = 65686358.0 / 487910083.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

  i = RK_stage * 12;
  *(va+i++) = 403863854.0 / 491063109.0;
  *(va+i++) = 0.0;
  *(va+i++) = 0.0;
  *(va+i++) = -5068492393.0 / 434740067.0;
  *(va+i++) = -411421997.0 / 543043805.0;
  *(va+i++) = 652783627.0 / 914296604.0;
  *(va+i++) = 11173962825.0 / 925320556.0;
  *(va+i++) = -13158990841.0 / 6184727034.0;
  *(va+i++) = 3936647629.0 / 1978049680.0;
  *(va+i++) = -160528059.0 / 685178525.0;
  *(va+i++) = 248638103.0 / 1413531060.0;
  *(va+i++) = 0.0;
  *(va+i  ) = 0.0;

/*
--------
*/

  vb7[ 0] = 13451932.0 / 455176623.0;
  vb7[ 1] = 0.0;
  vb7[ 2] = 0.0;
  vb7[ 3] = 0.0;
  vb7[ 4] = 0.0;
  vb7[ 5] = -808719846.0 / 976000145.0; 
  vb7[ 6] = 1757004468.0 / 5645159321.0;
  vb7[ 7] = 656045339.0 / 265891186.0;
  vb7[ 8] = -3867574721.0 / 1518517206.0; 
  vb7[ 9] = 465885868.0 / 322736535.0;
  vb7[10] = 53011238.0 / 667516719.0;
  vb7[11] = 2.0 / 45.0;
  vb7[12] = 0.0;

  vb8[ 0] = 14005451.0 / 335480064.0;
  vb8[ 1] = 0.0;
  vb8[ 2] = 0.0;
  vb8[ 3] = 0.0;
  vb8[ 4] = 0.0;
  vb8[ 5] = -59238493.0 / 1068277825.0;
  vb8[ 6] = 181606767.0 / 758867731.0; 
  vb8[ 7] = 561292985.0 / 797845732.0;
  vb8[ 8] = -1041891430.0 / 1371343529.0;
  vb8[ 9] = 760417239.0 / 1151165299.0; 
  vb8[10] = 118820643.0 / 751138087.0;
  vb8[11] = -528747749.0 / 2220607170.0;
  vb8[12] = 1.0 / 4.0;

  for (i=0; i<13; i++) {
    vb[i] = vb7[i];
  }

  vc[ 0] = 0.0;
  vc[ 1] = 1.0 / 18.0;
  vc[ 2] = 1.0 / 12.0;
  vc[ 3] = 1.0 / 8.0;
  vc[ 4] = 5.0 / 16.0;
  vc[ 5] = 3.0 / 8.0;
  vc[ 6] = 59.0 / 400.0;
  vc[ 7] = 93.0 / 200.0;
  vc[ 8] = 5490023248.0 / 9719169821.0;
  vc[ 9] = 13.0 / 20.0;
  vc[10] = 1201146811.0 / 1299019798.0;
  vc[11] = 1.0;
  vc[12] = 1.0;

  return;
}



void RK04_init(double *va, double *vb, double *vc)
{
  int   i, RK_stage=4;

  i = RK_stage *  0;
  *(va+i++) = 0.0; *(va+i++) = 0.0; *(va+i++) = 0.0; *(va+i  ) = 0.0;
  i = RK_stage *  1;
  *(va+i++) = 0.5; *(va+i++) = 0.0; *(va+i++) = 0.0; *(va+i  ) = 0.0;
  i = RK_stage *  2;
  *(va+i++) = 0.0; *(va+i++) = 0.5; *(va+i++) = 0.0; *(va+i  ) = 0.0;
  i = RK_stage *  3;
  *(va+i++) = 0.0; *(va+i++) = 0.0; *(va+i++) = 1.0; *(va+i  ) = 0.0;

  vb[0] = 1.0 / 6.0;
  vb[1] = 2.0 / 6.0;
  vb[2] = 2.0 / 6.0;
  vb[3] = 1.0 / 6.0;

  vc[0] = 0.0;
  vc[1] = 0.5;
  vc[2] = 0.5;
  vc[3] = 1.0;

  return;
}


void  RK_end(double *va, double *vb, double *vc)
{
  free (va);
  free (vb);
  free (vc);

  return;
}
