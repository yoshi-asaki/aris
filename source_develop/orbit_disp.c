#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <unistd.h>
#include <aris.h>

void  oe_plot_label(float , float , float , float ,
                    char *, char *, char *, char *, int *);

/****
#define __DESCRIPTION__
#define __ORBIT_ERROR__
****/

/****
#define __ORBD_ANIMATION__
****/
#define __ORBD_ANIMATION__
#ifdef __ORBD_ANIMATION__
/****
  #define __ANIMA_CLOCK__
  #define __TRACKING_STATUS__
****/
  #define __ANIMA_CLOCK__
  #define __TRACKING_STATUS__
  #define __TRT_STATUS__
#endif


/****
#define __COLOR__
****/
#define __COLOR__
#ifndef __COLOR__
  #define __BW__
#endif


/****
#define __PROPAGATION__
****/


void orbit_disp(int SRT_NUM, int  nobs,
                struct srt_orbit_parameter *srt,
                int *TimUTC, double UT1_UTC,
                double sep_angle_limit,
                struct srt_data_link *srt_link,
                double OBS_T,       double OBS_p,
                int    TRK_NUM,
                int    GRT_NUM,
                struct st_observable  *int_obs[],
                struct antenna_parameter *trk_pos,
                struct source_parameter src,
                struct source_parameter sun,
                char   *ascii_out)
{
  int    i, j, I, iant, ib, iobs, NOBS;
  int    iiant;
  int    mrk;
  _Bool  shadow_swt, color_swt, axis_swt;
  float  fai;
  float  pgxmin, pgxmax, pgymin, pgymax;
  float  vpxmin, vpxmax, vpymin, vpymax;
  float  *pgx, *pgy;
  float  *pg1, *pg2, *pg3;
  float  xmin, xmax, ymin, ymax;
  float  pgxx[10], pgyy[10];
  double R, L, P[3], E[3], V_xyz[3], Oe[3];
  double p[3];
  int    timUTC[6];
  double *srt_pos, dpos[3];
  double DPI;
  char   string[100];
  int    clrplt[10];
  double A, B, C, D;
  double ang1, ang2, ang3;
  double GST_theta;
  double Omega[SRTMAX], omega[SRTMAX], d_day[SRTMAX];
  double mean_anomaly[SRTMAX];
  double eccentric_anomaly[SRTMAX];
  double init_l[SRTMAX];
  int    RK_stage = 13;

  float  pg_day[2][1080];
  float  pg_ngt[2][1080];
  float  pg_shd[2][721];

  int    itrk;
  double AZ, EL, dAZdt, dELdt;
  double srt_AZ, srt_EL, srt_dAZdt, srt_dELdt;
  double ptrk[3];
  int    trk_status;
  int    onbd_status, gtrk_status;
  float  Re2, ri, ro;
  _Bool  mrk_swt[10];
  float  tpos[10][2];
  float  a, b, c, d, e, rim_x, rim_y;
  _Bool  draw_swt[SRTMAX];

  float  clkx[2][2], clky[2][2], hh, mm;
  float  hand[2], dial;
  float  dial_circ[2][361];
  float  pg_width;
  float  charx[4], chary[4];
  struct earth_shape earth_shape;

  FILE   *log_fp;

#ifdef __PROPAGATION__
  double *va, *vb, *vc;
#endif

/*
----------------------------------------------
*/

#ifdef __PROPAGATION__
  if ((va = (double *)calloc(RK_stage*RK_stage, sizeof(double))) == NULL) {
    printf("ERROR: va in RK_init: Stop.\n");
    exit (0);
  }
  if ((vb = (double *)calloc(RK_stage,          sizeof(double))) == NULL) {
    printf("ERROR: vb in RK_init: Stop.\n");
    free (va);
    exit (0);
  }
  if ((vc = (double *)calloc(RK_stage,          sizeof(double))) == NULL) {
    printf("ERROR: vc in RK_init: Stop.\n");
    free (va);
    free (vb);
    exit (0);
  }
  RK07_init(va, vb, vc);
#endif

/*
----------------------------------------------
*/

  DPI = dpi / 180.0;

  for (itrk=0; itrk<TRK_NUM; itrk++) {
    mrk_swt[itrk] = false;
  }

/*
----------------------------------------------
*/

  if (ascii_out[0] == '!') {
    if ((log_fp = fopen(ascii_out+1, "w")) == NULL) {
      printf("Warning; ORBIT_DISP: %s cannot be made.\n", ascii_out+1);
    } else {
      NOBS = nobs / TIME_STEP;
      if (NOBS < 360) {
        NOBS = 360;
      }
      for (i=0; i<NOBS; i++) {
        iobs = i * TIME_STEP;
        timUTC[5] = TimUTC[5] + iobs;

        for (iant=0; iant<SRT_NUM; iant++) {
          ib = 40 * iant;

#ifdef __PROPAGATION__
          orbit_propagation(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                        P, E, Oe, V_xyz, init_l+iant, RK_stage, va, vb, vc);
#else
          spacecraft_position(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                        P, E, Oe, V_xyz, init_l+iant);
#endif

          if (log_fp != NULL) {
            fprintf(log_fp,
                "%7.1lf,  %lf,  %lf,  %lf,  %lf,  %lf,  %lf,  ",
                (double)iobs,
                P[0], P[1], P[2], V_xyz[0], V_xyz[1], V_xyz[2]);
            fprintf(log_fp,
                "%lf,  %lf,  %lf,  %d\n", E[0], E[1], E[2],
                earth_eclipse(sun.s, P, &R, &L, sep_angle_limit));
          }

        }
      }
      if (log_fp != NULL) {
        fclose (log_fp);
      }
    }
    return;
  }

/*
----------------------------------------------
*/

#ifdef __COLOR__
  color_swt = true;
#elif defined __BW__
  color_swt = false;
#endif

  pg_width = 0.0;
  for (iant=0; iant<SRT_NUM; iant++) {
    if (srt[iant].apogee > pg_width) {
      pg_width = srt[iant].apogee;
    }
  }
  pg_width += earth_radius;
  pg_width *= 1.2e-3;
  axis_swt  = true;

#ifdef __ORBD_ANIMATION__
  axis_swt  = false;
  vpxmin =  0.08;
  vpxmax =  0.94;
  vpymin =  0.08;
  vpymax =  0.94;
  pgxmin = -pg_width;
  pgxmax =  pg_width;
  pgymin = -pg_width;
  pgymax =  pg_width;
  if (pgxmax - pgxmin > pgymax - pgymin) {
    pg_width = pgymax - pgymin;
  } else {
    pg_width = pgxmax - pgxmin;
  }
  clkx[0][0] = pgxmax - 0.15 * pg_width;
  clky[0][0] = pgymax - 0.15 * pg_width;
  clkx[1][0] = pgxmax - 0.15 * pg_width;
  clky[1][0] = pgymax - 0.15 * pg_width;
  hand[0]    = 0.04 * pg_width;
  hand[1]    = 0.06 * pg_width;
  dial       = 0.08 * pg_width;
  for (i=0; i<=360; i++) {
    dial_circ[0][i] = clkx[0][0] + 0.09 * pg_width * cos((float)i*DPI);
    dial_circ[1][i] = clky[0][0] + 0.09 * pg_width * sin((float)i*DPI);
  }
  charx[0] = clkx[0][0];
  chary[0] = clky[0][0] + 0.08 * pg_width - 0.010 * pg_width;
  charx[1] = clkx[0][0] + 0.08 * pg_width - 0.010 * pg_width;
  chary[1] = clky[0][0]                   - 0.005 * pg_width;
  charx[2] = clkx[0][0];
  chary[2] = clky[0][0] - 0.08 * pg_width;
  charx[3] = clkx[0][0] - 0.08 * pg_width - 0.005 * pg_width;
  chary[3] = clky[0][0];
#else
  axis_swt  = true;
  vpxmin =  0.08;
  vpxmax =  0.47;
  vpymin =  0.60;
  vpymax =  0.99;
  pgxmin = -3.7e4;
  pgxmax =  3.7e4;
  pgymin = -3.7e4;
  pgymax =  3.7e4;
/********
********/

/****
  vpxmin =  0.08;
  vpxmax =  0.47;
  vpymin =  0.67;
  vpymax =  0.98;

  vpxmin =  0.08;
  vpxmax =  0.86;
  vpymin =  0.36;
  vpymax =  0.98;
  pgxmin = -3.7e4;
  pgxmax =  3.7e4;
  pgymin = -2.0e4;
  pgymax =  3.7e4;
****/

#endif /* __ORBD_ANIMATION__ */

  ang1 = -src.RA - dpi/2.0;
  ang2 =  src.DC - dpi/2.0;

#ifdef __COLOR__
  for (i=0; i<10; i++) {
    clrplt[i] = i;
  }
#elif defined __BW__
  for (i=0; i<10; i++) {
    cpgscr(i+1, (float)i/10.0, (float)i/10.0, (float)i/10.0);
    clrplt[i] = i + 1;
  }
#endif

#ifdef __ORBD_ANIMATION__
  shadow_swt = false;
#else
  shadow_swt = true;
#endif /* __ORBD_ANIMATION__ */
/*XXXXXXXXX*/
  shadow_swt = true;

  Re2 = earth_radius * earth_radius * 1.0e-6;

  for (i=0; i<6; i++) {
    timUTC[i] = TimUTC[i];
  }

  NOBS = nobs / TIME_STEP;
  if (NOBS < 360) {
    NOBS = 360;
  }
  pgx = (float *)calloc(NOBS, sizeof(float));
  pgy = (float *)calloc(NOBS, sizeof(float));
  pg1 = (float *)calloc(NOBS, sizeof(float));
  pg2 = (float *)calloc(NOBS, sizeof(float));
  pg3 = (float *)calloc(NOBS, sizeof(float));
  NOBS = nobs / TIME_STEP;

  if ((srt_pos=(double *)calloc(2*3*SRT_NUM, sizeof(double))) == NULL) {
    printf("ERROR: orbit_disp: memory allocation of srt_pos\n");
    return;
  }

/****
  cpgpap(1.00*pgpap_prm, 1.0);
  cpgpap(0.50*pgpap_prm, 1.0);
  cpgpap(0.70*pgpap_prm, 1.0);
  cpgpap(0.75*pgpap_prm, 1.0);
  cpgpap(1.50*pgpap_prm, 1.0);
  cpgpap(0.75*pgpap_prm, 1.0);
****/
  cpgpap(1.00*pgpap_prm, 1.0);

  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);

/*
------------------------------------------
*/

  cpgsci(1);
  cpgsch(0.75);

  cpgsvp(vpxmin, vpxmax, vpymin, vpymax);
  cpgswin(pgxmin, pgxmax, pgymin, pgymax);
  cpgbox("BCNTS", 0, 4.0e+4, "BCNTS", 0, 4.0e+4);
/****xxxx
****/

/********YYYY
  cpglab("\\fiu \\fn[km]", "\\fiv \\fn[km]",
         "Projected orbit to (\\fiu\\fn, \\fiv\\fn)");
********/
  cpglab("\\fiu \\fn[km]", "\\fiv \\fn[km]", "");
/*xxxx
xxxx*/
  GST_theta = GST(TimUTC, 0.0, UT1_UTC);
  draw_earth(earth_radius, EPSIRON(ET(timUTC, UT1_UTC)), GST_theta,
             src, sun, &earth_shape,
             clrplt, shadow_swt, color_swt, axis_swt);

  for (iant=0; iant<SRT_NUM; iant++) {
    draw_swt[iant] = false;
  }

  for (iant=0; iant<SRT_NUM; iant++) {
    init_l[iant] = dpi;
  }

  for (iant=0; iant<SRT_NUM; iant++) {
    d_day[iant] = MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                      TimUTC[3], TimUTC[4], TimUTC[5], UT1_UTC)
                - srt[iant].t0;
    spacecraft_position(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                        P, E, Oe, V_xyz, init_l+iant);
    Omega[iant] = srt[iant].Omega + srt[iant].d_Omega * d_day[iant];
    omega[iant] = srt[iant].omega + srt[iant].d_omega * d_day[iant];

    Omega[iant] -= 2.0 * dpi * (double)floor(Omega[iant]  / 2.0 / dpi);
    omega[iant] -= 2.0 * dpi * (double)floor(omega[iant]  / 2.0 / dpi);
    eccentric_anomaly[iant] = init_l[iant]
                 - 2.0 * dpi * (double)floor(init_l[iant] / 2.0 / dpi);
    mean_anomaly[iant] = 
       eccentric_anomaly[iant] - srt[iant].e * sin(eccentric_anomaly[iant]);
  }

  for (i=0; i<NOBS; i++) {
    cpgbbuf();
    iobs = i * TIME_STEP;

#ifdef __ORBD_ANIMATION__
    usleep(50000);
    cpgsch(0.75);
    cpgbox("BCNTS", 0, 4.0e+4, "BCNTS", 0, 4.0e+4);
    cpgsci(0);
    for (iant=0; iant<SRT_NUM; iant++) {
      if (draw_swt[iant] == true) {
        ib = 40 * iant;

/*** for erase the orbit trajectries ****/
        cpgpt(1, pgx+ib, pgy+ib, mrk);

        cpgsch(0.50);
        cpgarro(pgx[ib+1], pgy[ib+1], pgx[ib+2], pgy[ib+2]);
        cpgsch(0.75);

        cpgsls(4);
        cpgline(2, pgx+ib+3, pgy+ib+3);
        cpgline(2, pgx+ib+5, pgy+ib+5);
        for (itrk=0; itrk<TRK_NUM; itrk++) {
          cpgline(2, pgx+ib+20+2*itrk, pgy+ib+21+2*itrk);
        }
      }
    }
    cpgsls(1);
  #ifdef __ANIMA_CLOCK__
    cpgline(2, clkx[0], clky[0]);
    cpgline(2, clkx[1], clky[1]);
  #endif /*__ANIMA_CLOCK__*/

    cpgsci(8);
    cpgpt(earth_shape.n_day, earth_shape.s_day[0], earth_shape.s_day[1], 1);
    cpgsci(4);
    cpgpt(earth_shape.n_ngt, earth_shape.s_ngt[0], earth_shape.s_ngt[1], 1);
#endif /* __ORBD_ANIMATION__ */

    timUTC[5] = TimUTC[5] + iobs;
    for (iant=0; iant<SRT_NUM; iant++) {
      ib = 40 * iant;

#ifdef __PROPAGATION__
      orbit_propagation(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                        P, E, Oe, V_xyz, init_l+iant, RK_stage, va, vb, vc);
#else
      spacecraft_position(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                        P, E, Oe, V_xyz, init_l+iant);
#endif

      srt_pos[3*iant+1] = P[1];
      srt_pos[3*iant+2] = P[2];
      drotate(srt_pos+3*iant, ang1, "z");
      drotate(srt_pos+3*iant, ang2, "x");
      pgx[ib] = (float)srt_pos[3*iant  ] * 1.0e-3;
      pgy[ib] = (float)srt_pos[3*iant+1] * 1.0e-3;


      draw_swt[iant] = false;
      if (! (vlen2(srt_pos+3*iant) <= earth_radius &&
             srt_pos[3*iant+2] <  0.0)) {
        draw_swt[iant] = true;
        if (earth_eclipse(sun.s, P, &R, &L, sep_angle_limit) == NIGHT) {
          cpgsci(clrplt[4]);
#ifndef __ORBD_ANIMATION__
          mrk = 2;
#endif
        } else {
          if (earth_eclipse(src.s, P, &R, &L, sep_angle_limit) == NIGHT) {
            cpgsci(clrplt[1]);
          } else {
            cpgsci(clrplt[iant+2]);
          }
#ifndef __ORBD_ANIMATION__
          mrk = 1;
#endif
        }
#ifdef __ORBD_ANIMATION__
        mrk = 17;
#endif /* __ORBD_ANIMATION__ */
        cpgpt(1, pgx+ib, pgy+ib, mrk);

#ifdef __ORBD_ANIMATION__
        cpgsci(1);
        pgx[ib+1] = pgx[ib];
        pgy[ib+1] = pgy[ib];

        dpos[0] = E[0];
        dpos[1] = E[1];
        dpos[2] = E[2];
        drotate(dpos, ang1, "z");
        drotate(dpos, ang2, "x");
        pgx[ib+2] = (float)dpos[0] * 0.5e6 + pgx[ib];
        pgy[ib+2] = (float)dpos[1] * 0.5e6 + pgy[ib];
        cpgsch(0.50);
        cpgarro(pgx[ib+1], pgy[ib+1], pgx[ib+2], pgy[ib+2]);
        cpgsch(0.75);

        A = -atan2(E[2], vlen2(E));
        B =  atan2(E[1],     E[0]);
        C = vlen3(E) * 0.5e6;
        cpgsls(4);

        pgx[ib+3] = pgx[ib+1];
        pgy[ib+3] = pgy[ib+1];
        dpos[0] = 0.0;
        dpos[1] = C;
        dpos[2] = 0.0;
        drotate(dpos, A, "y");
        drotate(dpos, B, "z");
        drotate(dpos, ang1, "z");
        drotate(dpos, ang2, "x");
        pgx[ib+4] = (float)dpos[0] + pgx[ib];
        pgy[ib+4] = (float)dpos[1] + pgy[ib];
        cpgline(2, pgx+ib+3, pgy+ib+3);

        pgx[ib+5] = pgx[ib+1];
        pgy[ib+5] = pgy[ib+1];
        dpos[0] = 0.0;
        dpos[1] = 0.0;
        dpos[2] = C;
        drotate(dpos, A, "y");
        drotate(dpos, B, "z");
        drotate(dpos, ang1, "z");
        drotate(dpos, ang2, "x");
        pgx[ib+6] = (float)dpos[0] + pgx[ib];
        pgy[ib+6] = (float)dpos[1] + pgy[ib];
        cpgline(2, pgx+ib+5, pgy+ib+5);

        cpgsls(1);

#endif /* __ORBD_ANIMATION__ */

/*
---- TRACKING CONDITION ----
*/

#ifdef __ORBD_ANIMATION__
        for (itrk=0; itrk<TRK_NUM; itrk++) {
          cpgsci(0);
          cpgline(2, pgx+ib+20+2*itrk, pgy+ib+20+2*itrk);
          if (mrk_swt[itrk] == true) {
            cpgpt(1, tpos[itrk], tpos[itrk]+1, 22);
          }
        }
#endif /* __ORBD_ANIMATION__ */

#ifdef __TRACKING_STATUS__
        for (itrk=0; itrk<TRK_NUM; itrk++) {
          cpgsci(itrk+2);
          onbd_status = __GO__;
          gtrk_status = __GO__;
          trk_status = tracking_condition(timUTC, UT1_UTC, trk_pos[itrk],
                                          OBS_T, OBS_p, false, P, V_xyz,
                                          &AZ, &EL, &dAZdt, &dELdt,
                                          srt_link, src, sun,
                                          &srt_AZ, &srt_EL,
                                          &srt_dAZdt, &srt_dELdt,
                                          ptrk, dpos, srt[iant].BODY_X_SUN);
          drotate(ptrk, ang1, "z");
          drotate(ptrk, ang2, "x");
          drotate(dpos, ang1, "z");
          drotate(dpos, ang2, "x");

          mrk_swt[itrk] = false;
          if (ptrk[2] > 0.0) {
            mrk_swt[itrk] = true;
            tpos[itrk][0] = (float)ptrk[0] * 1.0e-3;
            tpos[itrk][1] = (float)ptrk[1] * 1.0e-3;
            cpgpt(1, tpos[itrk], tpos[itrk]+1, 22);
          }

          pgx[ib+10] = 0.0;
          pgy[ib+10] = 0.0;
          pgx[ib+11] = 0.0;
          pgy[ib+11] = 0.0;

          if (trk_status == 0) {
            pgx[ib+10] = (float)ptrk[0] * 1.0e-3;
            pgy[ib+10] = (float)ptrk[1] * 1.0e-3;
            pgx[ib+11] = pgx[ib];
            pgy[ib+11] = pgy[ib];
          } else {
            if (trk_status / 4 == 1) {
              gtrk_status = __NG__;
            }
            trk_status %= 4;
            if (trk_status / 2 == 1) {
              gtrk_status = __NG__;
            }
            trk_status %= 2;
            if (trk_status == 1) {
              onbd_status = __NG__;
            }

            if (gtrk_status == __NG__ && onbd_status == __GO__) {
              pgx[ib+11] = pgx[ib];
              pgy[ib+11] = pgy[ib];
              pgx[ib+10] = pgx[ib+11] - 0.2e-3 * (float)dpos[0];
              pgy[ib+10] = pgy[ib+11] - 0.2e-3 * (float)dpos[1];
            } else if (gtrk_status == __GO__ && onbd_status == __NG__) {
              pgx[ib+10] = (float)ptrk[0] * 1.0e-3;
              pgy[ib+10] = (float)ptrk[1] * 1.0e-3;
              pgx[ib+11] = pgx[ib+10] + 0.2e-3 * (float)dpos[0];
              pgy[ib+11] = pgy[ib+10] + 0.2e-3 * (float)dpos[1];
            }
          }

          ri = pgx[ib+10]*pgx[ib+10] + pgy[ib+10]*pgy[ib+10];
          ro = pgx[ib+11]*pgx[ib+11] + pgy[ib+11]*pgy[ib+11];

/*XXXXXXXXXXXXXXXXXXXXXXXXXX*/
          if (ptrk[2] < 0.0) {
            if (ri <  Re2 && ro >= Re2) {
              if (pgx[ib+10] != pgx[ib+11]) {
                a = (pgy[ib+10] - pgy[ib+11]) / (pgx[ib+10] - pgx[ib+11]);
                b = (pgx[ib+10]*pgy[ib+11] - pgx[ib+11]*pgy[ib+10])
                  / (pgx[ib+10] - pgx[ib+11]);
                e = 1.0 + a * a;
                c = -a * b / e;
                d = sqrt(Re2 * e - b*b) / e;

                if (pgx[ib+10] < pgx[ib+11]) {
                  if (pgx[ib+10] <= c + d && pgx[ib+11] >= c + d) {
                    rim_x = c + d;
                  } else {
                    rim_x = c - d;
                  }
                } else {
                  if (pgx[ib+11] <= c + d && pgx[ib+10] >= c + d) {
                    rim_x = c + d;
                  } else {
                    rim_x = c - d;
                  }
                }
                rim_y = a * rim_x + b;
              } else {
                rim_x = pgx[ib+10];
                if (pgy[ib+11] >= 0.0) {
                  rim_y =  sqrt(Re2 - rim_x*rim_x);
                } else {
                  rim_y = -sqrt(Re2 - rim_x*rim_x);
                }
              }
              pgx[ib+10] = rim_x;
              pgy[ib+10] = rim_y;
              cpgline(2, pgx+ib+10, pgy+ib+10);

            } else if (ri >= Re2 && ro >= Re2) {
              cpgline(2, pgx+ib+10, pgy+ib+10);
            }
          } else {
            cpgsls(1);
            cpgline(2, pgx+ib+10, pgy+ib+10);
          }


#ifdef __ORBD_ANIMATION__
          pgx[ib+20+2*itrk] = pgx[ib+10];
          pgy[ib+20+2*itrk] = pgy[ib+10];
          pgx[ib+21+2*itrk] = pgx[ib+11];
          pgy[ib+21+2*itrk] = pgy[ib+11];
#endif /* __ORBD_ANIMATION__ */
        }
#endif /* __TRACKING_STATUS__ */
      }
      srt_pos[3*iant] = P[0];
    }

#ifdef __TRT_STATUS__
    for (iiant=0; iiant<GRT_NUM; iiant++) {
      I = iiant * nobs + iobs;
      if (int_obs[0][I].w > 0.0) {
        pgxx[0] = int_obs[0][I].u * 1.0e-3;
        pgyy[0] = int_obs[0][I].v * 1.0e-3;
        cpgpt(1, pgxx, pgyy, 15);
      }
    }
#endif /* __TRT_STATUS__ */

#ifdef __ORBD_ANIMATION__
  #ifdef __ANIMA_CLOCK__
    hh = (float)timUTC[3] / 12.0
       + (float)timUTC[4] / 60.0 / 12.0
       + (float)timUTC[5] / 3600.0 / 12.0 - 0.25;
    hh *= (float)(-2.0 * dpi);
    mm = (float)timUTC[4] / 60.0 + (float)timUTC[5] / 3600.0 - 0.25;
    mm *= (float)(-2.0 * dpi);

    clkx[0][1] = clkx[0][0] + hand[0] * cos(hh);
    clky[0][1] = clky[0][0] + hand[0] * sin(hh);
    clkx[1][1] = clkx[1][0] + hand[1] * cos(mm);
    clky[1][1] = clky[1][0] + hand[1] * sin(mm);

    cpgsci(1);
    cpgpt(1, clkx[0], clky[0], 17);
    cpgline(2, clkx[0], clky[0]);
    cpgline(2, clkx[1], clky[1]);
    cpgptxt(charx[0], chary[0], 0.0, 0.5, "0");
    cpgptxt(charx[1], chary[1], 0.0, 0.5, "3");
    cpgptxt(charx[2], chary[2], 0.0, 0.5, "6");
    cpgptxt(charx[3], chary[3], 0.0, 0.5, " 9");
    cpgline(361, dial_circ[0], dial_circ[1]);
  #endif /* __ANIMA_CLOCK__ */
#endif /* __ORBD_ANIMATION__ */
    cpgebuf();
  }

#ifdef __ORBD_ANIMATION__
/********
  draw_earth(earth_radius, EPSIRON(ET(timUTC, UT1_UTC)), GST_theta,
             src, sun, &earth_shape,
             clrplt, shadow_swt, color_swt, axis_swt);

  cpgscr(20, 0.05, 0.05, 0.05);
  cpgsci(20);

  cpgsvp(0.00, 0.50, 0.00, 0.42);
  cpgswin(0.0, 1.0, 0.0, 1.0);
  cpgrect(0.0, 1.0, 0.0, 1.0);

  cpgsvp(0.52, 1.00, 0.00, 0.66);
  cpgswin(0.0, 1.0, 0.0, 1.0);
  cpgrect(0.0, 1.0, 0.0, 1.0);
********/

/****
  cpgsvp(0.48, 1.00, 0.68, 1.00);
  cpgswin(0.0, 1.0, 0.0, 1.0);
  cpgrect(0.0, 1.0, 0.0, 1.0);
****/
#endif

/*
--------------------------------------
*/

  for (i=0; i<6; i++) {
    timUTC[i] = TimUTC[i];
  }

/*
--------------------------------------
*/

/********YYYY
  cpgsch(0.75);
  cpgsci(1);
  cpgsvp(0.08, 0.47, 0.10, 0.35);
  xmin = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               TimUTC[3], TimUTC[4], TimUTC[5], UT1_UTC)
                         - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,        0,        0,        UT1_UTC));
  xmax = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               TimUTC[3], TimUTC[4], TimUTC[5]+nobs, UT1_UTC)
                         - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,        0,        0,        UT1_UTC));
  ymin = 0.0;
  ymax = 1.2e-3 * srt[0].apogee;
  cpgswin(xmin, xmax, ymin, ymax);
  cpgtbox("BCNSTZH", 0, 0, "BCNTS", 0, 0);
  cpglab("Time (UT)", "Altitude [km]", "Spacecraft altitude");
  cpgsch(0.75);
  for (iant=0; iant<SRT_NUM; iant++) {
    init_l[iant] = dpi;
    for (i=0; i<NOBS; i++) {
      iobs = i * TIME_STEP;
      timUTC[5] = TimUTC[5] + iobs;
      spacecraft_position(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                          P, E, Oe, V_xyz, init+l+iant);
      pgx[0] = 86400.0 * (float)(MJD(timUTC[0], timUTC[1], timUTC[2],
                                     timUTC[3], timUTC[4], timUTC[5], UT1_UTC)
                               - MJD(timUTC[0], timUTC[1], timUTC[2],
                                     0,        0,        0,        UT1_UTC));
      pgy[0] = (float)(vlen3(P) - earth_radius) * 1.0e-3;

      dpos[0] = P[0];
      dpos[1] = P[1];
      dpos[2] = P[2];
      drotate(dpos, ang1, "z");
      drotate(dpos, ang2, "x");

      if (! (vlen2(dpos) <= earth_radius && dpos[2] <  0.0)) {
        if (earth_eclipse(sun.s, P, &R, &L, sep_angle_limit) == NIGHT) {
          cpgsci(clrplt[4]);
          mrk = 2;
        } else {
          if (earth_eclipse(src.s, P, &R, &L, sep_angle_limit) == NIGHT) {
            cpgsci(clrplt[1]);
          } else {
            cpgsci(clrplt[iant+2]);
          }
          mrk = 1;
        }
        cpgpt(1, pgx, pgy, mrk);
      }
    }
  }
YYYYY********/

/*
------------------------------------------
*/

#ifdef __ORBIT_ERROR__
  cpgsch(0.75);
  cpgsci(1);
  cpgsvp(0.58, 0.97, 0.10, 0.35);
  cpgsvp(0.08, 0.47, 0.07, 0.32);
  xmin = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               TimUTC[3], TimUTC[4], TimUTC[5], UT1_UTC)
                         - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,        0,        0,        UT1_UTC));
  xmax = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               TimUTC[3], TimUTC[4], TimUTC[5]+nobs, UT1_UTC)
                         - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,        0,        0,        UT1_UTC));
  if (srt[0].ODDA != 0.0) {
    ymin = -1.2 * srt[0].ODDA * 1.0e2;
    ymax =  1.9 * srt[0].ODDA * 1.0e2;
  } else {
    ymin = -1.0;
    ymax =  1.0;
  }
  cpgswin(xmin, xmax, ymin, ymax);
  cpgtbox("BCSTNZH", 0, 0, "BCNTS", 0, 0);
  cpglab("Time (UT)", "OD displacement [cm]", "OD Error");
  cpgsch(0.75);
  for (iant=0; iant<SRT_NUM; iant++) {
    I = 0;
    init_l[iant] = dpi;
    for (i=0; i<NOBS; i++) {
      iobs = i * TIME_STEP;
      timUTC[5] = TimUTC[5] + iobs;
      pgx[I] = 86400.0 * (float)(MJD(timUTC[0], timUTC[1], timUTC[2],
                                     timUTC[3], timUTC[4], timUTC[5], UT1_UTC)
                               - MJD(timUTC[0], timUTC[1], timUTC[2],
                                     0,        0,        0,        UT1_UTC));
      spacecraft_position(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                          P, E, Oe, V_xyz, init_l+iant);
      pg1[I] = (float)E[0] * 1.0e2;
      pg2[I] = (float)E[1] * 1.0e2;
      pg3[I] = (float)E[2] * 1.0e2;
      pgy[I] = (float)vlen3(E) * 1.0e2;
      I++;
    }
    cpgsci(clrplt[8]);
    cpgsls(2);
    cpgline(I, pgx, pg1);

    cpgsci(clrplt[3]);
/*xxxx*/
    cpgsci(clrplt[2]);
/*xxxx*/
    cpgsls(3);
    cpgline(I, pgx, pg2);

    cpgsci(clrplt[7]);
/*xxxx*/
    cpgsci(clrplt[4]);
/*xxxx*/
    cpgsls(4);
    cpgline(I, pgx, pg3);

    cpgsci(clrplt[1]);
    cpgsls(1);
    cpgline(I, pgx, pgy);
    cpgsls(1);
  }

  pgx[0] = xmin;
  pgx[1] = xmax;
  pgy[0] = 0.0;
  pgy[1] = 0.0;
  cpgsci(1);
  cpgsls(4);
  cpgline(2, pgx, pgy);
  cpgsls(1);
  oe_plot_label(xmin, xmax, ymin, ymax,
             " RSS (geocentric coordinate)", " x", " y", " z", clrplt);
#endif /* __ORBIT_ERROR__ */

/*
-------------------
*/

/*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*/
#ifndef __ORBD_ANIMATION__
  cpgsch(0.75);
  cpgsci(1);
/********YYYY
  cpgsvp(0.58, 0.97, 0.35, 0.60);
********/
  cpgsvp(0.58, 0.97, 0.07, 0.32);
  xmin = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               TimUTC[3], TimUTC[4], TimUTC[5], UT1_UTC)
                         - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,        0,        0,        UT1_UTC));
  xmax = 86400.0 * (float)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               TimUTC[3], TimUTC[4], TimUTC[5]+nobs, UT1_UTC)
                         - MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                               0,        0,        0,        UT1_UTC));
  if (srt[0].ODDA != 0.0) {
    ymin = -1.2 * srt[0].ODDA * 1.0e2;
    ymax =  1.9 * srt[0].ODDA * 1.0e2;
  } else {
    ymin = -1.0;
    ymax =  1.0;
  }
  cpgswin(xmin, xmax, ymin, ymax);
/********
  cpgtbox("BCSTZH", 0, 0, "BCNTS", 0, 0);
cpglab("", "OD displacement [cm]", "OD displacement");
********/
  cpgtbox("BCNSTZH", 0, 0, "BCNTS", 0, 0);
  cpglab("Time (UT)", "OD displacement [cm]", "OD Error");
  cpgsch(0.75);

/*xxxxxxxx*/
  iant = 0;
  spacecraft_position(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                      P, E, Oe, V_xyz, init_l+iant);
/*xxxxxxxx*/

  for (iant=0; iant<SRT_NUM; iant++) {
    I = 0;
    init_l[iant] = dpi;
    for (i=0; i<NOBS; i++) {
      iobs = i * TIME_STEP;
      timUTC[5] = TimUTC[5] + iobs;

#ifdef __PROPAGATION__
      orbit_propagation(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                        P, E, Oe, V_xyz, init_l+iant, RK_stage, va, vb, vc);
#else
      spacecraft_position(srt[iant], timUTC, UT1_UTC, (double)TIME_STEP,
                          P, E, Oe, V_xyz, init_l+iant);
#endif

      pgx[I] = 86400.0 * (float)(MJD(timUTC[0], timUTC[1], timUTC[2],
                                     timUTC[3], timUTC[4], timUTC[5], UT1_UTC)
                               - MJD(timUTC[0], timUTC[1], timUTC[2],
                                     0,        0,        0,        UT1_UTC));
      pg1[I] = (float)Oe[0] * 1.0e2;
      pg2[I] = (float)Oe[1] * 1.0e2;
      pg3[I] = (float)Oe[2] * 1.0e2;
      pgy[I] = (float)vlen3(Oe) * 1.0e2;
      I++;
    }
    cpgsci(clrplt[8]);
    cpgsls(2);
    cpgline(I, pgx, pg1);

    cpgsci(clrplt[3]);
/*xxxx*/
    cpgsci(clrplt[2]);
/*xxxx*/
    cpgsls(3);
    cpgline(I, pgx, pg2);

    cpgsci(clrplt[7]);
/*xxxx*/
    cpgsci(clrplt[4]);
/*xxxx*/
    cpgsls(4);
    cpgline(I, pgx, pg3);

    cpgsci(clrplt[1]);
    cpgsls(1);
    cpgline(I, pgx, pgy);
    cpgsls(1);
  }

  pgx[0] = xmin;
  pgx[1] = xmax;
  pgy[0] = 0.0;
  pgy[1] = 0.0;
  cpgsci(1);
  cpgsls(4);
  cpgline(2, pgx, pgy);
  cpgsls(1);
  oe_plot_label(xmin, xmax, ymin, ymax,
             " RSS (co-moving coordinate)",
             " Radial", " Along_Track", " Cross_Track", clrplt);
#endif /* __ORBD_ANIMATION__ */

/*
------------------------------------------
*/

#ifdef __DESCRIPTION__
  iant = 0;
  for (iant=0; iant<SRT_NUM; iant++) {
    cpgsvp(0.51, 0.97, 0.690-0.310*(float)iant, 0.990-0.310*(float)iant);
    cpgswin(0.0, 1.0, 0.0, 1.0);
    cpgbox("BC", 0.0, 0, "BC", 0.0, 0);
    cpgsfs(2);
    cpgrect(0.01, 0.99, 0.01, 0.99);

    sprintf(string, "Observation start time\0");
    cpgtext(0.03, 0.92, string);
    sprintf(string, ": %4d/%2d/%2d %2d:%2d:%2d\0",
          TimUTC[0], TimUTC[1], TimUTC[2], TimUTC[3], TimUTC[4], TimUTC[5]);
    cpgtext(0.50, 0.92, string);

    sprintf(string, "Observation duration\0");
    cpgtext(0.03, 0.84, string);
    sprintf(string, ": %4.1f hours\0", (float)nobs/3600.0);
    cpgtext(0.50, 0.84, string);

    sprintf(string, "Target (RA)\0");
    cpgtext(0.03, 0.76, string);
    sprintf(string, ": %4.1f [deg] (J2000)\0", src.RA2k / DPI);
    cpgtext(0.40, 0.76, string);

    sprintf(string, "Target (DEC)\0");
    cpgtext(0.03, 0.68, string);
    sprintf(string, ": %4.1f [deg] (J2000)\0", src.DC2k / DPI);
    cpgtext(0.40, 0.68, string);

    sprintf(string, "Orbit of the spacecraft\0");
    cpgtext(0.03, 0.60, string);

    sprintf(string, "\\fii\\fn\0");
    cpgtext(0.13, 0.36, string);
    sprintf(string, ": %6.1f [deg]\0", srt[iant].inclination / DPI);
    cpgtext(0.40, 0.36, string);

    sprintf(string, "\\gW\0");
    cpgtext(0.13, 0.52, string);
    sprintf(string, ": %8.3lf [deg]\0", Omega[iant] / DPI);
    cpgtext(0.40, 0.52, string);

    sprintf(string, "\\gw\0");
    cpgtext(0.13, 0.44, string);
    sprintf(string, ": %8.3lf [deg]\0", omega[iant] / DPI);
    cpgtext(0.40, 0.44, string);

    sprintf(string, "\\fie\\fn\0");
    cpgtext(0.13, 0.28, string);
    sprintf(string, ": %9.7lf\0", srt[iant].e);
    cpgtext(0.40, 0.28, string);

    sprintf(string, "apogee\0");
    cpgtext(0.13, 0.20, string);
    sprintf(string, ": %7.1lf [km]\0", (float)srt[iant].apogee * 1.0e-3);
    cpgtext(0.40, 0.20, string);

    sprintf(string, "M\0");
    cpgtext(0.13, 0.12, string);
    sprintf(string, ": %8.3lf [deg]\0", (float)mean_anomaly[iant] / DPI);
    cpgtext(0.40, 0.12, string);
  }
#endif /* __DESCRIPTION__ */

/*
------------------------------------------
*/

  free (pgx);
  free (pgy);
  free (pg1);
  free (pg2);
  free (pg3);

#ifdef __PROPAGATION__
  RK_end(va, vb, vc);
#endif

  return;
}


void   oe_plot_label(float xmin, float xmax, float ymin, float ymax,
                     char *lab0, char *lab1, char *lab2, char *lab3,
                     int *clrplt)
{
  int    i;
  float  pgx[2], pgy[2];

/*
--------
*/

  pgx[0] = 0.98 * xmin + 0.02 * xmax;
  pgx[1] = 0.87 * xmin + 0.13 * xmax;
  pgy[0] = 0.08 * ymin + 0.92 * ymax;
  pgy[1] = pgy[0];
  cpgsci(clrplt[1]);
  cpgsls(1);
  cpgline(2, pgx, pgy);
  pgy[1] = 0.10 * ymin + 0.90 * ymax;
  cpgsci(clrplt[1]);
  cpgsls(1);
  cpgtext(pgx[1], pgy[1], lab0);

  pgy[0] = 0.14 * ymin + 0.86 * ymax;
  pgy[1] = pgy[0];
  cpgsci(clrplt[8]);
  cpgsls(2);
  cpgline(2, pgx, pgy);
  pgy[1] = 0.16 * ymin + 0.84 * ymax;
  cpgsci(clrplt[1]);
  cpgsls(1);
  cpgtext(pgx[1], pgy[1], lab1);

  pgy[0] = 0.20 * ymin + 0.80 * ymax;
  pgy[1] = pgy[0];
  cpgsci(clrplt[3]);
/*xxxx*/
  cpgsci(clrplt[2]);
/*xxxx*/
  cpgsls(3);
  cpgline(2, pgx, pgy);
  pgy[1] = 0.22 * ymin + 0.78 * ymax;
  cpgsci(clrplt[1]);
  cpgsls(1);
  cpgtext(pgx[1], pgy[1], lab2);

  pgy[0] = 0.26 * ymin + 0.74 * ymax;
  pgy[1] = pgy[0];
  cpgsci(clrplt[7]);
/*xxxx*/
  cpgsci(clrplt[4]);
/*xxxx*/
  cpgsls(4);
  cpgline(2, pgx, pgy);
  pgy[1] = 0.28 * ymin + 0.72 * ymax;
  cpgsci(clrplt[1]);
  cpgsls(1);
  cpgtext(pgx[1], pgy[1], lab3);

  return;
}
