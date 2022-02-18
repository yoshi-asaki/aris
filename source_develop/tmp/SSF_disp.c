#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <aris.h>


void  SSF_disp(double *wdist,   int  emax,  int  dmax,
                                   double distance_on_phase_screen_X,
                                   double distance_on_phase_screen_Y,
                                   double pix_atm,
                                   double C1, double D1, int cpgid)
{
  static int    i, j, k, l, I, J, IREF, JREF;
  static int    II, JJ;
  static int    NOFF_X, NOFF_Y;
  static char   string[100];
  static float  **dist;
  static char   **title;
  static float  r, pmax, pmin, noise, err_x, err_y, delta_x, delta_y;
  static float  pgx[2], pgy[2];
  static float  *s_x_cnt, *s_y_cnt, *s_x_w, *s_y_w;

  static int    N_SSF, *NDIS;
  static float  *rms, *bln;

/*
-------------------------------
*/

  N_SSF = (int)rint( sqrt(2.0) * (float)(emax));
  rms   = (float *)calloc(N_SSF, sizeof(float));
  bln   = (float *)calloc(N_SSF, sizeof(float));
  NDIS  = (int   *)calloc(N_SSF, sizeof(int));

  dist     = (float **)calloc(1,         sizeof(float *));
  dist[0]  = (float  *)calloc(emax*emax, sizeof(float));
  title    = (char  **)calloc(1,         sizeof(char *));
  title[0] = (char   *)calloc(10,        sizeof(char *));
  sprintf(title[0], "Phase Scnreen\0");

/*
-------------------------------
*/

  NOFF_X = (int)rint(distance_on_phase_screen_X / pix_atm);
  NOFF_Y = (int)rint(distance_on_phase_screen_Y / pix_atm);
  IREF   = emax/2 + NOFF_X;
  JREF   = emax/2 + NOFF_Y;
  r = (float)*(wdist + IREF*emax + JREF);
  for (i=0; i<dmax; i++) {
    I = emax/2 - dmax/2 + i + NOFF_X;
    for (j=0; j<dmax; j++) {
      J = emax/2 - dmax/2 + j + NOFF_Y;
      *(dist[0] + i*dmax + j) = (float)*(wdist + I*emax + J);
    }
  }
  pmax = 0.0;
  for (i=0; i<dmax; i++) {
    for (j=0; j<dmax; j++) {
      if (*(dist[0] + i*dmax + j) > pmax) {
        pmax = *(dist[0] + i*dmax + j);
      }
    }
  }
  pmin = pmax;
  for (i=0; i<dmax; i++) {
    for (j=0; j<dmax; j++) {
      if (*(dist[0] + i*dmax + j) < pmin) {
        pmin = *(dist[0] + i*dmax + j);
      }
    }
  }

  for (i=0; i<dmax; i++)
    for (j=0; j<dmax; j++)
      *(dist[0] + i*dmax + j)
          = (*(dist[0] + i*dmax + j) - pmin) / (pmax - pmin);
  brightness_disp(1, dmax, dmax/2, dmax/2, pix_atm, pix_atm,
                  1024.0*pix_atm, 1024.0*pix_atm,
                  0.0, 0.0, OFF, s_x_cnt, s_y_cnt, s_x_w, s_y_w, &dist[0],
                  "[m]", title, ON, ON, OFF, OFF, OFF, 64, "clr",
                  &pmin, &pmax, &noise, &err_x, &err_y, &delta_x, &delta_y);

/*
-------------------------------
*/

  NOFF_X = (int)rint(distance_on_phase_screen_X / pix_atm);
  NOFF_Y = (int)rint(distance_on_phase_screen_Y / pix_atm);
  IREF   = emax/2 + NOFF_X;
  JREF   = emax/2 + NOFF_Y;
  r = (float)*(wdist + IREF*emax + JREF);
  for (i=0; i<dmax; i++) {
    I = emax/2 - dmax/2 + i + NOFF_X;
    for (j=0; j<dmax; j++) {
      J = emax/2 - dmax/2 + j + NOFF_Y;
      *(dist[0] + i*dmax + j) = (float)*(wdist + I*emax + J) - r;
    }
  }
  pmax = 0.0;
  for (i=0; i<dmax; i++) {
    for (j=0; j<dmax; j++) {
      if (fabs(*(dist[0] + i*dmax + j)) > pmax) {
        pmax = fabs(*(dist[0] + i*dmax + j));
      }
    }
  }

  for (i=0; i<dmax; i++)
    for (j=0; j<dmax; j++)
      *(dist[0] + i*dmax + j) = fabs(*(dist[0] + i*dmax + j)) / pmax;

  brightness_disp(1, dmax, dmax/2, dmax/2,
                  pix_atm, pix_atm,
                  1024.0*pix_atm, 1024.0*pix_atm,
                  0.0, 0.0, OFF, s_x_cnt, s_y_cnt, s_x_w, s_y_w, &dist[0],
                  "[m]", title, ON, ON, OFF, OFF, OFF, 64, "clr",
                  &pmin, &pmax, &noise, &err_x, &err_y, &delta_x, &delta_y);

/*
==========================================================
*/

  printf("Display of SSF? : ");
  fgets(string, sizeof(string), stdin);
  if (! (string[0] == 'y' || string[0] == 'Y')) {
    return;
  }

  cpgbeg(0, "?", 1, 1);
  cpgpap(pgpap_prm, 1.0);
  for (I=0; I<5; I++) {
    if (I < 4) {
      if (I == 0) {
        II = 0; JJ = 0;
      } else if (I == 1) {
        II = 0; JJ = 1;
      } else if (I == 2) {
        II = 1; JJ = 0;
      } else if (I == 3) {
        II = 1; JJ = 1;
      }
      N_SSF = (int)rint( sqrt(2.0) * (float)(emax));
      IREF  = II * (emax -1);
      JREF  = JJ * (emax -1);
    } else {
      N_SSF = (int)rint( sqrt(2.0) * (float)(emax/2));
      IREF = emax/2;
      JREF = emax/2;
    }

    for (i=0; i<N_SSF; i++) {
      rms[i]  = 0.0;
      NDIS[i] = 0;
    }

    if (I == 0) {
      r = (float)*(wdist + IREF*emax + JREF);
      for (i=0; i<emax; i++)
        for (j=0; j<emax; j++)
          *(dist[0] + i*emax + j) = (float)*(wdist + i*emax + j) - r;
    } else {
      r = (float)*(dist[0] + IREF*emax + JREF);
      for (i=0; i<emax; i++)
        for (j=0; j<emax; j++)
          *(dist[0] + i*emax + j) -= r;
    }

    for (i=0; i<emax; i++) {
      for (j=0; j<emax; j++) {
        J = (int)rint(
              sqrt((float)((i-IREF)*(i-IREF) + (j-JREF)*(j-JREF))));
        rms[J] += pow(*(dist[0] + i*emax +j), 2.0);
        NDIS[J]++;
      }
    }
    J = 0;
    for (i=1; i<N_SSF; i++) {
      if (NDIS[i] != 0) {
        bln[J] = log10((float)i * pix_atm);
        rms[J] = log10(sqrt(rms[i] / (float)NDIS[i]));
        J++;
      }
    }

    if (I < 4) {
      if (II == 0) {
        pgx[0] = 0.05;
      } else {
        pgx[0] = 0.67;
      }
      pgx[1] = pgx[0] + 0.28;
      if (JJ == 0) {
        pgy[0] = 0.05;
      } else {
        pgy[0] = 0.67;
      }
      pgy[1] = pgy[0] + 0.28;
    } else {
      pgx[0] = 0.36;
      pgx[1] = 0.64;
      pgy[0] = 0.36;
      pgy[1] = 0.64;
    }

    cpgsci(1);
    cpgsvp(pgx[0], pgx[1], pgy[0], pgy[1]);
    cpgswin(log10(2.0), log10(2.0e4), log10(1.0e-6), log10(1.0e-2));
    cpgbox("BCNLTS", 0.0, 0, "BCNLTS", 0.0, 0);
    cpgpt(J, bln, rms, 15);
    if (I == 4) {
      cpglab("Baseline [m]", "RMS [m]", "SSF");
    }
    pgx[0] = log10(   1.0*pix_atm);
    pgx[1] = log10(D1);
    pgy[1] = log10(C1 * pow(D1, 5.0/6.0));
    pgy[0] = pgy[1] - (pgx[1] - pgx[0]) * 5.0/6.0;
    cpgsci(2);
    cpgline(2, pgx, pgy);
    pgx[0] = log10(D1);
    pgx[1] = log10((float)emax * pix_atm);
    pgy[0] = log10(C1 * pow(D1, 5.0/6.0));
    pgy[1] = pgy[0] + (pgx[1] - pgx[0]) * 2.0/6.0;
    cpgsci(2);
    cpgline(2, pgx, pgy);
  }
  cpgend();

/*
------------------------------------
*/

  cpgbeg(1, "?", 1, 1);
  N_SSF = (int)rint( sqrt(2.0) * (float)(emax));
  for (i=0; i<N_SSF; i++) {
    rms[i]  = 0.0;
    NDIS[i] = 0;
  }

  for (i=0; i<emax; i+=2) {
    for (j=0; j<emax; j+=2) {
      for (k=i; k<emax; k+=2) {
        for (l=j; l<emax; l+=2) {
          if (! (i == k && l == j)) {
            J = (int)rint(sqrt((float)((i-k)*(i-k) + (j-l)*(j-l))));
            rms[J] += pow((*(dist + i*emax +j)-*(dist + k*emax +l)), 2.0);
            NDIS[J]++;
          }
        }
      }
    }
  }

  J = 0;
  for (i=1; i<N_SSF; i++) {
    if (NDIS[i] != 0) {
      bln[J] = log10((float)i * pix_atm);
      rms[J] = log10(sqrt(rms[i] / (float)NDIS[i]));
      J++;
    }
  }

  pgx[0] = 0.10;
  pgx[1] = pgx[0] + 0.80;
  pgy[0] = 0.10;
  pgy[1] = pgy[0] + 0.80;

  cpgsci(1);
  cpgsvp(pgx[0], pgx[1], pgy[0], pgy[1]);
  cpgswin(log10(2.0), log10(2.0e4), log10(1.0e-6), log10(1.0e-2));
  cpgbox("BCNLTS", 0.0, 0, "BCNLTS", 0.0, 0);
  cpgpt(J, bln, rms, 15);
  cpglab("Baseline [m]", "RMS [m]", "SSF");
  pgx[0] = log10(   1.0*pix_atm);
  pgx[1] = log10(D1);
  pgy[1] = log10(C1 * pow(D1, 5.0/6.0));
  pgy[0] = pgy[1] - (pgx[1] - pgx[0]) * 5.0/6.0;
  cpgsci(2);
  cpgline(2, pgx, pgy);
  pgx[0] = log10(D1);
  pgx[1] = log10((float)emax * pix_atm);
  pgy[0] = log10(C1 * pow(D1, 5.0/6.0));
  pgy[1] = pgy[0] + (pgx[1] - pgx[0]) * 2.0/6.0;
  cpgsci(2);
  cpgline(2, pgx, pgy);
  cpgend();

/*
------------------------------------
*/

  free (dist[0]);
  free (dist);
  free (title[0]);
  free (title);

  free (rms);
  free (bln);
  free (NDIS);
}
