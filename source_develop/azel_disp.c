#include <stdio.h>
#include <stdlib.h>
#include <cpgplot.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


void  azel_disp(int GRT_NUM,  int nobs,
                struct antenna_parameter  *ant_prm,
                struct st_observable *int_obs[],
                int *TimUT,   double elevation_limit,
                char *ascii_out)
{
  int    i, iant, iobs, ns;
  float  *pgx, *pgy;
  double DPI;
  FILE   *log_fp;
  double EL1, EL2, rho, Z, z, delta_z;

/*
----------
*/

  DPI = dpi / 180.0;

/*
----------
*/

  if (ascii_out[0] == '!') {
    if ((log_fp=fopen(ascii_out+1, "w")) == NULL) {
      printf("AZEL_DISP: ERROR.\n");
      return;
    } else {
      fprintf(log_fp, "%4d %2d %2d %2d %2d %2d %8d\n",
              TimUT[0], TimUT[1], TimUT[2], TimUT[3], TimUT[4], TimUT[5], nobs);
    }

    for (iant=0; iant<GRT_NUM; iant++) {
      for (iobs=0; iobs<nobs; iobs+=60) {
        i = iant * nobs + iobs;

        spherical_geometry(450.0e3, int_obs[0][i].el, &z,
                           &rho, &EL1, &Z, &delta_z);
        spherical_geometry(450.0e3, int_obs[1][i].el, &z,
                           &rho, &EL2, &Z, &delta_z);
        fprintf(log_fp,
                "%4d,%s,%8d,%7.2lf,%7.2lf,%7.2lf,%7.2lf,%7.2lf,%7.2lf\n",
                iant, ant_prm[iant].IDC, iobs,
                int_obs[0][i].az / DPI, int_obs[0][i].el / DPI,
                int_obs[1][i].az / DPI, int_obs[1][i].el / DPI,
                (int_obs[0][i].az - int_obs[1][i].az) / DPI,
                (int_obs[0][i].el - int_obs[1][i].el) / DPI);
/****/
        printf( "%4d,%s,%8d,%7.2lf,%7.2lf,%7.2lf,%7.2lf,%7.2lf,%7.2lf\n",
                iant, ant_prm[iant].IDC, iobs,
                int_obs[0][i].az / DPI, int_obs[0][i].el / DPI,
                int_obs[1][i].az / DPI, int_obs[1][i].el / DPI,
                (int_obs[0][i].az - int_obs[1][i].az) / DPI,
                (int_obs[0][i].el - int_obs[1][i].el) / DPI);
/****/
      }
    }
    fclose(log_fp);
    return;
  }

/*
----------
*/

  pgx = (float *)calloc(nobs, sizeof(float));
  pgy = (float *)calloc(nobs, sizeof(float));
  for (iobs=0; iobs<nobs; iobs++) {
    pgx[iobs] = (float)(TimUT[3]*3600 + TimUT[4]*60 + TimUT[5] + iobs);
  }

  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);

  cpgpap(pgpap_prm, 1.0);
  cpgsch(0.8);
  cpgsci(1);
  cpgsvp(0.2, 0.8, 0.6, 0.9);
  cpgswin(pgx[0], pgx[nobs-1], -360.0, 360.0);
  cpgtbox("BCSTNZHI", 0.0, 0, "BCNTS", 0.0, 0);
  cpglab("", "Azimuth [deg]", "Azimuth Plot");


  for (iant=0; iant<GRT_NUM; iant++) {
    for (ns=0; ns<SRC_NUM; ns++) {
      for (iobs=0; iobs<nobs; iobs++) {
        i = iant * nobs + iobs;
        pgy[iobs] = (float)(int_obs[ns][i].az / DPI);
      }
      cpgsci(ns+1);
      cpgsls(1+3*ns);
      cpgslw(2.0);
      cpgline(nobs, pgx, pgy);
    }
  }
  cpgslw(1.0);
  cpgsls(1);
  cpgsci(1);

/*
-------------------------------------------------------
*/

  cpgsvp(0.2, 0.8, 0.1, 0.4);
  cpgswin(pgx[0], pgx[nobs-1], 0.0, 90.0);
  cpgtbox("BCSTNZHI", 0.0, 0, "BCNTS", 0.0, 0);
  cpglab("time (UT)", "Elevation [deg]", "Elevation Plot");
  for (iant=0; iant<GRT_NUM; iant++) {
    for (ns=0; ns<SRC_NUM; ns++) {
      for (iobs=0; iobs<nobs; iobs++) {
        i = iant * nobs + iobs;
        pgy[iobs] = (float)(int_obs[ns][i].el / DPI);
      }
      cpgsci(ns+1);
      cpgsls(1+3*ns);
      cpgslw(2.0);
      cpgline(nobs, pgx, pgy);
    }
  }
  cpgslw(1.0);
  cpgsls(1);
  cpgsci(1);

/*
-------------------------------------------------------
*/

  pgx[1] = pgx[nobs-1];
  pgy[0] = (float)(elevation_limit / DPI);
  pgy[1] = (float)(elevation_limit / DPI);
  cpgsls(2);
  cpgline(2, pgx, pgy);
  cpgsls(1);

/*
-------------------------------------------------------
*/

  free (pgx);
  free (pgy);
  return;
}
