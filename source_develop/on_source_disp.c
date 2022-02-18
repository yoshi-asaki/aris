#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

int  on_source_disp(int ANT_NUM, int GRT_NUM, int nobs,
                 int  *TimUTC, double UT1_UTC, int nswt,
                 struct antenna_parameter *ant_prm,
                 struct st_observable *int_obs[],
                 double elevation_limit)
{
  int    i, itmp1, itmp2;
  int    iant, ns=0;
  int    timUTC[6];
  double ut1_utc;

  float  *pgx, *pgy, *pgz;
  int    iobs, nslewt;
  double azl=0.0, ell=0.0;
  float  sens_min, sens_max, s_tmp;

/*
------------------------------------------
*/

  nslewt = 0;

  if (nswt == 0) {
    return 1;
  }

  for (i=0; i<6; i++) {
    timUTC[i] = TimUTC[i];
  }
  ut1_utc = UT1_UTC;

  pgx = (float *)calloc(nobs, sizeof(float));
  pgy = (float *)calloc(nobs, sizeof(float));
  pgz = (float *)calloc(nobs, sizeof(float));
  for (iobs=0; iobs<nobs; iobs++) {
    pgx[iobs] = (float)(timUTC[3]*3600 + timUTC[4]*60 + timUTC[5] + iobs);
  }

/*
--------
*/

  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);
  cpgpap(pgpap_prm, 1.0);
  cpgsch(0.8);
  cpgsci(1);
  cpgsvp(0.2, 0.8, 0.6, 0.9);
  cpgswin(pgx[0], pgx[nobs-1], -1.0, 0.5 * (float)nswt);
  cpgtbox("BCSTNZHI", 0.0, 0, "BCNTS", 0.0, 0);
  cpglab("time (UT)", "ON source time [sec]", "");

  sens_min = -1.0;
  sens_max = -1.0;
  for (iant=0; iant<ANT_NUM; iant++) {
    if (ant_prm[iant].WID[0] > 0 || ant_prm[iant].WID[1] > 0) {
      if (iant < GRT_NUM) {
        azl = ant_prm[iant].AZSV * ant_prm[iant].AZSV / ant_prm[iant].AZSA;
        ell = ant_prm[iant].ELSV * ant_prm[iant].ELSV / ant_prm[iant].ELSA;
      } else if (iant >= GRT_NUM) {
        nslewt = (int)lrint(ant_prm[iant].slewt);
      }

      iobs = 0;
      for (iobs=0; iobs<nobs; iobs++) {
        pgz[iobs] = 1.0;
        i = iant * nobs + iobs;
        if (iant < GRT_NUM) {
          if (int_obs[0][i].el >= elevation_limit &&
              int_obs[1][i].el >= elevation_limit) {

/********
            if (iobs % 10 == 0) {
              printf("XXXX %f   %f   %f\n", pgx[iobs],
                      fabs(int_obs[0][i].az - int_obs[1][i].az) / 3.14159265*180.0,
                      fabs(int_obs[0][i].el - int_obs[1][i].el) / 3.14159265*180.0);
            }
********/

            nslewt = (int)lrint(slew_time(int_obs[0][i].az,
                                         int_obs[1][i].az, azl,
                             ant_prm[iant].AZSV, ant_prm[iant].AZSA,
                             int_obs[0][i].el, int_obs[1][i].el, ell,
                             ant_prm[iant].ELSV, ant_prm[iant].ELSA));
            pgy[iobs] = (float)((nswt - 2 * nslewt) / 2);
            if (pgy[iobs] > 0.0) {
              s_tmp = sqrt(pgy[iobs])
                       / SEFD(ant_prm[iant].Trx[ns],
                              ant_prm[iant].Tsky[ns] / sin(int_obs[ns][i].el),
                              ant_prm[iant].Dm[ns], ant_prm[iant].Ae[ns]);
              if (s_tmp > sens_max) {
                sens_max = s_tmp;
              } else if (s_tmp < sens_min || sens_min < 0.0) {
                sens_min = s_tmp;
              }
            }
          } else {
            pgy[iobs] = 0.0;
            pgz[iobs] = -1.0;
          }
        } else {
          pgy[iobs] = (float)((nswt - 2 * nslewt) / 2);
          if (pgy[iobs] > 0.0) {
            s_tmp = sqrt(pgy[iobs]) / SEFD(ant_prm[iant].Trx[ns],
                                    ant_prm[iant].Tsky[ns],
                                    ant_prm[iant].Dm[ns], ant_prm[iant].Ae[ns]);
            if (s_tmp > sens_max) {
              sens_max = s_tmp;
            } else if (s_tmp < sens_min || sens_min < 0.0) {
              sens_min = s_tmp;
            }
          }
        }
      }
      itmp1 = 0;
      cpgsci(iant+1);
      while (itmp1 <= nobs) {
        itmp2 = 0;
        for (iobs=itmp1; iobs<nobs; iobs++) {
          if (pgz[iobs] == -1.0) {
            break;
          } else {
            itmp2++;
          }
        }
        cpgline(itmp2, pgx+itmp1, pgy+itmp1);
        itmp1 += (itmp2 + 1);
      }
      cpgsci(1);
    }
  }

  sens_min *= 1.0e-23;
  sens_max *= 1.0e-23;

/*
--------
*/


  cpgsvp(0.2, 0.8, 0.1, 0.4);
  cpgswin(pgx[0], pgx[nobs-1], log10(0.8 * sens_min), log10(1.2 * sens_max));
  cpgtbox("BCSTNZHI", 0.0, 0, "BCLNTS", 0.0, 0);
  cpglab("time (UT)", "Antenna Sensitivity",
         "sqrt(ON source time) / SEFD / 10\\u23");

  for (iant=0; iant<ANT_NUM; iant++) {
    if (ant_prm[iant].WID[0] > 0 && ant_prm[iant].WID[1] > 0) {
      if (iant < GRT_NUM) {
        azl = ant_prm[iant].AZSV * ant_prm[iant].AZSV / ant_prm[iant].AZSA;
        ell = ant_prm[iant].ELSV * ant_prm[iant].ELSV / ant_prm[iant].ELSA;
      } else if (iant >= GRT_NUM) {
        nslewt = (int)lrint(ant_prm[iant].slewt);
      }

      iobs = 0;
      for (iobs=0; iobs<nobs; iobs++) {
        pgz[iobs] = 1.0;
        i = iant * nobs + iobs;
        if (iant < GRT_NUM) {
          if (int_obs[0][i].el >= elevation_limit &&
              int_obs[1][i].el >= elevation_limit) {
            nslewt = (int)lrint(slew_time(int_obs[0][i].az,
                                         int_obs[1][i].az, azl,
                             ant_prm[iant].AZSV, ant_prm[iant].AZSA,
                             int_obs[0][i].el, int_obs[1][i].el, ell,
                             ant_prm[iant].ELSV, ant_prm[iant].ELSA));
            if (nswt - 2 * nslewt > 0) {
              pgy[iobs] = sqrt((float)((nswt - 2 * nslewt) / 2))
                         / SEFD(ant_prm[iant].Trx[ns],
                                ant_prm[iant].Tsky[ns] / sin(int_obs[ns][i].el),
                                ant_prm[iant].Dm[ns], ant_prm[iant].Ae[ns]);
            } else {
              pgy[iobs] = sens_min;
            }
            pgy[iobs] = log10(pgy[iobs] * 1.0e-23);
          } else {
            pgy[iobs] = 0.0;
            pgz[iobs] = -1.0;
          }
        } else {
          if (nswt - 2 * nslewt > 0) {
            pgy[iobs] = sqrt((float)((nswt - 2 * nslewt) / 2))
                         / SEFD(ant_prm[iant].Trx[ns],
                                ant_prm[iant].Tsky[ns],
                                ant_prm[iant].Dm[ns], ant_prm[iant].Ae[ns]);
          } else {
            pgy[iobs] = sens_min;
          }
          pgy[iobs] = log10(pgy[iobs] * 1.0e-23);
        }
      }
      itmp1 = 0;
      cpgsci(iant+1);
      while (itmp1 <= nobs) {
        itmp2 = 0;
        for (iobs=itmp1; iobs<nobs; iobs++) {
          if (pgz[iobs] == -1.0) {
            break;
          } else {
            itmp2++;
          }
        }
        cpgline(itmp2, pgx+itmp1, pgy+itmp1);
        itmp1 += (itmp2 + 1);
      }
      cpgsci(1);
    }
  }

/*
--------
*/

  free (pgx);
  free (pgy);
}
