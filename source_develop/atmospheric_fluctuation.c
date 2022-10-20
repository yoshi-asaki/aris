#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdbool.h>
#include <cpgplot.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


/**
#define __DEBUG1__
#define __DEBUG2__
#define __DEBUG3__
#define __DEBUG4__
#define __DEBUG5__
**/


int    atmospheric_fluctuation(
         int NSRC,
         int NMTRX, int GRT_NUM, int nobs, double *ds[],
         int *TimUT, double UT1_UTC,
         struct antenna_parameter      *ant_prm,
         struct phase_screen_parameter *atm,
         struct source_parameter       *src,
         double OBS_T, double OBS_p, _Bool ATM_SWT,
         int    nseries,     int   SITE_NUM,   int NELEM,
         _Bool  AZEL_FIX,
         struct comment_param *cmnt, char   comment[][NCOMLEN],
         _Bool  TV_SWT)
{
  int    i, j, pscale, I, J, iobs, iobs_p1, ns, ne, iant;
  int    nc, NC, NS;
  int    PROC_SWT;
  int    ndat;
  int    NSP, NSCREEN, screen_limit;
  int    *IX, *IY, *SRCFLG, *OB;
  int    *IX_DUM=NULL, *IY_DUM=NULL;
  _Bool  CONT_SWT_S;
  double *DS;
  int    **NX, **NY;
  int    **screen;
  int    NY_MIN, NY_MAX;
  int    NOFFX, NOFFY;
  int    SOBS;
  int    *noffx, *sobs, *eobs;
  int    timUT[6];
  double **az, **el;
  double *seed_dist;
  double dAZdt, dELdt;
  double rho, theta_el, dz;
  int    corner_position[2][2];
  double soffx, soffy;
  int    COUNT_NOD, y_len;
  int    *npsnum, npsmax, NSPNEW, NSPOLD;
  int    Y_SHIFT=0;
  _Bool  Y_SHIFT_SWT=false;
  double z, Z, vtmp[3];
  char   string[100];

/*
===================================================================
*/

  corner_position[0][0] = 0;
  corner_position[0][1] = 0;
  corner_position[1][0] = 0;
  corner_position[1][1] = 0;

  for (i=0; i<6; i++) {
    timUT[i] = TimUT[i];
  }

  pscale = 0;
  for (iant=0; iant<GRT_NUM; iant++) {
    spherical_geometry(atm[iant].H_d, ant_prm[iant].ELLIM, &z,
                       &rho, &theta_el, &Z, &dz);
    j = (int)((2.0 * rho + (double)nobs * sqrt(pow(atm[iant].v[0], 2.0)
      + pow(atm[iant].v[1], 2.0))) / atm[iant].pixel) / NMTRX + 1;
    if (j > pscale) {
      pscale = j;
    }
  }
  pscale++;

  if ((seed_dist = (double *)calloc(NMTRX+1, sizeof(double))) == NULL) {
    printf("calloc fail in allocating memories for seed_dist.\n");
    exit (-1);
  }

/*
==============================================================
*/

  if ((az = (double **)calloc(nseries, sizeof(double *))) == NULL) {
    printf("calloc fail in allocating a pointer to az.\n");;
    exit (-1);
  }
  if ((el = (double **)calloc(nseries, sizeof(double *))) == NULL) {
    printf("calloc fail in allocating a pointer to el.\n");;
    exit (-1);
  }
  if ((NX = (int    **)calloc(nseries, sizeof(int *))) == NULL) {
    printf("calloc fail in allocating a pointer to el.\n");;
    exit (-1);
  }
  if ((NY = (int    **)calloc(nseries, sizeof(int *))) == NULL) {
    printf("calloc fail in allocating a pointer to el.\n");;
    exit (-1);
  }
  if ((screen = (int **)calloc(nseries, sizeof(int *))) == NULL) {
    printf("calloc fail in allocating a pointer to el.\n");;
    exit (-1);
  }
  for (ns=0; ns<nseries; ns++) {
    if ((az[ns]     = (double *)calloc(nobs,  sizeof(double))) == NULL) {
      printf("calloc fail in az: %d.\n", ns);
      exit (-1);
    }
    if ((el[ns]     = (double *)calloc(nobs,  sizeof(double))) == NULL) {
      printf("calloc fail in el: %d.\n", ns);
      exit (-1);
    }
    if ((NX[ns]     = (int    *)calloc(nobs,  sizeof(int))) == NULL) {
      printf("calloc fail in NX: %d.\n", ns);
      exit (-1);
    }
    if ((NY[ns]     = (int    *)calloc(nobs,  sizeof(int))) == NULL) {
      printf("calloc fail in NY: %d.\n", ns);
      exit (-1);
    }
    if ((screen[ns] = (int    *)calloc(nobs,  sizeof(int))) == NULL) {
      printf("calloc fail in screen: %d.\n", ns);
      exit (-1);
    }
  }

  if ((noffx = (int *)calloc(nseries, sizeof(int ))) == NULL) {
    printf("calloc fail in allocating a pointer to noffx.\n");;
    exit (-1);
  }
  if ((sobs  = (int *)calloc(nseries, sizeof(int ))) == NULL) {
    printf("calloc fail in allocating a pointer to sobs.\n");;
    exit (-1);
  }
  if ((eobs  = (int *)calloc(nseries, sizeof(int ))) == NULL) {
    printf("calloc fail in allocating a pointer to eobs.\n");;
    exit (-1);
  }

/*
--------------------------------------------
*/

  for (iant=0; iant<SITE_NUM; iant++) {

/*
--------------------------------------------
*/

    for (iobs=0; iobs<nobs; iobs++) {
      timUT[5] = TimUT[5] + iobs;

      for (ns=0; ns<NSRC; ns++) {
        for (ne=0; ne<NELEM; ne++) {

          i = ne * NSRC + ns;

          if (AZEL_FIX == false) {
            azel_position(timUT, UT1_UTC,
                ant_prm[iant].LLH[0],
                ant_prm[iant].LLH[1],
                ant_prm[iant].LLH[2],
                OBS_T, OBS_p, ATM_SWT, src[ns].RA, src[ns].DC,
                &az[i][iobs], &el[i][iobs], &dELdt, &dAZdt, 0.0, vtmp);
          } else if (AZEL_FIX == true) {
            az[i][iobs] =      src[ns].RA2k;
            el[i][iobs] = fabs(src[ns].DC2k);
          }

          position_on_screen((double)iobs,
                  ant_prm[ne].OFS, atm[iant], az[i][iobs], el[i][iobs],
                  &soffx, &soffy);

          NX[i][iobs] = soffx / atm[iant].pixel;
          NY[i][iobs] = soffy / atm[iant].pixel;
          screen[i][iobs] = -1;
        }
      }
    }

/*
--------------------------------------------
*/

    for (ns=0; ns<nseries; ns++) {
      sobs[ns]  = nobs - 1;
      noffx[ns] = NX[ns][sobs[ns]];
      for (iobs=0; iobs<nobs; iobs++) {
        if (el[ns][iobs] >= ant_prm[iant].ELLIM) {
          sobs[ns]  = iobs;
          noffx[ns] = NX[ns][iobs];
          break;
        }
      }
      eobs[ns] = sobs[ns];
      for (iobs=sobs[ns]+1; iobs<nobs; iobs++) {
        if (el[ns][iobs] >= ant_prm[iant].ELLIM) {
          eobs[ns] = iobs + 1;
        }
      }
    }

    NOFFX = noffx[0];
    for (ns=1; ns<nseries; ns++) {
      if (noffx[ns] < NOFFX) {
        NOFFX = noffx[ns];
      }
    }
/****
    NOFFX -= NMTRX / 2;
****/

    for (ns=0; ns<nseries; ns++) {
      for (iobs=0; iobs<nobs; iobs++) {
        NX[ns][iobs] -= NOFFX;
      }
    }

    SOBS = sobs[0];
    NS   = 0;
    for (ns=1; ns<nseries; ns++) {
      if (sobs[ns] < SOBS) {
        SOBS = sobs[ns];
        NS   = ns;
      }
    }

    for (ns=0; ns<nseries; ns++) {
      for (iobs=sobs[ns]; iobs<eobs[ns]; iobs++) {
        if (NX[ns][iobs] >= NMTRX) {
          nc = 0;
          while (NX[ns][iobs] >= NMTRX) {
            NX[ns][iobs] -= NMTRX;
            nc++;
          }
          NC = nc * NMTRX;
          for (j=iobs+1; j<eobs[ns]; j++) {
            NX[ns][j] -= NC;
          }

          if (iobs == sobs[ns]) {
            screen[ns][iobs] = nc;
          } else {
            screen[ns][iobs] = screen[ns][iobs-1] + nc;
          }
        } else {
          if (iobs == sobs[ns]) {
            screen[ns][iobs] = 0;
          } else {
            screen[ns][iobs] = screen[ns][iobs-1];
          }
        }
      }
    }

    NY_MIN = NY[NS][sobs[NS]];
    NY_MAX = NY[NS][sobs[NS]];
    for (ns=0; ns<nseries; ns++) {
      for (iobs=sobs[ns]; iobs<eobs[ns]; iobs++) {
        if (NY[ns][iobs] > NY_MAX) {
          NY_MAX = NY[ns][iobs];
        }
        if (NY[ns][iobs] < NY_MIN) {
          NY_MIN = NY[ns][iobs];
        }
      }
    }

#ifdef __DEBUG1__
    printf("__DEBUG1__  %d  %d  %d\n", NY_MIN, NY_MAX, NMTRX); fflush(stdout);
#endif

    if (NY_MAX - NY_MIN <= NMTRX) {
      Y_SHIFT_SWT = false;

      NOFFY = (int)lrint(0.5 * (float)(NY_MAX + NY_MIN));
      NOFFY -= NMTRX / 2;

      for (ns=0; ns<nseries; ns++) {
        for (iobs=sobs[ns]; iobs<eobs[ns]; iobs++) {
          NY[ns][iobs] -= NOFFY;
        }
      }
    } else {
      Y_SHIFT_SWT = true;

      NOFFY = 0;
      nc = 0;
      for (ns=0; ns<nseries; ns++) {
        for (iobs=sobs[ns]; iobs<eobs[ns]; iobs++) {
          if (screen[ns][iobs] == 0) {
            NOFFY += NY[ns][iobs];
            nc++;
          }
        }
      }
      NOFFY /= nc;
      NOFFY -= NMTRX / 2;

      for (ns=0; ns<nseries; ns++) {
        for (iobs=sobs[ns]; iobs<eobs[ns]; iobs++) {
          NY[ns][iobs] -= NOFFY;
        }
      }
    }

#ifdef __DEBUG2__
    for (iobs=0; iobs<nobs; iobs++) {
      for (ns=0; ns<nseries; ns++) {
        printf("__DEBUG2__ :  %4d  %4d  %4d:  %4d  %4d  %4d  %4d  %4d\n",
              iant, iobs, ns, screen[ns][iobs],
              NX[ns][iobs], NY[ns][iobs], sobs[ns], eobs[ns]);
        fflush (stdout);
      }
/**
      if (iobs%10 == 9) {
        getchar();
      }
**/
    }
#endif /*__DEBUG2__*/

/*
---------------------------------------------
*/

    sprintf(string, "Frozen Flow  SITE :%10s  ", ant_prm[iant].IDC);
    if (TV_SWT == false) {
      printf("%s", string);
      fflush(stdout);
    } else if (TV_SWT == true) {
      comment_disp(cmnt, comment, string, true);
    }

/*
---------------------------------------------
*/

    screen_limit = -1;
    for (iobs=0; iobs<nobs; iobs++) {
      for (ns=0; ns<nseries; ns++) {
        if (screen[ns][iobs] > screen_limit) {
          screen_limit = screen[ns][iobs];
        }
      }
    }
    if (screen_limit == -1) {
      printf("CAUTION [%s]: ", ant_prm[iant].IDC);
      printf("elevations are below the limit (%4.1f).\n",
              ant_prm[iant].ELLIM / dpi * 180.0);
      printf("CAUTION [%s]: this station is not used at all.\n",
              ant_prm[iant].IDC);
    } else {
      if ((npsnum = (int  *)calloc(screen_limit+1,  sizeof(int))) == NULL) {
        printf("calloc fail in npsnum.\n");
        exit (-1);
      }
      for (iobs=0; iobs<nobs; iobs++) {
        for (ns=0; ns<nseries; ns++) {
          if ((i=screen[ns][iobs]) >= 0) {
            npsnum[i]++;
          }
        }
      }
      npsmax = 0;
      for (i=0; i<=screen_limit; i++) {
        if (npsnum[i] > npsmax) {
          npsmax = npsnum[i];
        }
      }

      if ((IX = (int    *)calloc(npsmax,  sizeof(int))) == NULL) {
        printf("ERROR: in ATMOSPHERIC_FLUCTUATION.\n");
        printf("ERROR: fain in allocating memories for IX.\n");
        exit (-1);
      }
      if ((IY = (int    *)calloc(npsmax,  sizeof(int))) == NULL) {
        printf("ERROR: in ATMOSPHERIC_FLUCTUATION.\n");
        printf("ERROR: fain in allocating memories for IY.\n");
        exit (-1);
      }
      if ((OB = (int    *)calloc(npsmax,  sizeof(int))) == NULL) {
        printf("ERROR: in ATMOSPHERIC_FLUCTUATION.\n");
        printf("ERROR: fain in allocating memories for OB.\n");
        exit (-1);
      }
      if ((SRCFLG = (int    *)calloc(npsmax,  sizeof(int))) == NULL) {
        printf("ERROR: in ATMOSPHERIC_FLUCTUATION.\n");
        printf("ERROR: fain in allocating memories for SRCFLG.\n");
        exit (-1);
      }

/*
------------------------------
*/

      NSP        = 0;
      NSCREEN    = 0;
      CONT_SWT_S = false;
      while (NSCREEN <= screen_limit) {

        ndat = 0;
        for (iobs=0; iobs<nobs; iobs++) {
          for (ns=0; ns<nseries; ns++) {
            if (screen[ns][iobs] == NSCREEN) {
              IX[ndat]     = NX[ns][iobs];
              IY[ndat]     = NY[ns][iobs];
              OB[ndat]     = iobs;
              SRCFLG[ndat] = ns;
              ndat++;
            }
          }
        }

        corner_position[0][0] = 0;
        corner_position[1][0] = NMTRX - 1;

        corner_position[0][1] = NMTRX - 1;
        corner_position[1][1] = 0;
        for (i=0; i<ndat; i++) {
          if (IY[i] < corner_position[0][1]) {
            corner_position[0][1] = IY[i];
          }
          if (IY[i] > corner_position[1][1]) {
            corner_position[1][1] = IY[i];
          }
        }

#ifdef __DEBUG3__
        printf("__DEBUG3__  %d %d %d %d\n",
                 corner_position[0][0], corner_position[0][1],
                 corner_position[1][0], corner_position[1][1]);
        fflush(stdout);
#endif /* __DEBUG3__ */

        y_len = corner_position[1][1] - corner_position[0][1] + 1;
        if ((DS=(double *)calloc(NMTRX * y_len,  sizeof(double)) ) == NULL) {
          printf("ERROR: in ATMOSPHERIC_FLUCTUATION.\n");
          printf("ERROR: fail in allocating memories for DS. ");
          printf("(y_len=%d)\n", y_len);
          exit (-1);
        }

        if (Y_SHIFT_SWT == true) {
          if (corner_position[0][1] >=  (7 * NMTRX) / 10 &&
              corner_position[0][1] <=   NMTRX - 1       &&
              corner_position[1][1] >=  (7 * NMTRX) / 10 &&
              corner_position[1][1] <=   NMTRX - 1       &&
              Y_SHIFT >= 0) {
            Y_SHIFT =  NMTRX / 2;
          } else if (corner_position[0][1] <= (3 * NMTRX) / 10 &&
                     corner_position[0][1] >=                0 &&
                     corner_position[1][1] <= (3 * NMTRX) / 10 &&
                     corner_position[1][1] >=                0 &&
                     Y_SHIFT <= 0) {
            Y_SHIFT = -NMTRX / 2;
          } else if (corner_position[0][1] < 0         ||
                     corner_position[0][1] > NMTRX - 1 ||
                     corner_position[1][1] < 0         ||
                     corner_position[1][1] > NMTRX - 1) {
            printf("ERROR: in ATMOSPHERIC_TURBULENCE.\n");
            printf("ERROR: screen number : %d\n", NSCREEN);
            printf("ERROR: Y_SHIFT       : %d\n", Y_SHIFT);
            printf("ERROR: corner position (Y) : (%d, %d)\n",
                                      corner_position[0][1],
                                      corner_position[1][1]);
            printf("ERROR: phase screen cannot cover Y range.\n");
            exit (-1);
          } else {
            Y_SHIFT = 0;
          }
        } else {
          Y_SHIFT = 0;
        }

#ifdef __DEBUG4__
        printf("__DEBUG4__   %d  %d  %d  %d  %d\n",
                NSCREEN, corner_position[0][1], corner_position[1][1],
                y_len, Y_SHIFT);
        fflush(stdout);
#endif /* __DEBUG4__ */

        if ((COUNT_NOD = turbulent_phase_screen
                (NMTRX, 0, seed_dist, &atm[iant],
                 CONT_SWT_S,
                 0, IX_DUM, IY_DUM, DS, corner_position, Y_SHIFT)) == -1) {
          printf("ERROR: atmospheric_fluctuation: turbulent_phase_screen.\n");
          exit (-1);
        } else if (COUNT_NOD != NMTRX * y_len) {
          printf(
           "ERROR: atmospheric_fluctuation: wrong number of COUNT_NOD %d\n",
           COUNT_NOD);
        } else {
          CONT_SWT_S = true;
          for (i=0; i<ndat; i++) {
            for (ns=0; ns<NSRC; ns++) {
              for (ne=0; ne<NELEM; ne++) {
                j = ne * NSRC + ns;
                if  (SRCFLG[i] == j) {
                  I = (IX[i] - corner_position[0][0]) * y_len
                        + IY[i] - corner_position[0][1];
                  *(ds[ns] + (iant*NELEM + ne) * nobs + OB[i]) = DS[I];
                }
              }
            }
          }
        }
        free (DS);

#ifdef __DEBUG5__
        printf("__DEBUG5__   %d  %d  %d  %d  %d\n",
                NSCREEN, corner_position[0][1], corner_position[1][1],
                y_len, Y_SHIFT); fflush(stdout);
        fflush(stdout);
#endif /* __DEBUG5__ */

        if (Y_SHIFT != 0) {
          for (iobs=0; iobs<nobs; iobs++) {
            for (ns=0; ns<nseries; ns++) {
              if (screen[ns][iobs] > NSCREEN) {
                NY[ns][iobs] -= Y_SHIFT;
              }
            }
          }
        }

/*
---------------------------------------------------------
*/

        if (NSCREEN % 10 == 0) {
          sprintf(string, "%2d", (NSCREEN/10) * 10);
        } else if (NSCREEN % 5 == 0) {
          sprintf(string, "-");
        }
        if (NSCREEN % 5 == 0) {
          if (TV_SWT == false) {
            printf("%s", string);
            fflush(stdout);
          } else if (TV_SWT == true) {
            comment_disp(cmnt, comment, string, false);
          }
        }

        NSCREEN++;

      }

/*
---------------------------
*/

      if (TV_SWT == false) {
        printf("\n");
        fflush(stdout);
      }

      free (npsnum);
      free (IX);
      free (IY);
      free (OB);
      free (SRCFLG);
    }
  }

  free (seed_dist);
  for (ns=0; ns<nseries; ns++) {
    free (NX[ns]);
    free (NY[ns]);
    free (screen[ns]);
    free (az[ns]);
    free (el[ns]);
  }
  free (NX);
  free (NY);
  free (screen);
  free (az);
  free (el);

  free (noffx);
  free (sobs);
  free (eobs);

  return 1;
}
