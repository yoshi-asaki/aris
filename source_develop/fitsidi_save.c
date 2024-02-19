#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <fitsio.h>
#include <astrotools.h>
#include <aris.h>

/********
#define __FITS_SAVE_DEBUG__
********/

/***********************************************************

                        CAUTION !!!!

In ARIS, the coordinate is right-hand system to treat 
antenna position. But, in AIPS, the coordinate of the 
antenna position seems left-hand syste. Since data 
created with ARIS will properly be treated in AIPS, 
the sign of the y-axis of the antenna coordinate mush 
be in the left-hand system. 

***********************************************************/

int  fitsidi_save(char *fits_fname,             int nobs,
                  int ANT_NUM,   int GRT_NUM,   int SRT_NUM,
                  int BGN_ANT_I, int END_ANT_I,
                  int BGN_ANT_J, int END_ANT_J,
                  struct  antenna_parameter *ant_prm,
                  int *TimUTC, double UT1_UTC, double inttim,
                  int nfrq, double band_width, double nu, 
                  struct EOP_data EOP,
                  double *tim,
                  struct baseline_uvw *bluvw,
                  struct fringe *frng, float *fringe_weight,
                  int AN_ANT_NUM,
                  int AN_BGN_ANT_I, int AN_END_ANT_I,
                  int AN_BGN_ANT_J, int AN_END_ANT_J,
                  struct  antenna_parameter *AN_ant_prm,
                  struct  source_parameter    src,
                  struct  srt_orbit_parameter *srt)
{
  fitsfile *fptr;
  int      i, I, J, iobs, iant, jant, ifrq;
  int      gcount, pcount;
  int      gmem, tfield;
  int      status = 0;
  char     string[200];
  char     comment[80];
  float    *fits_uvdata;
  double   pscale[10], pzero[10];
  double   JD;
  long     irows, nrows, group, firstelem, nelements;
  char     *tunit[20], *ttype[20], *tform[20];
  double   srtpos[3], srtvel[3], errxyz[3], oe_drc[3];
  double   antxyz[3];
  float    ant_pol[2], antprm, srtprm[4];
  int      *antflg, ant_no;
  int      srt_num;
  int      itmp;
  float    ftmp;
  double   lftmp;
  char     **ch_data;
  double   init_l[SRTMAX];

/*
------------------------------------------------------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (01)\n"); fflush(stdout);
#endif

/*
------------------------------------------------------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (02)\n"); fflush(stdout);
#endif

/*
------------------------------------------------------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (03)\n"); fflush(stdout);
#endif

  if (fits_create_file(&fptr, fits_fname, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(comment, " ");

/*
------------------------------------------------------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (04)\n"); fflush(stdout);
#endif

  gcount = 0;
  pcount = 9;
  gmem   = pcount + 3 * nfrq;

  JD = MJD(TimUTC[0], TimUTC[1], TimUTC[2], 0, 0, 0, 0.0) + 2400000.5;

  pscale[0] = 4.49653284296e-11;
  pscale[1] = 4.49653284296e-11;
  pscale[2] = 4.49653284296e-11;
  pscale[3] = 1.0;
  pscale[4] = 1.0;
  pscale[5] = 1.0;
  pscale[6] = 1.0;
  pscale[7] = 1.0;
  pscale[8] = 1.0;
  for (i=0; i<9; i++) {
    pzero[i] = 0.0;
  }
  pzero[3] = JD;

  for (iobs=0; iobs<nobs; iobs++) {
    for (iant=0; iant<ANT_NUM; iant++) {
      for (jant=iant+1; jant<ANT_NUM; jant++) {
        if (iant >= BGN_ANT_I && iant < END_ANT_I &&
            jant >= BGN_ANT_J && jant < END_ANT_J) {
          I = baseline_number(ANT_NUM, iant, jant) * nobs + iobs;
          if (fringe_weight[I] > 0.0) {
            gcount++;
          }
        }
      }
    }
  }

/*
------------------------------------------------------------
*/

  FILE   *fp;
  int    SOBS, EOBS;
  double phs;

  SOBS = 74;
  SOBS = 80;
  EOBS = SOBS + 60;

  printf("save data..............\n");
  fp = fopen("phase_check.dat", "w");
  for (iobs=SOBS; iobs<EOBS; iobs++) {
    for (ifrq=0; ifrq<nfrq; ifrq++) {
      fprintf(fp, " 0,  0, %d, 0.0\n", iobs-SOBS);
    }
  }

  iant = 0;
  for (jant=iant+1; jant<ANT_NUM; jant++) {
    if (jant >= BGN_ANT_J && jant < END_ANT_J) {
      for (iobs=SOBS; iobs<EOBS; iobs++) {
        I = baseline_number(ANT_NUM, iant, jant) * nobs + iobs;
        for (ifrq=0; ifrq<nfrq; ifrq++) {
          J = I * nfrq + ifrq;
          phs = atan2((float)frng[J].im, (float)frng[J].rl) / 2.0 / dpi / nu * 1.0e12;
          fprintf(fp, "%2d, %2d, %d, %lf\n", jant, iant, iobs-SOBS, phs);
        }
      }
    }
  }
  fclose(fp);
  printf("Done!\n");

/*
------------------------------------------------------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (05)\n"); fflush(stdout);
#endif

  itmp = 1;
  /*   if FALSE, i = 0. */
  if (fits_write_key(fptr, TLOGICAL, "SIMPLE", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = -32;
  if (fits_write_key(fptr, TINT, "BITPIX", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 7;
  if (fits_write_key(fptr, TINT, "NAXIS", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 0;
  if (fits_write_key(fptr, TINT, "NAXIS1", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 3;
  if (fits_write_key(fptr, TINT, "NAXIS2", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_key(fptr, TINT, "NAXIS3", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = nfrq;
  if (fits_write_key(fptr, TINT, "NAXIS4", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_key(fptr, TINT, "NAXIS5", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_key(fptr, TINT, "NAXIS6", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_key(fptr, TINT, "NAXIS7", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  /*   if FALSE, i = 0. */
  if (fits_write_key(fptr, TLOGICAL, "EXTEND", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_key(fptr, TLOGICAL, "BLOCKED", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "%s", src.name);
  if (fits_write_key(fptr, TSTRING, "OBJECT", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "ARIS    ");
  if (fits_write_key(fptr, TSTRING, "TELESCOP", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "ARIS    ");
  if (fits_write_key(fptr, TSTRING, "INSTRUME", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "ARIS    ");
  if (fits_write_key(fptr, TSTRING, "OBSERVER", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  fits_time2str(TimUTC[0], TimUTC[1], TimUTC[2], TimUTC[3],
                TimUTC[4], (double)TimUTC[5], -1, string, &status);
  if (fits_write_key(fptr, TSTRING, "DATE-OBS", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  if (fits_write_key(fptr, TSTRING, "DATE-MAP", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "BSCALE  =   %18.11E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  sprintf(string, "BZERO   =   %18.11E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "UNCALIB ");
  if (fits_write_key(fptr, TSTRING, "BUNIT", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = (double)EPOCH;
  sprintf(string, "EPOCH   =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 257;
  if (fits_write_key(fptr, TINT, "VELREF", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "OBSRA   =   %18.11E /", src.RA2k * 180.0 / dpi);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "OBSDEC  =   %18.11E /", src.DC2k * 180.0 / dpi);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
Fringe Data (Complex Number)
*/

  sprintf(string, "COMPLEX ");
  if (fits_write_key(fptr, TSTRING, "CTYPE2", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CRVAL2  =   %18.11E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CDELT2  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CRPIX2  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  sprintf(string, "CROTA2  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
Polarization Section
*/

  sprintf(string, "STOKES  ");
  if (fits_write_key(fptr, TSTRING, "CTYPE3", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = -2.0;
  sprintf(string, "CRVAL3  =   %18.11E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = -1.0;
  sprintf(string, "CDELT3  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CRPIX3  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  sprintf(string, "CROTA3  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
Frequency
*/

  sprintf(string, "FREQ    ");
  if (fits_write_key(fptr, TSTRING, "CTYPE4", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "CRVAL4  =   %18.11E /", nu);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "CDELT4  =     %16.9E /", band_width/(double)nfrq);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CRPIX4  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  sprintf(string, "CROTA4  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "IF      ");
  if (fits_write_key(fptr, TSTRING, "CTYPE5", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CRVAL5  =   %18.11E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CDELT5  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CRPIX5  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  sprintf(string, "CROTA5  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
Right Ascension
*/

  sprintf(string, "RA      ");
  if (fits_write_key(fptr, TSTRING, "CTYPE6", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "CRVAL6  =   %18.11E /", src.RA2k * 180.0 / dpi);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CDELT6  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CRPIX6  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  sprintf(string, "CROTA6  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
Declination
*/

  sprintf(string, "DEC     ");
  if (fits_write_key(fptr, TSTRING, "CTYPE7", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "CRVAL7  =   %18.11E /", src.DC2k * 180.0 / dpi);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CDELT7  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 1.0;
  sprintf(string, "CRPIX7  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  sprintf(string, "CROTA7  =     %16.9E /", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "T");
  if (fits_write_key(fptr, TLOGICAL, "GROUPS", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = gcount;
  if (fits_write_key(fptr, TINT, "GCOUNT", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = pcount;
  if (fits_write_key(fptr, TINT, "PCOUNT", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "UU---SIN");
  if (fits_write_key(fptr, TSTRING, "PTYPE1", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL1  =   %18.11E /", pscale[0]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO1  =   %18.11E /", pzero[0]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "VV---SIN");
  if (fits_write_key(fptr, TSTRING, "PTYPE2", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL2  =   %18.11E /", pscale[1]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO2  =   %18.11E /", pzero[1]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "WW---SIN");
  if (fits_write_key(fptr, TSTRING, "PTYPE3", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL3  =   %18.11E /", pscale[2]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO3  =   %18.11E /", pzero[2]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "DATE    ");
  if (fits_write_key(fptr, TSTRING, "PTYPE4", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL4  =   %18.11E /", pscale[3]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO4  =   %18.11E /", pzero[3]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "DATE    ");
  if (fits_write_key(fptr, TSTRING, "PTYPE5", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL5  =   %18.11E /", pscale[4]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO5  =   %18.11E /", pzero[4]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "BASELINE");
  if (fits_write_key(fptr, TSTRING, "PTYPE6", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL6  =   %18.11E /", pscale[5]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO6  =   %18.11E /", pzero[5]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "INTTIM  ");
  if (fits_write_key(fptr, TSTRING, "PTYPE7", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL7  =   %18.11E /", pscale[6]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO7  =   %18.11E /", pzero[6]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "GATEID  ");
  if (fits_write_key(fptr, TSTRING, "PTYPE8", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL8  =   %18.11E /", pscale[7]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO8  =   %18.11E /", pzero[7]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

  sprintf(string, "CORR-ID ");
  if (fits_write_key(fptr, TSTRING, "PTYPE9", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PSCAL9  =   %18.11E /", pscale[8]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "PZERO9  =   %18.11E /", pzero[8]);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
---- Additional Comments
*/

/*----
The next comment is needed to let AIPS identify TB sorted 
FITS files. 
----*/

  sprintf(string, "HISTORY AIPS   SORT ORDER = 'TB'");
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "              / Where T means TIME (IAT)");
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "              / Where B means BASELINE NUM");
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
----
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (06)\n"); fflush(stdout);
#endif

  gcount = 0;
  group = 1;
  nelements = gmem;

  if ((fits_uvdata = (float *)calloc(gmem, sizeof(float))) == NULL) {
    printf("FITSIDI_SAVE: ERROR: memory allocation: fits_uvdata.\n");
    free (fits_uvdata);
    return (-1);
  }

  for (iobs=0; iobs<nobs; iobs++) {
    for (iant=0; iant<ANT_NUM; iant++) {
      for (jant=iant+1; jant<ANT_NUM; jant++) {
        if (iant >= BGN_ANT_I && iant < END_ANT_I &&
            jant >= BGN_ANT_J && jant < END_ANT_J) {
          I = baseline_number(ANT_NUM, iant, jant) * nobs + iobs;
          if (fringe_weight[I] > 0.0) {


/*
-------------------
*/

            fits_uvdata[0] = (float)
              ((bluvw[I].u - pzero[0]) / speed_of_light  / pscale[0]);
            fits_uvdata[1] = (float)
              ((bluvw[I].v - pzero[1]) / speed_of_light  / pscale[1]);
            fits_uvdata[2] = (float)
              ((bluvw[I].w - pzero[2]) / speed_of_light  / pscale[2]);
            fits_uvdata[3] = (float)(tim[iobs] / 86400.0 / pscale[3]);
            fits_uvdata[4] = (float)0.0;
            fits_uvdata[5] = (float)
                  ((256.0*(double)(ant_prm[iant].NOSTA)
                        + (double)(ant_prm[jant].NOSTA) - pzero[5])
                                                            / pscale[5]);
            fits_uvdata[6] = (float)((inttim - pzero[6])    / pscale[6]);
            fits_uvdata[7] = (float)0.0;
            fits_uvdata[8] = (float)1.0;
/* AIPS COORDINATE CHANGE OPTION */
/****xxxx
            fits_uvdata[1] *= -1.0;
****/

/*
-------------------
*/

            for (ifrq=0; ifrq<nfrq; ifrq++) {
              J = I * nfrq + ifrq;
              fits_uvdata[pcount + 3*ifrq + 0]  = (float)frng[J].rl;
              fits_uvdata[pcount + 3*ifrq + 1]  = (float)frng[J].im;
              fits_uvdata[pcount + 3*ifrq + 2]  = (float)1.0;
              fits_uvdata[pcount + 3*ifrq + 2]  = (float)0.2;
            }

            firstelem = gcount * gmem + 1;
            if (fits_write_grppar_flt(fptr, group, firstelem, nelements,
                            fits_uvdata, &status)) {
              fits_report_error(stderr, status);
              free (fits_uvdata);
              return (-1);
            }
            gcount++;
          }
        }
      }
    }
  }

  free (fits_uvdata);

/*
----------- FITS EXTENSION TABLE : AIPS FQ -------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (07)\n"); fflush(stdout);
#endif

  nrows  = 1;
  tfield = 5;
  for (i=0; i<tfield; i++) {
    ttype[i] = (char *)calloc(20, sizeof(char));
    tform[i] = (char *)calloc(20, sizeof(char));
    tunit[i] = (char *)calloc(20, sizeof(char));
  }
  sprintf(ttype[0], "FRQSEL          ");
  sprintf(ttype[1], "IF FREQ         ");
  sprintf(ttype[2], "CH WIDTH        ");
  sprintf(ttype[3], "TOTAL BANDWIDTH ");
  sprintf(ttype[4], "SIDEBAND        ");

  sprintf(tform[0], "1J      ");
  sprintf(tform[1], "1D      ");
  sprintf(tform[2], "1E      ");
  sprintf(tform[3], "1E      ");
  sprintf(tform[4], "1J      ");

  sprintf(tunit[0], "        ");
  sprintf(tunit[1], "HZ      ");
  sprintf(tunit[2], "HZ      ");
  sprintf(tunit[3], "HZ      ");
  sprintf(tunit[4], "        ");
  if (fits_create_tbl(fptr, BINARY_TBL, nrows, tfield,
                      ttype, tform, tunit, "AIPS FQ ", &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_key(fptr, TINT, "EXTVER", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_key(fptr, TINT, "NO_IF", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_col(fptr, TINT, 1, 1, 1, 1, &itmp, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  if (fits_write_col(fptr, TDOUBLE, 2, 1, 1, 1, &lftmp, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  ftmp = (float)band_width / (float)nfrq;
  if (fits_write_col(fptr, TFLOAT, 3, 1, 1, 1, &ftmp, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  ftmp = (float)band_width;
  if (fits_write_col(fptr, TFLOAT, 4, 1, 1, 1, &ftmp, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_col(fptr, TINT, 5, 1, 1, 1, &itmp, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  for (i=0; i<tfield; i++) {
    free (ttype[i]);
    free (tform[i]);
    free (tunit[i]);
  }

/*
----------- FITS EXTENSION TABLE : AIPS AN -------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (08)\n"); fflush(stdout);
#endif

  antflg = (int  *)calloc(AN_ANT_NUM, sizeof(int));
  nrows = 0;
  for (iant=0; iant<AN_ANT_NUM; iant++) {
    if ((iant >= AN_BGN_ANT_I && iant < AN_END_ANT_I) ||
        (iant >= AN_BGN_ANT_J && iant < AN_END_ANT_J)) {
      antflg[iant] = iant;
      nrows++;
    } else {
      antflg[iant] = -1;
    }
  }

  tfield = 12;
  for (i=0; i<tfield; i++) {
    ttype[i] = (char *)calloc(20, sizeof(char));
    tform[i] = (char *)calloc(20, sizeof(char));
    tunit[i] = (char *)calloc(20, sizeof(char));
  }
  sprintf(ttype[ 0], "ANNAME          ");
  sprintf(ttype[ 1], "STABXYZ         ");
  sprintf(ttype[ 2], "ORBPARM         ");
  sprintf(ttype[ 3], "NOSTA           ");
  sprintf(ttype[ 4], "MNTSTA          ");
  sprintf(ttype[ 5], "STAXOF          ");
  sprintf(ttype[ 6], "POLTYA          ");
  sprintf(ttype[ 7], "POLAA           ");
  sprintf(ttype[ 8], "POLCALA         ");
  sprintf(ttype[ 9], "POLTYB          ");
  sprintf(ttype[10], "POLAB           ");
  sprintf(ttype[11], "POLCALB         ");

  sprintf(tform[ 0], "8A      ");
  sprintf(tform[ 1], "3D      ");
  sprintf(tform[ 2], "0D      ");
  sprintf(tform[ 3], "1J      ");
  sprintf(tform[ 4], "1J      ");
  sprintf(tform[ 5], "1E      ");
  sprintf(tform[ 6], "1A      ");
  sprintf(tform[ 7], "1E      ");
  sprintf(tform[ 8], "2E      ");
  sprintf(tform[ 9], "1A      ");
  sprintf(tform[10], "1E      ");
  sprintf(tform[11], "2E      ");

  sprintf(tunit[ 0], "        ");
  sprintf(tunit[ 1], "METERS  ");
  sprintf(tunit[ 2], "        ");
  sprintf(tunit[ 3], "        ");
  sprintf(tunit[ 4], "        ");
  sprintf(tunit[ 5], "METERS  ");
  sprintf(tunit[ 6], "        ");
  sprintf(tunit[ 7], "DEGREES ");
  sprintf(tunit[ 8], "        ");
  sprintf(tunit[ 9], "        ");
  sprintf(tunit[10], "DEGREES ");
  sprintf(tunit[11], "        ");
  if (fits_create_tbl(fptr, BINARY_TBL, nrows, tfield,
                      ttype, tform, tunit, "AIPS AN ", &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  itmp = 1;
  if (fits_write_key(fptr, TINT, "EXTVER", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 0.0;
  sprintf(string, "ARRAYX  =  %24.17E", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  lftmp = 0.0;
  sprintf(string, "ARRAYY  =  %24.17E", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  lftmp = 0.0;
  sprintf(string, "ARRAYZ  =  %24.17E", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  lftmp = 180.0 / dpi * GST(TimUTC, 0.0, UT1_UTC);
  lftmp -= 360.0 * floor(lftmp / 360.0);
  if (lftmp < 0.0) {
    lftmp += 360.0;
  }
  sprintf(string, "GSTIA0  =  %24.17E", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  lftmp = 360.98564497330005;
  sprintf(string, "DEGPDY  =  %24.17E", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "FREQ    =  %24.17E", nu);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  fits_time2str(TimUTC[0], TimUTC[1], TimUTC[2], TimUTC[3],
                TimUTC[4], (double)TimUTC[5], -1, string, &status);
  if (fits_write_key(fptr, TSTRING, "RDATE", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "POLARX  =  %24.17E", EOP.WX);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "POLARY  =  %24.17E", EOP.WY);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

  sprintf(string, "UT1UTC  =  %24.17E", UT1_UTC);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  lftmp = 0.0;
  sprintf(string, "DATUTC  =  %24.17E", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "UTC     ");
  if (fits_write_key(fptr, TSTRING, "TIMSYS", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "VLBA    ");
  if (fits_write_key(fptr, TSTRING, "ARRNAM", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "RIGHT   ");
  if (fits_write_key(fptr, TSTRING, "XYZHAND", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  sprintf(string, "ITRF    ");
  if (fits_write_key(fptr, TSTRING, "FRAME", string, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  itmp = 0;
  if (fits_write_key(fptr, TINT, "NUMORB", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  itmp = 0;
  if (fits_write_key(fptr, TINT, "NOPCAL", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  itmp = -1;
  if (fits_write_key(fptr, TINT, "FREQID", &itmp, comment, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }
  lftmp = 0.0;
  sprintf(string, "IATUTC  =  %24.17E", lftmp);
  if (fits_write_record(fptr, string, &status)) {
    fits_report_error(stderr, status);
    return (-1);
  }

/*
---------------------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (09)\n"); fflush(stdout);
#endif

  irows = 1;
  ch_data = (char **)malloc(sizeof(char *));
  ch_data[0] = (char *)malloc(8);

  for (iant=0; iant<AN_ANT_NUM; iant++) {
    if (antflg[iant] >= 0) {
      sprintf(ch_data[0], "%s", AN_ant_prm[iant].IDC);
      if (fits_write_col(fptr, TSTRING,  1, irows, 1, 1, ch_data, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

/* AIPS COORDINATE CHANGE OPTION */

/****xxxx
      if (fits_write_col(fptr, TDOUBLE,  2, irows, 1, 3,
                         AN_ant_prm[iant].XYZ, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }
****/

      antxyz[0] =  AN_ant_prm[iant].XYZ[0];
      antxyz[1] = -AN_ant_prm[iant].XYZ[1];
      antxyz[2] =  AN_ant_prm[iant].XYZ[2];
      if (fits_write_col(fptr, TDOUBLE,  2, irows, 1, 3, antxyz, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      ant_no = AN_ant_prm[iant].NOSTA;
      if (fits_write_col(fptr, TINT,     4, irows, 1, 1, &ant_no, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      ant_no = AN_ant_prm[iant].MNTSTA;
      if (fits_write_col(fptr, TINT,     5, irows, 1, 1, &ant_no, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      antprm = AN_ant_prm[iant].STAXOF;
      if (fits_write_col(fptr, TFLOAT,   6, irows, 1, 1, &antprm, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      sprintf(ch_data[0], "R");
      if (fits_write_col(fptr, TSTRING,  7, irows, 1, 1, ch_data, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      antprm = 0.0;
      if (fits_write_col(fptr, TFLOAT,   8, irows, 1, 1, &antprm, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      ant_pol[0] = 0.0;
      ant_pol[1] = 0.0;
      if (fits_write_col(fptr, TFLOAT,   9, irows, 1, 1, ant_pol, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      sprintf(ch_data[0], "L");
      if (fits_write_col(fptr, TSTRING, 10, irows, 1, 1, ch_data, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      antprm = 0.0;
      if (fits_write_col(fptr, TFLOAT,  11, irows, 1, 1, &antprm, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      ant_pol[0] = 0.0;
      ant_pol[1] = 0.0;
      if (fits_write_col(fptr, TFLOAT,  12, irows, 1, 1, ant_pol, &status)) {
        fits_report_error(stderr, status);
        return (-1);
      }

      irows++;
    }
  }

  free (ch_data[0]);
  free (ch_data);

  for (i=0; i<tfield; i++) {
    free (ttype[i]);
    free (tform[i]);
    free (tunit[i]);
  }
  free (antflg);

/*
----------- FITS EXTENSION TABLE : AIPS OB -------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (10)\n"); fflush(stdout);
#endif

  antflg = (int  *)calloc(ANT_NUM, sizeof(int));
  jant = 0;
  for (iant=GRT_NUM; iant<ANT_NUM; iant++) {
    if (((iant >= BGN_ANT_I && iant < END_ANT_I) ||
         (iant >= BGN_ANT_J && iant < END_ANT_J)) && iant >= GRT_NUM) {
      antflg[jant] = iant;
      jant++;
    }
  }

  if (jant == 0) {
    free (antflg);
  } else {
    srt_num = jant;
    nrows = srt_num * nobs;

    tfield =  8;
    for (i=0; i<tfield; i++) {
      ttype[i] = (char *)calloc(20, sizeof(char));
      tform[i] = (char *)calloc(20, sizeof(char));
      tunit[i] = (char *)calloc(20, sizeof(char));
    }
    sprintf(ttype[0],  "ANTENNA_NO      ");
    sprintf(ttype[1],  "SUBARRAY        ");
    sprintf(ttype[2],  "TIME            ");
    sprintf(ttype[3],  "ORBXYZ          ");
    sprintf(ttype[4],  "VELXYZ          ");
    sprintf(ttype[5],  "SUN_ANGLE       ");
    sprintf(ttype[6],  "ECLIPSE         ");
    sprintf(ttype[7],  "ORIENTATION     ");

    sprintf(tform[0], "1J      ");
    sprintf(tform[1], "1J      ");
    sprintf(tform[2], "1D      ");
    sprintf(tform[3], "3D      ");
    sprintf(tform[4], "3D      ");
    sprintf(tform[5], "3E      ");
    sprintf(tform[6], "4E      ");
    sprintf(tform[7], "1E      ");

    sprintf(tunit[0], "        ");
    sprintf(tunit[1], "        ");
    sprintf(tunit[2], "DAYS    ");
    sprintf(tunit[3], "METERS  ");
    sprintf(tunit[4], "M/SEC   ");
    sprintf(tunit[5], "DEGREES ");
    sprintf(tunit[6], "DAYS    ");
    sprintf(tunit[7], "DEGREES ");
    if (fits_create_tbl(fptr, BINARY_TBL, nrows, tfield,
                        ttype, tform, tunit, "AIPS OB ", &status)) {
      fits_report_error(stderr, status);
      return (-1);
    }

    itmp = 1;
    if (fits_write_key(fptr, TINT, "EXTVER", &itmp, comment, &status)) {
      fits_report_error(stderr, status);
      return (-1);
    }

    itmp = 2;
    if (fits_write_key(fptr, TINT, "TABREV", &itmp, comment, &status)) {
      fits_report_error(stderr, status);
      return (-1);
    }

/*
---------------------------
*/

    for (iant=0; iant<srt_num; iant++) {
      init_l[iant] = dpi;
    }
    irows = 1;
    for (iobs=0; iobs<nobs; iobs++) {

      for (iant=0; iant<srt_num; iant++) {

        itmp = ant_prm[antflg[iant]].NOSTA;
        if (fits_write_col(fptr, TINT, 1, irows, 1, 1, &itmp, &status)) {
          fits_report_error(stderr, status);
          return (-1);
        }

        itmp = 1;
        if (fits_write_col(fptr, TINT, 2, irows, 1, 1, &itmp, &status)) {
          fits_report_error(stderr, status);
          return (-1);
        }

        lftmp = (double)(TimUTC[3] * 3600 + TimUTC[4] * 60 + TimUTC[5]);
        lftmp /= (double)86400.0;
        if (fits_write_col(fptr, TDOUBLE, 3, irows, 1, 1, &lftmp, &status)) {
          fits_report_error(stderr, status);
          return (-1);
        }

        itmp = antflg[iant] - GRT_NUM;
        spacecraft_position(srt[itmp], TimUTC, UT1_UTC, (double)1,
                            srtpos, errxyz, oe_drc, srtvel, init_l+iant);
/* AIPS COORDINATE CHANGE OPTION */
/* I am not sure whether this option should be added or not. */
/* (2007.11.16)                                              */
/********xxxx
        srtpos[1] *= -1.0;
        srtvel[1] *= -1.0;
********/

        if (fits_write_col(fptr, TDOUBLE, 4, irows, 1, 3, srtpos,  &status)) {
          fits_report_error(stderr, status);
          return (-1);
        }
        if (fits_write_col(fptr, TDOUBLE, 5, irows, 1, 3, srtvel,  &status)) {
          fits_report_error(stderr, status);
          return (-1);
        }

        srtprm[0] = 0.0;
        srtprm[1] = 0.0;
        srtprm[2] = 0.0;
        if (fits_write_col(fptr, TFLOAT, 6, irows, 1, 3, srtprm,  &status)) {
          fits_report_error(stderr, status);
          return (-1);
        }

        srtprm[0] = 0.0;
        srtprm[1] = 0.0;
        srtprm[2] = 0.0;
        srtprm[3] = 0.0;
        if (fits_write_col(fptr, TFLOAT, 7, irows, 1, 4, srtprm,  &status)) {
          fits_report_error(stderr, status);
          return (-1);
        }

        ftmp = 0.0;
        if (fits_write_col(fptr, TFLOAT, 8, irows, 1, 1, &ftmp,  &status)) {
          fits_report_error(stderr, status);
          return (-1);
        }

        irows++;
      }
    }

    free (antflg);
    for (i=0; i<tfield; i++) {
      free (ttype[i]);
      free (tform[i]);
      free (tunit[i]);
    }
  }

/*
---- FITS FILE CLOSE -----------------------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (11)\n"); fflush(stdout);
#endif

  fits_close_file(fptr, &status);

/*
---- FITS I/O END --------------------------------
*/

#ifdef __FITS_SAVE_DEBUG__
  printf("#### __DEBUG__ : fitsidi_save : (12)\n"); fflush(stdout);
#endif

  return (1);

}
