#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <aris.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define N_HOUR  25

/****
#define __DEBUG__
****/


int   TEC_fluctuation_amp(double Ci[][86400], int    *TimUTC)
{
  int    i, imonth, ihour, I, Ihem;
  FILE   *fp;
  double C[N_HOUR], sec[N_HOUR];
  double lfdum1, lfdum2;
  double CI, fm;
  double ci[12][24];
  char   string[100];

  gsl_interp_accel *acc    = gsl_interp_accel_alloc();
  const gsl_interp_type *t = gsl_interp_cspline_periodic;
  gsl_spline *spline       = gsl_spline_alloc(t, N_HOUR);

/*
============================================
*/

  CI = 0.0;
  if ((fp=fopen("aris_input/tecfluc.dat", "r")) == NULL) {
    printf("ERROR: TEC_fluctuation_amp: cannot open tecfluc.dat.\n");
    return (__NG__);
  } else {
    for (imonth=0; imonth<12; imonth++) {
      for (ihour=0; ihour<24; ihour++) {
        if (fgets(string, sizeof(string), fp) == NULL) {
          printf("ERROR: TEC_FLUCTUATION_AMP: Invalid input.\n");
          fclose (fp);
          return (__NG__);
        }
        sscanf(string, "%lf %lf %lf", &lfdum1, &lfdum2, &ci[imonth][ihour]);
        ci[imonth][ihour] *= 1.0e10;
        if (ci[imonth][ihour] > CI) {
          CI = ci[imonth][ihour];
        }
      }
      if (fgets(string, sizeof(string), fp) == NULL) {
        printf("ERROR: TEC_FLUCTUATION_AMP: Invalid input.\n");
        fclose (fp);
        return (__NG__);
      }
    }
    fclose (fp);
  }

  for (imonth=0; imonth<12; imonth++) {
    for (ihour=0; ihour<24; ihour++) {
      ci[imonth][ihour] /= CI;
    }
  }

#ifdef __DEBUG__
  if ((fp=fopen("tecfluc.dat", "w")) == NULL) {
    printf("ERROR: aris: cannot open tecfluc.dat.\n");
    return (__NG__);
  } else {
    for (imonth=0; imonth<12; imonth++) {
      for (ihour=0; ihour<24; ihour++) {
        fprintf(fp, "%2d %2d %lf\n", imonth+1, ihour+1, ci[imonth][ihour]);
      }
      fprintf(fp, "\n");
    }
    fclose (fp);
  }
#endif /*__DEBUG__*/

/*
============================================
*/

  for (Ihem=0; Ihem<2; Ihem++) {
    for (ihour=0; ihour<24; ihour++) {
      fm = (double)TimUTC[2] / 30.0;
      if (Ihem == 0) {
        I = TimUTC[1] - 1;
      } else {
        I = TimUTC[1] - 1 + 6;
        if (I > 11) {
          I -= 12;
        }
      }
      if (I == 11) {
        C[ihour] = (1.0 - fm) * ci[11][ihour] + fm * ci[  0][ihour];
      } else {
        C[ihour] = (1.0 - fm) * ci[ I][ihour] + fm * ci[I+1][ihour];
      }
      sec[ihour] = (float)ihour * 3600.0;
    }
    C[24] = C[0];
    sec[24] = 86400.0;

    gsl_spline_init(spline, sec, C, 25);
    for (i=0; i<86400; i++) {
      Ci[Ihem][i] = (double)gsl_spline_eval(spline, (double)i, acc);
    }
  }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return (__GO__);
}
