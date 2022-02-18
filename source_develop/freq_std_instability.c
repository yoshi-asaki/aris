#include <stdio.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <aris.h>

void    asd_time_series
        (
          double    *,    int        ,
          int        ,    int        ,
          _Bool
        );


int   freq_std_instability(int nobs, double *fqs_ds, int FRQSTD)
{

  int    i, j, n_tmp, n_seg;
  double *ds, asdf;

/*
--------
*/

  n_tmp = 1;
  while (n_tmp <= nobs) {
    n_tmp *= 2;
  }

/*
--------
*/

  if ((ds = (double *)calloc(n_tmp+1, sizeof(double))) == NULL) {
    printf("#### Clock error: memory allocation failed.\n");
    return -1;
  }

/*
--------
*/

  asdf = sqrt(3.0/2.0)
           * modified_asd_model((double)n_tmp, FRQSTD) * (double)n_tmp;
  ds[    0] = asdf * random_val1();
  ds[n_tmp] = asdf * random_val1();

  n_tmp /= 2;
  n_seg = n_tmp;
  asd_time_series(ds, n_tmp, n_seg, FRQSTD, false);

/*
--------
*/

  if (random_val1() >= 0.0) {
    j = (int)((double)(2*n_tmp - nobs) * random_val0());
    for (i=0; i<nobs; i++) {
      fqs_ds[i] = ds[j++];
    }
  } else {
    j = 2*n_tmp - (int)((double)(2*n_tmp - nobs) * random_val0());
    for (i=0; i<nobs; i++) {
      fqs_ds[i] = ds[j--];
    }
  }
  free (ds);

/*
--------
*/

  return 1;
}



void  asd_time_series(double *ds, int  n_tmp, int n_seg,
                      int FRQSTD, _Bool ZERO_SWT)
{
  int    i, nstep;
  int    n0, n1;
  double asdf, fac;

/*
--------
*/

  if (ZERO_SWT == true) {
    ds[      0] = 0.0;
    ds[2*n_seg] = 0.0;
  }

  fac = sqrt(3.0/2.0);
  nstep = 1;

  while (1) {
    asdf = fac * modified_asd_model((double)n_tmp, FRQSTD) * (double)n_tmp;
    for (i=0; i<nstep; i++) {
      n0 = 2 *  i      * n_seg;
      n1 = 2 * (i + 1) * n_seg;
      ds[n0 + n_seg] = 0.5 * (ds[n0] + ds[n1]) + asdf * random_val1();
    }
    nstep *= 2;
    n_seg /= 2;
    n_tmp /= 2;
    if (n_seg == 0) {
      break;
    }
  }

  return;
}
