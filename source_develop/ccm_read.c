#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <aris.h>


double residu(double , double );

int   ccm_read(char *cc_file, float *dist, int dmax,  int fmax,
               double *pix_uvl, double *pix_mas, double wave_length,
               double uv_max,   double *uv_factor)
{
  int    i, j, icnt, ncnt, ix, iy, idum;
  int    dref;
  FILE   *ccfp;
  char   string[1000];
  double delta_x, delta_y;
  double cflux, lfdum;
  double log_factor, pow_10;
  double d_pix, d_pix_tmp, *d_x, *d_y;

/*
------------------------------------
*/

  dref = dmax / 2;

/*
------------------------------------
*/

  if ((ccfp=fopen(cc_file, "r")) == NULL) {
    printf("ERROR: CCM_READ: CC FILE [%s] does not exit.\n", cc_file);
    return (__NG__);
  }

/*
--------
*/

  pixel_calc(pix_uvl, pix_mas, wave_length, uv_max, uv_factor, fmax, 1);
  d_pix = 1000.0 * *pix_mas;

/*
--------
*/

  while (1) {
    if (fgets(string, sizeof(string), ccfp) == NULL) {
      break;
    }
    if (strncmp(string, "    Comp", 8) == 0) {
      break;
    }
  }
  if (fgets(string, sizeof(string), ccfp) == NULL) {
    printf("ERROR: CCM_READ: Invalid input\n");
    fclose (ccfp);
    return (__NG__);
  }
  if (fgets(string, sizeof(string), ccfp) == NULL) {
    printf("ERROR: CCM_READ: Invalid input\n");
    fclose (ccfp);
    return (__NG__);
  }

/*
--------
*/

  ncnt = 1;
  while (1) {
    if (fgets(string, sizeof(string), ccfp) == NULL) {
      break;
    }

    if (string[0] != 12) {
      if (strncmp(string, "    Comp", 8) == 0) {
        if (fgets(string, sizeof(string), ccfp) == NULL) {
          break;
        }
        if (fgets(string, sizeof(string), ccfp) == NULL) {
          break;
        }
      } else {
        sscanf(string, "%d %lf %lf %lf %lf",
               &idum, &delta_x, &delta_y, &cflux, &lfdum);
        if (idum == ncnt) {
          ncnt++;
        }
      }
    }
  }
  fclose (ccfp);
  ncnt--;

/*
====================================
*/

  if ((ccfp=fopen(cc_file, "r")) == NULL) {
    printf("ERROR: CCM_READ: CC FILE [%s] does not exit.\n", cc_file);
    return (__NG__);
  }

/*
--------
*/

  while (1) {
    if (fgets(string, sizeof(string), ccfp) == NULL) {
      break;
    }
    if (strncmp(string, "    Comp", 8) == 0) {
      break;
    }
  }
  if (fgets(string, sizeof(string), ccfp) == NULL) {
    printf("ERROR: CCM_READ: Invalid input\n");
    fclose (ccfp);
    return (__NG__);
  }
  if (fgets(string, sizeof(string), ccfp) == NULL) {
    printf("ERROR: CCM_READ: Invalid input\n");
    fclose (ccfp);
    return (__NG__);
  }

/*
--------
*/

  if ((d_x = (double *)calloc(ncnt, sizeof(double))) == NULL) {
    printf("ERROR: CCM_READ: cannot allocate memory for d_x.\n");
    return (__NG__);
  }
  if ((d_y = (double *)calloc(ncnt, sizeof(double))) == NULL) {
    printf("ERROR: CCM_READ: cannot allocate memory for d_y.\n");
    free (d_x);
    return (__NG__);
  }

  icnt = 0;
  while (1) {
    if (fgets(string, sizeof(string), ccfp) == NULL) {
      break;
    }

    if (string[0] != 12) {
      if (strncmp(string, "    Comp", 8) == 0) {
        if (fgets(string, sizeof(string), ccfp) == NULL) {
          break;
        }
        if (fgets(string, sizeof(string), ccfp) == NULL) {
          break;
        }
      } else {
        sscanf(string, "%d %lf %lf %lf %lf",
               &idum, &delta_x, &delta_y, &cflux, &lfdum);
        if (cflux >= 0.0) {
          d_x[icnt] = delta_x;
          d_y[icnt] = delta_y;
          icnt++;
        }
        if (idum >= ncnt) {
          break;
        }
      }
    }
  }
  fclose (ccfp);
  ncnt = icnt;

/*
--------
*/

  d_pix = 1.0e9;
  if (d_x[0] != 0.0) {
    d_pix = fabs(d_x[0]);
  }
  for (i=0; i<ncnt; i++) {
    if (d_pix > fabs(d_x[i]) && d_x[i] != 0.0) {
      d_pix = fabs(d_x[i]);
    }
    if (d_pix > fabs(d_y[i]) && d_y[i] != 0.0) {
      d_pix = fabs(d_y[i]);
    }
  }

  for (i=0; i<ncnt; i++) {
    for (j=i+1; j<ncnt; j++) {
      d_pix_tmp = fabs(d_x[i] - d_x[j]);
      if (d_pix_tmp > 0.0) {
        if (d_pix_tmp < 0.75 * d_pix || d_pix < 0.0) {
          d_pix = d_pix_tmp;
        }
      }
      d_pix_tmp = fabs(d_y[i] - d_y[j]);
      if (d_pix_tmp > 0.0) {
        if (d_pix_tmp < 0.75 * d_pix || d_pix < 0.0) {
          d_pix = d_pix_tmp;
        }
      }
    }
  }

  free (d_x);
  free (d_y);

/*
--------
*/

  d_pix /= 1000.0;
  if (d_pix < *pix_mas) {
    *pix_mas = d_pix;
    pixel_calc(pix_uvl, pix_mas, wave_length, uv_max, uv_factor, fmax, -1);
  } else {

    log_factor = log10(*pix_mas);
    if (log_factor > 0.0) {
      log_factor = floor(log_factor);
    } else if (log_factor < 0.0) {
      log_factor = ceil(log_factor);
    }

    pow_10 = pow(10.0, log_factor);
    d_pix_tmp = floor(*pix_mas / pow_10);
    if (d_pix_tmp >= 7.0 && d_pix_tmp < 10.0) {
      d_pix = 10.0;
      if (residu(d_pix, d_pix_tmp) != 0.0) {
        d_pix = 5.0;
        if (residu(d_pix, d_pix_tmp) != 0.0) {
          d_pix = 2.0;
          if (residu(d_pix, d_pix_tmp) != 0.0) {
            d_pix = 1.0;
          }
        }
      }
    } else if (d_pix_tmp >= 4.0 && d_pix_tmp < 7.0) {
      d_pix = 5.0;
      if (residu(d_pix, d_pix_tmp) != 0.0) {
        d_pix = 2.0;
        if (residu(d_pix, d_pix_tmp) != 0.0) {
          d_pix = 1.0;
        }
      }
    } else if (d_pix_tmp >= 2.0 && d_pix_tmp < 4.0) {
      d_pix = 2.0;
      if (residu(d_pix, d_pix_tmp) != 0.0) {
        d_pix = 1.0;
      }
    } else {
      d_pix = 1.0;
    }
    *pix_mas = d_pix * pow_10;
    pixel_calc(pix_uvl, pix_mas, wave_length, uv_max, uv_factor, fmax, -1);
  }

/*
====================================
*/

  if ((ccfp=fopen(cc_file, "r")) == NULL) {
    printf("ERROR: CCM_READ: CC FILE [%s] does not exit.\n", cc_file);
    return (__NG__);
  }

/*
--------
*/

  while (1) {
    if (fgets(string, sizeof(string), ccfp) == NULL) {
      break;
    }
    if (strncmp(string, "    Comp", 8) == 0) {
      break;
    }
  }
  if (fgets(string, sizeof(string), ccfp) == NULL) {
    printf("ERROR: CCM_READ: Invalid input\n");
    fclose (ccfp);
    return (__NG__);
  }
  if (fgets(string, sizeof(string), ccfp) == NULL) {
    printf("ERROR: CCM_READ: Invalid input\n");
    fclose (ccfp);
    return (__NG__);
  }

/*
--------
*/

  ncnt = 0;
  while (1) {
    if (fgets(string, sizeof(string), ccfp) == NULL) {
      break;
    }

    if (string[0] != 12) {
      if (strncmp(string, "    Comp", 8) == 0) {
        if (fgets(string, sizeof(string), ccfp) == NULL) {
          break;
        }
        if (fgets(string, sizeof(string), ccfp) == NULL) {
          break;
        }
      } else {
        sscanf(string, "%d %lf %lf %lf %lf",
               &idum, &delta_x, &delta_y, &cflux, &lfdum);

        delta_x /= 1000.0;
        delta_y /= 1000.0;
        ix = (int)lrint(delta_x / *pix_mas) + dref;
        iy = (int)lrint(delta_y / *pix_mas) + dref;

        if (ix >= 0 && ix < dmax &&
            iy >= 0 && iy < dmax && cflux >= 0.0) {
          dist[dmax *  ix + iy ] += cflux;
        ncnt++;
        }
      }
    }
  }
  fclose (ccfp);

/*
--------
*/

  if (ncnt == 0) {
    printf("ERROR: CCM_READ: NO CC component exists.\n");
    return (__NG__);
  }

  return (__GO__);
}



double residu(double u, double t)
{
  return (u / t - u * lrint(u / t));
}
