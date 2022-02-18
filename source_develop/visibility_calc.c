#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <aris.h>

#define NPNT    10

/****
#define __VISIBILITY_DISP__
****/

int  visibility_calc(float *mapr, float *mapi, int fmax,
                     double *pix_uvl, double *pix_mas,
                     double wave_length,
                     double uv_max, double *uv_factor,
                     struct source_parameter src,
                     struct morphology_file_name  *ch_file,
                     float xmin, float xmax, float ymin, float ymax,
                     _Bool TV_SWT, int  NS, int pgid)
{
  int    i, j, dmax, dmax2, dmaxh, I;
  char   string[100], comp[100];
  float  *dist;
  double theta;
  float  pmin[2], pmax[2], noise[2];
  float  err_x[2], err_y[2], delta_x[2], delta_y[2];
  float  s_x_cnt[2], s_y_cnt[2], s_x_w[2], s_y_w[2];
  float  t_flux, bias;
  float  *DIST[2];
  char   *TITLE[2];
  FILE   *fp;
  struct comment_param cmnt;
  char   comment[10][NCOMLEN];

  int    nxmax, nymax, dwin;
  int    nxw_min, nxw_max, nyw_min, nyw_max;
  int    nx_cen, ny_cen;
  float  x_cen, y_cen, d_width;
  int    nxwin, nywin;
  float  *pgx, *pgy;
  float  *mskr, *mski;
  float  R, X, Y;
  float  a, b;
  double T;
  int    LEN, IFREF, JFREF, IDREF, JDREF;
  mskr = calloc(fmax*fmax, sizeof(float));
  mski = calloc(fmax*fmax, sizeof(float));


/*
-----------------------------------
*/

  cmnt.ncol  = 4;
  cmnt.xmin  = xmin;
  cmnt.xmax  = xmax;
  cmnt.ymin  = 0.515;
  cmnt.ymax  = 0.595;
  cmnt.pitch = 0.03;
/*
  comment_init(&cmnt, comment, true);
*/

/*
-----------------------------------
*/

  dmax = fmax / 4;
  dmax2 = dmax * dmax;
  dmaxh = dmax / 2;

  if ((dist = (float *)calloc(dmax2, sizeof(float))) == NULL) {
    printf("ERROR: visibility_calc: calloc failure in visibility_calc.\n");
    return (__NG__);
  }

  IFREF = fmax / 2;
  JFREF = fmax / 2;
  IDREF = dmax / 2;
  JDREF = dmax / 2;

/*
-----------------------------------
*/

  if ((pgx = (float *)calloc(fmax, sizeof(float))) == NULL) {
    printf("ERROR: visibility_calc: calloc failure in visibility_calc.\n");
    free (dist);
    return (__NG__);
  }

  if ((pgy = (float *)calloc(fmax, sizeof(float))) == NULL) {
    printf("ERROR: visibility_calc: calloc failure in visibility_calc.\n");
    free (pgx);
    free (dist);
    return (__NG__);
  }

/*
-----------------------------------
*/

  if (src.morphology == SRC_POINT) {
    dist[dmaxh*dmax + dmaxh] = (float)src.flux;
    for (i=0; i<fmax; i++) {
      for (j=0; j<fmax; j++) {
        *(mapr + fmax*i + j) = (float)src.flux;
        *(mapi + fmax*i + j) = (float)0.0;
      }
    }

/*
-----------------------------------
*/

  } else {
    if        (src.morphology == SRC_DISK_JETCJET) {
      theta = (45.0 - 90.0) / 180.0 * dpi;
      jet_cjet  (dist, dmax, theta, *pix_mas, 10.0 * src.flux);
      ADAF_disk (dist, dmax, theta, *pix_mas,  1.0 * src.flux);
    } else if (src.morphology == SRC_DISK_VSOP2JET) {
      theta = (45.0 - 90.0) / 180.0 * dpi;
      VSOP2_logo(dist, dmax, theta, *pix_mas, 10.0 * src.flux);
      ADAF_disk (dist, dmax, theta, *pix_mas,  1.0 * src.flux);
    } else if (src.morphology == SRC_CC_COMP) {
      if (ccm_read(ch_file->cct, dist, dmax, fmax, pix_uvl, pix_mas,
                   wave_length, uv_max, uv_factor) == __NG__) {
        printf("ERROR: visibility_calc: ccm_read returned ERROR.\n");
        free (dist);
        free (pgx);
        free (pgy);
        return (__NG__);
      }
    } else if (src.morphology == SRC_MULTI_COMP) {
      if (mcm_read(ch_file->mcm, dist, dmax, *pix_mas) == __NG__) {
        printf("ERROR: visibility_calc: mcm_read returned ERROR.\n");
        free (dist);
        free (pgx);
        free (pgy);
        return (__NG__);
      }

    } else if (src.morphology == SRC_BHS_MOD) {
      if (bhs_model(ch_file->bhs, dist, dmax, fmax, pix_uvl, pix_mas,
                   wave_length, uv_max, uv_factor) == __NG__) {
        printf("ERROR: visibility_calc: bhs_model returned ERROR.\n");
        free (dist);
        free (pgx);
        free (pgy);
        return (__NG__);
      }
    }

    if (src.morphology != SRC_CC_COMP &&
        src.morphology != SRC_MULTI_COMP &&
        src.morphology != SRC_BHS_MOD &&
        src.morphology != SRC_POINT) {
      t_flux = 0.0;
      for (i=0; i<dmax2; i++) {
        t_flux += dist[i];
      }
      t_flux /= src.flux;
      for (i=0; i<dmax2; i++) {
        dist[i] /= t_flux;
      }
    }
    get_vis(mapr, mapi, fmax, dist, dmax);
  }

/**
  for (i=0; i<fmax*fmax; i++) {
    mskr[i] = 1.0;
  }
  for (i=0; i<fmax; i++) {
    for (j=0; j<fmax; j++) {
      R = (dpi / 4.0) * exp(-0.005*0.005/2.0/log(2.0)*pow(*pix_uvl, 2.0) 
        *(float)((i-IFREF)*(i-IFREF) + (j-JFREF)*(j-JFREF)));
      R = 1.0;
      mapr[fmax*i+j] *= R;
      mapi[fmax*i+j] *= R;
    }
  }
  wfft2d(mapr, mapi, mskr, fmax, 1);
  LEN = sizeof(float) * dmax;
  for (i=0; i<dmax; i++) {
    I = IFREF - IDREF + i;
    memcpy(dist + i*dmax, mapr + I*fmax + JFREF - JDREF, LEN);
  }
**/
  free (mskr);
  free (mski);

/*
--------
*/

  if (pgid != -1) {
    if (NS == 0) {
      sprintf(string, "Source Image (Target)");
    } else if (NS == 1) {
      sprintf(string, "Source Image (Reference)");
    }

    pmax[0] = 0.0;
    nxmax = 0;
    nymax = 0;
    nxw_min = dmax - 1;
    nyw_min = dmax - 1;
    nxw_max = 0;
    nyw_max = 0;
    for (i=0; i<dmax; i++) {
      I = dmax * i;
      for (j=0; j<dmax; j++) {
        if (dist[I] > pmax[0]) {
          pmax[0] = dist[I];
          nxmax = i;
          nymax = j;
        }
        if (dist[I] != 0.0) {
          if (i < nxw_min) {
            nxw_min = i;
          }
          if (i > nxw_max) {
            nxw_max = i;
          }
          if (j < nyw_min) {
            nyw_min = j;
          }
          if (j > nyw_max) {
            nyw_max = j;
          }
        }
        I++;
      }
    }
    delta_x[0] = *pix_mas * (nxmax - dmaxh);
    delta_y[0] = *pix_mas * (nymax - dmaxh);

    if (nxw_max - nxw_min > nyw_max - nyw_min) {
      dwin = nxw_max - nxw_min;
    } else {
      dwin = nyw_max - nyw_min;
    }
    nx_cen = (nxw_max + nxw_min) / 2;
    ny_cen = (nyw_max + nyw_min) / 2;
    x_cen = *pix_mas * (float)(nx_cen - dmaxh);
    y_cen = *pix_mas * (float)(ny_cen - dmaxh);
    d_width = *pix_mas * (float)(dwin + dmax / 16);
    d_width = *pix_mas * (float)(dwin + dmax / 2);

    t_flux = 0.0;
    bias = dist[0];
    for (i=0; i<dmax; i++) {
      I = dmax * i;
      for (j=0; j<dmax; j++) {
        if (dist[I] > t_flux) {
          t_flux = dist[I];
        }
        if (dist[I] < bias) {
          bias = dist[I];
        }
        I++;
      }
    }

    pg_color_map(dmax, nx_cen, ny_cen, *pix_mas, *pix_mas,
               bias, t_flux, d_width, d_width, 0.0, 0.0, dist,
               xmin, xmax, ymin, ymax,
               "[mas]", "[mas]", string, false, true, "clr");
/********
    pgcont_map(dmax, nx_cen, ny_cen, *pix_mas, *pix_mas,
               bias, t_flux, d_width, d_width, 0.0, 0.0, dist,
               xmin, xmax, ymin, ymax,
               "[mas]", "[mas]", string, false, true, "clr");
********/

    sprintf(string, "Maximum Peak: %7.2E\n", pmax[0]);
    if (TV_SWT == true) {
      cpgsvp(0.0, 1.0, 0.0, 1.0);
      cpgswin(0.0, 1.0, 0.0, 1.0);
      cpgsch(0.70);
      comment_init(&cmnt, comment, true);
      comment_disp(&cmnt, comment, string, true);
    } else {
      printf("%s", string);
    }
    sprintf(string,  "(X, Y) @ Peak : (%f, %f)\n", delta_x[0], delta_y[0]);
    if (TV_SWT == true) {
      comment_disp(&cmnt, comment, string, true);
    } else {
      printf("%s", string);
    }
  }

/*
----------
*/

#ifdef __VISIBILITY_DISP__
  if (pgid != -1) {
    cpgsvp(0.0, 1.0, 0.0, 1.0);
    cpgswin(0.0, 1.0, 0.0, 1.1);
    cpgsci(0);
    cpgrect(0.0, 1.0, 0.3, 1.1);
    cpgsci(1);

    string[0] = 0;
    DIST[0] = mapr;
    DIST[1] = mapi;
    TITLE[0] = string;
    TITLE[1] = string;
    s_x_cnt[0] = 0.0;
    s_y_cnt[0] = 0.0;
    s_x_w[0]   = 0.0;
    s_y_w[0]   = 0.0;
    s_x_cnt[1] = 0.0;
    s_y_cnt[1] = 0.0;
    s_x_w[1]   = 0.0;
    s_y_w[1]   = 0.0;
    cpgslct(pgid);
    brightness_disp(2, fmax, fmax/2+1, fmax/2+1,
                    (float)(*pix_uvl), (float)(*pix_uvl),
                    (float)fmax * (float)(*pix_uvl),
                    (float)fmax * (float)(*pix_uvl),
                    0.0, 0.0,
                    false, s_x_cnt, s_y_cnt, s_x_w, s_y_w, DIST,
                    "", TITLE,
                    TV_SWT, true, true, true, true, 128, "clr",
                    pmin, pmax, noise, err_x, err_y, delta_x, delta_y);
    printf("Hit RETURN : "); getchar();


    cpgsvp(0.0, 1.0, 0.0, 1.0);
    cpgswin(0.0, 1.0, 0.0, 1.1);
    cpgsci(0);
    cpgrect(0.0, 1.0, 0.3, 1.1);
    cpgsci(1);

    string[0] = 0;
    for (i=0; i<fmax; i++) {
      for (j=0; j<fmax; j++) {
        *(mapr + fmax*i + j) = sqrt(pow(*(mapr + fmax*i + j), 2.0)
                                 +  pow(*(mapi + fmax*i + j), 2.0));
      }
    }
    DIST[0] = mapr;
    TITLE[0] = string;
    s_x_cnt[0] = 0.0;
    s_y_cnt[0] = 0.0;
    s_x_w[0]   = 0.0;
    s_y_w[0]   = 0.0;
    cpgslct(pgid);
    brightness_disp(1, fmax, fmax/2+1, fmax/2+1,
                    (float)(*pix_uvl), (float)(*pix_uvl),
                    0.9 * (float)fmax * (float)(*pix_uvl),
                    0.9 * (float)fmax * (float)(*pix_uvl),
                    0.0, 0.0,
                    false, s_x_cnt, s_y_cnt, s_x_w, s_y_w, DIST,
                    "", TITLE,
                    TV_SWT, true, false, false, true, 128, "clr",
                    pmin, pmax, noise, err_x, err_y, delta_x, delta_y);
    printf("Hit RETURN : "); getchar();

    for (i=0; i<fmax; i++) {
      for (j=0; j<fmax; j++) {
        *(mapr + fmax*i + j) = sqrt(pow(*(mapr + fmax*i + j), 2.0)
                                 +  pow(*(mapi + fmax*i + j), 2.0));
      }
    }
    cpgsvp(0.0, 1.0, 0.0, 1.0);
    cpgswin(0.0, 1.0, 0.0, 1.1);
    cpgsci(0);
    cpgrect(0.0, 1.0, 0.3, 1.1);
    cpgsci(1);
    cpgsvp(0.15, 0.90, 0.38, 0.93);
    j = 0;
    for (i=fmax/2; i<fmax; i++) {
      pgx[j] = (float)(*pix_uvl) * (float)j;
      pgy[j] = *(mapr + fmax*i + fmax/2);
      j++;
    }
/****
    cpgsch(1.3);
    cpgswin(pgx[0], pgx[j-1], 0.0, 1.1 * pgy[0]);
    cpgswin(pgx[0],    5.0e7, 0.0, 1.2 * pgy[0]);
****/
    cpgswin(pgx[0], pgx[j-1], 0.0, 1.1 * pgy[0]);
    cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
    cpgline(j, pgx, pgy);

    j = 0;
    for (i=fmax/2; i<fmax; i++) {
      pgx[j] = (float)(*pix_uvl) * (float)j;
      pgy[j] = *(mapr + fmax*fmax/2 + i);
      j++;
    }
    cpgsci(2);
    cpgline(j, pgx, pgy);
    cpgsci(7);

    cpglab("\\fr(u, v) distance", "\\frVisibility amplitude [Jy]", "");
    printf("Hit RETURN : "); getchar();
  }
#endif /* __VISIBILITY_DISP__ */

/*
----------
*/

  free (dist);
  free (pgx);
  free (pgy);

  return (__GO__);
}
