#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


int  ant_pos_err_disp(int GRT_NUM,
                 struct antenna_parameter *ant_prm,
                 struct antenna_error_parameter *ant_err,
                 struct atmospheric_zenith_error *dz)
{
  int    i, j, I;
  int    iant, ns;
  char   string[100];
  float  **xyz, *path_err;
  float  ftmp, errtmp, pgmax;
  float  x[50], y[50];
  double *tim, obs_start_time_utc;
  double xyz_tmp[3];

/*
------------------------------------------
*/

  xyz = (float **)calloc(GRT_NUM, sizeof(float *));
  for (iant=0; iant<GRT_NUM; iant++) {
    xyz[iant] = (float *)calloc(3, sizeof(float));
  }
  path_err = (float  *)calloc(GRT_NUM, sizeof(float *));

/*
------------------------------------------
*/

  errtmp = 0.0;
  pgmax  = 0.0;

  for (iant=0; iant<GRT_NUM; iant++) {
    if (ant_prm[iant].UFL == true) {
      xyz_tmp[0] = ant_err[iant].ERR[0];
      xyz_tmp[1] = ant_err[iant].ERR[1];
      xyz_tmp[2] = ant_err[iant].ERR[2];
      drotate(xyz_tmp, -ant_prm[iant].LLH[0]-0.5*dpi, "z");
      drotate(xyz_tmp,  ant_prm[iant].LLH[1]-0.5*dpi, "x");
      xyz[iant][0] = xyz_tmp[0];
      xyz[iant][1] = xyz_tmp[1];
      xyz[iant][2] = xyz_tmp[2];

      ftmp = xyz[iant][0]*xyz[iant][0]
           + xyz[iant][1]*xyz[iant][1]
           + xyz[iant][2]*xyz[iant][2];

      if (ftmp > errtmp) {
        errtmp = ftmp;
      }

      if (fabs(ant_prm[iant].OFS[0]) > pgmax) {
        pgmax = fabs(ant_prm[iant].OFS[0]);
      }
      if (fabs(ant_prm[iant].OFS[1]) > pgmax) {
        pgmax = fabs(ant_prm[iant].OFS[1]);
      }
    }
  }

  errtmp = sqrt(errtmp);

  for (iant=0; iant<GRT_NUM; iant++) {
    if (ant_prm[iant].UFL == true) {
      path_err[iant] = dz[iant].trp * speed_of_light;
      if (fabs(path_err[iant]) > errtmp) {
        errtmp = path_err[iant];
      }
    }
  }

  cpgscr(0, 1.0, 1.0, 1.0);
  cpgscr(1, 0.0, 0.0, 0.0);

  for (I=0; I<4; I++) {
    if (I == 0) {
      cpgsvp(0.15, 0.45, 0.60, 0.90);
    } else if (I == 1) {
      cpgsvp(0.55, 0.85, 0.60, 0.90);
    } else if (I == 2) {
      cpgsvp(0.15, 0.45, 0.15, 0.45);
    } else if (I == 3) {
      cpgsvp(0.55, 0.85, 0.15, 0.45);
    }

    cpgswin(-1.2*pgmax, 1.2*pgmax, -1.2*pgmax, 1.2*pgmax);
    cpgsci(1);
    cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
    if (I == 0) {
      cpglab("", "y (m)", "E-W error");
    } else if (I == 1) {
      cpglab("", "", "N-S error");
    } else if (I == 2) {
      cpglab("x (m)", "y (m)", "Zenith error");
    } else if (I == 3) {
      cpglab("x (m)", "", "Delta Z");
    }

    for (iant=0; iant<GRT_NUM; iant++) {
      if (ant_prm[iant].UFL == true) {
        if (I < 3) {
          ftmp = fabs(xyz[iant][I])   / errtmp * 0.2 * pgmax;
        } else {
          ftmp = fabs(path_err[iant]) / errtmp * 0.2 * pgmax;
        }
/**
        ftmp = 0.05 * pgmax;
**/
        j = 0;
        for (i=0; i<=360; i+=30) {
          x[j] = ant_prm[iant].OFS[0] + ftmp * cos((float)i/180.0*dpi);
          y[j] = ant_prm[iant].OFS[1] + ftmp * sin((float)i/180.0*dpi);
          j++;
        }
        if (xyz[iant][I] < 0.0) {
          cpgsci(2);
        } else {
          cpgsci(4);
        }
/**
        cpgsci(1);
**/
        cpgsfs(2);
        cpgpoly(12, x, y);
        cpgsfs(1);
        cpgsci(1);
        x[0] = ant_prm[iant].OFS[0];
        y[0] = ant_prm[iant].OFS[1];
        cpgpt(1, x, y, 1);
      }
    }
  }

  cpgsvp(0.70, 1.00, 0.45, 0.75);
  cpgswin(-1.2*pgmax, 1.2*pgmax, -1.2*pgmax, 1.2*pgmax);
  ftmp = 0.2 * pgmax * 0.5;
  j = 0;
  for (i=0; i<=360; i+=30) {
    x[j] =  0.00 * pgmax + ftmp * cos((float)i/180.0*dpi);
    y[j] = -0.30 * pgmax + ftmp * sin((float)i/180.0*dpi);
    j++;
  }
  cpgsfs(2);
  cpgsci(2);
  cpgpoly(12, x, y);
  j = 0;
  for (i=0; i<=360; i+=30) {
    x[j] =  0.00 * pgmax + ftmp * cos((float)i/180.0*dpi);
    y[j] = -0.60 * pgmax + ftmp * sin((float)i/180.0*dpi);
    j++;
  }
  cpgsci(4);
  cpgpoly(12, x, y);
  cpgsci(1);
  sprintf(string, "%3.1f mm", -errtmp * 0.5e3);
  cpgptxt(0.20 * pgmax,  -0.35 * pgmax, 0.0, 0.0, string);
  sprintf(string, "%3.1f mm",  errtmp * 0.5e3);
  cpgptxt(0.20 * pgmax,  -0.65 * pgmax, 0.0, 0.0, string);


  for (iant=0; iant<GRT_NUM; iant++) {
    free (xyz[iant]);
  }
  free (xyz);
  free (path_err);

  return 1;
}
