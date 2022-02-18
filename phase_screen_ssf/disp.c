#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cpgplot.h>

#define NDAT  3400000
#define NSMP     1000

int  main()
{
  int    i, j, k, m, n, I, NEL, ndat;
  int    irow, nrow=32, ncol=0;
  int    NSWT;
  float  pgx[NDAT], pgy[NDAT];
  FILE   *fp;
  char   string[1000];
  char   fname[1000];
  float  ftmp1, ftmp2;
  float  sx[NSMP], sy[NSMP];
  int    nsmp,   nskip=10;
  float  scale1, scale2;
  float  expon1, expon2;
  float  coeff1, coeff2, coeff3;
  float  amp0, scale0;
  float  pixel;

/*
------
*/

  scale1 =  1638.4;
  scale2 = 13107.2;
  expon1 =     1.2;
  expon2 =     0.6;

  scale0 =  100.0;
  amp0    =   0.2;
  scale0 =  200000.0;
  amp0    =   3.0;
  scale0 =  20000.0;
  amp0    =   4.0;

  if (scale0 <= scale1) {
    coeff1   = amp0   / pow(scale0, 0.5*expon1);
    coeff2   = coeff1 * pow(scale1, 0.5*(expon1 - expon2));
    coeff3   = coeff2 * pow(scale2, 0.5*expon2);
  } else if (scale0 > scale1 && scale0 <= scale2) {
    coeff2   = amp0   / pow(scale0, 0.5*expon2);
    coeff1   = coeff2 / pow(scale1, 0.5*(expon1 - expon2));
    coeff3   = coeff2 * pow(scale2, 0.5*expon2);
  } else {
    coeff3   = amp0;
    coeff2   = coeff3 / pow(scale2, 0.5*(expon2));
    coeff1   = coeff2 / pow(scale1, 0.5*(expon1 - expon2));
  }

/*
------
*/

  cpgbeg(1, "/xs", 1, 1);
  cpgpap(8.0, 1.0);
  cpgscr(0, 1, 1, 1);
  cpgscr(1, 0, 0, 0);

/*
------
*/

  NSWT = 0;
  for (irow=0; irow<nrow; irow+=200) {
    n = NDAT;
    sprintf(fname, "sample00.dat");
    if ((fp = fopen("sample00.dat", "r")) == NULL) {
      break;
    }
    while (1) {
      fgets(string, sizeof(string), fp);
      if (strncmp(string, "PIXEL SIZE       ", 17) == 0) {
        sscanf(string+18, "%f", &pixel);
      } else if (strncmp(string, "SCREEN WIDTH (NX)", 17) == 0) {
        sscanf(string+18, "%d", &NEL);
      } else if (strncmp(string, "SCREEN WIDTH (NY)", 17) == 0) {
        sscanf(string+18, "%d", &ncol);
      } else if (strncmp(string, "# END OF HEADER", 15) == 0) {
        break;
      }
    }
    for (i=0; i<NEL; i++) {
      if (i == NDAT) {
        break;
      }
      for (j=0; j<ncol; j++) {
        fgets(string, sizeof(string), fp);
        if (j == irow) {
          sscanf(string, "%f,%f,%f", &pgx[i], &ftmp1, &pgy[i]);
          pgx[i] *= pixel;
        }
      }
    }
    fclose(fp);

    n = i;
    if (n == 0) {
      break;
    }

    cpgsvp(0.15, 0.90, 0.15, 0.40);
    cpgswin(pgx[0], pgx[n-1], -10.0, 10.0);
    if (NSWT == 0) {
      cpgbox("BCNTS", 0, 0, "BCNTS", 0, 0);
      cpglab("x", "y", "");
    }
    cpgline(n, pgx, pgy);

/*
-------------------
*/

    m = nskip;
    for (i=0; i<NSMP; i++) {
      if (m > n - 1) {
        break;
      }
      sx[i] = log10((float)(pgx[m]));
      sy[i] = 0.0;
      k = 0;
      for (j=0; j+m<n; j++) {
        sy[i] += pow((pgy[j+m] - pgy[j]), 2.0);
        k++;
      }
      sy[i] /= (float)k;
      sy[i] = 0.5 * log10(sy[i]);
      m = (int)rint((float)m * 1.2);
    }
    nsmp = i;

    cpgsvp(0.15, 0.90, 0.45, 1.00);
    cpgswin(sx[0], sx[nsmp-1], log10(0.01), log10(10.0));
    if (NSWT == 0) {
      cpgbox("BCNLTS", 0, 0, "BCNLTS", 0, 0);
    }
    cpgsci(1);
    cpgpt(nsmp, sx, sy, 1);
    cpgline(nsmp, sx, sy);
    cpgsci(1);
    printf("%f\n", pow(10.0, sy[0]));

    NSWT++;
  }
  printf("\n\n");

/*
-------------------
*/

  m = nskip;
  for (i=0; i<NSMP; i++) {
    if (m > n - 1) {
      break;
    }
    sx[i] = log10((float)(pgx[m]));
    if (pgx[m] < scale1) {
      sy[i] = log10(coeff1) + 0.5 * expon1 * sx[i];
    } else if (pgx[m] >= scale1 && pgx[m] < scale2) {
      sy[i] = log10(coeff2) + 0.5 * expon2 * sx[i];
    } else if (pgx[m] >= scale2) {
      sy[i] = log10(coeff3);
    }
    m = (int)rint((float)m * 1.2);
  }
  nsmp = i;

  cpgsvp(0.15, 0.90, 0.45, 1.00);
  cpgswin(sx[0], sx[nsmp-1], log10(0.01), log10(10.0));
  cpgbox("BCNLTS", 0, 0, "BCNLTS", 0, 0);
  cpgsci(4);
  cpgpt(nsmp, sx, sy, 1);
  cpgline(nsmp, sx, sy);
  cpgsci(1);
  printf("%f\n", pow(10.0, sy[0]));
  printf("\n\n");

/*
-------------------
*/

  NSWT = 0;
  for (irow=0; irow<nrow; irow+=200) {
    n = 0;
    I = 0;
    while (1) {
      sprintf(fname, "sample00.dat-%3d", I);
      if (fname[12] == ' ') {
        fname[12] = '0';
      }
      if (fname[13] == ' ') {
        fname[13] = '0';
      }
      if (fname[14] == ' ') {
        fname[14] = '0';
      }
      if ((fp = fopen(fname, "r")) == NULL) {
        break;
      }
      while (1) {
        fgets(string, sizeof(string), fp);
        if (strncmp(string, "PIXEL SIZE       ", 17) == 0) {
          sscanf(string+18, "%f", &pixel);
        } else if (strncmp(string, "SCREEN WIDTH (NX)", 17) == 0) {
          sscanf(string+18, "%d", &NEL);
        } else if (strncmp(string, "SCREEN WIDTH (NY)", 17) == 0) {
          sscanf(string+18, "%d", &ncol);
        } else if (strncmp(string, "# END OF HEADER", 15) == 0) {
          break;
         }
      }
      for (i=0; i<NEL; i++) {
        for (j=0; j<ncol; j++) {
          fgets(string, sizeof(string), fp);
          if (j == irow) {
            sscanf(string, "%f,%f,%f", &ftmp1, &ftmp2, &pgy[n]);
            pgx[n] = (float)n;
            pgx[n] *= pixel;
#ifdef __DEBUG__
            if (i == 0 || i == NEL-1) {
              printf("%s", string);
            } 
#endif
            n++;
            if (n == NDAT) {
              fclose(fp);
              break;
            }
          }
        }
        if (n == NDAT) {
          fclose(fp);
          break;
        }
      }
      if (n == NDAT) {
        fclose(fp);
        break;
      }
      fclose(fp);
      I++;
    }
    if (n == 0) {
      break;
    }

    cpgsvp(0.15, 0.90, 0.15, 0.40);
    cpgswin(pgx[0], pgx[n-1], -10.0, 10.0);
    cpgsci(2);
    cpgline(n, pgx, pgy);
    cpgsci(1);

/*
----------
*/

    m = nskip;
    for (i=0; i<NSMP; i++) {
      if (m > n - 1) {
        break;
      }
      sx[i] = log10((float)(pgx[m]));
      sy[i] = 0.0;
      k = 0;
      for (j=0; j+m<n; j++) {
        sy[i] += pow((pgy[j+m] - pgy[j]), 2.0);
        k++;
      }
      sy[i] /= (float)k;
      sy[i] = 0.5 * log10(sy[i]);
      m = (int)rint((float)m * 1.2);
    }
    nsmp = i;

    cpgsvp(0.15, 0.90, 0.45, 1.00);
    cpgswin(sx[0], sx[nsmp-1], log10(0.01), log10(10.0));
    cpgsci(2);
    cpgpt(nsmp, sx, sy, 1);
    cpgline(nsmp, sx, sy);
    cpgsci(1);
    printf("%f\n", pow(10.0, sy[0]));

    NSWT++;

  }

/*
--------
*/

  cpgend(); 

  return 1;
}
