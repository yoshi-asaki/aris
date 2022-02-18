/* LeaST SQuaRe fitting                          */
/*  0: least square fitting                      */
/*  1: average                                   */
/*  2: offset from 0                             */
/*                                               */

#include "mathtools.h"

void   lstsqr(int n, float *x, float *y, float *a, float *b,
              float *variance, int mode)
{
    int    i;
    double xm, ym, sx, sy, u;

    *a = 0.0;
    *b = 0.0;
    xm = 0.0;
    ym = 0.0;
    sx = 0.0;
    sy = 0.0;
    *variance = 0;

    if (mode == 0) {
        for (i=0; i<n; i++) {
            xm += x[i];
            ym += y[i];
            sx += x[i] * x[i];
            sy += x[i] * y[i];
        }
        xm /= (double)n;
        ym /= (double)n;
        sx /= (double)n;
        sy /= (double)n;
        u = sx - xm*xm;
        *a = (sy -xm*ym) / u;
        *b = (ym*sx -xm*sy) / u;
        for (i=0; i<n; i++) {
            *variance += (y[i] - *a*x[i] - *b) * (y[i] - *a*x[i] - *b);
        }
    }

    else if (mode == 1) {
        for (i=0; i<n; i++) {
            ym += y[i];
        }
        ym /= (double)n;
        for (i=0; i<n; i++) {
            *variance += (y[i] -ym) * (y[i] -ym);
        }
        *a = 0.0;
        *b = ym;
    }

    else if (mode == 2) {
        for (i=0; i<n; i++) {
            *variance += y[i] * y[i];
        }
        *a = 0.0;
        *b = 0.0;
    }
    else
    {
        printf("lstsqr: invalid process mode <%d>\n", mode);
    }

    if (n != 1) {
        *variance /= (double)(n -1);
    } else {
        *variance = 0.0;
    }


    return;
}
