#include "mathtools.h"

/******************************************************************/
/*                                                                */
/*   SUBROUTINE: 1D FAST FOURIER TARANSFORM                       */
/*                  (SSLIB)                                       */
/*                                                                */
/******************************************************************/


int    fftrfm(float *xr, float *xi, int *n, int *iter, int *nflag)
{
    /* Initialized data */

    float         pi=(float)3.141592653589793238462643;


    /* Local variables */
    int           i, j, k;
    int           j1, j2, it;
    int           nxp, nxp2;
    float         ti, wi, tr, wr, di1, di2, dr1, dr2, arg, sgn, w;


    /* Function Body */

    if (*n < 2)
    {
        return 999;
    }

    if (*iter <= 0)
    {
        *iter = 0;
        i = *n;
        while ((i>>=1) != 0)
        {
            (*iter)++;
        }
    }

    j = 1;
    for (i=0; i<*iter; i++) {
        j <<= 1;
    }

    if (*n != j)
    {
        return 0;
    }

    if (*nflag == 1) {
        sgn = (float)1.0;
    }
    else {
        sgn = (float) -1.0;
    }

    nxp2 = *n;
    for (it=0; it<*iter; it++)
    {
        nxp = nxp2;
        nxp2 >>= 1;
        w = pi / (float)nxp2;
        for (k=0; k<nxp2; k++)
        {
            arg = (float)k * w;
            wr =      cos(arg);
            wi = sgn *sin(arg);
            i = k - nxp;
            for (j=nxp; j<*n+1; j+=nxp)
            {
                j1 = j + i;
                j2 = j1 + nxp2;
                dr1 = xr[j1];
                dr2 = xr[j2];
                di1 = xi[j1];
                di2 = xi[j2];
                tr = dr1 - dr2;
                ti = di1 - di2;
                xr[j1] = dr1 + dr2;
                xi[j1] = di1 + di2;
                xr[j2] = tr * wr - ti * wi;
                xi[j2] = ti * wr + tr * wi;
            }
        }
    }

    j1 = *n >> 1;
    j2 = *n - 2;
    j = 0;
    for (i=0; i<j2; i++)
    {
        if (i < j)
        {
            tr = xr[j];
            ti = xi[j];
            xr[j] = xr[i];
            xi[j] = xi[i];
            xr[i] = tr;
            xi[i] = ti;
        }
        k = j1;

        while (k <= j)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    if (*nflag == 1)
    {
        w = (float) (*n);
        for (i=0; i<*n; i++)
        {
            xr[i] /= w;
            xi[i] /= w;
        }
    }
    return 1;
}
