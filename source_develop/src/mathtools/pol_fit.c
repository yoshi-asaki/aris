#include "mathtools.h"

int    pol_fit(int n, float *x, float *y, int mp0, float *c)
{
    /* Builtin functions */

    /* Local variables */
    float a[420]	/* was [20][21] */;
    int   i, j, k, l, mp1;
    int   m2;
    float pivot, w1, w2, w3;
    int   mp2;
    float aik;
    float *w;


    /* Parameter adjustments */
    --y;
    --x;
    --c;

    /* Function Body */
    mp1 = mp0 +1;
    if (mp1 > n || mp1 < 2 || mp1 > 20) {
	return 0;
    }
    else
    {
	mp2 = mp1 + 1;
	m2 = mp0 << 1;
/******
        if ((a=(float *)calloc(mp1*mp2, sizeof(float))) == NULL ||
            (w=(float *)calloc(m2,      sizeof(float))) == NULL)
******/
        if ((w=(float *)calloc(m2,      sizeof(float))) == NULL)
        {
            printf("fiting: mamory allocation error.\n");
            return 0;
        }
    }

    for (i = 1; i <= m2; ++i)
    {
        w1 = (float)0.;
        for (j = 1; j <= n; ++j)
        {
            w2 = x[j];
            w1 += pow(w2, i);
        }
        w[i - 1] = w1;
    }
    for (i = 1; i <= mp1; ++i)
    {
        for (j = 1; j <= mp1; ++j)
        {
            l = i + j - 2;
            if (l != 0)
            {
                a[i + j * 20 - 21] = w[l - 1];
            }
        }
    }
    a[0] = (float)n;
    w1 = (float)0.;
    for (i = 1; i <= n; ++i)
    {
        w1 += y[i];
    }
    a[mp2 * 20 - 20] = w1;
    for (i = 1; i <= mp0; ++i) {
        w1 = (float)0.;
        for (j = 1; j <= n; ++j) {
            w2 = y[j];
            w3 = x[j];
            w1 += w2 * pow(w3, i);
        }
        a[i + 1 + mp2 * 20 - 21] = w1;
    }
    for (k = 1; k <= mp1; ++k) {
        pivot = a[k + k * 20 - 21];
        for (j = k; j <= mp2; ++j) {
            a[k + j * 20 - 21] /= pivot;
        }
        for (i = 1; i <= mp1; ++i) {
            if (i != k) {
               aik = a[i + k * 20 - 21];
               for (j = k; j <= mp2; ++j) {
                   a[i + j * 20 - 21] -= aik * a[k + j * 20 - 21];
                }
            }
        }
    }
    for (i = 1; i <= mp1; ++i)
    {
        c[i] = a[i + mp2 * 20 - 21];
    }

    free (w);
    return 1;
}
