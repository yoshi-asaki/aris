#include <stdio.h>
#include <math.h>
#include <aris.h>

/****
Neill Mapping Function (1996)
A. E. Niell, Journal of Geophysical Research, 101, B2, pp.3227-3246, 1996.
****/

double marini_mfc(double , double , double , double );



int  nmf20(double *hmf,   double *wmf,
           double *h,     double *h_hcr, double *w,
           double H,      double el)
{
  double sin_el, mfh, mf_hcr;

  sin_el  = sin(el);

  *wmf    = marini_mfc(sin_el,     w[0],     w[1],     w[2]);
  *wmf   /= marini_mfc(   1.0,     w[0],     w[1],     w[2]);

  *hmf    = marini_mfc(sin_el,     h[0],     h[1],     h[2]);
  *hmf   /= marini_mfc(   1.0,     h[0],     h[1],     h[2]);
  mf_hcr  = marini_mfc(sin_el, h_hcr[0], h_hcr[1], h_hcr[2]);
  mf_hcr /= marini_mfc(   1.0, h_hcr[0], h_hcr[1], h_hcr[2]);
  *hmf   += H*1.0e-3 * (1.0/sin_el - mf_hcr);

  return 1;
}



double marini_mfc(double sin_el, double a, double b, double c)
{
  double marini_mf;

  marini_mf = 1.0 / (  sin_el + a/(  sin_el + b/(  sin_el + c)));

  return marini_mf;
}
