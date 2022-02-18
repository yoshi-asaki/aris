#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"


void  CC_to_GC(double RA, double DC, double *s)
{
  double e1[3], e2[3], phi1, phi2;

  GC_CC(e1, e2, &phi1, &phi2);

/*
-------------------------------------------
*/

  s[0] = cos(DC) * cos(RA);
  s[1] = cos(DC) * sin(RA);
  s[2] = sin(DC);
  vector_rotation(s, e1, phi1);
  vector_rotation(s, e2, phi2);

  return;

}



void  GC_to_CC(double GL, double GB, double *s)
{
  double e1[3], e2[3], phi1, phi2;

  GC_CC(e1, e2, &phi1, &phi2);

/*
-------------------------------------------
*/

  s[0] = cos(GB) * cos(GL);
  s[1] = cos(GB) * sin(GL);
  s[2] = sin(GB);
  vector_rotation(s, e2, -phi2);
  vector_rotation(s, e1, -phi1);

  return;

}



void  vector_CC_to_GC(double *s)
{
  double e1[3], e2[3], phi1, phi2;
  double p[3];

  GC_CC(e1, e2, &phi1, &phi2);
  vector_rotation(s, e1, phi1);
  vector_rotation(s, e2, phi2);

  return;

}



void  vector_GC_to_CC(double *s)
{
  double e1[3], e2[3], phi1, phi2;

  GC_CC(e1, e2, &phi1, &phi2);
  vector_rotation(s, e2, -phi2);
  vector_rotation(s, e1, -phi1);

  return;

}



void  GC_CC(double *e1, double *e2, double *phi1, double *phi2)
/**** Galactic Coordinate <----> Celestial Coordinate ****/

{
  double GCRA, GCDC, GPRA, GPDC;
  double dpi=3.141592653589793238462643;
  double DPI;
  double x[3], y[3], z[3];
  double u[3], v[3], w[3];

/*
----------------------------
*/

  DPI = dpi / 180.0;


/**** Galactic Center ****/

  GCRA =  ( 17.0 + 45.0/60.0 + 37.224/3600.0) * 15.0 * DPI;
  GCDC = -( 28.0 + 56.0/60.0 + 10.230/3600.0)        * DPI;
  x[0] = cos(GCDC) * cos(GCRA);
  x[1] = cos(GCDC) * sin(GCRA);
  x[2] = sin(GCDC);

/**** Galactic North Pole ****/

  GPRA = ( 12.0 + 51.0/60.0 + 26.282/3600.0) * 15.0 * DPI;
  GPDC = ( 27.0 +  7.0/60.0 + 42.010/3600.0)        * DPI;
  z[0] = cos(GPDC) * cos(GPRA);
  z[1] = cos(GPDC) * sin(GPRA);
  z[2] = sin(GPDC);

  y[0] = z[1]*x[2] - z[2]*x[1];
  y[1] = z[2]*x[0] - z[0]*x[2];
  y[2] = z[0]*x[1] - z[1]*x[0];

  u[0] = 1.0;
  u[1] = 0.0;
  u[2] = 0.0;

  v[0] = 0.0;
  v[1] = 1.0;
  v[2] = 0.0;

  w[0] = 0.0;
  w[1] = 0.0;
  w[2] = 1.0;

  Euler_parameter(z,  w, e1, phi1);
  vector_rotation(x,     e1, *phi1);
  Euler_parameter(x,  u, e2, phi2);

/*
----------------------------
*/

  return;
}
