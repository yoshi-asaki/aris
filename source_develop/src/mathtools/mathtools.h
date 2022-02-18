#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void   allanv(
           int     , int     , float  *, float  *, float  *, float   );
void   lstsqr(
           int     , float  *, float  *, float  *, float  *, float  *, int);
void   lstsqr2(
           int     , float  *, float  *, float  *, float  *,
           float  *, float  *, float  *);
int    frotate(
           float  *, float   , char   *);
int    drotate(
           double *, double  , char   *);
double dvec2scalar(
           double  *, int     );
float  vec2scalar(
           float   *, int     );
int    pol_fit(
           int     , float  *, float  *, int     , float  *);
int    pol_fit_d(
           int     , double *, double *, int     , double *);
int    fftrfm(
           float  *, float  *, int    *, int    *, int    *);
void   fmaxmin(
                     int     , float  *, float  *, float  *,
                     int    *, int    *);
void   dmaxmin(
                     int     , double *, double *, double *,
                     int    *, int    *);
int    drotation_matrix_set(
           double *, double  , char   *);
void   dvector_calc(double [][3], double *);
double inner_product3(double *, double *);
void   cross_product3(double *, double *, double *);
