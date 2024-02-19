#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


int    UTC_minus_TAI
       (
          int   *
       );
double MJD
       (
          int   , int    , int   , int   , int   , int   , double 
       );
int    MJD2date
       (
          double , int   * , int   * , int    *,
                         int   * , int   * , int    * 
       );
void   coordinate_rotation
       (
         double *, double , double *, double *,
                           double *, double *
       );
void   vector_rotation
       (
         double *, double *, double 
       );
void   coordinate_rotation_parameter
       (
         double *, double , double *, double *
       );
void   Q_parameter
       (
         double *, double *, double *,  double *, double *
       );
void   Q_param
       (
         double *, double *,  double 
       );
void   Euler_parameter
       (
         double *, double *, double *,  double *
       );
void   xyz2q
       (
         double *, double *, double *, double *, double *
       );
void   xyz2Euler_param
       (
         double *, double *, double *, double *, double *
       );
void   xyz2radec
       (
         double *, int *, int *, double *, int *, int *, double *
       );
void   xyz2radec_rad
       (
         double *, double *, double *
       );
void   matrix_rotation
       (
         double *, double *, double *
       );
void   q2xyz
       (
         double *, double *, double *, double *
       );
void   q2matrix
       (
         double *, double *
       );
void   radec_rad2xyz
       (
         double *, double *, double *
       );
double vlen2
       (
         double *
       );
double vlen3
       (
         double *
       );
double sepang
       (
         double *, double *
       );
void   heirocentric_equatorial_rectangular_coordinates_earth_position
       (
         int *, double , double *
       );
void   sun_position
       (
         int  *, double  ,  double *, double *
       );
void   azel_position
       (
         int  *, double  , double , double , double ,
         double ,  double, _Bool  ,
         double , double , double *, double *,
         double *, double *,  double  , double *
       );
void   azel_rot_speed
       (
         double *, double *, double , double  , double  , double *, double *
       );
double LST
       (
         int  *,  double  , double  , double 
       );
double GST
       (
         int  *,  double  , double  
       );
double  mean_obliquity_of_the_ecliptic
        (
          int       *,    double
        );
void    nutation_calc
        (
          int       *,    double     ,    double    *,    double    *,    
          double    *
        );
double atmospheric_effect
       (
         double , double, double
       );
void   luna_position
       (
         int  *,  double , double  , double  , double  , 
                   double *, double *
       );
void   nutation
       (
         int  *,  double , double *, int 
       );
void   nutation_matrix_set
       (
         char   *, char   *, int   , double  *
       );
void   precession
       (
         int  *,  double , double *, int 
       );
void   precession_calc
       (
         int  *,  double , double *, double *, double *
       );
void   topocentric_equatorial_rectangular_coordinate
       (
         int    *, double  , double , double , double , double *
       );
void   J_system_geocentric_equatorial_rectangular_coordinate
       (
         double ,
         double ,
         double ,
         double *
       );
void   J_system_geocentric_equatorial_rectangular_coordinate2llh
       (
         double *, double *, double *, double *
       );
void   earth_rotation
       (
         int    *, double , double *,   double *
       );
double EPSIRON
       (
         double 
       );
double ET
       (
         int *, double 
       );
double mean_obiliquity_of_the_equipric
       (
         int *, double 
       );
int    minor_shift_refpos(
         double *, double *, double *, double *,
         double *,    double *
       );
void   vector_CC_to_GC(
         double *
       );
void   vector_GC_to_CC(
         double *
       );
void   GC_to_CC(
         double, double, double *
       );
void   CC_to_GC(
         double, double, double *
       );
void   GC_CC(
         double *, double *, double *, double *
       );
int    baseline_number
       (
         int     , int     , int    
       );
int    baseline2antenna_number
       (
         int     , int     , int    *, int    *
       );
double baseline_base2antenna_base_solution
       (
         size_t  , int     , double *, double *, double *, double *, int , int 
       );
void   B1950toJ2000(
         double  , double , double *, double *
       );
double mean_obliquity_of_the_eliptic(
         int    *, double
       );
int    solar_motion_J2000(
         double *, double *, int    , char   *
       );
void  annual_parallax(
         double *, double  , double ,
         int    *, double 
       );
