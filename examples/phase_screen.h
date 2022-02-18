#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

/*
---------------------------------------------------------
      DIFINITIONS
---------------------------------------------------------
*/

#define ON                1
#define OFF               0

#define OK                1
#define NG               -1

/*
---------------------------------------------------------
      STRUCTURES
---------------------------------------------------------
*/

struct phase_screen_parameter
        {
          double pixel;
          double H_d;
          double H_s;
          double i_coeffi;
          double o_coeffi;
          double c_coeffi;
          double i_expon;
          double o_expon;
          double i_scale[2];
          double o_scale[2];
          double v[2];
        };

/*
---------------------------------------------------------
      PARAMETERS' DEFINITION
---------------------------------------------------------
*/

#define __RANDOM_SEED__

static double dpi            = 3.141592653589793238462643;

/*
---------------------------------------------------------
      FUNCTIONS
---------------------------------------------------------
*/

void    seed_random
        (
          int
        );
double  random_val1
        (
          void
        );
void    char_copy
        (
          char      *,    char      *
        );
double  gauss_dev
        (
          void
        );
int     turbulent_phase_screen
        (
          int        ,    int        ,    double    *,    double    *,
          double    *,
          struct phase_screen_parameter             *,
          int        ,    int        ,    int        ,
          int       *,    int       *,    double    *,
          int    [][],    int
        );
/*
---------------------------------------------------------
      END
---------------------------------------------------------
*/
