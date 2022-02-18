#include <stdio.h>
#include <math.h>
#include "../mathtools/mathtools.h"
#include "astrotools.h"


int  solar_motion_J2000(double *vs, double *d_vs,
                        int  DEFINITION_SWT, char *FRAME_SWT)
{
  int    i;
  double DPI, dpi=3.141592653589793238462643;
  double dv[3][3];
  DPI = dpi / 180.0;

/****
  DEFINITION_SWT: 
      definition of Solar motion
        1. Belhaye (1965) in Galactic Structure ed. Blaauw and Schmidt, 61
        2. AIPS Solar Motion
            The definition of the Solar Motion in AIPS is provided by
            Gordon (1975). In the definition, the speed of the solar 
            motion is 20.0 km/s towards (18h, +30deg) at epoch 1900.0.
            According to an NRAO site the above solar motion with J2000 
            Vector from Barycenter is described as 20 km/s towards 
            (18h03m50.29s, +30deg00'16.8").
            See more detail from the following NRAO site:
               - http://www.gb.nrao.edu/~fghigo/gbtdoc/doppler.html
        3. IAU Standard Solar Motion
        4. Kerr and Lynden-Bell (1986)
        5. Dehnen and Binney (1998) [HIPPARCOS]
        6. Francis and Anderson (2008) [HIPPARCOS]
        7. Schonrich, Binney, and Dehnen (2010) [HIPPARCOS]

  FRAME_SWT:
      definition of frame
        G : Galactic Coordinate (UVW)
        C : Celestial Coordinate (RA-Dec)
****/

  if        (DEFINITION_SWT == 1) {
    vs[0]   = 16.5;
    vs[1]   =  0.0;
    vs[2]   =  0.0;
    drotate(vs, -25.0*DPI, "y");
    drotate(vs,  53.0*DPI, "z");
    if (*FRAME_SWT == 'C' || *FRAME_SWT == 'c') {
      vector_GC_to_CC(vs);
    }
    d_vs[0] =  0.0;
    d_vs[1] =  0.0;
    d_vs[2] =  0.0;

  } else if (DEFINITION_SWT == 2) {
    vs[0]   = 20.0;
    vs[1]   =  0.0;
    vs[2]   =  0.0;
    drotate(vs, -(30.0 +  0.0/60.0 + 16.83 /3600.0)        * DPI, "y");
    drotate(vs,  (18.0 +  3.0/60.0 + 50.280/3600.0) * 15.0 * DPI, "z");
/**
    drotate(vs, -(30.0 +  0.0/60.0 + 16.8 /3600.0)        * DPI, "y");
    drotate(vs,  (18.0 +  3.0/60.0 + 50.29/3600.0) * 15.0 * DPI, "z");
**/
    if (*FRAME_SWT == 'G' || *FRAME_SWT == 'g') {
      vector_CC_to_GC(vs);
    }
    d_vs[0] =  0.0;
    d_vs[1] =  0.0;
    d_vs[2] =  0.0;

  } else if (DEFINITION_SWT == 3) {
    vs[0]   = 19.5;
    vs[1]   =  0.0;
    vs[2]   =  0.0;
    drotate(vs, -22.0*DPI, "y");
    drotate(vs,  55.0*DPI, "z");
    if (*FRAME_SWT == 'C' || *FRAME_SWT == 'c') {
      vector_GC_to_CC(vs);
    }
    d_vs[0] =  0.0;
    d_vs[1] =  0.0;
    d_vs[2] =  0.0;

  } else if (DEFINITION_SWT == 4) {
    vs[0]   = 10.0;
    vs[1]   = 15.4;
    vs[2]   =  7.8;
    if (*FRAME_SWT == 'C' || *FRAME_SWT == 'c') {
      vector_GC_to_CC(vs);
    }
    d_vs[0] =  0.0;
    d_vs[1] =  0.0;
    d_vs[2] =  0.0;

  } else if (DEFINITION_SWT == 5) {
    vs[0]   =  10.00;
    vs[1]   =   5.25;
    vs[2]   =   7.17;
    d_vs[0] =   0.36;
    d_vs[1] =   0.62;
    d_vs[2] =   0.38;
    if (*FRAME_SWT == 'C' || *FRAME_SWT == 'c') {
      vector_GC_to_CC(vs);
      for (i=0; i<3; i++) {
        dv[i][0] = 0.0;
        dv[i][1] = 0.0;
        dv[i][2] = 0.0;
        dv[i][i] = d_vs[i];
        vector_GC_to_CC(dv[i]);
      }
      for (i=0; i<3; i++) {
        d_vs[i] = sqrt(dv[0][i]*dv[0][i]
                     + dv[1][i]*dv[1][i]
                     + dv[2][i]*dv[2][i]);
      }
    }

  } else if (DEFINITION_SWT == 6) {
    vs[0]   =   7.5;
    vs[1]   =  13.5;
    vs[2]   =   6.8;
    d_vs[0] =   1.0;
    d_vs[1] =   0.3;
    d_vs[2] =   0.1;
    if (*FRAME_SWT == 'C' || *FRAME_SWT == 'c') {
      vector_GC_to_CC(vs);
    } else if (*FRAME_SWT == 'G' || *FRAME_SWT == 'g') {
      for (i=0; i<3; i++) {
        dv[i][0] = 0.0;
        dv[i][1] = 0.0;
        dv[i][2] = 0.0;
        dv[i][i] = d_vs[i];
        vector_GC_to_CC(dv[i]);
      }
      for (i=0; i<3; i++) {
        d_vs[i] = sqrt(dv[0][i]*dv[0][i]
                     + dv[1][i]*dv[1][i]
                     + dv[2][i]*dv[2][i]);
      }
    }

  } else if (DEFINITION_SWT == 7) {
    vs[0]   =  11.10;
    vs[1]   =  12.24;
    vs[2]   =   7.25;
    d_vs[0] =   0.75;
    d_vs[1] =   0.47;
    d_vs[2] =   0.37;
    if (*FRAME_SWT == 'C' || *FRAME_SWT == 'c') {
      vector_GC_to_CC(vs);
      for (i=0; i<3; i++) {
        dv[i][0] = 0.0;
        dv[i][1] = 0.0;
        dv[i][2] = 0.0;
        dv[i][i] = d_vs[i];
        vector_GC_to_CC(dv[i]);
      }
      for (i=0; i<3; i++) {
        d_vs[i] = sqrt(dv[0][i]*dv[0][i]
                     + dv[1][i]*dv[1][i]
                     + dv[2][i]*dv[2][i]);
      }
    }
  }

  return 1;
}
