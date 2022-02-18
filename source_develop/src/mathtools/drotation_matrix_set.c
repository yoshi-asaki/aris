#include <math.h>
#include "mathtools.h"

int   drotation_matrix_set(double *mtrx, double e, char *axis)
{
  if (*axis == 'x' || *axis == 'X') {
    *(mtrx + 3 * 0 + 0) =     1.0;
    *(mtrx + 3 * 1 + 0) =     0.0;
    *(mtrx + 3 * 2 + 0) =     0.0;
    *(mtrx + 3 * 0 + 1) =     0.0;
    *(mtrx + 3 * 1 + 1) =  cos(e);
    *(mtrx + 3 * 2 + 1) = -sin(e);
    *(mtrx + 3 * 0 + 2) =     0.0;
    *(mtrx + 3 * 1 + 2) =  sin(e);
    *(mtrx + 3 * 2 + 2) =  cos(e);
  } else if (*axis == 'y' || *axis == 'Y') {
    *(mtrx + 3 * 0 + 0) =  cos(e);
    *(mtrx + 3 * 1 + 0) =     0.0;
    *(mtrx + 3 * 2 + 0) =  sin(e);
    *(mtrx + 3 * 0 + 1) =     0.0;
    *(mtrx + 3 * 1 + 1) =     1.0;
    *(mtrx + 3 * 2 + 1) =     0.0;
    *(mtrx + 3 * 0 + 2) = -sin(e);
    *(mtrx + 3 * 1 + 2) =     0.0;
    *(mtrx + 3 * 2 + 2) =  cos(e);
  } else if (*axis == 'z' || *axis == 'Z') {
    *(mtrx + 3 * 0 + 0) =  cos(e);
    *(mtrx + 3 * 1 + 0) = -sin(e);
    *(mtrx + 3 * 2 + 0) =     0.0;
    *(mtrx + 3 * 0 + 1) =  sin(e);
    *(mtrx + 3 * 1 + 1) =  cos(e);
    *(mtrx + 3 * 2 + 1) =     0.0;
    *(mtrx + 3 * 0 + 2) =     0.0;
    *(mtrx + 3 * 1 + 2) =     0.0;
    *(mtrx + 3 * 2 + 2) =     1.0;
  }
}
