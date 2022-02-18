#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


int  ref_pos_shift(
                   struct source_parameter *tgt_src,
                   struct source_parameter *ref_src,
                   double DS,
                   struct comment_param *cmnt,
                   char  comment[][NCOMLEN], _Bool TV_SWT)
/********
    tgt_src : TARGET SOURCE
    ref_src : REFERENCE SOURCE
********/
{
  char   string[100];
  double theta1, theta2, p[3];
  double s_e[3];
  double dra, ddc;
  double RA_e, DC_e;
  double DPI_sec;
  double q[4], A[9];

/*
===================================================================
*/

  DPI_sec = dpi / 180.0 / 3600.0;

  theta1 = 2.0 * dpi * random_val0();
  theta2 = 2.0 * dpi * random_val0();
  p[0] = cos(theta1) * cos(theta2);
  p[1] = cos(theta1) * sin(theta2);
  p[2] = sin(theta1);
  vector_rotation(ref_src->s_e, p, DS);

  ref_src->RA_e  = atan2(ref_src->s_e[1], ref_src->s_e[0]);
  ref_src->DC_e  = atan2(ref_src->s_e[2], vlen2(ref_src->s_e));

  dra = (ref_src->RA_e - ref_src->RA2k) / DPI_sec * cos(ref_src->DC_e);
  ddc = (ref_src->DC_e - ref_src->DC2k) / DPI_sec;
  sprintf(string,
       "REF Position Shift [sec] : (RA, Dec)=(%12.9lf, %12.9lf)", dra, ddc);
  if (TV_SWT == true) {
    comment_disp(cmnt, comment, string, true);
  } else {
    printf("%s\n", string);
  }

/*
----------------
*/

  minor_shift_refpos(ref_src->s2k, ref_src->s_e, tgt_src->s2k, s_e, q, A);

  RA_e = atan2(s_e[1], s_e[0]);
  DC_e = atan2(s_e[2], vlen2(s_e));

  dra = (RA_e - tgt_src->RA2k) / DPI_sec * cos(tgt_src->DC2k);
  ddc = (DC_e - tgt_src->DC2k) / DPI_sec;
  sprintf(string,
       "TGT Position Shift [sec] : (RA, Dec)=(%12.9lf, %12.9lf)", dra, ddc);
  if (TV_SWT == true) {
    comment_disp(cmnt, comment, string, true);
  } else {
    printf("%s\n", string);
  }

  return 1;
}
