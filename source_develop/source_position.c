#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

/****
#define __DEBUG__
****/

int  source_position(struct source_parameter *src,
                     struct pair_src_info *pair_src,
                     struct char_src_info *ch_src,
                     int SEP_MODE, int POS_MODE)
{
  int    i, ns;
  double DPI;
  double ra, dc, vl;
  double dra, ddc;
  double sa, pa, mid_s[3], s[3];

/*
--------
*/

  DPI = dpi / 180.0;

  if (SEP_MODE == SRC__RA__DEC) {
    input_star_position(ch_src->tgt_ra, ch_src->tgt_dc,
                        &src[0].RA2k, &src[0].DC2k);
    input_star_position(ch_src->ref_ra, ch_src->ref_dc,
                        &src[1].RA2k, &src[1].DC2k);
    for (ns=0; ns<SRC_NUM; ns++) {
      src[ns].s2k[0] = cos(src[ns].DC2k) * cos(src[ns].RA2k);
      src[ns].s2k[1] = cos(src[ns].DC2k) * sin(src[ns].RA2k);
      src[ns].s2k[2] = sin(src[ns].DC2k);
    }

    mid_s[0] = src[0].s2k[0] + src[1].s2k[0];
    mid_s[1] = src[0].s2k[1] + src[1].s2k[1];
    mid_s[2] = src[0].s2k[2] + src[1].s2k[2];
    pair_src->mid_ra = atan2(mid_s[1], mid_s[0]);
    pair_src->mid_dc = atan2(mid_s[2], vlen2(mid_s));
    mid_s[0] = cos(pair_src->mid_dc) * cos(pair_src->mid_ra);
    mid_s[1] = cos(pair_src->mid_dc) * sin(pair_src->mid_ra);
    mid_s[2] = sin(pair_src->mid_dc);

    s[0] = src[0].s2k[0];
    s[1] = src[0].s2k[1];
    s[2] = src[0].s2k[2];
    drotate(s, -pair_src->mid_ra, "z");
    drotate(s,  pair_src->mid_dc, "y");
    pair_src->dlt_ra = -2.0 * atan2(s[1], s[0])     / DPI;
    pair_src->dlt_dc = -2.0 * atan2(s[2], vlen2(s)) / DPI;
    pair_src->sepang = sepang(src[0].s2k, src[1].s2k)            / DPI;
    pair_src->posang = atan2(pair_src->dlt_ra, pair_src->dlt_dc) / DPI;

/*
----
*/

  } else if (SEP_MODE == SRC_dRA_dDEC) {
    if (       POS_MODE == SRC_POS1) {
      input_star_position(ch_src->tgt_ra, ch_src->tgt_dc,
                          &src[0].RA2k, &src[0].DC2k);
      ra = src[0].RA2k;
      dc = src[0].DC2k;
      src[0].s2k[0] = cos(dc) * cos(ra);
      src[0].s2k[1] = cos(dc) * sin(ra);
      src[0].s2k[2] = sin(dc);

      dra = 0.5 * pair_src->dlt_ra * DPI;
      ddc = 0.5 * pair_src->dlt_dc * DPI;
      src[1].s2k[0] = cos(ddc) * cos(dra);
      src[1].s2k[1] = cos(ddc) * sin(dra);
      src[1].s2k[2] = sin(ddc);
      drotate(src[1].s2k,  dra, "z");
      drotate(src[1].s2k, -ddc, "y");
      drotate(src[1].s2k,  -dc, "y");
      drotate(src[1].s2k,   ra, "z");
      vl = vlen3(src[1].s2k);
      src[1].s2k[0] /= vl;
      src[1].s2k[1] /= vl;
      src[1].s2k[2] /= vl;
      src[1].RA2k = atan2(src[1].s2k[1], src[1].s2k[0]);
      src[1].DC2k = atan2(src[1].s2k[2], vlen2(src[1].s2k));

      mid_s[0] = src[0].s2k[0] + src[1].s2k[0];
      mid_s[1] = src[0].s2k[1] + src[1].s2k[1];
      mid_s[2] = src[0].s2k[2] + src[1].s2k[2];
      pair_src->mid_ra = atan2(mid_s[1], mid_s[0]);
      pair_src->mid_dc = atan2(mid_s[2], vlen2(mid_s));
      mid_s[0] = cos(pair_src->mid_dc) * cos(pair_src->mid_ra);
      mid_s[1] = cos(pair_src->mid_dc) * sin(pair_src->mid_ra);
      mid_s[2] = sin(pair_src->mid_dc);

    } else if (POS_MODE == SRC_POS2) {
      input_star_position(ch_src->mid_ra,   ch_src->mid_dc,
                          &(pair_src->mid_ra), &(pair_src->mid_dc));
      ra = pair_src->mid_ra;
      dc = pair_src->mid_dc;
      mid_s[0] = cos(dc) * cos(ra);
      mid_s[1] = cos(dc) * sin(ra);
      mid_s[2] = sin(dc);

      dra = -0.5 * pair_src->dlt_ra * DPI;
      ddc = -0.5 * pair_src->dlt_dc * DPI;
      src[0].s2k[0] = cos(ddc) * cos(dra);
      src[0].s2k[1] = cos(ddc) * sin(dra);
      src[0].s2k[2] = sin(ddc);
      drotate(src[0].s2k, -dc, "y");
      drotate(src[0].s2k,  ra, "z");
      vl = vlen3(src[0].s2k);
      src[0].s2k[0] /= vl;
      src[0].s2k[1] /= vl;
      src[0].s2k[2] /= vl;
      src[0].RA2k = atan2(src[0].s2k[1], src[0].s2k[0]);
      src[0].DC2k = atan2(src[0].s2k[2], vlen2(src[0].s2k));

      mid_s[0] = src[0].s2k[0] + src[1].s2k[0];
      mid_s[1] = src[0].s2k[1] + src[1].s2k[1];
      mid_s[2] = src[0].s2k[2] + src[1].s2k[2];
      pair_src->mid_ra = atan2(mid_s[1], mid_s[0]);
      pair_src->mid_dc = atan2(mid_s[2], vlen2(mid_s));
      mid_s[0] = cos(pair_src->mid_dc) * cos(pair_src->mid_ra);
      mid_s[1] = cos(pair_src->mid_dc) * sin(pair_src->mid_ra);
      mid_s[2] = sin(pair_src->mid_dc);

      dra =  0.5 * pair_src->dlt_ra * DPI;
      ddc =  0.5 * pair_src->dlt_dc * DPI;
      src[1].s2k[0] = cos(ddc) * cos(dra);
      src[1].s2k[1] = cos(ddc) * sin(dra);
      src[1].s2k[2] = sin(ddc);
      drotate(src[1].s2k, -dc, "y");
      drotate(src[1].s2k,  ra, "z");
      vl = vlen3(src[1].s2k);
      src[1].s2k[0] /= vl;
      src[1].s2k[1] /= vl;
      src[1].s2k[2] /= vl;
      src[1].RA2k = atan2(src[1].s2k[1], src[1].s2k[0]);
      src[1].DC2k = atan2(src[1].s2k[2], vlen2(src[1].s2k));
    }

    pair_src->sepang  = sepang(src[0].s2k, src[1].s2k)            / DPI;
    pair_src->posang  = atan2(pair_src->dlt_ra, pair_src->dlt_dc) / DPI;

/*
----
*/

  } else if (SEP_MODE == SRC_SEP_POSA) {
    sa = pair_src->sepang * DPI;
    pa = pair_src->posang * DPI;

    if (       POS_MODE == SRC_POS1) {
/******
      input_star_position(ch_src->tgt_ra, ch_src->tgt_dc,
                          &src[0].RA2k, &src[0].DC2k);
      ra = src[0].RA2k;
      dc = src[0].DC2k;
      src[0].s2k[0] = cos(dc) * cos(ra);
      src[0].s2k[1] = cos(dc) * sin(ra);
      src[0].s2k[2] = sin(dc);

      dra =  0.0;
      ddc =  0.5 * sa;
      src[1].s2k[0] = cos(ddc) * cos(dra);
      src[1].s2k[1] = cos(ddc) * sin(dra);
      src[1].s2k[2] = sin(ddc);
      drotate(src[1].s2k, -pa, "x");

      dra = atan2(src[1].s2k[1], src[1].s2k[0]);
      ddc = atan2(src[1].s2k[2], vlen2(src[1].s2k));
      drotate(src[1].s2k,  dra, "z");
      drotate(src[1].s2k, -ddc, "y");
      drotate(src[1].s2k,  -dc, "y");
      drotate(src[1].s2k,   ra, "z");
      vl = vlen3(src[1].s2k);
      src[1].s2k[0] /= vl;
      src[1].s2k[1] /= vl;
      src[1].s2k[2] /= vl;
      src[1].RA2k = atan2(src[1].s2k[1], src[1].s2k[0]);
      src[1].DC2k = atan2(src[1].s2k[2], vlen2(src[1].s2k));

      pair_src->dlt_ra = 2.0 * dra / DPI;
      pair_src->dlt_dc = 2.0 * ddc / DPI;

      mid_s[0] = src[0].s2k[0] + src[1].s2k[0];
      mid_s[1] = src[0].s2k[1] + src[1].s2k[1];
      mid_s[2] = src[0].s2k[2] + src[1].s2k[2];
      pair_src->mid_ra = atan2(mid_s[1], mid_s[0]);
      pair_src->mid_dc = atan2(mid_s[2], vlen2(mid_s));
      mid_s[0] = cos(pair_src->mid_dc) * cos(pair_src->mid_ra);
      mid_s[1] = cos(pair_src->mid_dc) * sin(pair_src->mid_ra);
      mid_s[2] = sin(pair_src->mid_dc);
******/


      input_star_position(ch_src->tgt_ra, ch_src->tgt_dc,
                          &src[0].RA2k, &src[0].DC2k);
      ra = src[0].RA2k;
      dc = src[0].DC2k;
      src[0].s2k[0] = cos(dc) * cos(ra);
      src[0].s2k[1] = cos(dc) * sin(ra);
      src[0].s2k[2] = sin(dc);

      dra =  0.0;
      ddc =  sa;
      src[1].s2k[0] = cos(ddc) * cos(dra);
      src[1].s2k[1] = cos(ddc) * sin(dra);
      src[1].s2k[2] = sin(ddc);
      drotate(src[1].s2k, -pa, "x");
      drotate(src[1].s2k, -dc, "y");
      drotate(src[1].s2k,  ra, "z");

      vl = vlen3(src[1].s2k);
      src[1].s2k[0] /= vl;
      src[1].s2k[1] /= vl;
      src[1].s2k[2] /= vl;
      src[1].RA2k = atan2(src[1].s2k[1], src[1].s2k[0]);
      src[1].DC2k = atan2(src[1].s2k[2], vlen2(src[1].s2k));

      pair_src->dlt_ra = (src[1].RA2k - src[0].RA2k) / DPI;
      if (fabs(pair_src->dlt_ra) >  180.0) {
        pair_src->dlt_ra -= pair_src->dlt_ra / fabs(pair_src->dlt_ra) * 360.0;
      }
      pair_src->dlt_dc = (src[1].DC2k - src[0].DC2k) / DPI;

      mid_s[0] = src[0].s2k[0] + src[1].s2k[0];
      mid_s[1] = src[0].s2k[1] + src[1].s2k[1];
      mid_s[2] = src[0].s2k[2] + src[1].s2k[2];
      pair_src->mid_ra = atan2(mid_s[1], mid_s[0]);
      pair_src->mid_dc = atan2(mid_s[2], vlen2(mid_s));
      mid_s[0] = cos(pair_src->mid_dc) * cos(pair_src->mid_ra);
      mid_s[1] = cos(pair_src->mid_dc) * sin(pair_src->mid_ra);
      mid_s[2] = sin(pair_src->mid_dc);

      output_star_position(ch_src->ref_ra, ch_src->ref_dc, (src+1)->s2k);
      output_star_position(ch_src->mid_ra, ch_src->mid_dc, mid_s);
      sprintf(ch_src->dlt_ra, "%lf", pair_src->dlt_ra);
      sprintf(ch_src->dlt_dc, "%lf", pair_src->dlt_dc);

    } else if (POS_MODE == SRC_POS2) {
      input_star_position(ch_src->mid_ra, ch_src->mid_dc,
                          &(pair_src->mid_ra), &(pair_src->mid_dc));
      ra = pair_src->mid_ra;
      dc = pair_src->mid_dc;
      mid_s[0] = cos(dc) * cos(ra);
      mid_s[1] = cos(dc) * sin(ra);
      mid_s[2] = sin(dc);

      dra =  0.0;
      ddc = -0.5 * sa;
      src[0].s2k[0] = cos(ddc) * cos(dra);
      src[0].s2k[1] = cos(ddc) * sin(dra);
      src[0].s2k[2] = sin(ddc);
      drotate(src[0].s2k,  pa, "x");
      drotate(src[0].s2k, -dc, "y");
      drotate(src[0].s2k,  ra, "z");
      vl = vlen3(src[0].s2k);
      src[0].s2k[0] /= vl;
      src[0].s2k[1] /= vl;
      src[0].s2k[2] /= vl;
      src[0].RA2k = atan2(src[0].s2k[1], src[0].s2k[0]);
      src[0].DC2k = atan2(src[0].s2k[2], vlen2(src[0].s2k));

      dra =  0.0;
      ddc =  0.5 * sa;
      src[1].s2k[0] = cos(ddc) * cos(dra);
      src[1].s2k[1] = cos(ddc) * sin(dra);
      src[1].s2k[2] = sin(ddc);
      drotate(src[1].s2k, -pa, "x");
      dra = atan2(src[1].s2k[1], src[1].s2k[0]);
      ddc = atan2(src[1].s2k[2], vlen2(src[1].s2k));
      drotate(src[1].s2k, -dc, "y");
      drotate(src[1].s2k,  ra, "z");
      vl = vlen3(src[1].s2k);
      src[1].s2k[0] /= vl;
      src[1].s2k[1] /= vl;
      src[1].s2k[2] /= vl;
      src[1].RA2k = atan2(src[1].s2k[1], src[1].s2k[0]);
      src[1].DC2k = atan2(src[1].s2k[2], vlen2(src[1].s2k));

      pair_src->dlt_ra = 2.0 * dra / DPI;
      pair_src->dlt_dc = 2.0 * ddc / DPI;

      mid_s[0] = src[0].s2k[0] + src[1].s2k[0];
      mid_s[1] = src[0].s2k[1] + src[1].s2k[1];
      mid_s[2] = src[0].s2k[2] + src[1].s2k[2];
      pair_src->mid_ra = atan2(mid_s[1], mid_s[0]);
      pair_src->mid_dc = atan2(mid_s[2], vlen2(mid_s));
      mid_s[0] = cos(pair_src->mid_dc) * cos(pair_src->mid_ra);
      mid_s[1] = cos(pair_src->mid_dc) * sin(pair_src->mid_ra);
      mid_s[2] = sin(pair_src->mid_dc);

    }
  }

#ifdef __DEBUG__
  printf("__POS_DEBUG_1_TGT__     : %s   %s\n", ch_src->tgt_ra,  ch_src->tgt_dc);
  printf("__POS_DEBUG_1_REF__     : %s   %s\n", ch_src->ref_ra,  ch_src->ref_dc);
  printf("__POS_DEBUG_1_MID__     : %s   %s\n", ch_src->mid_ra,  ch_src->mid_dc);
  printf("__POS_DEBUG_1_DLTRADC__ : %s   %s\n", ch_src->dlt_ra,  ch_src->dlt_dc);
  printf("__POS_DEBUG_1_SA_PA____ : %s   %s\n", ch_src->sepang,  ch_src->posang);
#endif /*__DEBUG__*/

/*
---------------------------------------------------------
*/

  if (       SEP_MODE == SRC__RA__DEC) {
    output_star_position(ch_src->mid_ra, ch_src->mid_dc, mid_s);

    sprintf(ch_src->sepang, "%15.10lf", pair_src->sepang);
    sprintf(ch_src->posang, "%15.10lf", pair_src->posang);
    sprintf(ch_src->dlt_ra,    "%15.10lf", pair_src->dlt_ra);
    sprintf(ch_src->dlt_dc,    "%15.10lf", pair_src->dlt_dc);

  } else if (SEP_MODE == SRC_dRA_dDEC) {
    if (       POS_MODE == SRC_POS1) {
      output_star_position(ch_src->ref_ra, ch_src->ref_dc, src[1].s2k);
    } else if (POS_MODE == SRC_POS2) {
      output_star_position(ch_src->mid_ra, ch_src->mid_dc, mid_s);
    }
    sprintf(ch_src->sepang, "%15.10lf", pair_src->sepang);
    sprintf(ch_src->posang, "%15.10lf", pair_src->posang);
  } else if (SEP_MODE == SRC_SEP_POSA) {
    if (       POS_MODE == SRC_POS1) {
      output_star_position(ch_src->ref_ra, ch_src->ref_dc, src[1].s2k);
    } else if (POS_MODE == SRC_POS2) {
      output_star_position(ch_src->mid_ra, ch_src->mid_dc, mid_s);
    }
    sprintf(ch_src->dlt_ra,    "%15.10lf", pair_src->dlt_ra);
    sprintf(ch_src->dlt_dc,    "%15.10lf", pair_src->dlt_dc);
  }

  while ((i=strlen(ch_src->sepang) - 1) >= 6) {
    if (*(ch_src->sepang+i) == '0') {
      *(ch_src->sepang+i) = '\0';
    } else {
      break;
    }
  }
  while ((i=strlen(ch_src->posang) - 1) >= 6) {
    if (*(ch_src->posang+i) == '0') {
      *(ch_src->posang+i) = '\0';
    } else {
      break;
    }
  }
  while ((i=strlen(ch_src->dlt_ra) - 1) >= 6) {
    if (*(ch_src->dlt_ra+i) == '0') {
      *(ch_src->dlt_ra+i) = '\0';
    } else {
      break;
    }
  }
  while ((i=strlen(ch_src->dlt_dc) - 1) >= 6) {
    if (*(ch_src->dlt_dc+i) == '0') {
      *(ch_src->dlt_dc+i) = '\0';
    } else {
      break;
    }
  }

#ifdef __DEBUG__
  printf("__POS_DEBUG_2_TGT__     : %s   %s\n", ch_src->tgt_ra,  ch_src->tgt_dc);
  printf("__POS_DEBUG_2_REF__     : %s   %s\n", ch_src->ref_ra,  ch_src->ref_dc);
  printf("__POS_DEBUG_2_MID__     : %s   %s\n", ch_src->mid_ra,  ch_src->mid_dc);
  printf("__POS_DEBUG_1_DLTRADC__ : %s   %s\n", ch_src->dlt_ra,  ch_src->dlt_dc);
  printf("__POS_DEBUG_1_SA_PA____ : %s   %s\n", ch_src->sepang,  ch_src->posang);
#endif /*__DEBUG__*/

/*
---------------------------------------------------------
*/

  for (ns=0; ns<2; ns++) {
    src[ns].RA2k = atan2(src[ns].s2k[1], src[ns].s2k[0]);
    src[ns].DC2k = atan2(src[ns].s2k[2], vlen2(src[ns].s2k));

    src[ns].RA   = src[ns].RA2k;
    src[ns].DC   = src[ns].DC2k;
    src[ns].RA_e = src[ns].RA2k;
    src[ns].DC_e = src[ns].DC2k;
    for (i=0; i<3; i++) {
      src[ns].s[i]   = src[ns].s2k[i];
      src[ns].s_e[i] = src[ns].s2k[i];
    }
  }

/*
---------------------------------------------------------
*/

  return 0;
}
