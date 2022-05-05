#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

#define ERROR_SECTION    10
#define SRCPROC_SECTION  30
#define SOURCE_SECTION   40
#define TIME_SECTION     50
#define GRT_SECTION      60
#define ARRAY_SECTION    70
#define SRT_SECTION      90
#define ANTLST_SECTION  100

#define CHAR_LEN         30
#define SECT_DAT        120
#define ARRAY_MAX        20
#define SEPANG_LIM     40.0

int  ant_list_chk(char *, float *,
                  struct comment_param      *,
                  char [][NCOMLEN],  _Bool  );


int  obs_param_input(_Bool *ERROR_FLAG,
                     struct array_parameter *array,
                     int  *wave_id,
                     int   *ANT_NUM,    int  *GRT_NUM,  int  *SRT_NUM,
                     struct srt_orbit_parameter *srt,
                     double *grt_elevation_limit,
                     double sep_angle_limit_from_earth_limb,
                     int  *TimUTC, double *UT1_UTC, double *obs_duration,
                     struct source_parameter *src,
                     struct source_parameter *sun,
                     char   *antenna_list_file,
                     struct antenna_parameter  *ant_prm,
                     char   ant_code[][10],
                     struct comment_param      *cmnt,
                     char comment[][NCOMLEN],
                     _Bool  TV_SWT, float *cursor_pos, int  *pgid)
{
  int     i, j, k, I, J;
  int     idum, nstr, iant, iarray;
  char    string[500];
  char    antenna_list_file_tmp[500];
  int     timUTC[6];
  int     ERROR_MENU_NUM;
  int     ALL_ERROR_FLAG;
  _Bool   SEPANG_UPDATE;
  int     ANT_NUM_tmp, GRT_NUM_tmp, grt_num;
  struct  char_obs_time  ch_obs_t;

  char    error_source[ERROR_NUM_P2][CHAR_LEN], srtaer[CHAR_LEN];
  float   y_pos, source_y_pos;
  float   bttn_box[SECT_DAT][4];
  float   ant_bttn_box[ANTMAX][4];

  int     ARRAY_NUM;
  struct  array_id {
            _Bool  flag;
            int    id;
            int    type;
            char   name[10];
          } array_id[ARRAY_MAX];
  struct  array_type {
            _Bool  flag;
            char   name[20];
          } array_type[2];
  _Bool   SRT_SWT, START_FLAG;

  int     BASIC_SEP_MODE=3, SRCPROC_MODE;
  _Bool   SRCPROC_SWT;
  int     POS_MODE, SEP_MODE;
  _Bool   srcproc_flag[NSRCPROC];
  char    srcproc_name[NSRCPROC][20];
  struct  char_src_info  ch_src;
  struct  pair_src_info pair_src;
  struct  antenna_parameter ant_prm_tmp;
  FILE    *fp;

  char    ch_grt_el_lim[20];
  struct  char_srt_info  ch_srt[SRTMAX];
  float   pitch=0.03;

/*
-----------------------------------------------------
*/

  START_FLAG     = false;
  SRCPROC_MODE   = SRCPROC1;
  ERROR_MENU_NUM = ERROR_NUM - 1;
  nstr           = sizeof(string);

/*
-----------------------------------------------------
*/

  sprintf((src+0)->name, "SOURCE1");
  sprintf((src+1)->name, "SOURCE2");
  sprintf((src+2)->name, "SOURCE1");

/*
-----------------------------------------------------
*/

  array_type[0].flag = true;
  array_type[1].flag = false;

  obs_param_file_io(ERROR_FLAG, antenna_list_file,
                    array, ANT_NUM, GRT_NUM, SRT_NUM,
                    srt, grt_elevation_limit,
                    sep_angle_limit_from_earth_limb,
                    TimUTC, UT1_UTC, obs_duration,
                    &SRCPROC_MODE, src, sun, ant_code,
                    ch_grt_el_lim, &pair_src, &ch_src,
                    ch_srt, &ch_obs_t, 0);

  if (       array->TYPE == _VLBI_ARRAY_) {
    array_type[0].flag = true;
    array_type[1].flag = false;
  } else if (array->TYPE == __CONNECTED_) {
    array_type[0].flag = false;
    array_type[1].flag = true;
  }

  if (ant_list_chk(antenna_list_file, bttn_box[I=ANTLST_SECTION+2],
                   cmnt, comment, false) == 1) {
    ANT_NUM_tmp = array_config(ALL_ANT, wave_id, *SRT_NUM, &GRT_NUM_tmp,
                               ant_prm,   "",
                               antenna_list_file,  false,   true);
    for (iant=0; iant<ANT_NUM_tmp; iant++) {
      ant_prm[iant].UFL = false;
    }
    for (j=0; j<*GRT_NUM; j++) {
      for (iant=0; iant<GRT_NUM_tmp; iant++) {
        if (strncmp(ant_code[j], (ant_prm+iant)->IDC,
                    strlen(ant_code[j])) == 0) {
          ant_prm[iant].UFL = true;
          break;
        }
      }
    }
  }

/*
--------
*/

  in__src_proc(SRCPROC_MODE, &SEP_MODE, &POS_MODE);
  if (SEP_MODE == SRC__RA__DEC) {
    if (strlen(ch_src.tgt_ra) == 0 || strlen(ch_src.tgt_dc) == 0 ||
        strlen(ch_src.ref_ra) == 0 || strlen(ch_src.ref_dc) == 0) {
      if (strlen(ch_src.tgt_ra) != 0 && strlen(ch_src.tgt_dc) != 0) {
        if (strlen(ch_src.dlt_ra) != 0 && strlen(ch_src.dlt_dc) != 0) {
          SEP_MODE = SRC_dRA_dDEC;
          POS_MODE = SRC_POS1;
        } else if (strlen(ch_src.sepang) != 0 && strlen(ch_src.posang) != 0) {
          SEP_MODE = SRC_dRA_dDEC;
          POS_MODE = SRC_POS2;
        }
      } else if (strlen(ch_src.mid_ra) != 0 && strlen(ch_src.mid_dc) != 0) {
        if (strlen(ch_src.dlt_ra) != 0 && strlen(ch_src.dlt_dc) != 0) {
          SEP_MODE = SRC_SEP_POSA;
          POS_MODE = SRC_POS1;
        } else if (strlen(ch_src.sepang) != 0 && strlen(ch_src.posang) != 0) {
          SEP_MODE = SRC_SEP_POSA;
          POS_MODE = SRC_POS2;
        }
      }
    }
  } else if (SEP_MODE == SRC_dRA_dDEC) {
    if (strlen(ch_src.tgt_ra) == 0 || strlen(ch_src.tgt_dc) == 0) {
      if (strlen(ch_src.mid_ra) != 0 && strlen(ch_src.mid_dc) != 0) {
        if (strlen(ch_src.dlt_ra) != 0 && strlen(ch_src.dlt_dc) != 0) {
          SEP_MODE = SRC_SEP_POSA;
          POS_MODE = SRC_POS1;
        } else if (strlen(ch_src.sepang) != 0 && strlen(ch_src.posang) != 0) {
          SEP_MODE = SRC_SEP_POSA;
          POS_MODE = SRC_POS2;
        }
      }
    }
  } else if (SEP_MODE == SRC_SEP_POSA) {
    if (strlen(ch_src.mid_ra) == 0 || strlen(ch_src.mid_dc) == 0) {
      if (strlen(ch_src.tgt_ra) != 0 && strlen(ch_src.tgt_dc) != 0) {
        if (strlen(ch_src.ref_ra) != 0 && strlen(ch_src.ref_dc) != 0) {
          SEP_MODE = SRC__RA__DEC;
        } else if (strlen(ch_src.dlt_ra) != 0 && strlen(ch_src.dlt_dc) != 0) {
          SEP_MODE = SRC_dRA_dDEC;
          POS_MODE = SRC_POS1;
        } else if (strlen(ch_src.posang) != 0 && strlen(ch_src.sepang) != 0) {
          SEP_MODE = SRC_dRA_dDEC;
          POS_MODE = SRC_POS2;
        }
      }
    }
  }

  SRCPROC_MODE = out_src_proc(SEP_MODE, POS_MODE);
  source_position(src, &pair_src, &ch_src, SEP_MODE, POS_MODE);

  for (i=0; i<NSRCPROC; i++) {
    srcproc_flag[i] = false;
  }
  srcproc_flag[SRCPROC_MODE] = true;

/*
-----------------------------------------------------
*/

  sprintf(error_source[APOSER], "Antenna Position Error");
  sprintf(error_source[TPOSER], "Source-1 Position Error");
  sprintf(error_source[RPOSER], "Source-2 Position Error");
  sprintf(error_source[EOPERR], "EOP Error");
  sprintf(error_source[TDSECZ], "Trospospheric Error");
  sprintf(error_source[IDSECZ], "Ionospheric Error");
  sprintf(error_source[TWVTRB], "Water Vapor Turbulence");
  sprintf(error_source[DRYTRB], "Dry Air Turbulence");
  sprintf(error_source[IONTRB], "TEC Turbulence");
  sprintf(error_source[THRMNS], "Thermal Noise Error");
  sprintf(error_source[FQSERR], "Frequency Standard Error");
  sprintf(error_source[LOPOFS], "LO Phase Offset");
  sprintf(error_source[LOPJMP], "LO Phase Jump");
  sprintf(error_source[AMPERR], "Antenna Gain Bias Error");
  sprintf(srtaer,               "SRT Attitude Error");

  sprintf(error_source[ERROR_MENU_NUM],   "ALL ON");
  sprintf(error_source[ERROR_MENU_NUM+1], "ALL OFF");

/*
-------------------------------------------
*/

  ARRAY_NUM = 12;

  array_id[ 0].id = NO_ANT;
  array_id[ 1].id = VLBA;
  array_id[ 2].id = EVN;
  array_id[ 3].id = HSA;
  array_id[ 4].id = VERA;
  array_id[ 5].id = JVN;
  array_id[ 6].id = KVN;
  array_id[ 7].id = LBA;
  array_id[ 8].id = KAVA;
  array_id[ 9].id = EALMA;
  array_id[10].id = ALMA;
  array_id[11].id = ACA;

  array_id[ 0].type = NO_DEF_ARRAY;
  array_id[ 1].type = _VLBI_ARRAY_;
  array_id[ 2].type = _VLBI_ARRAY_;
  array_id[ 3].type = _VLBI_ARRAY_;
  array_id[ 4].type = _VLBI_ARRAY_;
  array_id[ 5].type = _VLBI_ARRAY_;
  array_id[ 6].type = _VLBI_ARRAY_;
  array_id[ 7].type = _VLBI_ARRAY_;
  array_id[ 8].type = _VLBI_ARRAY_;
  array_id[ 9].type = _VLBI_ARRAY_;
  array_id[10].type = __CONNECTED_;
  array_id[11].type = __CONNECTED_;

  sprintf(array_id[ 0].name, "RESET");
  sprintf(array_id[ 1].name, "VLBA");
  sprintf(array_id[ 2].name, "EVN");
  sprintf(array_id[ 3].name, "HSA");
  sprintf(array_id[ 4].name, "VERA");
  sprintf(array_id[ 5].name, "JVN");
  sprintf(array_id[ 6].name, "KVN");
  sprintf(array_id[ 7].name, "LBA");
  sprintf(array_id[ 8].name, "KaVA");
  sprintf(array_id[ 9].name, "EALMA");
  sprintf(array_id[10].name, "ALMA");
  sprintf(array_id[11].name, "ACA");

  for (iarray=0; iarray<ARRAY_NUM; iarray++) {
    array_id[iarray].flag = false;
  }

  sprintf(array_type[0].name, "VLBI");
  sprintf(array_type[1].name, "Connected Array");

/*
-------------------------------------------
*/

  sprintf(srcproc_name[0], "R.A.-DEC.");
  sprintf(srcproc_name[1], "\\gDR.A.-\\gDDEC.");
  sprintf(srcproc_name[2], "\\gD\\gh-P.A.");

  sprintf(srcproc_name[3], "Source-1 Position");
  sprintf(srcproc_name[4], "MID Point");

/*
-------------------------------------------
*/

  if (TV_SWT == false) {
    START_FLAG = true;

    printf("---- ERROR PROC MODE ----\n");
    for (i=0; i<ERROR_MENU_NUM; i++) {
      if (ERROR_FLAG[i] == true) {
        printf("<%2d> %s\n", i+1, error_source[i]);
      } else {
        printf(" %2d  %s\n", i+1, error_source[i]);
      }
    }
    printf("---- Current ERROR FLAG Status : ");
    for (i=0; i<ERROR_MENU_NUM; i++) {
      printf("%1d", ERROR_FLAG[i]);
    }
    printf("\n");
    printf("Input ERROR FLAG: ");
    if (fgets(string, nstr, stdin) == NULL) {
      printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
      return (__NG__);
    }
    for (i=0; i<ERROR_MENU_NUM; i++) {
      if (string[i] == '1') {
        ERROR_FLAG[i] = true;
      } else if (string[i] == '0') {
        ERROR_FLAG[i] = false;
      } else if (string[i] == '-') {
        ERROR_FLAG[i] = ERROR_FLAG[i];
      } else if (string[i] == '\n') {
        break;
      }
    }
    printf("---- New ERROR FLAG Status : ");
    for (i=0; i<ERROR_MENU_NUM; i++) {
      printf("%1d", ERROR_FLAG[i]);
    }
    printf("\n");

/*
------------------------------
*/

    printf("Which input type for the source positions?\n");
    in__src_proc(SRCPROC_MODE, &SEP_MODE, &POS_MODE);
    printf("1. R.A.       -   DEC.\n");
    printf("2. Delta_R.A. -   Delta_DEC.\n");
    printf("3. Sep_ang    -   P.A.\n");
    printf("(CR->%d)\n", SEP_MODE + 1);
    while (1) {
      printf("Input number: ");
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        break;
      } else {
        sscanf(string, "%d", &i);
        if        (i == 1) {
          SEP_MODE = SRC__RA__DEC;
          break;
        } else if (i == 2) {
          SEP_MODE = SRC_dRA_dDEC;
          break;
        } else if (i == 3) {
          SEP_MODE = SRC_SEP_POSA;
          break;
        }
      }
    }

    if (SEP_MODE == SRC_dRA_dDEC || SEP_MODE == SRC_SEP_POSA) {
      printf("Which position for the input?\n");
      printf("1. Source-1 Position\n");
      printf("2. MID Point\n");
      printf("(CR->%d)\n", POS_MODE + 1);
      while (1) {
        printf("Input number: ");
        if (fgets(string, nstr, stdin) == NULL) {
          printf("ERROR: OBS_PARAM_INPUT : Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] == '\n') {
          break;
        } else {
          sscanf(string, "%d", &i);
          if        (i == 1) {
            POS_MODE = SRC_POS1;
            break;
          } else if (i == 2) {
            POS_MODE = SRC_POS2;
            break;
          } else {
            printf("Invalid number as an input type. Select again: ");
          }
        }
      }
    }
    SRCPROC_MODE = out_src_proc(SEP_MODE, POS_MODE);

/*
------------------------------
*/

    if (SRCPROC_MODE == SRCPROC1 || SRCPROC_MODE == SRCPROC2) {
      printf("**** Source-1 (J2000) ****\n");
      printf("Input RA (CR->%s) : ", ch_src.tgt_ra);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.tgt_ra, string);
      }
      printf("Input DEC(CR->%s) : ", ch_src.tgt_dc);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.tgt_dc, string);
      }
      input_star_position(ch_src.tgt_ra, ch_src.tgt_dc,
                          &src[0].RA2k, &src[0].DC2k);

      printf("**** Source-2 (J2000) ****\n");
      printf("Input RA (CR->%s) : ", ch_src.ref_ra);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.ref_ra, string);
      }
      printf("Input DEC(CR->%s) : ", ch_src.ref_dc);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.ref_dc, string);
      }
      input_star_position(ch_src.ref_ra, ch_src.ref_dc,
                          &src[1].RA2k, &src[1].DC2k);

    } else if (SRCPROC_MODE == SRCPROC3) {
      printf("**** Source-1 (J2000) ****\n");
      printf("Input RA (CR->%s) : ", ch_src.tgt_ra);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.tgt_ra, string);
      }
      printf("Input DEC(CR->%s) : ", ch_src.tgt_dc);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.tgt_dc, string);
      }
      input_star_position(ch_src.tgt_ra, ch_src.tgt_dc,
                          &src[0].RA2k, &src[0].DC2k);

      printf("Separation Angle (RA) [deg] (CR->%lf) : ", pair_src.dlt_ra);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        sprintf(ch_src.dlt_ra, "%lf", pair_src.dlt_ra);
      } else {
        char_copy(ch_src.dlt_ra, string);
        sscanf(ch_src.dlt_ra, "%lf", &(pair_src.dlt_ra));
      }

      printf("Separation Angle (DC) [deg] (CR->%lf) : ", pair_src.dlt_dc);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        sprintf(ch_src.dlt_dc, "%lf", pair_src.dlt_dc);
      } else {
        char_copy(ch_src.dlt_dc, string);
        sscanf(ch_src.dlt_dc, "%lf", &(pair_src.dlt_dc));
      }

    } else if (SRCPROC_MODE == SRCPROC4) {
      printf("**** MID-POTINT OF SOURCES (J2000) ****\n");
      printf("Input RA (CR->%s) : ", ch_src.mid_ra);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.mid_ra, string);
      }

      printf("Input DEC (CR->%s) : ", ch_src.mid_dc);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.mid_dc, string);
      }

      input_star_position(ch_src.mid_ra, ch_src.mid_dc,
                          &pair_src.mid_ra, &pair_src.mid_dc);

      printf("Separation Angle (RA) [deg] (CR->%lf) : ", pair_src.dlt_ra);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        sprintf(ch_src.dlt_ra, "%lf", pair_src.dlt_ra);
      } else {
        char_copy(ch_src.dlt_ra, string);
        sscanf(ch_src.dlt_ra, "%lf", &(pair_src.dlt_ra));
      }

      printf("Separation Angle (DC) [deg] (CR->%lf) : ", pair_src.dlt_dc);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        sprintf(ch_src.dlt_dc, "%lf", pair_src.dlt_dc);
      } else {
        char_copy(ch_src.dlt_dc, string);
        sscanf(ch_src.dlt_dc, "%lf", &(pair_src.dlt_dc));
      }

    } else if (SRCPROC_MODE == SRCPROC5) {
      printf("**** Source-1 (J2000) ****\n");
      printf("Input RA (CR->%s) : ", ch_src.tgt_ra);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.tgt_ra, string);
      }
      printf("Input DEC(CR->%s) : ", ch_src.tgt_dc);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.tgt_dc, string);
      }
      input_star_position(ch_src.tgt_ra, ch_src.tgt_dc,
                          &src[0].RA2k, &src[0].DC2k);

      printf("Separation Angle [deg] (CR->%lf) : ", pair_src.sepang);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        sprintf(ch_src.sepang, "%lf", pair_src.sepang);
      } else {
        char_copy(ch_src.sepang, string);
        sscanf(ch_src.sepang, "%lf", &(pair_src.sepang));
      }

      printf("Position Angle   [deg] (CR->%lf) : ", pair_src.posang);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        sprintf(ch_src.posang, "%lf", pair_src.posang);
      } else {
        char_copy(ch_src.posang, string);
        sscanf(ch_src.posang, "%lf", &(pair_src.posang));
      }

    } else if (SRCPROC_MODE == SRCPROC6) {
      printf("**** MID-POTINT OF SOURCES (J2000) ****\n");
      printf("Input RA (CR->%s) : ", ch_src.mid_ra);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.mid_ra, string);
      }

      printf("Input DEC (CR->%s) : ", ch_src.mid_dc);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(ch_src.mid_dc, string);
      }

      input_star_position(ch_src.mid_ra, ch_src.mid_dc,
                          &pair_src.mid_ra, &pair_src.mid_dc);

      printf("Separation Angle [deg] (CR->%lf) : ", pair_src.sepang);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        sprintf(ch_src.sepang, "%lf", pair_src.sepang);
      } else {
        char_copy(ch_src.sepang, string);
        sscanf(ch_src.sepang, "%lf", &(pair_src.sepang));
      }

      printf("Position Angle   [deg] (CR->%lf) : ", pair_src.posang);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        sprintf(ch_src.posang, "%lf", pair_src.posang);
      } else {
        char_copy(ch_src.posang, string);
        sscanf(ch_src.posang, "%lf", &(pair_src.posang));
      }
    }
    source_position(src, &pair_src, &ch_src, SEP_MODE, POS_MODE);
    printf("Separation Angle [deg]: %s\n", ch_src.sepang);

/*
------------------------------
*/

    printf("Input Observation Date (YYYYMMDD) (CR->%4s%2s%2s) : ",
            ch_obs_t.start_t[0], ch_obs_t.start_t[1], ch_obs_t.start_t[2]);
    if (fgets(string, nstr, stdin) == NULL) {
      printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
      return (__NG__);
    }
    if (string[0] != '\n') {
      char_ncopy(ch_obs_t.start_t[0], string,   4);
      char_ncopy(ch_obs_t.start_t[1], string+4, 2);
      char_ncopy(ch_obs_t.start_t[2], string+6, 2);
    }

    printf("Input Start UTC (hhmmss) (CR->%s%s%s) : ",
            ch_obs_t.start_t[3], ch_obs_t.start_t[4], ch_obs_t.start_t[5]);
    if (fgets(string, nstr, stdin) == NULL) {
      printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
      return (__NG__);
    }
    if (string[0] != '\n') {
      char_ncopy(ch_obs_t.start_t[3], string,   2);
      char_ncopy(ch_obs_t.start_t[4], string+2, 2);
      char_ncopy(ch_obs_t.start_t[5], string+4, 2);
    }

    printf("Observing Duration [hour] (CR->%s) : ", ch_obs_t.obsd);
    if (fgets(string, nstr, stdin) == NULL) {
      printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
      return (__NG__);
    }
    if (string[0] != '\n') {
      char_copy(ch_obs_t.obsd, string);
      sscanf(ch_obs_t.obsd, "%lf", obs_duration);
    }

/*
------------------------------
*/

    while (1) {
      printf("GRT MINIMUM elevation (CR->%s [deg]) : ", ch_grt_el_lim);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        sscanf(string, "%lf", grt_elevation_limit);
        char_copy(ch_grt_el_lim, string);
        if (*grt_elevation_limit < 5.0) {
          printf("WARNING: The lowest elevation angle must be above 5 deg.\n");
          printf("WARNING: Limit elevation is set to 5 deg.\n");
          *grt_elevation_limit = 5.0;
          sprintf(ch_grt_el_lim, "5.0");
          break;
        } else if (*grt_elevation_limit >= 90.0) {
          printf("WARNING: Wrong number for the limit elevation.\n");
        } else {
          break;
        }
      } else {
        break;
      }
    }

/*
------------------------------
*/

    while (1) {
      printf("Which array type?\n");
      printf("1. VLBI (indepndent atmosphere and frequency standard)\n");
      printf("2. Conneced array\n");
      printf("(CR->%d) : ", array->TYPE);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == '\n') {
        break;
      } else {
        sscanf(string, "%d", &i);
        if        (i == 1) {
          array->TYPE = _VLBI_ARRAY_;
          printf("Array type: VLBI\n");
          break;
        } else if (i == 2) {
          array->TYPE = __CONNECTED_;
          printf("Array type: Connected Array\n");
          break;
        } else {
          printf("WARNING: Wrong number for the array type.\n");
        }
      }
    }

/*
------------------------------
*/

    while (1) {
      printf("Antenna list file     (CR->%s) : ", antenna_list_file);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT :");
        printf("There is something wrong in standard input.\n");
        return (__NG__);
      }
      if (string[0] != '\n') {
        char_copy(antenna_list_file_tmp, string);
        if (ant_list_chk(antenna_list_file_tmp, bttn_box[I=ANTLST_SECTION+2],
                         cmnt, comment, false) == 1) {
          if (strcmp(antenna_list_file, antenna_list_file_tmp) != 0) {
            char_copy(antenna_list_file, antenna_list_file_tmp);
            ANT_NUM_tmp = array_config(ALL_ANT, wave_id, *SRT_NUM, &GRT_NUM_tmp,
                                       ant_prm,   "",
                                       antenna_list_file, false,    true);
            for (iant=0; iant<ANT_NUM_tmp; iant++) {
              ant_prm[iant].UFL = false;
            }

            for (j=0; j<*GRT_NUM; j++) {
              for (iant=0; iant<GRT_NUM_tmp; iant++) {
                if (strncmp(ant_code[j], (ant_prm+iant)->IDC,
                            strlen(ant_code[j])) == 0) {
                  ant_prm[iant].UFL = true;
                  break;
                }
              }
            }
          }
        }
      }
      break;
    }

/*
------------------------------
*/

    printf("\n");
    while (1) {
      printf("SELECTED Antennas:\n");
      for (i=0; i<GRT_NUM_tmp; i++) {
        if (ant_prm[i].UFL == true) {
          printf("%3d. %s", i+1, ant_prm[i].IDC);
          k = strlen(ant_prm[i].IDC);
          if (k > 10) {
            k = 10;
          }
          for (j=0; j<11-k; j++) {
            printf(" ");
          }
          if (i % 5 == 4) {
            printf("\n");
          }
        } else {
          printf("%3d.            ", i+1);
          if (i % 5 == 4) {
            printf("\n");
          }
        }
      }
      printf("\n");

      printf("Array List:\n");
      for (iarray=0; iarray<ARRAY_NUM; iarray++) {
        string[0] = ' ';
        string[1] = ' ';
        string[2] = 'a' + iarray;
        string[3] = 0;

        printf("%s. %s", string, array_id[iarray].name);
        k = strlen(array_id[iarray].name);
        if (k > 10) {
          k = 10;
        }
        for (j=0; j<11-k; j++) {
          printf(" ");
        }
        if (iarray % 5 == 4) {
          printf("\n");
        }
      }
      printf("\n");

      printf("Antenna List:\n");
      for (i=0; i<GRT_NUM_tmp; i++) {
        if (ant_prm[i].UFL == false) {
          printf("%3d. %s", i+1, ant_prm[i].IDC);
          k = strlen(ant_prm[i].IDC);
          if (k > 10) {
            k = 10;
          }
          for (j=0; j<11-k; j++) {
            printf(" ");
          }
          if (i % 5 == 4) {
            printf("\n");
          }
        } else {
          printf("%3d.            ", i+1);
          if (i % 5 == 4) {
            printf("\n");
          }
        }
      }
      printf("\n");

      printf("Input array character or antenna number (0 or CR->Exit) : ");
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standard Input may have a problem.\n");
        return (__NG__);
      }

      if (string[0] == '0' || string[0] == '\n') {
        break;
      } else if (string[0] >= 'a' && string[0] <= 'a' + ARRAY_NUM) {
        array->ID = string[0] - 'a';
        if (array->ID == 0) {
          for (i=0; i<GRT_NUM_tmp; i++) {
            ant_prm[i].UFL = false;
          }
        } else {
          for (i=0; i<GRT_NUM_tmp; i++) {
            idum = array_config(array->ID, wave_id, 0, &grt_num,
                                &ant_prm_tmp, ant_prm[i].IDC,
                                antenna_list_file,  true, true);
            if (idum == 1 &&
                strncmp(ant_prm_tmp.IDC,
                        ant_prm[i].IDC, strlen(ant_prm[i].IDC)) == 0) {
              ant_prm[i].UFL = true;
            }
          }
        }
      } else {
        sscanf(string, "%d", &k);
        k--;
        if (ant_prm[k].UFL == false) {
          ant_prm[k].UFL = true;
        } else if (ant_prm[k].UFL == true) {
          ant_prm[k].UFL = false;
        }
      }
    }

/*
------------------------------
*/

    if (array->ID != ACA) {
      idum = *SRT_NUM;
      while (1) {
        printf("How many space telescopes [0-%d] (CR->%d) : ",
               SRTMAX, *SRT_NUM);
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          sscanf(string, "%d", &idum);
        }
        if (idum < 0 || idum > SRTMAX) {
          printf("Input the number btween 0 and %d: ", SRTMAX);
        } else {
          break;
        }
      }
      *SRT_NUM = idum;
    }

    if (*SRT_NUM > 0) {
      if (ERROR_FLAG[SRTAER] == true) {
        printf("SRT pointing error (y/n; CR->y) : ");
      } else {
        printf("SRT pointing error (y/n; CR->n) : ");
      }
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
        return (__NG__);
      }
      if (string[0] == 'y') {
        ERROR_FLAG[SRTAER] = true;
      } else if (string[0] == 'n') {
        ERROR_FLAG[SRTAER] = false;
      }

      for (i=0; i<*SRT_NUM; i++) {
        printf("Input orbit parameters for SRT %d\n", i + 1);
        printf("Apogee altitude  [km] (CR->%s km) : ", ch_srt[i].apo);
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          char_copy(ch_srt[i].apo, string);
        }

        printf("Perigee altitude [km] (CR->%s km) : ", ch_srt[i].per);
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          char_copy(ch_srt[i].per, string);
        }

        printf("Inclination     [deg] (CR->%s deg) : ", ch_srt[i].inc);
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          char_copy(ch_srt[i].inc, string);
        }

        printf("Omega           [deg] (CR->%s deg) : ", ch_srt[i].OMG);
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          char_copy(ch_srt[i].OMG, string);
        }

        printf("omega           [deg] (CR->%s deg) : ", ch_srt[i].omg);
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          char_copy(ch_srt[i].omg, string);
        }

        printf("time of the passage at perigee (YYYYMMDDhhmmss)\n");
        printf("(CR->%s) : ", ch_srt[i].t_0); 
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          char_copy(ch_srt[i].t_0, string);
        }

        printf("d_Omega/dt[deg/yr] (CR->%s deg/yr) : ", ch_srt[i].d_OMG);
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          char_copy(ch_srt[i].d_OMG, string);
        }

        printf("d_omega/dt[deg/yr] (CR->%s deg/yr) : ", ch_srt[i].d_omg);
        if (fgets(string, nstr, stdin) == NULL) {
          printf(
            "ERROR: OBS_PARAM_INPUT : Standar Input may have a problem.\n");
          return (__NG__);
        }
        if (string[0] != '\n') {
          char_copy(ch_srt[i].d_omg, string);
        }
      }
    }

/*
============================================================
*/

  } else if (TV_SWT == true) {

/*****
    cpgbbuf();

    cpgslct(pgid[0]);
    cpgpap(2.00*pgpap_prm, 1.0/1.4);
    cpgsch(1.5*pgpap_prm/13.0);
    cpgsvp(0.0,  1.0, 0.0, 1.0);
    cpgswin(0.0, 1.4, 0.0, 1.0);
    comment_init(cmnt, comment, true);

    sprintf(string,
     "ARIS: Observation Parameter Setup --GRAPHICAL USER INTERFACE MODE--");
    comment_disp(cmnt, comment, string, true);
    sprintf(string, "Welcome!! Please set the observing parameters.");
    comment_disp(cmnt, comment, string, true);
****/

/*
-----
*/

    y_pos = 0.190;

    I = 0;
    bttn_box[I][0] = 0.02;
    bttn_box[I][1] = 0.20;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    off_button(&idum, "START\0", bttn_box[I]);

    I = 1;
    bttn_box[I][0] = 0.26;
    bttn_box[I][1] = 0.44;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    off_button(&idum, "ANTENNA VISIBILITY\0", bttn_box[I]);

    I = 2;
    bttn_box[I][0] = 0.45;
    bttn_box[I][1] = 0.63;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    off_button(&idum, "PHASE SCREEN\0", bttn_box[I]);

    I = 3;
    bttn_box[I][0] = 0.80;
    bttn_box[I][1] = 0.98;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    off_button(&idum, "EXIT\0", bttn_box[I]);

/*
-----
*/

    y_pos = 0.910;

    cpgsci(1);
    cpgsfs(2);
    cpgrect(0.010, 1.400, y_pos-0.040, y_pos+0.080);
    cpgsfs(1);

    y_pos = 0.940;

    for (i=0; i<ERROR_MENU_NUM; i++) {
      I = ERROR_SECTION + i;
      bttn_box[I][0] = 0.015 + 0.230 * (float)((i+3) % 6);
      bttn_box[I][1] = bttn_box[I][0] + 0.230;
      bttn_box[I][2] = y_pos - 0.033 * (float)((i+3) / 6);
      bttn_box[I][3] = bttn_box[I][2] + pitch;
    }
    y_pos = 0.910;
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM  ][0] = 0.015;
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM  ][1] = 0.280;
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM  ][2] = y_pos + 0.040;
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM  ][3] = y_pos + 0.040 + pitch;

    bttn_box[ERROR_SECTION+ERROR_MENU_NUM+1][0] =
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM  ][0] + 0.27;
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM+1][1] =
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM  ][1] + 0.27;
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM+1][2] =
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM  ][2];
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM+1][3] =
    bttn_box[ERROR_SECTION+ERROR_MENU_NUM  ][3];

    for (i=0; i<ERROR_MENU_NUM; i++) {
      I = ERROR_SECTION + i;
      if (ERROR_FLAG[i] == true) {
        on_button( &ALL_ERROR_FLAG, error_source[i], bttn_box[I]);
      } else if (ERROR_FLAG[i] == false) {
        off_button(&ALL_ERROR_FLAG, error_source[i], bttn_box[I]);
      }
    }
    for (i=ERROR_MENU_NUM; i<ERROR_MENU_NUM+2; i++) {
      I = ERROR_SECTION + i;
      off_button(&ALL_ERROR_FLAG, error_source[i], bttn_box[I]);
    }

/*
----
*/

    y_pos -= ((float)(ERROR_MENU_NUM / 6) * 0.030 + 0.008);

    cpgsci(1);
    cpgsfs(2);
    cpgrect(0.010, 1.400, y_pos-0.111, y_pos+0.020);
    cpgsfs(1);

    in__src_proc(SRCPROC_MODE, &SEP_MODE, &POS_MODE);
    for (i=0; i<BASIC_SEP_MODE; i++) {
      I = SRCPROC_SECTION + i;
      bttn_box[I][0] = 0.02 + 0.2 * (float)i;
      bttn_box[I][1] = 0.03 + 0.2 * (float)i;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = y_pos + 0.01;
      if (i == SEP_MODE) {
        _on_button( &srcproc_flag[i], "", bttn_box[I]);
      } else {
        _off_button(&srcproc_flag[i], "", bttn_box[I]);
      }
      cpgsci(1);
      cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[i]);
    }
    y_pos -= 0.020;
    for (i=0; i<2; i++) {
      I = SRCPROC_SECTION + BASIC_SEP_MODE + i;
      bttn_box[I][0] = 0.02 + 0.2 * (float)i;
      bttn_box[I][1] = 0.03 + 0.2 * (float)i;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = y_pos + 0.01;
    }

/****
  SRCPROC_MODE = 0;
  SEP_MODE = 0;
  POS_MODE = 0;
****/

    if (SEP_MODE == SRC_dRA_dDEC || SEP_MODE == SRC_SEP_POSA) {
      if (       POS_MODE == SRC_POS1) {
        I = SRCPROC_SECTION + BASIC_SEP_MODE;
        _on_button (&srcproc_flag[BASIC_SEP_MODE], "", bttn_box[I]);
        cpgsci(1);
        cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[BASIC_SEP_MODE]);
        I++;
        _off_button(&srcproc_flag[BASIC_SEP_MODE+1], "", bttn_box[I]);
        cpgsci(1);
        cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[BASIC_SEP_MODE+1]);
      } else if (POS_MODE == SRC_POS2) {
        I = SRCPROC_SECTION + BASIC_SEP_MODE;
        _off_button(&srcproc_flag[BASIC_SEP_MODE], "", bttn_box[I]);
        cpgsci(1);
        cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[BASIC_SEP_MODE]);
        I++;
        _on_button (&srcproc_flag[BASIC_SEP_MODE+1], "", bttn_box[I]);
        cpgsci(1);
        cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[BASIC_SEP_MODE+1]);
      }
    }

/*
--------
*/

    source_y_pos = y_pos - 0.045;
    sprintf(string, "Separation Angle [deg]: %5.2lf",
            sepang(src[0].s2k, src[1].s2k) * 180.0 / dpi);
    source_info_disp(SEP_MODE, POS_MODE, pitch,  &bttn_box[SOURCE_SECTION],
                     source_y_pos,   &ch_src);

/*
----
*/

    y_pos = source_y_pos - 0.095;

    cpgsci(1);
    cpgsfs(2);
    cpgrect(0.010, 1.400, y_pos-0.005, y_pos+0.040);
    cpgsfs(1);

    I = TIME_SECTION;
    bttn_box[I][0] = 0.21;
    bttn_box[I][1] = 0.28;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    cpgptxt(bttn_box[I][0], text_bottom(bttn_box[I][2], bttn_box[I][3]),
            0.0, 1.05, "Obs Date (YYYYMMDD)");
    tv_button_disp(bttn_box[I], ch_obs_t.start_t[0]);
    I++;

    bttn_box[I][0] = 0.29;
    bttn_box[I][1] = 0.34;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    tv_button_disp(bttn_box[I], ch_obs_t.start_t[1]);
    I++;

    bttn_box[I][0] = 0.35;
    bttn_box[I][1] = 0.40;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    tv_button_disp(bttn_box[I], ch_obs_t.start_t[2]);
    I++;

    bttn_box[I][0] = 0.60;
    bttn_box[I][1] = 0.65;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    cpgptxt(bttn_box[I][0], text_bottom(bttn_box[I][2], bttn_box[I][3]),
            0.0, 1.05, "Start UTC (hhmmss)");
    tv_button_disp(bttn_box[I],  ch_obs_t.start_t[3]);
    cpgsci(1);
    I++;

    bttn_box[I][0] = 0.66;
    bttn_box[I][1] = 0.71;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    tv_button_disp(bttn_box[I],  ch_obs_t.start_t[4]);
    cpgsci(1);
    I++;

    bttn_box[I][0] = 0.72;
    bttn_box[I][1] = 0.77;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    tv_button_disp(bttn_box[I],  ch_obs_t.start_t[5]);
    cpgsci(1);
    I++;

    bttn_box[I][0] = 0.96;
    bttn_box[I][1] = 1.03;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    cpgptxt(bttn_box[I][0], text_bottom(bttn_box[I][2], bttn_box[I][3]),
            0.0, 1.05, "Obs Duration [h]");
    tv_button_disp(bttn_box[I], ch_obs_t.obsd);
    cpgsci(1);
    I++;

/*
----
*/

    I = GRT_SECTION;
    bttn_box[I][0] = 1.29;
    bttn_box[I][1] = 1.37;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    cpgptxt(bttn_box[I][0], text_bottom(bttn_box[I][2], bttn_box[I][3]),
            0.0, 1.05, "GRT Elevation Limit [deg]");
    tv_button_disp(bttn_box[I], ch_grt_el_lim);
    cpgsci(1);

/*
----
*/

    I = ANTLST_SECTION;
    y_pos -= 0.043;
    cpgsci(1);
    for (i=0; i<2; i++) {
      I = ANTLST_SECTION + i;
      bttn_box[I][0] = 0.02 + 0.10 * (float)i;
      bttn_box[I][1] = 0.03 + 0.10 * (float)i;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = y_pos + 0.01;

      if (array_type[i].flag == true) {
        _on_button( &array_type[i].flag, "", bttn_box[I]);
      } else if (array_type[i].flag == false) {
        _off_button(&array_type[i].flag, "", bttn_box[I]);
      }
      cpgsci(1);
      cpgptxt(bttn_box[I][1]+0.01, text_bottom(bttn_box[I][2], bttn_box[I][3]),
            0.0, 0.0, array_type[i].name);
    }

/*
----
*/

    I = ANTLST_SECTION + 2;
    y_pos -= 0.006;
    bttn_box[I][0] = 0.67;
    bttn_box[I][1] = 1.07;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    cpgptxt(bttn_box[I][0], text_bottom(bttn_box[I][2], bttn_box[I][3]),
            0.0, 1.05, "Antenna list file: ");
    tv_button_disp(bttn_box[I], antenna_list_file);
    cpgsci(1);
    ant_list_chk(antenna_list_file, bttn_box[I], cmnt, comment, TV_SWT);

/*
----
*/

    SRT_SWT = false;
    I = ARRAY_SECTION;
    y_pos -= 0.035;
    for (i=0; i<ARRAY_NUM; i++) {
      I = ARRAY_SECTION + i;
      bttn_box[I][0] = 0.020 + 0.075 * (float)i;
      bttn_box[I][1] = 0.090 + 0.075 * (float)i;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
    }

    for (iarray=0; iarray<ARRAY_NUM; iarray++) {
      array_id[iarray].flag = false;
      I = ARRAY_SECTION + iarray;
      _off_button(&array_id[iarray].flag, array_id[iarray].name, bttn_box[I]);
    }

/*
----
*/

    y_pos -= 0.035;
    for (i=0; i<GRT_NUM_tmp; i++) {
      ant_bttn_box[i][0] = 0.023 + 0.097 * (float)(i % 14);
      ant_bttn_box[i][1] = 0.120 + 0.097 * (float)(i % 14);
      ant_bttn_box[i][2] = y_pos - pitch * (float)(i / 14);
      ant_bttn_box[i][3] = ant_bttn_box[i][2] + pitch;
    }
    for (i=0; i<GRT_NUM_tmp; i++) {
      if (ant_prm[i].UFL == true) {
        _on_button (&ant_prm[i].UFL, &ant_prm[i].IDC[0], ant_bttn_box[i]);
      } else {
        _off_button(&ant_prm[i].UFL, &ant_prm[i].IDC[0], ant_bttn_box[i]);
      }
    }

    cpgsci(1);
    cpgsfs(2);
    cpgrect(0.010, 1.400, bttn_box[ARRAY_SECTION][3]+0.040, 0.345);
    cpgsci(1);
    cpgsfs(1);

/*
----
*/

    if (array->ID == ACA || array->ID == ALMA) {
      SRT_SWT = false;
    } else {
      SRT_SWT = true;
    }

/*
----
*/

    y_pos = 0.303;
    I = SRT_SECTION;
    bttn_box[I][0] = 0.195;
    bttn_box[I][1] = 0.255;
    bttn_box[I][2] = ant_bttn_box[GRT_NUM_tmp - 1][2] - 0.05;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    if (SRT_SWT == true) {
      sprintf(string, "SRT Number (0-%d)", SRTMAX);
      cpgsci(1);
      cpgtext(0.02, text_bottom(bttn_box[I][2], bttn_box[I][3]), string);
      sprintf(string, "%d", *SRT_NUM);
      tv_button_disp(bttn_box[I], string);
    }

/*
----
*/

    srt_info_disp(*SRT_NUM,  pitch, bttn_box+SRT_SECTION+1,
                  bttn_box[SRT_SECTION][2],
                  srtaer, &ERROR_FLAG[SRTAER], ch_srt);

/*
----
*/

    cursor_pos[0] = 0.5 * (bttn_box[0][0] + bttn_box[0][1]);
    cursor_pos[1] = 0.5 * (bttn_box[0][2] + bttn_box[0][3]);

    while (1) {

      cpgcurs(cursor_pos, cursor_pos+1, string);

/*
---------------------------------
*/

      if (_button_chk(cursor_pos, bttn_box[0]) == true) {
        I = 0;
        on_button(&idum, "START\0", bttn_box[I]);
        START_FLAG = true;
        break;
      }

      if (_button_chk(cursor_pos, bttn_box[1]) == true) {
        I = 1;
        on_button(&idum, "ANTENNA VISIBILITY\0", bttn_box[I]);
        START_FLAG = false;

        if (obs_param_set(ERROR_FLAG, *SRT_NUM,
                          ch_srt, srt,
                          ch_grt_el_lim, grt_elevation_limit,
                          sep_angle_limit_from_earth_limb,
                          &ch_obs_t, TimUTC, UT1_UTC,
                          obs_duration, SRCPROC_MODE,
                          &ch_src, &pair_src, src, sun) == -1) {
          printf("WARNING: OBS_PARAM_INPUT: SUN ANGLE CHECK.\n");
        }
        *GRT_NUM = 0;
        for (iant=0; iant<GRT_NUM_tmp; iant++) {
          if (ant_prm[iant].UFL == true) {
            sprintf(ant_code[*GRT_NUM], "%s", (ant_prm+iant)->IDC);
            (*GRT_NUM)++;
          }
        }
        obs_param_file_io(ERROR_FLAG, antenna_list_file,
                          array, ANT_NUM, GRT_NUM, SRT_NUM,
                          srt, grt_elevation_limit,
                          sep_angle_limit_from_earth_limb,
                          TimUTC, UT1_UTC, obs_duration,
                          &SRCPROC_MODE, src, sun, ant_code,
                          ch_grt_el_lim, &pair_src, &ch_src,
                          ch_srt, &ch_obs_t, 1);
        return (ANT_VIS);
      }

      if (_button_chk(cursor_pos, bttn_box[2]) == true) {
        I = 2;
        on_button(&idum, "PHASE SCREEN\0", bttn_box[I]);
        START_FLAG = false;
        return (PHS_SCR);
      }

      if (_button_chk(cursor_pos, bttn_box[3]) == true) {
        I = 3;
        on_button(&idum, "EXIT\0", bttn_box[I]);
        START_FLAG = false;
        break;
      }

/*
---------------------------------
*/

      for (i=0; i<ERROR_MENU_NUM+2; i++) {
        I = ERROR_SECTION + i;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          if (i < ERROR_MENU_NUM) {
            if (ERROR_FLAG[i] == true) {
              idum = 1;
            } else {
              idum = 0;
            }
            toggle_button(&idum, error_source[i], bttn_box[I]);
            if (idum == 1) {
              ERROR_FLAG[i] = true;
            } else {
              ERROR_FLAG[i] = false;
            }
            break;
          } else if (i == ERROR_MENU_NUM) {
            for (j=0; j<ERROR_MENU_NUM; j++) {
              J = ERROR_SECTION + j;
              on_button(&idum, error_source[j], bttn_box[J]);
              ERROR_FLAG[j] = true;
            }
            break;
          } else if (i == ERROR_MENU_NUM+1) {
            for (j=0; j<ERROR_MENU_NUM; j++) {
              J = ERROR_SECTION + j;
              off_button(&idum, error_source[j], bttn_box[J]);
              ERROR_FLAG[j] = false;
            }
            break;
          }
        }
      }

/*
-----------------------------------------
*/

      SRCPROC_SWT = false;
      for (i=0; i<3; i++) {
        I = SRCPROC_SECTION + i;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          SEP_MODE = i;
          SRCPROC_SWT = true;
          break;
        }
      }

      if (SEP_MODE == SRC_dRA_dDEC || SEP_MODE == SRC_SEP_POSA) {
        for (i=0; i<2; i++) {
          I = SRCPROC_SECTION + BASIC_SEP_MODE + i;
          if (_button_chk(cursor_pos, bttn_box[I]) == true) {
            POS_MODE = i;
            SRCPROC_SWT = true;
          }
        }
      }

      if (SRCPROC_SWT == true) {
        for (i=0; i<3; i++) {
          I = SRCPROC_SECTION + i;
          if (i == SEP_MODE) {
            _on_button( &srcproc_flag[i], "", bttn_box[I]);
          } else {
            _off_button(&srcproc_flag[i], "", bttn_box[I]);
          }
        }

        y_pos = source_y_pos + 0.045;
        if (SEP_MODE == SRC__RA__DEC) {
          TV_menu_hatch(0.011, 0.550, y_pos, y_pos+0.012, 0, 1);
        } else {
          if (       POS_MODE == SRC_POS1) {
            I = SRCPROC_SECTION + BASIC_SEP_MODE;
            _on_button (&srcproc_flag[BASIC_SEP_MODE], "", bttn_box[I]);
            cpgsci(1);
            cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[BASIC_SEP_MODE]);
            I++;
            _off_button(&srcproc_flag[BASIC_SEP_MODE+1], "", bttn_box[I]);
            cpgsci(1);
            cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[BASIC_SEP_MODE+1]);
          } else if (POS_MODE == SRC_POS2) {
            I = SRCPROC_SECTION + BASIC_SEP_MODE;
            _off_button(&srcproc_flag[BASIC_SEP_MODE], "", bttn_box[I]);
            cpgsci(1);
            cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[BASIC_SEP_MODE]);
            I++;
            _on_button (&srcproc_flag[BASIC_SEP_MODE+1], "", bttn_box[I]);
            cpgsci(1);
            cpgtext(bttn_box[I][0]+0.02, y_pos, srcproc_name[BASIC_SEP_MODE+1]);
          }
        }

        SRCPROC_MODE = out_src_proc(SEP_MODE, POS_MODE);
        TV_menu_hatch(0.011, 1.00, bttn_box[SOURCE_SECTION+1][2],
                                  bttn_box[SOURCE_SECTION+2][3], 0, 1);
        source_info_disp(SEP_MODE, POS_MODE, pitch, &bttn_box[SOURCE_SECTION],
                         source_y_pos, &ch_src);
        SRCPROC_SWT = false;
      }

/*
-----------------------------------------
*/

      SEPANG_UPDATE = false;
      if (_button_chk(cursor_pos, bttn_box[I=SOURCE_SECTION]) == true) {
        if (SEP_MODE == SRC__RA__DEC ||
           (SEP_MODE != SRC__RA__DEC && POS_MODE == SRC_POS1)) {
          char_copy(string, ch_src.tgt_ra);
          tv_get_param("char", cursor_pos, bttn_box[I],
                       pitch, string, 0.0, 0.0);
          char_copy(ch_src.tgt_ra, string);
          str_init(string, nstr);
          input_star_position(ch_src.tgt_ra, ch_src.tgt_dc,
                              &src[0].RA2k, &src[0].DC2k);
        } else {
          char_copy(string, ch_src.mid_ra);
          tv_get_param("char", cursor_pos, bttn_box[I],
                       pitch, string, 0.0, 0.0);
          char_copy(ch_src.mid_ra, string);
          str_init(string, nstr);
          input_star_position(ch_src.mid_ra, ch_src.mid_dc,
                              &pair_src.mid_ra, &pair_src.mid_dc);
        }
        SEPANG_UPDATE = true;
      }

      if (_button_chk(cursor_pos, bttn_box[I=SOURCE_SECTION+1]) == true) {
        if (SEP_MODE == SRC__RA__DEC ||
           (SEP_MODE != SRC__RA__DEC && POS_MODE == SRC_POS1)) {
          char_copy(string, ch_src.tgt_dc);
          tv_get_param("char", cursor_pos, bttn_box[I],
                       pitch, string, 0.0, 0.0);
          char_copy(ch_src.tgt_dc, string);
          str_init(string, nstr);
          input_star_position(ch_src.tgt_ra, ch_src.tgt_dc,
                              &src[0].RA2k, &src[0].DC2k);
        } else {
          char_copy(string, ch_src.mid_dc);
          tv_get_param("char", cursor_pos, bttn_box[I],
                       pitch, string, 0.0, 0.0);
          char_copy(ch_src.mid_dc, string);
          str_init(string, nstr);
          input_star_position(ch_src.mid_ra, ch_src.mid_dc,
                              &pair_src.mid_ra, &pair_src.mid_dc);
        }
        SEPANG_UPDATE = true;
      }

      if (_button_chk(cursor_pos, bttn_box[I=SOURCE_SECTION+2]) == true) {
        if (SEP_MODE == SRC__RA__DEC) {
          char_copy(string, ch_src.ref_ra);
          tv_get_param("char", cursor_pos, bttn_box[I],
                       pitch, string, 0.0, 0.0);
          char_copy(ch_src.ref_ra, string);
          str_init(string, nstr);
          input_star_position(ch_src.ref_ra, ch_src.ref_dc,
                              &src[1].RA2k, &src[1].DC2k);
        } else if (SEP_MODE == SRC_dRA_dDEC) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_src.dlt_ra,
                      -(double)SEPANG_LIM, (double)SEPANG_LIM);
          sscanf(ch_src.dlt_ra, "%lf", &pair_src.dlt_ra);
          sprintf(string, "Delta RA is set to %lf [deg].", pair_src.dlt_ra);
          comment_disp(cmnt, comment, string, true);
        } else if (SEP_MODE == SRC_SEP_POSA) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_src.sepang, 0.0, (double)SEPANG_LIM);
          sscanf(ch_src.sepang, "%lf", &pair_src.sepang);
          sprintf(string, "Separation angle is set to %lf [deg].",
                  pair_src.sepang);
          comment_disp(cmnt, comment, string, true);
        }
        SEPANG_UPDATE = true;
      }

      if (_button_chk(cursor_pos, bttn_box[I=SOURCE_SECTION+3]) == true) {
        if (SEP_MODE == SRC__RA__DEC) {
          char_copy(string, ch_src.ref_dc);
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, string, 0.0, 0.0);
          char_copy(ch_src.ref_dc, string);
          str_init(string, nstr);
          input_star_position(ch_src.ref_ra, ch_src.ref_dc,
                              &src[1].RA2k, &src[1].DC2k);
        } else if (SEP_MODE == SRC_dRA_dDEC) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_src.dlt_dc, -10.0, 10.0);
          sscanf(ch_src.dlt_dc, "%lf", &pair_src.dlt_dc);
          sprintf(string, "Delta DEC is set to %lf [deg].", pair_src.dlt_dc);
          comment_disp(cmnt, comment, string, true);
        } else if (SEP_MODE == SRC_SEP_POSA) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_src.posang, 0.0, 0.0);
          sscanf(ch_src.posang, "%lf", &pair_src.posang);
          sprintf(string, "Position angle is set to %lf [deg].",
                  pair_src.posang);
          comment_disp(cmnt, comment, string, true);
        }
        SEPANG_UPDATE = true;
      }

      if (SEPANG_UPDATE == true) {
        source_position(src, &pair_src, &ch_src, SEP_MODE, POS_MODE);
        source_info_disp(SEP_MODE,  POS_MODE,  pitch,
                         &bttn_box[SOURCE_SECTION], source_y_pos, &ch_src);
        SEPANG_UPDATE = false;
      }

/*
-----------------------------------------
*/

      if (_button_chk(cursor_pos, bttn_box[I=GRT_SECTION]) == true) {
        tv_get_param("double", cursor_pos, bttn_box[I],
                     pitch, ch_grt_el_lim, 0.0, 0.0);
        sscanf(ch_grt_el_lim, "%lf", grt_elevation_limit);
        sprintf(string, "GRT minimum elevation is set to %lf [deg].",
                *grt_elevation_limit);
        comment_disp(cmnt, comment, string, true);
      }

/*
-----------------------------------------
*/

      for (i=ANTLST_SECTION; i<ANTLST_SECTION+2; i++) {
        if (_button_chk(cursor_pos, bttn_box[i]) == true) {
          I = i - ANTLST_SECTION;
          if (array->TYPE == _VLBI_ARRAY_ && I == 1) {
            _toggle_button(&array_type[0].flag, "", bttn_box[ANTLST_SECTION  ]);
            _toggle_button(&array_type[1].flag, "", bttn_box[ANTLST_SECTION+1]);
            array->TYPE = __CONNECTED_;
          } else if (array->TYPE == __CONNECTED_ && I == 0) {
            _toggle_button(&array_type[0].flag, "", bttn_box[ANTLST_SECTION  ]);
            _toggle_button(&array_type[1].flag, "", bttn_box[ANTLST_SECTION+1]);
            array->TYPE = _VLBI_ARRAY_;
          }
          obs_param_file_io(ERROR_FLAG, antenna_list_file,
                    array, ANT_NUM, GRT_NUM, SRT_NUM,
                    srt, grt_elevation_limit,
                    sep_angle_limit_from_earth_limb,
                    TimUTC, UT1_UTC, obs_duration,
                    &SRCPROC_MODE, src, sun, ant_code,
                    ch_grt_el_lim, &pair_src, &ch_src,
                    ch_srt, &ch_obs_t, 1);
          return (RE_LOAD);
        }
      }

/*
-----------------------------------------
*/

      if (_button_chk(cursor_pos, bttn_box[I=ANTLST_SECTION+2]) == true) {
        sprintf(string, "%s", antenna_list_file);
        tv_get_param("char", cursor_pos, bttn_box[I],
                     pitch, string, 0.0, 0.0);
        tv_button_disp(bttn_box[I], string);
        char_copy(antenna_list_file_tmp, string);
        if (ant_list_chk(antenna_list_file_tmp, bttn_box[I=ANTLST_SECTION+2],
                         cmnt, comment, true) == 1) {

          if (strcmp(antenna_list_file, antenna_list_file_tmp) != 0) {
            char_copy(antenna_list_file, antenna_list_file_tmp);

            cpgsci(0);
            cpgsfs(1);
            cpgrect(0.009, 1.402, bttn_box[ARRAY_SECTION][3]+0.001, 0.344);
            cpgsci(1);
            cpgsfs(1);

            cpgsci(1);
            cpgsfs(2);
            cpgrect(0.010, 1.400, bttn_box[ARRAY_SECTION][3]+0.040, 0.345);
            cpgsci(1);
            cpgsfs(1);

            ANT_NUM_tmp = array_config(ALL_ANT, wave_id, *SRT_NUM, &GRT_NUM_tmp,
                                       ant_prm,   "",
                                       antenna_list_file, false,  true);
            for (iant=0; iant<ANT_NUM_tmp; iant++) {
              ant_prm[iant].UFL = false;
            }

            for (j=0; j<*GRT_NUM; j++) {
              for (iant=0; iant<GRT_NUM_tmp; iant++) {
                if (strncmp(ant_code[j], (ant_prm+iant)->IDC,
                            strlen(ant_code[j])) == 0) {
                  ant_prm[iant].UFL = true;
                  break;
                }
              }
            }
            for (iarray=0; iarray<ARRAY_NUM; iarray++) {
              array_id[iarray].flag = false;
              I = ARRAY_SECTION + iarray;
              _off_button(&array_id[iarray].flag, array_id[iarray].name,
                          bttn_box[I]);
            }
            array_config(ALL_ANT, wave_id, *SRT_NUM, GRT_NUM,
                         ant_prm,   "", antenna_list_file, false,  true);
            for (iant=0; iant<*GRT_NUM; iant++) {
              strncpy(ant_code[iant], (ant_prm+iant)->IDC,
                      strlen((ant_prm+iant)->IDC));
              ant_prm[iant].UFL = false;
              _off_button(&ant_prm[iant].UFL, &ant_prm[iant].IDC[0],
                          ant_bttn_box[iant]);
            }
          }
        } else {
          cpgsci(0);
          cpgsfs(1);
          cpgrect(0.009, 1.402, bttn_box[ARRAY_SECTION][3]+0.001, 0.344);
          cpgsci(1);
          cpgsfs(1);

          cpgsci(1);
          cpgsfs(2);
          cpgrect(0.010, 1.400, bttn_box[ARRAY_SECTION][3]+0.040, 0.345);
          cpgsci(1);
          cpgsfs(1);
        }
        char_copy(antenna_list_file, antenna_list_file_tmp);
        obs_param_file_io(ERROR_FLAG, antenna_list_file,
                    array, ANT_NUM, GRT_NUM, SRT_NUM,
                    srt, grt_elevation_limit,
                    sep_angle_limit_from_earth_limb,
                    TimUTC, UT1_UTC, obs_duration,
                    &SRCPROC_MODE, src, sun, ant_code,
                    ch_grt_el_lim, &pair_src, &ch_src,
                    ch_srt, &ch_obs_t, 1);
        return (RE_LOAD);
      }

/*
-----------------------------------------
*/

      for (i=ARRAY_SECTION; i<ARRAY_SECTION+ARRAY_NUM; i++) {
        if (_button_chk(cursor_pos, bttn_box[i]) == true) {
          I = i - ARRAY_SECTION;
          array->ID = array_id[I].id;
          _toggle_button(&array_id[I].flag, array_id[I].name, bttn_box[i]);
          if (array->ID == NO_ANT) {
            for (j=0; j<ARRAY_NUM; j++) {
              J = ARRAY_SECTION + j;
              _off_button(&array_id[j].flag, array_id[j].name, bttn_box[J]);
            }
            for (j=0; j<GRT_NUM_tmp; j++) {
              _off_button(&ant_prm[j].UFL, &ant_prm[j].IDC[0], ant_bttn_box[j]);
            }
          } else {
            for (j=0; j<GRT_NUM_tmp; j++) {
              idum = array_config(array->ID, wave_id, 0, &grt_num,
                                  &ant_prm_tmp, ant_prm[j].IDC,
                                  antenna_list_file,   true, true);
              if (idum == 1 &&
                  strncmp(ant_prm_tmp.IDC,
                          ant_prm[j].IDC, strlen(ant_prm[j].IDC)) == 0) {
                if (array_id[I].flag == true) {
                  _on_button (&ant_prm[j].UFL,
                              &ant_prm[j].IDC[0], ant_bttn_box[j]);
                } else if (array_id[I].flag == false) {
                  _off_button(&ant_prm[j].UFL,
                              &ant_prm[j].IDC[0], ant_bttn_box[j]);
                }
              }
            }
          }

          if (array->ID == ACA) {
            cpgsci(0);
          } else {
            cpgsci(1);
          }
          I = SRT_SECTION;
          sprintf(string, "SRT Number (0-%d)", SRTMAX);
          cpgsci(1);
          cpgtext(0.02, text_bottom(bttn_box[I][2], bttn_box[I][3]), string);
          cpgrect(bttn_box[I][0], bttn_box[I][1],
                  bttn_box[I][2], bttn_box[I][3]);
          if (array->ID == ACA) {
            *SRT_NUM = 0;
            SRT_SWT = false;
          } else {
            sprintf(string, "%d", *SRT_NUM);
            cpgsci(0);
            cpgptxt(bttn_box[I][1] - 0.015,
                    text_bottom(bttn_box[I][2], bttn_box[I][3]),
                    0.0, 1.0, string);
            cpgsci(1);
            SRT_SWT = true;
          }
        }
      }

      for (i=0; i<GRT_NUM_tmp; i++) {
        if (_button_chk(cursor_pos, ant_bttn_box[i]) == true) {
          _toggle_button(&ant_prm[i].UFL, &ant_prm[i].IDC[0], ant_bttn_box[i]);
        }
      }

/*
--------
*/

      if (_button_chk(cursor_pos, bttn_box[SRT_SECTION]) == true &&
          SRT_SWT == true) {
        str_init(string, nstr);
        I = SRT_SECTION;
        sprintf(string, "%d", *SRT_NUM);
        tv_get_param("int", cursor_pos, bttn_box[I], pitch, string,
                     0.0, (double)SRTMAX);
        sscanf(string, "%d", SRT_NUM);
        str_init(string, nstr);
        TV_menu_hatch(0.250, 1.400, bttn_box[0][3],
                  bttn_box[SRT_SECTION][2] + 0.030, 0, 1);
        TV_menu_hatch(0.010, 0.250, bttn_box[0][3],
                  bttn_box[SRT_SECTION][2] + 0.000, 0, 1);
        if (*SRT_NUM != 0) {
          srt_info_disp(*SRT_NUM,  pitch, bttn_box+SRT_SECTION+1,
                        bttn_box[SRT_SECTION][2],
                        srtaer, &ERROR_FLAG[SRTAER], ch_srt);
        }
      }

/*
-----------------------------------------
*/

      if (_button_chk(cursor_pos, bttn_box[I=TIME_SECTION  ]) == true) {
        sprintf(string, "%s", ch_obs_t.start_t[0]);
        tv_get_param("char", cursor_pos, bttn_box[I],
                     pitch, string, 0.0, 0.0);
        sprintf(ch_obs_t.start_t[0], "%s", string);
        char_ncopy(ch_obs_t.start_t[0], string, 4);
        str_init(string, nstr);
      }

      if (_button_chk(cursor_pos, bttn_box[I=TIME_SECTION+1]) == true) {
        sprintf(string, "%s", ch_obs_t.start_t[1]);
        tv_get_param("char", cursor_pos, bttn_box[I],
                     pitch, string, 0.0, 0.0);
        sprintf(ch_obs_t.start_t[1], "%s", string);
        char_ncopy(ch_obs_t.start_t[1], string, 2);
        str_init(string, nstr);
      }

      if (_button_chk(cursor_pos, bttn_box[I=TIME_SECTION+2]) == true) {
        sprintf(string, "%s", ch_obs_t.start_t[2]);
        tv_get_param("char", cursor_pos, bttn_box[I],
                     pitch, string, 0.0, 0.0);
        sprintf(ch_obs_t.start_t[2], "%s", string);
        char_ncopy(ch_obs_t.start_t[2], string, 2);
        str_init(string, nstr);
      }

      if (_button_chk(cursor_pos, bttn_box[I=TIME_SECTION+3]) == true) {
        sprintf(string, "%s", ch_obs_t.start_t[3]);
        tv_get_param("char", cursor_pos, bttn_box[I],
                     pitch, string, 0.0, 0.0);
        sprintf(ch_obs_t.start_t[3], "%s",  string);
        char_ncopy(ch_obs_t.start_t[3], string, 2);
        str_init(string, nstr);
      }

      if (_button_chk(cursor_pos, bttn_box[I=TIME_SECTION+4]) == true) {
        sprintf(string, "%s", ch_obs_t.start_t[4]);
        tv_get_param("char", cursor_pos, bttn_box[I],
                     pitch, string, 0.0, 0.0);
        sprintf(ch_obs_t.start_t[4], "%s",  string);
        char_ncopy(ch_obs_t.start_t[4], string, 2);
        str_init(string, nstr);
      }

      if (_button_chk(cursor_pos, bttn_box[I=TIME_SECTION+5]) == true) {
        sprintf(string, "%s", ch_obs_t.start_t[5]);
        tv_get_param("char", cursor_pos, bttn_box[I],
                     pitch, string, 0.0, 0.0);
        sprintf(ch_obs_t.start_t[5], "%s",  string);
        char_ncopy(ch_obs_t.start_t[5], string, 2);
        str_init(string, nstr);
      }

      if (_button_chk(cursor_pos, bttn_box[I=TIME_SECTION+6]) == true) {
        tv_get_param("double", cursor_pos, bttn_box[I],
                     pitch, ch_obs_t.obsd, 0.0, 0.0);
        sscanf(ch_obs_t.obsd, "%lf", obs_duration);
        sprintf(string, "Oberving duration is set to %lf [h].",
                *obs_duration);
        comment_disp(cmnt, comment, string, true);
      }

/*
-----------------------------------------
*/

      if (*SRT_NUM != 0) {
        I = SRT_SECTION + 1;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          if (ERROR_FLAG[SRTAER] == true) {
            idum = 1;
          } else {
            idum = 0;
          }
          toggle_button(&idum, srtaer, bttn_box[I]);
          if (idum == 1) {
            ERROR_FLAG[SRTAER] = true;
          } else {
            ERROR_FLAG[SRTAER] = false;
          }
        }

        for (i=0; i<*SRT_NUM; i++) {
          for (j=0; j<8; j++) {
            I++;
            if (_button_chk(cursor_pos, bttn_box[I]) == true) {
              if (j == 0) {
                tv_get_param("char", cursor_pos, bttn_box[I],
                             pitch, ch_srt[i].apo, 0.0, 0.0);
              } else if (j == 1) {
                tv_get_param("char", cursor_pos, bttn_box[I],
                             pitch, ch_srt[i].per, 0.0, 0.0);
              } else if (j == 2) {
                tv_get_param("char", cursor_pos, bttn_box[I],
                             pitch, ch_srt[i].inc, 0.0, 0.0);
              } else if (j == 3) {
                tv_get_param("char", cursor_pos, bttn_box[I],
                             pitch, ch_srt[i].OMG, 0.0, 0.0);
              } else if (j == 4) {
                tv_get_param("char", cursor_pos, bttn_box[I],
                             pitch, ch_srt[i].omg, 0.0, 0.0);
              } else if (j == 5) {
                tv_get_param("char", cursor_pos, bttn_box[I],
                             pitch, ch_srt[i].t_0, 0.0, 0.0);
              } else if (j == 6) {
                tv_get_param("char", cursor_pos, bttn_box[I],
                             pitch, ch_srt[i].d_OMG, 0.0, 0.0);
              } else if (j == 7) {
                tv_get_param("char", cursor_pos, bttn_box[I],
                             pitch, ch_srt[i].d_omg, 0.0, 0.0);
              }
            }
          }
        }
      }

/*
----
*/

    }
    cpgebuf();
  }

/*
============================================================
*/

  if (TV_SWT == false) {
    printf("ARIS starts....\n");
  } else if (TV_SWT == true) {
    if (START_FLAG == true) {
      comment_disp(cmnt, comment, "ARIS starts...\0", true);
    } else if (START_FLAG == false) {
      comment_disp(cmnt, comment, "EXIT.\0", true);
    }
  }

/*
============================================================
*/

  if (obs_param_set(ERROR_FLAG, *SRT_NUM,
                    ch_srt, srt,
                    ch_grt_el_lim, grt_elevation_limit,
                    sep_angle_limit_from_earth_limb,
                    &ch_obs_t, TimUTC, UT1_UTC,
                    obs_duration, SRCPROC_MODE,
                    &ch_src, &pair_src, src, sun) == -1) {
    printf("ERROR: OBS_PARAM_INPUT: SUN ANGLE CHECK.\n");
    return (__NG__);
  }

/*
============================================================
*/

  *GRT_NUM = 0;
  for (iant=0; iant<GRT_NUM_tmp; iant++) {
    if (ant_prm[iant].UFL == true) {
      sprintf(ant_code[*GRT_NUM], "%s", (ant_prm+iant)->IDC);
      (*GRT_NUM)++;
    }
  }
  obs_param_file_io(ERROR_FLAG, antenna_list_file,
                    array, ANT_NUM, GRT_NUM, SRT_NUM,
                    srt, grt_elevation_limit,
                    sep_angle_limit_from_earth_limb,
                    TimUTC, UT1_UTC, obs_duration,
                    &SRCPROC_MODE, src, sun, ant_code,
                    ch_grt_el_lim, &pair_src, &ch_src,
                    ch_srt, &ch_obs_t, 1);

/*
-----------------------------------------------------
*/

  if (START_FLAG == true) {
    j = 0;
    for (i=0; i<GRT_NUM_tmp; i++) {
      if (ant_prm[i].UFL == true) {
        j++;
      }
    }
    *GRT_NUM = j;
    if (j + *SRT_NUM == 0) {
      sprintf(string, "ERROR: OBS_PARAM_INPUT: NO ANTENNA selected.");
      printf("%s\n", string);
      if (TV_SWT == true) {
        comment_disp(cmnt, comment, string, true);
      }
      return (__NG__);
    } else if (j + *SRT_NUM == 1) {
      sprintf(string, "WARNING: OBS_PARAM_INPUT: NO BASELINE selected.");
      printf("%s\n", string);
      if (TV_SWT == true) {
        comment_disp(cmnt, comment, string, true);
      }
      return (__GO__);
    } else {
      return (__GO__);
    }
  } else if (START_FLAG == false) {
    return (_EXIT_);
  }
}


int  ant_list_chk(char *antenna_list_file, float *bttn_box,
                  struct comment_param      *cmnt,
                  char comment[][NCOMLEN], _Bool  TV_SWT)
{
  FILE   *fp;
  char   string[500];

  if ((fp=fopen(antenna_list_file, "r")) == NULL) {
    sprintf(string, "ERROR: ANT_LIST_CHK: ");
    sprintf(string, "antenna list file does not exit or cannot be read: %s",
         antenna_list_file);

    if (TV_SWT == true) {
      cpgsci(2);
      cpgptxt(bttn_box[1], text_bottom(bttn_box[2], bttn_box[3]),
              0.0, -0.05, string);
      cpgsci(1);
      comment_disp(cmnt, comment, string, true);
    } else {
      printf("%s\n", string);
    }

/********
    sprintf(string, "Input antenna list file.");
    if (TV_SWT == true) {
      cpgsci(2);
      cpgptxt(bttn_box[1], text_bottom(bttn_box[2], bttn_box[3]),
              0.0, -0.05, string);
      cpgsci(1);
      comment_disp(cmnt, comment, string, true);
    } else {
      printf("%s\n", string);

      printf("Antenna list file     (CR->%s) : ", antenna_list_file);
      if (fgets(string, nstr, stdin) == NULL) {
        printf("ERROR: OBS_PARAM_INPUT :");
        printf("There is something wrong in standard input.\n");
        return (__NG__);
      }
    }
********/

    return -1;
  } else {
/**
    if (TV_SWT == true) {
      cpgsci(0);
      cpgrect(bttn_box[1]+0.02, 1.200, bttn_box[2], bttn_box[3]);
    }
**/
    return  1;
  }
}
