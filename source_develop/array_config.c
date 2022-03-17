#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
/****
#include <astrotools.h>
****/
#include <aris.h>

/****
#define __DEBUG__
****/

int  array_config(
                  int    ARRAY_ID,
                  int    *wave_id,
                  int    SRT_NUM,
                  int    *GRT_NUM,
                  struct antenna_parameter *ant_prm,
                  char   *antenna_code,
                  char   *antenna_list_file,
                  _Bool  INDIVIDUAL_PICK_UP_SWT,
                  _Bool  ERR_RESET_SWT)
{
  int    i, j, k, nline, array_id;
  int    iant, ANT_NUM, ns;
  int    dd, mm;
  double ss, DPI;
  double XYZ[3];
  char   string[1000];
  FILE   *fp;
  _Bool  REC_FLAG, ARRAY_REC_FLAG, STATION_REC_FLAG;
  struct antenna_parameter ant_tmp, srt_tmp;

/*
----------------
*/

  DPI = dpi / 180.0;
  *GRT_NUM = 0;

/*
----------------
*/

  if ((fp = fopen(antenna_list_file, "r")) == NULL) {
    printf("ERROR: array_config: aris_input/antenna.prm not found.\n");
    return -1;
  }

  nline    = 0;
  ANT_NUM  = 0;
  while (1) {
    if (fgets(string, sizeof(string), fp) == NULL) {
      break;
    } else {
      nline++;
    }

    if (strncmp(string, "BEGIN_ANT", 9) == 0) {
#ifdef __DEBUG_1__
      printf("__DEBUG_1__ : BEGIN_ANT\n");
#endif
      REC_FLAG         = false;
      STATION_REC_FLAG = false;
      ARRAY_REC_FLAG   = false;
      if (ARRAY_ID == ALL_ANT) {
        REC_FLAG = true;
      }

      for (i=0; i<10; i++) {
        ant_tmp.IDC[i] = 0;
      }
      ant_tmp.ARRAY    = 0;
      ant_tmp.UFL      = false;
      ant_tmp.FRQSTD   = NONE;
      for (i=0; i<3; i++) {
        ant_tmp.LLH[i] = 0.0;
        ant_tmp.XYZ[i] = 0.0;
        ant_tmp.ERR[i] = 0.0;
      }
      for (ns=0; ns<SRC_NUM; ns++) {
        if (wave_id[ns] == ALLBAND) {
          ant_tmp.WID[ns]  =  0;
        } else {
          ant_tmp.WID[ns]  = -1;
        }
        ant_tmp.Dm[ns]     = 0.0;
        ant_tmp.Ae[ns]     = 0.0;
        ant_tmp.Trx[ns]    = 0.0;
        ant_tmp.Tsky[ns]   = 0.0;
      }
      while (1) {
        if (fgets(string, sizeof(string), fp) == NULL) {
          printf("WARNING: something wrong with %s.\n", antenna_list_file);
          printf("Check line %d in %s.\n", nline, antenna_list_file);
          break;
        } else {
          nline++;
          if (strncmp(string, "END_ANT", 7) == 0) {
#ifdef __DEBUG_1__
            printf("__DEBUG_1__ : BEGIN_ANT\n");
#endif
            break;
          } else {

/*
--------
*/

            if (strncmp(string, "STATION", 7) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == '\"') {
                  break;
                }
              }
              j = 0;
              while (j != 8) {
                if (string[i+j] == '\"') {
                  break;
                }
                j++;
              }
              if (j == 8 && string[i+j] != '\"') {
                k = 1;
                while (1) {
                  if (string[i+j+k] == '\"') {
                    break;
                  }
                  k++;
                }
                string[i+j+k] = '\0';
                printf("WARNING: array_config: ");
                printf("Length of STATION (%s) must be shorter\n", string+i);
                printf("         than nine characters.\n");
              }

              strncpy(ant_tmp.IDC, string+i, j);
              ant_tmp.IDC[j] = '\0';
              if (strlen(antenna_code) == j &&
                  strncmp(antenna_code, ant_tmp.IDC, j) == 0) {
                STATION_REC_FLAG = true;
                if (ARRAY_ID == -1) {
                  REC_FLAG = true;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "ARRAY", 5) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == '\"') {
                  break;
                }
              }
              j = 0;
              while (1) {
                if (string[i+j] == '\"') {
                  break;
                }
                j++;
              }
              string[i+j] = 0;
              if (strncmp(string+i, "VLBA",                    j) == 0) {
                array_id = VLBA;
              } else if (strncmp(string+i, "EVN",              j) == 0) {
                array_id = EVN;
              } else if (strncmp(string+i, "HSA",              j) == 0) {
                array_id = HSA;
              } else if (strncmp(string+i, "VERA",             j) == 0) {
                array_id = VERA;
              } else if (strncmp(string+i, "JVN",              j) == 0) {
                array_id = JVN;
              } else if (strncmp(string+i, "KVN",              j) == 0) {
                array_id = KVN;
              } else if (strncmp(string+i, "LBA",              j) == 0) {
                array_id = LBA;
              } else if (strncmp(string+i, "TRACKING_NETWORK", j) == 0) {
                array_id = TRACKING_NETWORK;
              } else if (strncmp(string+i, "ALMA",             j) == 0) {
                array_id = ALMA;
              } else if (strncmp(string+i, "ACA",              j) == 0) {
                array_id = ACA;
              } else if (strncmp(string+i, "EALMA",            j) == 0) {
                array_id = EALMA;
              } else if (strncmp(string+i, "STAND_ALONE",      j) == 0) {
                array_id = STAND_ALONE;
              } else if (strncmp(string+i, "ORBITING",         j) == 0) {
                array_id = ORBITING;
              } else if (strncmp(string+i, "SLR_NETWORK",      j) == 0) {
                array_id = SLR_NETWORK;
              } else if (strncmp(string+i, "KAVA",             j) == 0) {
                array_id = KAVA;
              } else {
                array_id = -100;
              }

              ant_tmp.ARRAY = array_id;
              if (array_id == ARRAY_ID) {
                REC_FLAG       = true;
                ARRAY_REC_FLAG = true;
              }
            }

/*
--------
*/

            if (strncmp(string, "MOUNT", 5) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == '\"') {
                  break;
                }
              }
              j = 0;
              while (1) {
                if (string[i+j] == '\"') {
                  break;
                }
                j++;
              }
              if (strncmp(string+i, "ALAZ", j) == 0) {
                ant_tmp.MNTSTA = ALAZ;
              } else if (strncmp(string+i, "EQUA", j) == 0) {
                ant_tmp.MNTSTA = EQUA;
              } else if (strncmp(string+i, "ORBI", j) == 0) {
                ant_tmp.MNTSTA = ORBI;
              }
            }

/*
--------
*/

            if (strncmp(string, "FREQ_STANDARD", 13) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == '\"') {
                  break;
                }
              }
              j = 0;
              while (1) {
                if (string[i+j] == '\"') {
                  break;
                }
                j++;
              }
              if (strncmp(string+i, "H_M", j) == 0) {
                ant_tmp.FRQSTD = H_M;
              } else if (strncmp(string+i, "CSO_10",  j) == 0) {
                ant_tmp.FRQSTD = CSO_10;
              } else if (strncmp(string+i, "CSO_100", j) == 0) {
                ant_tmp.FRQSTD = CSO_100;
              } else if (strncmp(string+i, "NONE",    j) == 0) {
                ant_tmp.FRQSTD = NONE;
              }
            }

/*
--------
*/

            if (strncmp(string, "POSITION", 8) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf,%lf,%lf",
                     &ant_tmp.XYZ[0], &ant_tmp.XYZ[1], &ant_tmp.XYZ[2]);
            }

/*
--------
*/

            if (strncmp(string, "POS_ERROR", 9) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf,%lf,%lf",
                     &ant_tmp.AAE[0], &ant_tmp.AAE[1], &ant_tmp.AAE[2]);
            }

/*
--------
*/

            if (strncmp(string, "POS_OFFSET", 10) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf,%lf,%lf",
                     &ant_tmp.OFS[0], &ant_tmp.OFS[1], &ant_tmp.OFS[2]);
            }

/*
--------
*/

            if (strncmp(string, "STAXOF", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf", &ant_tmp.STAXOF);
            }

/*
--------
*/

            if (strncmp(string, "LONGITUDE", 9) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%d %d %lf", &dd, &mm, &ss);
              if (dd >= 0) {
                ant_tmp.LLH[0] = (double)dd
                                        + (double)mm/60.0 + ss/3600.0;
              } else {
                ant_tmp.LLH[0] = fabs((double)dd)
                                        + (double)mm/60.0 + ss/3600.0;
                ant_tmp.LLH[0] *= -1.0;
              }
              ant_tmp.LLH[0] *= DPI;
            }

/*
--------
*/

            if (strncmp(string, "LATITUDE", 8) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%d %d %lf", &dd, &mm, &ss);
              if (dd >= 0) {
                ant_tmp.LLH[1] = (double)dd + (double)mm/60.0 + ss/3600.0;
              } else {
                ant_tmp.LLH[1] = fabs((double)dd)
                                            + (double)mm/60.0 + ss/3600.0;
                ant_tmp.LLH[1] *= -1.0;
              }
              ant_tmp.LLH[1] *= DPI;
            }

/*
--------
*/

            if (strncmp(string, "ALTITUDE", 8) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf", &ant_tmp.LLH[2]);
            }

/*
--------
*/

            if (strncmp(string, "AZ_RATE_ACCEL", 13) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf,%lf", &ant_tmp.AZSV, &ant_tmp.AZSA);
              ant_tmp.AZSV *= DPI;
              ant_tmp.AZSA *= DPI;
            }

/*
--------
*/

            if (strncmp(string, "EL_RATE_ACCEL", 13) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf,%lf", &ant_tmp.ELSV, &ant_tmp.ELSA);
              ant_tmp.ELSV *= DPI;
              ant_tmp.ELSA *= DPI;
            }

/*
--------
*/

            if (strncmp(string, "EL_LIMIT", 8) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf", &ant_tmp.ELLIM);
              ant_tmp.ELLIM *= DPI;
            }

/*
--------
*/

            if (strncmp(string, "FREQ_SWT_LAG", 12) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf", &ant_tmp.FQSST);
            }

/*
--------
*/

            if (strncmp(string, "SLEW_TIME", 9) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf", &ant_tmp.slewt);
            }

/*
--------
*/

            if (strncmp(string, "LO_PHASE_CHANGE", 15) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              sscanf(string+i, "%lf,%lf", &ant_tmp.lo_phs_jmp_val, &ant_tmp.lo_phs_jmp_tim);
            }

/*
--------
*/

            if (strncmp(string, "L_BAND", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == L_BAND) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = L_BAND;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "S_BAND", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == S_BAND) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = S_BAND;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "C_BAND", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == C_BAND) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = C_BAND;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "X_BAND", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == X_BAND) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = X_BAND;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "KU_BAND", 7) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == KU_BAND) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = KU_BAND;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "K_BAND", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == K_BAND) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = K_BAND;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "Q_BAND", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == Q_BAND) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = Q_BAND;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "W_BAND", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == W_BAND) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = W_BAND;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND01", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND01) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND01;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND02", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND02) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND02;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND03", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND03) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND03;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND04", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND04) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND04;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND05", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND05) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND05;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND06", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND06) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND06;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND07", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND07) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND07;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND08", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND08) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND08;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND09", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND09) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND09;
                }
              }
            }

/*
--------
*/

            if (strncmp(string, "BAND10", 6) == 0) {
              i = 0;
              while (1) {
                if (string[i++] == ':') {
                  break;
                }
              }
              for (ns=0; ns<SRC_NUM; ns++) {
                if (wave_id[ns] == BAND10) {
                  sscanf(string+i, "%lf,%lf,%lf,%lf",
                         &ant_tmp.Dm[ns],  &ant_tmp.Ae[ns],
                         &ant_tmp.Trx[ns], &ant_tmp.Tsky[ns]);
                  ant_tmp.WID[ns] = BAND10;
                }
              }
            }

/*
--------
*/

          }
        }
      }

/*
----------------------------------------
*/

      if ((INDIVIDUAL_PICK_UP_SWT == false && REC_FLAG == true) ||
          (INDIVIDUAL_PICK_UP_SWT == true  &&
           STATION_REC_FLAG == true && ARRAY_REC_FLAG == true)) {

        ant_tmp.UFL = true;
        if (ant_tmp.XYZ[0] == 0.0 &&
            ant_tmp.XYZ[1] == 0.0 &&
            ant_tmp.XYZ[2] == 0.0 &&
              (ant_tmp.LLH[0] != 0.0 ||
               ant_tmp.LLH[1] != 0.0 ||
               ant_tmp.LLH[2] != 0.0)) {
          J_system_geocentric_equatorial_rectangular_coordinate(
            ant_tmp.LLH[0], ant_tmp.LLH[1], ant_tmp.LLH[2], ant_tmp.XYZ);
        }

        if (ant_tmp.LLH[0] == 0.0 &&
            ant_tmp.LLH[1] == 0.0 &&
            ant_tmp.LLH[2] == 0.0 &&
              (ant_tmp.XYZ[0] != 0.0 ||
               ant_tmp.XYZ[1] != 0.0 ||
               ant_tmp.XYZ[2] != 0.0)) {
          J_system_geocentric_equatorial_rectangular_coordinate2llh(
            &ant_tmp.LLH[0], &ant_tmp.LLH[1], &ant_tmp.LLH[2], ant_tmp.XYZ);
        }

        if (ant_tmp.OFS[0] != 0.0 ||
            ant_tmp.OFS[1] != 0.0 ||
            ant_tmp.OFS[2] != 0.0) {
          XYZ[0] = ant_tmp.OFS[0];
          XYZ[1] = ant_tmp.OFS[1];
          XYZ[2] = ant_tmp.OFS[2];
          drotate(XYZ, 0.5*dpi-ant_tmp.LLH[1], "x");
          drotate(XYZ, ant_tmp.LLH[0]+0.5*dpi, "z");
          ant_tmp.XYZ[0] += XYZ[0];
          ant_tmp.XYZ[1] += XYZ[1];
          ant_tmp.XYZ[2] += XYZ[2];
          J_system_geocentric_equatorial_rectangular_coordinate2llh(
            &ant_tmp.LLH[0], &ant_tmp.LLH[1], &ant_tmp.LLH[2], ant_tmp.XYZ);
        }

        if (ERR_RESET_SWT == true) {
          ant_tmp.ERR[0] =  ant_tmp.XYZ[0];
          ant_tmp.ERR[1] =  ant_tmp.XYZ[1];
          ant_tmp.ERR[2] =  ant_tmp.XYZ[2];
        }

        if (ant_tmp.ARRAY == ORBITING) {
          srt_tmp          = ant_tmp;
        } else if (ant_tmp.ARRAY != SLR_NETWORK) {
          ant_prm[ANT_NUM] = ant_tmp;
          ANT_NUM++;
          if (ANT_NUM == ANTMAX) {
            printf("WARNING: ");
            printf("the number of antennas reasches ");
            printf("its maximum number: %d\n", ANT_NUM);
            fclose (fp);
            return (ANT_NUM);
          }
        }
      }
    }
  }
  fclose (fp);
  *GRT_NUM = ANT_NUM;

/*
--------------------
*/

  if (SRT_NUM >= 1) {
    for (iant=0; iant<SRT_NUM; iant++) {
      ant_prm[ANT_NUM] = srt_tmp;
      if (SRT_NUM == 1) {
        sprintf(ant_prm[ANT_NUM].IDC, "SRT");
      } else {
        if (iant+1 < 10) {
          sprintf(ant_prm[ANT_NUM].IDC, "SRT_0%1d", iant+1);
        } else if (iant+1 >= 10 && iant+1 < 100) {
          sprintf(ant_prm[ANT_NUM].IDC, "SRT_%2d",  iant+1);
        }
      }
      ANT_NUM++;
      if (ANT_NUM == ANTMAX) {
        printf("WARNING: ");
        printf("the number of antennas reasches ");
        printf("its maximum number: %d\n", ANT_NUM);
        fclose (fp);
        return (ANT_NUM);
      }
    }
  }

/*
--------------------
*/

#ifdef __DEBUG__
  for (iant=0; iant<ANT_NUM; iant++) {
    printf("__DEBUG__  %d  %s  %d  %d  %d\n",
      iant, ant_prm[iant].IDC, ant_prm[iant].ARRAY, ant_prm[iant].UFL,
      ant_prm[iant].FRQSTD);
    printf("__DEBUG__  (XYZ) %lf  %lf  %lf\n", 
      ant_prm[iant].XYZ[0], ant_prm[iant].XYZ[1], ant_prm[iant].XYZ[2]);
    printf("__DEBUG__  (LLH) %lf  %lf  %lf\n", 
      ant_prm[iant].LLH[0], ant_prm[iant].LLH[1], ant_prm[iant].LLH[2]);
    printf("__DEBUG__  (ERR) %lf  %lf  %lf\n", 
      ant_prm[iant].ERR[0], ant_prm[iant].ERR[1], ant_prm[iant].ERR[2]);
    printf("__DEBUG__  (OFS) %lf  %lf  %lf\n", 
      ant_prm[iant].OFS[0], ant_prm[iant].OFS[1], ant_prm[iant].OFS[2]);
    printf("__DEBUG__        %lf  %lf  %lf  %lf  %lf\n",
            ant_prm[iant].AZSV, ant_prm[iant].AZSA,
            ant_prm[iant].ELSV, ant_prm[iant].ELSA, ant_prm[iant].FQSST);
    for (ns=0; ns<SRC_NUM; ns++) {
      printf("__DEBUG__  %lf  %lf  %lf  %lf\n", 
        ant_prm[iant].Dm[ns],  ant_prm[iant].Ae[ns],
        ant_prm[iant].Trx[ns], ant_prm[iant].Tsky[ns]);
    }
    printf("__DEBUG__\n");
  }
#endif

  return ANT_NUM;
}
