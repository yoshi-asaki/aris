#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>

#define BOTTUN_NUM     200

#define TROPOS_SECTION  10
#define IONOS_SECTION   20
#define ORBIT_SECTION   30
#define SRTATT_SECTION  40
#define OBSTIM_SECTION  50
#define T_WAVE_SECTION  60
#define R_WAVE_SECTION  80
#define BW_SECTION     100
#define FRCHAN_SECTION 110
#define CFACT_SECTION  120
#define T_FLUX_SECTION 130
#define R_FLUX_SECTION 150
#define SWTCYC_SECTION 170
#define TRK_SECTION    180

#define ALLBL            0
#define GG_BL            1
#define GS_BL            2
#define SS_BL            3

double  band_width_set(int  );
double  cfact_set(int  );
int     channel_num_set(int  );

/****
#define __DEBUG__
****/

#define TRP_NUM   8

int   err_parameter_set(int ANT_NUM,  int GRT_NUM,  int SRT_NUM,
                int *BGN_ANT_I, int *END_ANT_I,
                int *BGN_ANT_J, int *END_ANT_J,
                struct array_parameter *array,
                int *TRP_CONDITION,  int *ION_CONDITION,
                double *Cw,     double *CI,
                struct phase_screen_parameter   *wvc,
                struct atmospheric_zenith_error *dz,
                struct source_parameter  *src,
                struct morphology_file_name *ch_file,
                int *wave_id, double *wave_length, double *nu,
                double *band_width, int *FRCHAN_NUM,
                double *cfact,
                int    *reference_phase_process_mode,
                _Bool   *ERROR_FLAG,
                struct srt_orbit_parameter *srt,
                struct data_number *data_num,
                int  *TimUTC, double *UT1_UTC,
                int *nswt,    double *apparent_tgt_on,
                struct antenna_parameter *ant_prm,
                int    *TRK_NUM,
                struct antenna_parameter  *trk_pos,
                struct comment_param *cmnt,
                char   comment[][NCOMLEN],
                _Bool  TV_SWT, float *cursor_pos, int cpgid1, int  cpgid2)
{
  int    i, j, k, iant, ns;
  int    ncol;
  int    idum, jdum, I, J;
  int    BL_SELECT;
  char   *c=0, string[100];
  _Bool  Bdum, QUIT_SWT;
  float  bttn_box[BOTTUN_NUM][4], y_pos, srtatt_y_pos, x_pos;
  float  ctrl_bttn_box[3][4], blsel_bttn_box[4][4];
  float  ant_bttn_box[ANTMAX][4];
  float  ftmp1, ftmp2;
  double swt_cyc_time;
  double orbit_error_ap, orbit_error_pe;
  double tdscz, idscz;
  double wl, fr;
  FILE   *fp;
  int    IANT, ITRP;

  int    blsel_num=0;
  struct blsel_status {
    _Bool  code;
    char   name[20];
  } blsel_status[4];

  int    trp_num, trp_status_num;
  int    trp_tune;
  struct trp_status {
    _Bool  code;
    char   name[20];
  } trp_status[8];

  int    ion_num;
  struct ion_status {
    _Bool  code;
    char   name[10];
  } ion_status[3];

  int    itrk, trk_num;
  int    TRK_GEN_CMD=2;
  char   trk_name[TRKMAX][10];
  int    trk_priority[TRKMAX], priority[TRKMAX];

  int    wave_num, wave_code[N_WAVE];
  _Bool  wave_swt[SRC_NUM][N_WAVE];
  char   wave_name__on[N_WAVE][10], wave_name_off[N_WAVE][10];

  int    srtatt_num, srtatt_code[SRTMAX];
  char   srtatt_name[10][10];
  float  pitch;

  char   ch_tdscz[20], ch_idscz[20];
  char   ch_tflux[20], ch_rflux[20], ch_swtim[20], ch_appon[20];
  char   ch_orbae[20], ch_orbpe[20];
  char   ch_tmorp[20], ch_rmorp[20];
  char   ch_tpost[20], ch_rpost[20];
  char   ch_timer[2][20];

  int    bw_id, bw_num, cfact_id, cfact_num;
  _Bool  bw_code[10], cfact_code[10];
  char   bw_name[10][10], cfact_name[10][10];

  int    frchan_id, frchan_num;
  _Bool  frchan_code[10];
  char   frchan_name[10][10];

  char   g_name[5][10];

  int    tgt_morpho_num, ref_morpho_num;
  char   tgt_morphology[SRC_MORPH_CANDIDATE][50];
  char   ref_morphology[SRC_MORPH_CANDIDATE][50];
  _Bool  tgt_proc_flag[SRC_MORPH_CANDIDATE];
  _Bool  ref_proc_flag[SRC_MORPH_CANDIDATE];
  int    FLX_SECT=2, MRP_SECT=3;

  int    nday, timUTC[6];
  int    id, ih, im, is;
  double ut1_utc, mjd;

  struct ant_trp_tune {
    char   IDC[10];
    int    code;
  } *ant_trp, *ant_trp_io;
  char   ant_tmp[10];

  float  trp_bttn_box[ANTMAX][7][4];

/*
---------------------------------------------------------
*/

  pitch   = 0.03;

/*
---------------------------------------------------------
*/

  sprintf((blsel_status  )->name, "Full");
  sprintf((blsel_status+1)->name, "Ground-Ground");
  sprintf((blsel_status+2)->name, "Ground-Space");
  sprintf((blsel_status+3)->name, "Space-Space");
  if (SRT_NUM == 0) {
    blsel_num = 0;
  } else if (SRT_NUM == 1) {
    blsel_num = 3;
  } else if (SRT_NUM >= 2) {
    blsel_num = 4;
  }
  for (i=0; i<blsel_num; i++) {
    blsel_status[i].code = false;
  }

/*
---------------------------------------------------------
*/

  trp_num  = 7;
  trp_tune = -1;

  trp_status_num = 7;
  if (       array->TYPE == __CONNECTED_) {
    trp_num  = 7;
    trp_tune = -1;
  } else if (array->TYPE == _VLBI_ARRAY_) {
    trp_num  = 8;
    trp_tune = 7;
  }

  sprintf((trp_status  )->name, "ALMA: 0.5 rad");
  sprintf((trp_status+1)->name, "ALMA: 0.75 rad");
  sprintf((trp_status+2)->name, "ALMA: 1 rad");
  sprintf((trp_status+3)->name, "Very Good");
  sprintf((trp_status+4)->name, "Good");
  sprintf((trp_status+5)->name, "Typical");
  sprintf((trp_status+6)->name, "Poor");
  sprintf((trp_status+7)->name, "Tuning");
  for (i=0; i<trp_num; i++) {
    trp_status[i].code = false;
  }

/*
---------------------------------------------------------
*/

  ion_num = 3;
  sprintf((ion_status  )->name, "Half");
  sprintf((ion_status+1)->name, "Nominal");
  sprintf((ion_status+2)->name, "Double");
  for (i=0; i<ion_num; i++) {
    ion_status[i].code = false;
  }

/*
---------------------------------------------------------
*/

  srtatt_num = 3;
  sprintf(srtatt_name[0], "+X SUN");
  sprintf(srtatt_name[1], "-X SUN");
  sprintf(srtatt_name[2], "NO MATTER");
  for (iant=0; iant<SRT_NUM; iant++) {
    srtatt_code[iant] = 2;
  }

/*
---------------------------------------------------------
*/

  wave_num = 16;
  wave_code[ 0] = L_BAND;
  wave_code[ 1] = S_BAND;
  wave_code[ 2] = C_BAND;
  wave_code[ 3] = X_BAND;
  wave_code[ 4] = KU_BAND;
  wave_code[ 5] = K_BAND;
  wave_code[ 6] = Q_BAND;
  wave_code[ 7] = W_BAND;
  wave_code[ 8] = BAND03;
  wave_code[ 9] = BAND04;
  wave_code[10] = BAND05;
  wave_code[11] = BAND06;
  wave_code[12] = BAND07;
  wave_code[13] = BAND08;
  wave_code[14] = BAND09;
  wave_code[15] = BAND10;
  sprintf(wave_name_off[ 0], "L Band");
  sprintf(wave_name_off[ 1], "S Band");
  sprintf(wave_name_off[ 2], "C Band");
  sprintf(wave_name_off[ 3], "X Band");
  sprintf(wave_name_off[ 4], "Ku Band");
  sprintf(wave_name_off[ 5], "K Band");
  sprintf(wave_name_off[ 6], "Q Band");
  sprintf(wave_name_off[ 7], "W Band");
  sprintf(wave_name_off[ 8], "BAND 3");
  sprintf(wave_name_off[ 9], "BAND 4");
  sprintf(wave_name_off[10], "BAND 5");
  sprintf(wave_name_off[11], "BAND 6");
  sprintf(wave_name_off[12], "BAND 7");
  sprintf(wave_name_off[13], "BAND 8");
  sprintf(wave_name_off[14], "BAND 9");
  sprintf(wave_name_off[15], "BAND10");
  for (i=0; i<wave_num; i++) {
    wave_select(wave_code[i], &wl, &fr);
    sprintf(wave_name__on[i], "%5.1fGHz", (float)fr*1.0e-9);
  }
  for (ns=0; ns<SRC_NUM; ns++) {
    for (i=0; i<wave_num; i++) {
      wave_swt[ns][i] = false;
    }
  }

/*
---------------------------------------------------------
*/

  bw_num = 10;
  sprintf(bw_name[0], "1.95kHz");
  sprintf(bw_name[1], "4MHz");
  sprintf(bw_name[2], "8MHz");
  sprintf(bw_name[3], "16MHz");
  sprintf(bw_name[4], "32MHz");
  sprintf(bw_name[5], "128MHz");
  sprintf(bw_name[6], "256MHz");
  sprintf(bw_name[7], "512MHz");
  sprintf(bw_name[8], "1024MHz");
  sprintf(bw_name[9], "2048MHz");
  for (i=0; i<bw_num; i++) {
    bw_code[i] = false;
  }

/*
---------------------------------------------------------
*/

  frchan_num = 5;
  sprintf(frchan_name[0], "1CH");
  sprintf(frchan_name[1], "4CH");
  sprintf(frchan_name[2], "8CH");
  sprintf(frchan_name[3], "16CH");
  sprintf(frchan_name[4], "32CH");
  for (i=0; i<frchan_num; i++) {
    frchan_code[i] = false;
  }

/*
---------------------------------------------------------
*/

  cfact_num = 3;
  sprintf(cfact_name[0], "1 bit");
  sprintf(cfact_name[1], "2 bit");
  sprintf(cfact_name[2], "3 bit");
  for (i=0; i<cfact_num; i++) {
    cfact_code[i] = false;
  }

/*
---------------------------------------------------------
*/

  tgt_morpho_num = SRC_MORPH_CANDIDATE;
  for (i=0; i<tgt_morpho_num; i++) {
    if (i == 0) {
      sprintf(tgt_morphology[0], "Single Point");
    } else if (i == 1) {
      sprintf(tgt_morphology[1], "Multi-Component");
    } else if (i == 2) {
      sprintf(tgt_morphology[2], "Disk & Jet-CJet");
    } else if (i == 3) {
      sprintf(tgt_morphology[3], "Disk & VSOP2-Jet");
    } else if (i == 4) {
      sprintf(tgt_morphology[4], "AIPS CC Table");
    } else if (i == 5) {
      sprintf(tgt_morphology[5], "Blackhole Shadow");
    }
  }

  ref_morpho_num = SRC_MORPH_CANDIDATE;
  for (i=0; i<ref_morpho_num; i++) {
    if (i == 0) {
      sprintf(ref_morphology[0], "Single Point");
    } else if (i == 1) {
      sprintf(ref_morphology[1], "Multi-Component");
    } else if (i == 2) {
      sprintf(ref_morphology[2], "Disk & Jet-CJet");
    } else if (i == 3) {
      sprintf(ref_morphology[3], "Disk & VSOP2-Jet");
    } else if (i == 4) {
      sprintf(ref_morphology[4], "AIPS CC Table");
    } else if (i == 5) {
      sprintf(ref_morphology[5], "Blackhole Shadow");
    }
  }

/*
---------------------------------------------------------
*/

  sprintf(g_name[0], "Lo(high)");
  sprintf(g_name[1], "GRT err");

/*
---------------------------------------------------------
*/

  *CI                = 0.0;
  orbit_error_ap     = 0.0;
  orbit_error_pe     = 0.0;
  trk_num            = 0;
  wave_id[0]         = 5;
  wave_id[1]         = 5;
  bw_id              = 0;
  frchan_id          = 0;
  cfact_id           = 0;
  src[0].flux        = 1.0;
  src[1].flux        = 1.0;
  src[0].morphology  = SRC_POINT;
  src[1].morphology  = SRC_POINT;
  src[0].positionID  = 0;
  src[1].positionID  = 1;
  swt_cyc_time       = 60.0;
  *apparent_tgt_on   = 50.0;
  BL_SELECT          = ALLBL;
  *TRP_CONDITION     = 5;
  *ION_CONDITION     = NOMINAL;
  *reference_phase_process_mode = 0;

  for (iant=0; iant<GRT_NUM; iant++) {
    Cw[iant]         = 0.0;
  }

  for (itrk=0; itrk<TRKMAX; itrk++) {
    trk_priority[itrk] = 0;
    priority[itrk]     = 0;
  }

  for (iant=0; i<SRT_NUM; iant++) {
    srtatt_code[iant] = 2;
  }

  if ((fp = fopen("aris_input/sim.prm", "r")) != NULL) {
    while (1) {
      if (fgets(string, sizeof(string), fp) == NULL) {
        break;
      } else {
        string[strlen(string)-1] = '\0';
      }

      if (       strncmp(string, "BASELINE SELECT     ", 20) == 0) {
        sscanf(string+20, "%d", &BL_SELECT);
      } else if (strncmp(string, "TROPOS CONDITION    ", 20) == 0) {
        sscanf(string+20, "%d", TRP_CONDITION);
      } else if (strncmp(string, "IONOS CONDITION     ", 20) == 0) {
        sscanf(string+20, "%d", ION_CONDITION);
      } else if (strncmp(string, "TROPOS ZENITH ERROR ", 20) == 0) {
        sscanf(string+20, "%lf", &tdscz);
        if (strlen(string+20) == 0) {
          sprintf(ch_tdscz, "0");
        } else {
          char_copy(ch_tdscz, string+20);
        }
      } else if (strncmp(string, "IONOS ZENITH ERROR  ", 20) == 0) {
        sscanf(string+20, "%lf", &idscz);
        if (strlen(string+20) == 0) {
          sprintf(ch_idscz, "0");
        } else {
          char_copy(ch_idscz, string+20);
        }
      } else if (strncmp(string, "ORBIT ERROR APOGEE  ", 20) == 0) {
        sscanf(string+20, "%lf", &orbit_error_ap);
        if (strlen(string+20) == 0) {
          sprintf(ch_orbae, "0");
        } else {
          char_copy(ch_orbae, string+20);
        }
        orbit_error_ap *= 1.0e-2;
      } else if (strncmp(string, "ORBIT ERROR PERIGEE ", 20) == 0) {
        sscanf(string+20, "%lf", &orbit_error_pe);
        if (strlen(string+20) == 0) {
          sprintf(ch_orbpe, "0");
        } else {
          char_copy(ch_orbpe, string+20);
        }
        orbit_error_pe *= 1.0e-2;
      } else if (strncmp(string, "X_AXIS_SUN_ANGLE    ", 20) == 0) {
        sscanf(string+20, "%d,%d", &idum, &jdum);
        idum--;
        if (jdum > 0) {
          srtatt_code[idum] = 0;
        } else if (jdum < 0) {
          srtatt_code[idum] = 1;
        } else if (jdum == 0) {
          srtatt_code[idum] = 2;
        }
      } else if (strncmp(string, "TRACKING NETWORK    ", 20) == 0) {
        trk_num = tracking_station_name_read(trk_num,  string+20,
                                             trk_name, priority);
      } else if (strncmp(string, "START TIME (UTC)    ", 20) == 0) {
        char_copy(ch_timer[0], string+20);
      } else if (strncmp(string, "STOP TIME  (UTC)    ", 20) == 0) {
        char_copy(ch_timer[1], string+20);
      } else if (strncmp(string, "TGT WAVE ID         ", 20) == 0) {
        sscanf(string+20, "%d", &wave_id[0]);
      } else if (strncmp(string, "REF WAVE ID         ", 20) == 0) {
        sscanf(string+20, "%d", &wave_id[1]);
      } else if (strncmp(string, "BAND WIDTH          ", 20) == 0) {
        sscanf(string+20, "%d", &bw_id);
      } else if (strncmp(string, "FREQUENCY CHANNEL   ", 20) == 0) {
        sscanf(string+20, "%d", &frchan_id);
      } else if (strncmp(string, "COHERENCE FACTOR    ", 20) == 0) {
        sscanf(string+20, "%d", &cfact_id);
      } else if (strncmp(string, "SWITCHING CYC TIME  ", 20) == 0) {
        sscanf(string+20, "%lf", &swt_cyc_time);
        if (strlen(string+20) == 0) {
          sprintf(ch_swtim, "0");
        } else {
          char_copy(ch_swtim, string+20);
        }
      } else if (strncmp(string, "APPARENT TARGET ON  ", 20) == 0) {
        sscanf(string+20, "%lf", apparent_tgt_on);
        if (strlen(string+20) == 0) {
          sprintf(ch_appon, "0");
        } else {
          char_copy(ch_appon, string+20);
        }
      } else if (strncmp(string, "TGT FLUX DENSITY    ", 20) == 0) {
        sscanf(string+20, "%lf", &src[0].flux);
        if (strlen(string+20) == 0) {
          sprintf(ch_tflux, "0");
        } else {
          char_copy(ch_tflux, string+20);
        }
      } else if (strncmp(string, "TGT POSITION        ", 20) == 0) {
        sscanf(string+20, "%d", &src[0].positionID);
      } else if (strncmp(string, "TGT MORPHOLOGY      ", 20) == 0) {
        sscanf(string+20, "%d", &src[0].morphology);
      } else if (strncmp(string, "TGT MULTI COMPONENT ", 20) == 0) {
        sscanf(string+20, "%s", ch_file[0].mcm);
      } else if (strncmp(string, "TGT CC TABLE        ", 20) == 0) {
        sscanf(string+20, "%s", ch_file[0].cct);
      } else if (strncmp(string, "TGT BHS MODEL       ", 20) == 0) {
        sscanf(string+20, "%s", ch_file[0].bhs);
      } else if (strncmp(string, "REF FLUX DENSITY    ", 20) == 0) {
        sscanf(string+20, "%lf", &src[1].flux);
        if (strlen(string+20) == 0) {
          sprintf(ch_rflux, "0");
        } else {
          char_copy(ch_rflux, string+20);
        }
      } else if (strncmp(string, "REF POSITION        ", 20) == 0) {
        sscanf(string+20, "%d", &src[1].positionID);
      } else if (strncmp(string, "REF MORPHOLOGY      ", 20) == 0) {
        sscanf(string+20, "%d", &src[1].morphology);
      } else if (strncmp(string, "REF MULTI COMPONENT ", 20) == 0) {
        sscanf(string+20, "%s", ch_file[1].mcm);
      } else if (strncmp(string, "REF CC TABLE        ", 20) == 0) {
        sscanf(string+20, "%s", ch_file[1].cct);
      } else if (strncmp(string, "REF BHS MODEL       ", 20) == 0) {
        sscanf(string+20, "%s", ch_file[1].bhs);
      } else if (strncmp(string, "REF PHASE PROCESS   ", 20) == 0) {
        sscanf(string+20, "%d", reference_phase_process_mode);
      }
    }
    fclose (fp);
  }

  blsel_status[BL_SELECT].code    = true;
  ion_status[*ION_CONDITION].code = true;
  trp_status[*TRP_CONDITION].code = true;


/****
  sscanf(ch_timer[0], "%d/%d:%d:%d", &sday, &shh, &smm, &sss);
  sscanf(ch_timer[1], "%d/%d:%d:%d", &eday, &ehh, &emm, &ess);
  smjd = MJD(TimUTC[0], TimUTC[1], TimUTC[2]+sday,
             shh,       smm,       sss, *UT1_UTC);
  emjd = MJD(TimUTC[0], TimUTC[1], TimUTC[2]+eday,
             ehh,       emm,       ess, *UT1_UTC);
  mjd  = MJD(TimUTC[0], TimUTC[1], TimUTC[2],
             TimUTC[3], TimUTC[4], TimUTC[5], *UT1_UTC);

  if (smjd < mjd) {
    sday = 0;
    shh  = TimUTC[3];
    smm  = TimUTC[4];
    sss  = TimUTC[5];
  } else if (emjd > mjd + (double)data_num.nobs/86400.0) {
    eday = (int)MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                    TimUTC[3], TimUTC[4], TimUTC[5] + data_num->nobs,
                    *UT1_UTC);
         - (int)MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                    TimUTC[3], TimUTC[4], TimUTC[5], *UT1_UTC);
    MJD2date(mjd, &timUTC[0], &timUTC[1], &timUTC[2],
                  &timUTC[3], &timUTC[4], &timUTC[5]);
    ehh  = timUTC[3];
    emm  = timUTC[4];
    ess  = timUTC[5];
  }
****/

  wave_swt[0][wave_select(wave_id[0], &wave_length[0], &nu[0])]  = true;
  wave_swt[1][wave_select(wave_id[1], &wave_length[1], &nu[1])]  = true;
  bw_code[bw_id]           = true;
  frchan_code[frchan_id]   = true;
  *FRCHAN_NUM = channel_num_set(frchan_id);
  cfact_code[cfact_id]     = true;
  if (SRT_NUM >= 1) {
    *TRK_NUM = tracking_init(trk_num, priority, trk_priority, trk_name, trk_pos);
  }

/*
---------------------------------------------------------
*/

  if (TV_SWT == false) {

/*
----------------------------------------------------
*/

    if (       array->TYPE == __CONNECTED_) {
      BL_SELECT == ALLBL;
    } else if (array->TYPE == _VLBI_ARRAY_) {
      while (1) {
        if (SRT_NUM >= 1) {
          printf("1. Full  2. Ground-Ground  3.Ground-Space  ");
          if (SRT_NUM >= 2) {
            printf("4. Space-Space    ");
          }
          printf("0. EXIT (CR->%d) : ", BL_SELECT+1);
          if (fgets(string, sizeof(string), stdin) == NULL) {
            printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
            return (-1);
          }
          if (string[0] == '\n') {
            break;
          } else {
            sscanf(string, "%d", &idum);
            if (idum == 0) {
              return (-1);
            } else if ((idum >= 1 && idum <= 3) ||
                       (idum == 4 && SRT_NUM >= 2)) {
              BL_SELECT = idum - 1;
              break;
            }
          }
        } else if (SRT_NUM == 0) {
          printf("1. Ground-Ground    0. EXIT (CR->%d) : ", BL_SELECT+1);
          if (fgets(string, sizeof(string), stdin) == NULL) {
            printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
            return (-1);
          }
          if (string[0] == '\n') {
            break;
          } else {
            sscanf(string, "%d", &idum);
            if (idum == 0) {
              return (-1);
            } else if (idum == 1) {
              BL_SELECT = ALLBL;
              break;
            }
          }
        }
      }
    }

/*
----------------------------------------------------
*/

    if (ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true) {
      while (1) {
        printf("Tropospheric condition:\n");
        for (i=0; i<trp_num; i++) {
          printf("%d. %s\n", i+1, (trp_status+i)->name);
        }
        printf("(CR->%d) : ", *TRP_CONDITION+1);
        if (fgets(string, sizeof(string), stdin) == NULL) {
          printf("ERROR: ERR_PARAMETER SET: Invalid input: %s", string);
          return (-1);
        }

        if (string[0] == '\n') {
          for (iant=0; iant<GRT_NUM; iant++) {
            ant_prm[iant].WVturb = *TRP_CONDITION;
          }
          break;
        } else {
          sscanf(string, "%d", &idum);
          if (idum >= 1 && idum < trp_num) {
            *TRP_CONDITION = idum - 1;
            for (iant=0; iant<GRT_NUM; iant++) {
              ant_prm[iant].WVturb = *TRP_CONDITION;
            }
            break;
          } else if (idum == trp_num) {
            if ((ant_trp = (struct ant_trp_tune *)
              calloc(GRT_NUM, sizeof(struct ant_trp_tune))) == NULL) {
              printf("ERROR: err_parameter: ");
              printf("tropospheric tuning parameter memry allocation.\n");
              exit (-1);
            } else {
              for (iant=0; iant<GRT_NUM; iant++) {
                strcpy((ant_trp+iant)->IDC, (ant_prm+iant)->IDC);
              }
            }

            if ((fp = fopen( "aris_input/ant_trp_tuning.prm", "r")) != NULL) {
              while (1) {
                if (fgets(string, sizeof(string), fp) == NULL) {
                  break;
                } else {
                  sscanf(string, "%s %d", &ant_tmp, &idum);
                  for (iant=0; iant<GRT_NUM; iant++) {
                    if (strncmp((ant_trp+iant)->IDC, ant_tmp,
                                strlen(ant_tmp)) == 0) {
                      ant_trp[iant].code = idum;
                      break;
                    }
                  }
                }
              }
              fclose (fp);
            } else {
              for (iant=0; iant<GRT_NUM; iant++) {
                ant_trp[iant].code = 0;
              }
            }

            printf("Antenna List:\n");
            for (iant=0; iant<GRT_NUM; iant++) {
              printf("%3d. %s", iant+1, (ant_trp+iant)->IDC);
              k = strlen(ant_trp[i].IDC);
              if (k > 10) {
                k = 10;
              }
              for (j=0; j<11-k; j++) {
                printf(" ");
              }
              printf("[%d]   ", ant_trp[iant].code+1);
              if (iant % 5 == 4) {
                printf("\n");
              }
            }
            if ((iant-1) % 5 != 4) {
              printf("\n");
            }

            while (1) {
              printf("Please select antenna (CR->quit): ");
              if (fgets(string, sizeof(string), stdin) == NULL) {
                printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
                return (-1);
              }
              if (string[0] == '\n') {
                break;
              }
              sscanf(string, "%d", &IANT);
              IANT--;
              printf("Please select the tropospheric option for %s (1->%d, CR->%d): ",
                     (ant_trp+IANT)->IDC, trp_status_num, ant_trp[IANT].code+1);
              if (fgets(string, sizeof(string), stdin) == NULL) {
                printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
                return (-1);
              }
              if (string[0] != '\n') {
                sscanf(string, "%d", &ITRP);
                ant_trp[IANT].code = ITRP-1;
              }
            }

            if ((fp = fopen(
                  "aris_input/ant_trp_tuning.prm", "w")) != NULL) {
              for (iant=0; iant<GRT_NUM; iant++) {
                fprintf(fp, "%s %d\n",
                     (ant_trp+iant)->IDC, ant_trp[iant].code);
              }
              fclose (fp);
            }

            for (iant=0; iant<GRT_NUM; iant++) {
              ant_prm[iant].WVturb = ant_trp[iant].code;
            }
            free (ant_trp);
            break;
          }
        }
      }
    }

/*
----------------------------------------------------
*/

    if (ERROR_FLAG[IONTRB] == true) {
      printf("Ionospheric Noise Power:\n");
      printf("1. Half Condition  ");
      printf("2. Nominal Codition  ");
      printf("3. Double Condition (CR->%d) : ", *ION_CONDITION+1);

      while (1) {
        if (fgets(string, sizeof(string), stdin) == NULL) {
          printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
          return (-1);
        }
        if (string[0] == '\n') {
          break;
        } else {
          sscanf(string, "%d", &idum);
          if (idum >= 1 && idum <= ion_num) {
            *ION_CONDITION = idum - 1;
            break;
          }
        }
      }
    }

/*
----------------------------------------------------
*/

    if (ERROR_FLAG[TDSECZ] == true && GRT_NUM > 0
     && array->TYPE == _VLBI_ARRAY_) {
      printf("Input Tropospheric Zenith Error [mm] (CR->%lf[mm]) : ", tdscz);

      if (fgets(ch_tdscz, sizeof(ch_tdscz), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (ch_tdscz[0] == '\n') {
        sprintf(ch_tdscz, "%lf", tdscz);
      } else {
        sscanf(ch_tdscz, "%lf", &tdscz);
        ch_tdscz[strlen(ch_tdscz)-1] = '\0';
      }
    }

/*
----------------------------------------------------
*/

    if (ERROR_FLAG[IDSECZ] == true && GRT_NUM > 0) {
      printf("Input Ionospheric Zenith Error [TECU] (CR->%lf[TECU]) : ", idscz);
      if (fgets(ch_idscz, sizeof(ch_idscz), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (ch_idscz[0] == '\n') {
        sprintf(ch_idscz, "%lf", idscz);
      } else {
        sscanf(ch_idscz, "%lf", &idscz);
        ch_idscz[strlen(ch_idscz)-1] = '\0';
      }
    }

/*
----------------------------------------------------
*/

    if (ERROR_FLAG[APOSER] == true && SRT_NUM > 0) {
      printf("Input the Orbit Error at Apogee  [cm] (CR->%lf[cm]) : ",
             orbit_error_ap * 1.0e2);
      if (fgets(ch_orbae, sizeof(ch_orbae), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (ch_orbae[0] == '\n') {
        sprintf(ch_orbae, "%lf", orbit_error_ap * 1.0e+2);
      } else {
        sscanf(ch_orbae, "%lf", &orbit_error_ap);
        orbit_error_ap *= 1.0e-2;
        ch_orbae[strlen(ch_orbae)-1] = '\0';
      }

      printf("Input the Orbit Error at Perigee [cm] (CR->%lf[cm]) : ",
             orbit_error_pe * 1.0e2);
      if (fgets(ch_orbpe, sizeof(ch_orbpe), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (ch_orbpe[0] == '\n') {
        sprintf(ch_orbpe, "%lf", orbit_error_pe * 1.0e+2);
      } else {
        sscanf(ch_orbpe, "%lf", &orbit_error_pe);
        orbit_error_pe *= 1.0e-2;
        ch_orbpe[strlen(ch_orbpe)-1] = '\0';
      }
    }

/*
----------------------------------------------------
*/

    if (SRT_NUM >= 1) {
      printf("Tracking Network\n");
      for (itrk=0; itrk<*TRK_NUM; itrk++) {
        printf("%d. %s  ", itrk+1, trk_name[itrk]);
        if ((itrk+1)%7 == 0) {
          printf("\n");
        }
      }
      if (*TRK_NUM%7 != 0) {
        printf("\n");
      }
      printf("Select tracking stations : ");
      printf("(CR->");
      for (itrk=0; itrk<*TRK_NUM; itrk++) {
        if (trk_priority[itrk] >= 1 && trk_priority[itrk] <= 9) {
          printf("%1d", trk_priority[itrk]);
        } else if (trk_priority[itrk] >= 10) {
          *c = trk_priority[itrk] - 10 + 'a';
          printf("%1s", c);
        } else if (trk_priority[itrk] == 0) {
          printf("0");
        }
      }
      printf(") : ");
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] != '\n') {
        for (itrk=0; itrk<*TRK_NUM; itrk++) {
          if (string[itrk] >= '1' && string[itrk] <= '9') {
            trk_priority[itrk] = string[itrk] - '0';
          } else if (string[itrk] >= 'a' && string[itrk] <= 'z') {
            trk_priority[itrk] = string[itrk] - 'a' + 10;
          } else {
            trk_priority[itrk] = 0;
          }
        }
      }

      trk_num = 0;
      for (itrk=0; itrk<*TRK_NUM; itrk++) {
        if (trk_priority[itrk] >= 1 && trk_priority[itrk] <= 9) {
          printf("%1d", trk_priority[itrk]);
          trk_num++;
        } else if (trk_priority[itrk] >= 10) {
          *c = trk_priority[itrk] - 10 + 'a';
          printf("%1s", c);
          trk_num++;
        } else if (trk_priority[itrk] == 0) {
          printf("0");
        }
      }
      printf("\n");
    }

/*
----------------------------------------------------
*/

    if (SRT_NUM >= 1 && trk_num >= 1) {
      for (iant=0; iant<SRT_NUM; iant++) {
        printf("SRT [%d]  Attitude (CR->%d) :\n", iant, srtatt_code[iant]+1);
        printf("( 1. PLUS  X-axis -> SUN )\n");
        printf("( 2. MINUS X-axis -> SUN )\n");
        printf("( 3. does not matter.    ) : ");
        while (1) {
          if (fgets(string, sizeof(string), stdin) == NULL) {
            printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
            return (-1);
          }
          if (string[0] == '\n') {
            break;
          }
          sscanf(string, "%d", srtatt_code+iant);
          if (srtatt_code[iant] >= 1 && srtatt_code[iant] <= 3) {
            srtatt_code[iant]--;
            break;
          } else {
            printf("Invalid input: again : ");
          }
        }
      }
    }

/*
----------------------------------------------------
*/

    printf("Wave Band (Target)\n");
    ncol = 5;

    I = 0;
    for (i=0; i<=wave_num/ncol; i++) {
      for (j=0; j<ncol; j++) {
        printf("%2d. %s     ", I+1, wave_name_off[I]);
        I++;
        if (I == wave_num) {
          break;
        }
      }
      printf("\n");
      if (I == wave_num) {
        break;
      }
    }
    printf("(CR->%d) : ",
           1 + wave_select(wave_id[0], &wave_length[0], &nu[0]));
    while (1) {
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] == '\n') {
        break;
      } else {
        sscanf(string, "%d", &idum);
        if (idum >= 1 && idum <= wave_num) {
          wave_id[0] = wave_code[idum-1];
          break;
        }
      }
    }

    printf("Wave Band (Reference)\n");
    I = 0;
    for (i=0; i<=wave_num/ncol; i++) {
      for (j=0; j<ncol; j++) {
        printf("%2d. %s     ", I+1, wave_name_off[I]);
        I++;
        if (I == wave_num) {
          break;
        }
      }
      printf("\n");
      if (I == wave_num) {
        break;
      }
    }
    printf("(CR->%d) : ",
           1 + wave_select(wave_id[1], &wave_length[1], &nu[1]));
    while (1) {
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] == '\n') {
        break;
      } else {
        sscanf(string, "%d", &idum);
        if (idum >= 1 && idum <= wave_num) {
          wave_id[1] = wave_code[idum-1];
          break;
        }
      }
    }

/*
----------------------------------------------------
*/

    printf("Observing Band Width\n");
    ncol = 4;
    k = 0;
    for (i=0; i<bw_num/ncol+1; i++) {
      for (j=0; j<ncol; j++) {
        sprintf(string, "         ");
        strcpy(string, bw_name[k]);
        string[strlen(bw_name[k])] = ' ';
        printf("%2d. %10s  ", k+1, string);
        k++;
        if (k == bw_num) {
          break;
        }
      }
      printf("\n");
    }
    printf("(CR->%d) : ", bw_id+1);
    while (1) {
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] == '\n') {
        break;
      } else {
        sscanf(string, "%d", &bw_id);
        if (bw_id >= 1 && bw_id <= bw_num) {
          bw_id--;
          break;
        } else {
          printf("Invalid. Enter again: ");
        }
      }
    }
    *band_width = band_width_set(bw_id);

/*
----------------------------------------------------
*/

    printf("Frequency Channel\n");
    printf("[1. 1CH   2. 4CH   3. 8CH   4. 16CH  5. 32CH ]\n");
    printf("(CR->%d) : ", frchan_id+1);
    while (1) {
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] == '\n') {
        break;
      } else {
        if (string[0]-'0' >= 1 && string[0]-'0' <= frchan_num) {
          frchan_id = string[0] - '0' - 1;
          break;
        } else {
          printf("Invalid. Enter again: ");
        }
      }
    }
    *FRCHAN_NUM = channel_num_set(frchan_id);

/*
----------------------------------------------------
*/

    printf("Sampling Bit Number [1.1-bit  2.2-bit  3.3-bit] (CR->%d) : ", cfact_id+1);
    while (1) {
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] == '\n') {
        break;
      } else {
        if (string[0]-'0' >= 1 && string[0]-'0' <= 2) {
          cfact_id = string[0] - '0' - 1;
          break;
        } else {
          printf("Invalid. Enter again: ");
        }
      }
    }
    *cfact = cfact_set(cfact_id);

/*
----------------------------------------------------
*/

    printf("Switching Cycle [sec] (CR->%lf) : ", swt_cyc_time);
    if (fgets(string, sizeof(string), stdin) == NULL) {
      printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
      return (-1);
    }
    if (string[0] != '\n') {
      sscanf(string, "%lf", &swt_cyc_time);
    }
    *nswt    = (int)lrint(swt_cyc_time);
    if (*nswt % 2 == 1) {
      (*nswt)++;
      printf("Switching cycle is changed to %d [sec].\n", *nswt);
    }
    swt_cyc_time = (double)(*nswt);
    sprintf(ch_swtim, "%d", *nswt);

/*
----------------------------------------------------
*/

    printf("Apparent Target ON Ratio (\%) (CR->%lf) : ", *apparent_tgt_on);
    if (fgets(string, sizeof(string), stdin) == NULL) {
      printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
      return (-1);
    }
    if (string[0] != '\n') {
      sscanf(string, "%lf", apparent_tgt_on);
    }
    sprintf(ch_appon, "%lf", *apparent_tgt_on);

/*
----------------------------------------------------
*/

    printf("Target Position    : 1. Source-1    2. Source-2\n");
    printf("Select number (CR->%d) : ", src[0].positionID + 1);
    if (fgets(ch_tpost, sizeof(ch_tpost), stdin) == NULL) {
      printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
      return (-1);
    }
    if (ch_tpost[0] == '\n') {
      sprintf(ch_tpost, "%d", src[0].positionID);
    } else {
      sscanf(ch_tpost, "%d", &src[0].positionID);
      ch_tpost[strlen(ch_tpost)-1] = '\0';
      src[0].positionID--;
    }

    printf("Target Morphology :\n");
    for (i=0; i<tgt_morpho_num; i++) {
      printf("%d. %s\n", i+1, tgt_morphology[i]);
    }
    printf("Select number (CR->%d) : ", src[0].morphology + 1);
    if (fgets(ch_tmorp, sizeof(ch_tmorp), stdin) == NULL) {
      printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
      return (-1);
    }
    if (ch_tmorp[0] == '\n') {
      sprintf(ch_tmorp, "%d", src[0].morphology);
    } else {
      sscanf(ch_tmorp, "%d", &src[0].morphology);
      ch_tmorp[strlen(ch_tmorp)-1] = '\0';
      src[0].morphology--;
    }

    if (src[0].morphology == SRC_MULTI_COMP) {
      printf("Multi-Component ASCII Data File : \n");
      if (strlen(ch_file[0].mcm) != 0) {
        printf("CR->%s : ", ch_file[0].mcm);
      }
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] != '\n') {
        sscanf(string, "%s", ch_file[0].mcm);
      }
    } else if (src[0].morphology == SRC_CC_COMP) {
      printf("CC Table File Name : \n");
      if (strlen(ch_file[0].cct) != 0) {
        printf("CR->%s : ", ch_file[0].cct);
      }
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] != '\n') {
        sscanf(string, "%s", ch_file[0].cct);
      }
    } else if (src[0].morphology == SRC_BHS_MOD) {
      printf("BHS Model Name : \n");
      if (strlen(ch_file[0].bhs) != 0) {
        printf("CR->%s : ", ch_file[0].bhs);
      }
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] != '\n') {
        sscanf(string, "%s", ch_file[0].bhs);
      }
    }

    if (src[0].morphology == SRC_POINT         ||
        src[0].morphology == SRC_DISK_JETCJET  ||
        src[0].morphology == SRC_DISK_VSOP2JET) {
      printf("Target Total Flux [Jy]   (CR->%lf) : ", src[0].flux);
      if (fgets(ch_tflux, sizeof(ch_tflux), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (ch_tflux[0] == '\n') {
        sprintf(ch_tflux, "%lf", src[0].flux);
      } else {
        sscanf(ch_tflux, "%lf", &src[0].flux);
        ch_tflux[strlen(ch_tflux)-1] = '\0';
      }
    }

/*
--------
*/

    printf("Reference Position : 1. Source-1    2. Source-2\n");
    printf("Select number (CR->%d) : ", src[1].positionID + 1);
    if (fgets(ch_rpost, sizeof(ch_rpost), stdin) == NULL) {
      printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
      return (-1);
    }
    if (ch_rpost[0] == '\n') {
      sprintf(ch_rpost, "%d", src[1].positionID);
    } else {
      sscanf(ch_rpost, "%d", &src[1].positionID);
      ch_rpost[strlen(ch_rpost)-1] = '\0';
      src[1].positionID--;
    }

    printf("Reference Morphology :\n");
    for (i=0; i<ref_morpho_num; i++) {
      printf("%d. %s\n", i+1, ref_morphology[i]);
    }
    printf("Select number (CR->%d) : ", src[1].morphology + 1);
    if (fgets(ch_rmorp, sizeof(ch_rmorp), stdin) == NULL) {
      printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
      return (-1);
    }
    if (ch_rmorp[0] == '\n') {
      sprintf(ch_rmorp, "%d", src[1].morphology);
    } else {
      sscanf(ch_rmorp, "%d", &src[1].morphology);
      ch_rmorp[strlen(ch_rmorp)-1] = '\0';
      src[1].morphology--;
    }

    if (src[1].morphology == SRC_MULTI_COMP) {
      printf("Multi-Component ASCII Data File : \n");
      if (strlen(ch_file[1].mcm) != 0) {
        printf("CR->%s : ", ch_file[1].mcm);
      }
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] != '\n') {
        sscanf(string, "%s", ch_file[1].mcm);
      }
    } else if (src[1].morphology == SRC_CC_COMP) {
      printf("CC Table File Name : \n");
      if (strlen(ch_file[1].cct) != 0) {
        printf("CR->%s : ", ch_file[1].cct);
      }
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] != '\n') {
        sscanf(string, "%s", ch_file[1].cct);
      }
    } else if (src[1].morphology == SRC_BHS_MOD) {
      printf("BHS Model Name : \n");
      if (strlen(ch_file[1].bhs) != 0) {
        printf("CR->%s : ", ch_file[1].bhs);
      }
      if (fgets(string, sizeof(string), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (string[0] != '\n') {
        sscanf(string, "%s", ch_file[1].bhs);
      }
    }

    if (src[1].morphology == SRC_POINT         ||
        src[1].morphology == SRC_DISK_JETCJET  ||
        src[1].morphology == SRC_DISK_VSOP2JET) {
      printf("Reference Total Flux [Jy]   (CR->%lf) : ", src[1].flux);
      if (fgets(ch_rflux, sizeof(ch_rflux), stdin) == NULL) {
        printf("ERROR: ERR_PARAMETER SET: Invalid input.\n");
        return (-1);
      }
      if (ch_rflux[0] == '\n') {
        sprintf(ch_rflux, "%lf", src[1].flux);
      } else {
        sscanf(ch_rflux, "%lf", &src[1].flux);
        ch_rflux[strlen(ch_rflux)-1] = '\0';
      }
    }

/*
--------
*/

    printf("Refrence phase: 1. Baseline-base  2. Antenna-base (CR->%d) : ",
            *reference_phase_process_mode+1);
    if (fgets(string, sizeof(string), stdin) == NULL) {
      printf("ERROR: ERR_PARAMETER_SET: Reference phase process mode.\n");
      return (-1);
    }
    if (string[0] != '\n') {
      sscanf(string, "%d", reference_phase_process_mode);
      *reference_phase_process_mode -= 1;
    } else {
      if (string[0] == '1') {
        *reference_phase_process_mode = 0;
      } else if (string[0] == '2') {
        *reference_phase_process_mode = 1;
      }
    }
/*
----------------------------------------------------
*/

  } else if (TV_SWT == true) {

/*
----------------------------------------------------
*/

    cpgbbuf();
    cpgslct(cpgid1);

    cpgpap(2.00*pgpap_prm, 1.0/1.4);
    cpgsch(1.5*pgpap_prm/13.0);
    cpgsvp(0.0,  1.0, 0.0, 1.0);
    cpgswin(0.0, 1.4, 0.0, 1.0);

/*
----------------------
*/

    comment_init(cmnt, comment, true);
    *nswt    = (int)lrint(swt_cyc_time);
    if (*nswt % 2 == 1) {
      (*nswt)++;
      printf("Switching cycle is changed to %d [sec].\n", *nswt);
    }
    swt_cyc_time = (double)(*nswt);

/*
----------------------
*/

    y_pos = 0.165;

    I = 0;
    ctrl_bttn_box[I][0] = 0.020;
    ctrl_bttn_box[I][1] = 0.250;
    ctrl_bttn_box[I][2] = y_pos;
    ctrl_bttn_box[I][3] = ctrl_bttn_box[I][2] + pitch;
    _off_button(&Bdum, "GO (baseline_base)\0", ctrl_bttn_box[I]);
    I = 1;
    ctrl_bttn_box[I][0] = 0.270;
    ctrl_bttn_box[I][1] = 0.500;
    ctrl_bttn_box[I][2] = y_pos;
    ctrl_bttn_box[I][3] = ctrl_bttn_box[I][2] + pitch;
    _off_button(&Bdum, "GO (antenna_base)\0",  ctrl_bttn_box[I]);
    I = 2;
    ctrl_bttn_box[I][0] = 0.650;
    ctrl_bttn_box[I][1] = 0.980;
    ctrl_bttn_box[I][2] = y_pos;
    ctrl_bttn_box[I][3] = ctrl_bttn_box[I][2] + pitch;
    _off_button(&Bdum, "EXIT (without saving parameters)\0",
               ctrl_bttn_box[I]);

    cursor_pos[0] = 0.5 * (ctrl_bttn_box[0][0] + ctrl_bttn_box[0][1]);
    cursor_pos[1] = y_pos + 0.5 * pitch;

/*
----------------------
*/

    y_pos = 1.065;

    if (SRT_NUM >= 1) {
      cpgsfs(2);
      cpgsci(1);
      cpgrect(0.020, 0.980, y_pos-0.002, y_pos+0.034);
      cpgsfs(1);

      cpgtext(0.035, y_pos + 0.3 * pitch, "Baseline Select\0");

      for (i=0; i<3; i++) {
        blsel_bttn_box[i][0] = 0.200 + 0.180 * (float)i;
        blsel_bttn_box[i][1] = blsel_bttn_box[i][0] + 0.16;
        blsel_bttn_box[i][2] = y_pos;
        blsel_bttn_box[i][3] = blsel_bttn_box[i][2] + pitch;
      }
      if (SRT_NUM >= 2) {
        i = 3;
        blsel_bttn_box[i][0] = 0.200 + 0.180 * (float)i;
        blsel_bttn_box[i][1] = blsel_bttn_box[i][0] + 0.16;
        blsel_bttn_box[i][2] = y_pos;
        blsel_bttn_box[i][3] = blsel_bttn_box[i][2] + pitch;
      }

      for (i=0; i<blsel_num; i++) {
        if (blsel_status[i].code == true) {
          _on_button (&blsel_status[i].code, (blsel_status+i)->name,
                      blsel_bttn_box[i]);
        } else if (blsel_status[i].code == false) {
          _off_button(&blsel_status[i].code, (blsel_status+i)->name,
                      blsel_bttn_box[i]);
        }
      }
    } else {
      BL_SELECT = ALLBL;
    }

/*
-------------
*/

    y_pos = 0.940;

    if (ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true) {
      cpgsfs(2);
      cpgsci(1);
      cpgrect(0.010, 1.400, y_pos-0.037, y_pos+0.035);
      cpgsfs(1);

      cpgtext(0.035, y_pos + 0.35 * pitch, "Tropospheric Condition\0");

      I = TROPOS_SECTION;
      for (i=0; i<trp_num; i++) {
        if (i < 3) {
          x_pos = 0.225;
        } else if (i == 7) {
          x_pos = 0.290;
        } else {
          x_pos = 0.255;
        }
        bttn_box[I][0] = x_pos + 0.128 * (float)i;
        bttn_box[I][1] = bttn_box[I][0] + 0.125;
        bttn_box[I][2] = y_pos;
        bttn_box[I][3] = bttn_box[I][2] + pitch;
        I++;
      }

      for (i=0; i<trp_num; i++) {
        I = TROPOS_SECTION + i;
        if (trp_status[i].code == true) {
          _on_button(&trp_status[i].code, (trp_status+i)->name, bttn_box[I]);
        } else {
          _off_button(&trp_status[i].code, (trp_status+i)->name, bttn_box[I]);
        }
      }
      if (i == trp_num - 1 && array->TYPE == _VLBI_ARRAY_) {
        _off_button(&trp_status[i].code, (trp_status+i)->name, bttn_box[I]);
      }
    }

/*
-------------
*/

    y_pos = 0.905;

    if (ERROR_FLAG[IONTRB] == true) {
      cpgsci(1);
      cpgtext(0.035, y_pos + 0.3 * pitch, "Ionospheric Disturbance Condition\0");

      I = IONOS_SECTION;
      for (i=0; i<ion_num; i++) {
        bttn_box[I][0] = 0.225 + 0.140 * (float)(i + 1);
        bttn_box[I][1] = bttn_box[I][0] + 0.130;
        bttn_box[I][2] = y_pos;
        bttn_box[I][3] = bttn_box[I][2] + pitch;
        I++;
      }

      for (i=0; i<ion_num; i++) {
        I = IONOS_SECTION + i;
        if (ion_status[i].code == true) {
          _on_button(&ion_status[i].code,  (ion_status+i)->name, bttn_box[I]);
        } else {
          _off_button(&ion_status[i].code, (ion_status+i)->name, bttn_box[I]);
        }
      }
    }

/*
--------------------
*/

    y_pos = 0.825;

    if (ERROR_FLAG[TDSECZ] == true && GRT_NUM > 0
     && array->TYPE ==  _VLBI_ARRAY_) {
      cpgsfs(2);
      cpgsci(1);
      cpgrect(0.010, 0.485, y_pos-0.039, y_pos+0.035);
      cpgsfs(1);
    }

    if (ERROR_FLAG[TDSECZ] == true && GRT_NUM > 0
     && array->TYPE ==  _VLBI_ARRAY_) {
      I = TROPOS_SECTION + trp_num;
      cpgsci(1);
      cpgtext(0.035, y_pos + 0.3 * pitch,
              "Tropospheric Zenith Delay Error [mm]\0");
      bttn_box[I][0] = 0.36;
      bttn_box[I][1] = bttn_box[I][0] + 0.11;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
              0.0, 1.0, ch_tdscz);
      cpgsci(1);
    }

/*
--------------------
*/

    y_pos -= 0.035;
    if (ERROR_FLAG[IDSECZ] == true && GRT_NUM > 0) {
      I = IONOS_SECTION + ion_num;
      cpgsci(1);
      cpgtext(0.035, y_pos + 0.3 * pitch,
              "Ionospheric Zenith Delay Error [TECU]\0");
      bttn_box[I][0] = 0.36;
      bttn_box[I][1] = bttn_box[I][0] + 0.11;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
              0.0, 1.0, ch_idscz);
      cpgsci(1);
    }

/*
-------------
*/

    if (ERROR_FLAG[APOSER] == true && SRT_NUM > 0) {
      I = ORBIT_SECTION;

      y_pos = 0.825;

      cpgsfs(2);
      cpgsci(1);
      cpgrect(0.515, 0.980, y_pos-0.039, y_pos+0.035);
      cpgsfs(1);

      bttn_box[I][0] = 0.860;
      bttn_box[I][1] = bttn_box[I][0] + 0.11;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      cpgsci(1);
      cpgtext(0.535, y_pos + 0.3 * pitch,
              "SRT displacement at Apogee [cm]\0");
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
              0.0, 1.0, ch_orbae);
      cpgsci(1);

      I = ORBIT_SECTION + 1;

      y_pos -= 0.035;
      bttn_box[I][0] = 0.860;
      bttn_box[I][1] = bttn_box[I][0] + 0.11;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      cpgsci(1);
      cpgtext(0.535, y_pos + 0.3 * pitch,
              "SRT displacement at Perigee [cm]\0");
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
              0.0, 1.0, ch_orbpe);
      cpgsci(1);
    }

/*
--------------------
*/

    if (SRT_NUM >= 1) {
      y_pos = 0.740;
      tracking_button_disp(TRK_SECTION, y_pos, pitch,
                           bttn_box, *TRK_NUM, trk_priority,
                           trk_name, trk_pos);
      srtatt_y_pos = y_pos - 0.005;
      srtatt_button_disp  (SRT_NUM, srtatt_num, srtatt_code, trk_num,
                           srtatt_y_pos, pitch, srtatt_name,
                           bttn_box+SRTATT_SECTION);
    }

/*
--------------------
*/

    y_pos = 0.630;

    cpgsfs(2);
    cpgsci(1);
    cpgrect(0.010, 1.400, y_pos-0.005, y_pos+0.034);
    cpgsfs(1);

    mjd = MJD(TimUTC[0], TimUTC[1], TimUTC[2],
              TimUTC[3], TimUTC[4], TimUTC[5] + data_num->nobs, *UT1_UTC);
    nday = (int)mjd
         - (int)MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                    TimUTC[3], TimUTC[4], TimUTC[5], *UT1_UTC);
    MJD2date(mjd, &timUTC[0], &timUTC[1], &timUTC[2],
                  &timUTC[3], &timUTC[4], &timUTC[5]);
    sprintf(string, "Start Time (UTC) 0/%2d:%2d:%2d",
                     TimUTC[3], TimUTC[4], TimUTC[5]);
    cpgtext(0.035, y_pos + 0.35 * pitch, string);
    sprintf(string, "Stop Time (UTC) %d/%2d:%2d:%2d",
            nday, timUTC[3], timUTC[4], timUTC[5]);
    cpgtext(0.535, y_pos + 0.35 * pitch, string);

    I = OBSTIM_SECTION;
    for (i=0; i<2; i++) {
      bttn_box[I][0] = 0.285 + 0.500 * (float)i;
      bttn_box[I][1] = bttn_box[I][0] + 0.150;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
              0.0, 1.0, ch_timer[i]);
      I++;
    }
    cpgsci(1);
    TV_menu_hatch(0.020, 0.980, y_pos-0.005, y_pos+0.034, 7, 4);

/*
--------------------
*/

    ncol = 16;
    y_pos = 0.558;

    cpgsfs(2);
    cpgsci(1);
    cpgrect(0.010, 1.400, y_pos-0.070, y_pos+0.060);
    cpgsfs(1);

    I = T_WAVE_SECTION;
    for (i=0; i<=wave_num/ncol; i++) {
      for (j=0; j<ncol; j++) {
        bttn_box[I][0] = 0.029 + (float)j * 0.085;
        bttn_box[I][1] = bttn_box[I][0] + 0.083;
        bttn_box[I][2] = y_pos - 0.035 * (float)i;
        bttn_box[I][3] = bttn_box[I][2] + pitch;
        I++;
        if (I == wave_num) {
          break;
        }
      }
    }
    cpgsci(1);
    cpgtext(0.035, bttn_box[T_WAVE_SECTION][3]+0.01,
            "Wave Band (Target)\0");
    cpgsci(8);
    cpgtext(0.235, bttn_box[T_WAVE_SECTION][3]+0.01,
        "[Wavelength setting: ./aris_input/fixed_parameter.prm]\0");
    cpgsci(1);

    for (i=0; i<wave_num; i++) {
      I = T_WAVE_SECTION + i;
      if (wave_swt[0][i] == true) {
        _on_button(&wave_swt[0][i], wave_name__on[i], bttn_box[I]);
      } else {
        _off_button(&wave_swt[0][i], wave_name_off[i], bttn_box[I]);
      }
    }

/*
-----------------
*/

    y_pos -= 0.061;
    I = R_WAVE_SECTION;
    for (i=0; i<=wave_num/ncol; i++) {
      for (j=0; j<ncol; j++) {
        bttn_box[I][0] = 0.029 + (float)j * 0.085;
        bttn_box[I][1] = bttn_box[I][0] + 0.083;
        bttn_box[I][2] = y_pos - 0.035 * (float)i;
        bttn_box[I][3] = bttn_box[I][2] + pitch;
        I++;
        if (I == wave_num) {
          break;
        }
      }
    }
    cpgsci(1);
    cpgtext(0.035, bttn_box[R_WAVE_SECTION][3]+0.01,
            "Wave Band (Reference)\0");
    cpgsci(8);
    cpgtext(0.235, bttn_box[R_WAVE_SECTION][3]+0.01,
        "[Wavelength setting: ./aris_input/fixed_parameter.prm]\0");
    cpgsci(1);

    for (i=0; i<wave_num; i++) {
      I = R_WAVE_SECTION + i;
      if (wave_swt[1][i] == true) {
        _on_button(&wave_swt[1][i], wave_name__on[i], bttn_box[I]);
      } else {
        _off_button(&wave_swt[1][i], wave_name_off[i], bttn_box[I]);
      }
    }

/*
-----------------
*/

    y_pos -= 0.074;

    cpgsfs(2);
    cpgsci(1);
    cpgrect(0.010, 1.400, y_pos-0.055, y_pos+0.060);
    cpgsfs(1);

    I = BW_SECTION;
    for (i=0; i<bw_num; i++) {
      bttn_box[I][0] = 0.029 + (float)i * 0.083;
      bttn_box[I][1] = bttn_box[I][0] + 0.080;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      I++;
    }
    cpgsci(1);
    cpgtext(0.035, bttn_box[BW_SECTION][3]+0.01, "Observing Band Width\0");

    for (i=0; i<bw_num; i++) {
      I = BW_SECTION + i;
      if (bw_code[i] == true) {
        _on_button(&bw_code[i], bw_name[i], bttn_box[I]);
      } else {
        _off_button(&bw_code[i], bw_name[i], bttn_box[I]);
      }
    }

/*
-----------------
*/

    I = FRCHAN_SECTION;
    for (i=0; i<frchan_num; i++) {
      bttn_box[I][0] = 0.950 + (float)i * 0.080;
      bttn_box[I][1] = bttn_box[I][0] + 0.075;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      I++;
    }
    cpgsci(1);
    cpgtext(0.950, bttn_box[FRCHAN_SECTION][3]+0.01, "Frequency Channel\0");

    for (i=0; i<frchan_num; i++) {
      I = FRCHAN_SECTION + i;
      if (frchan_code[i] == true) {
        _on_button(&frchan_code[i], frchan_name[i], bttn_box[I]);
      } else {
        _off_button(&frchan_code[i], frchan_name[i], bttn_box[I]);
      }
    }

/*
-----------------
*/

    y_pos -= 0.042;

    I = CFACT_SECTION;
    cpgsci(1);
    cpgtext(0.035, y_pos + 0.3 * pitch, "Sampling Bit Number\0");
    for (i=0; i<cfact_num; i++) {
      bttn_box[I][0] = 0.210 + (float)i * 0.060;
      bttn_box[I][1] = bttn_box[I][0]   + 0.055;
      bttn_box[I][2] = y_pos;
      bttn_box[I][3] = bttn_box[I][2] + pitch;
      I++;
    }

    for (i=0; i<cfact_num; i++) {
      I = CFACT_SECTION + i;
      if (cfact_code[i] == true) {
        _on_button(&cfact_code[i], cfact_name[i], bttn_box[I]);
      } else {
        _off_button(&cfact_code[i], cfact_name[i], bttn_box[I]);
      }
    }

/*
-----------------
*/

    I = SWTCYC_SECTION;
    cpgsci(1);
    bttn_box[I][0] = 0.600;
    bttn_box[I][1] = bttn_box[I][0]+0.070;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    cpgptxt(bttn_box[I][0]-0.003, y_pos + 0.3 * pitch, 0.0, 1.0,
           "Switching Cycle [sec]\0");
    cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
    cpgsci(0);
    cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
            0.0, 1.0, ch_swtim);

    I++;
    cpgsci(1);
    bttn_box[I][0] = 0.880;
    bttn_box[I][1] = bttn_box[I][0]+0.070;
    bttn_box[I][2] = y_pos;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    cpgsci(1);
    cpgptxt(bttn_box[I][0]-0.003, y_pos + 0.3 * pitch, 0.0, 1.0,
            "Target ON ratio (\%)\0");
    cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
    cpgsci(0);
    cpgptxt(bttn_box[I][1]-0.015, 0.6*bttn_box[I][2]+0.4*bttn_box[I][3],
            0.0, 1.0, ch_appon);

/*
--------------------------------------
*/

    y_pos = 0.330;

    cpgsci(1);
    cpgsfs(2);
    cpgsci(1);
    cpgrect(0.010, 0.685, y_pos-0.125, y_pos+0.031);
    cpgsfs(1);

    I = T_FLUX_SECTION;
    bttn_box[I][0] = 0.230;
    bttn_box[I][1] = 0.240;
    bttn_box[I][2] = y_pos + 0.3 * pitch;
    bttn_box[I][3] = bttn_box[I][2] + 0.01;
    if (src[0].positionID == 0) {
      _on_button (&Bdum, "", bttn_box[I]);
    } else if (src[0].positionID == 1) {
      _off_button(&Bdum, "", bttn_box[I]);
    }
    cpgsci(1);
    cpgtext(bttn_box[I][0]+0.02, bttn_box[I][2], "Source-1\0");
    I++;
    bttn_box[I][0] = 0.350;
    bttn_box[I][1] = 0.360;
    bttn_box[I][2] = y_pos + 0.3 * pitch;
    bttn_box[I][3] = bttn_box[I][2] + 0.01;
    if (src[0].positionID == 0) {
      _off_button(&Bdum, "", bttn_box[I]);
    } else if (src[0].positionID == 1) {
      _on_button (&Bdum, "", bttn_box[I]);
    }
    cpgsci(1);
    cpgtext(bttn_box[I][0]+0.02, bttn_box[I][2], "Source-2\0");
    cpgtext(0.035, bttn_box[I][2], "Target Position\0");

    I = T_FLUX_SECTION + FLX_SECT;
    bttn_box[I][0] = 0.260;
    bttn_box[I][1] = 0.330;
    bttn_box[I][2] = bttn_box[T_FLUX_SECTION][2] - 1.3 * pitch;
    bttn_box[I][3] = bttn_box[I][2] + pitch;

    cpgsci(1);
    cpgtext(0.035, bttn_box[I][2] + 0.3 * pitch, "Target Total Flux [Jy]\0");
    cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
    cpgsci(0);
    cpgptxt(bttn_box[I][1]-0.015, text_bottom(bttn_box[I][2], bttn_box[I][3]),
            0.0, 1.0, ch_tflux);
    cpgsci(1);

    for (i=0; i<tgt_morpho_num; i++) {
      I = T_FLUX_SECTION + MRP_SECT + i;
      bttn_box[I][0] = 0.040 + 0.2 * (float)(i%3);
      bttn_box[I][1] = 0.050 + 0.2 * (float)(i%3);
      bttn_box[I][2] = bttn_box[T_FLUX_SECTION][2] - 2.0 * pitch
                       - pitch * (float)(i/3);
      bttn_box[I][3] = bttn_box[I][2] + 0.01;
      if (i == src[0].morphology) {
        _on_button( &tgt_proc_flag[i], "", bttn_box[I]);
      } else {
        _off_button(&tgt_proc_flag[i], "", bttn_box[I]);
      }
      cpgsci(1);
      cpgtext(bttn_box[I][0]+0.02, bttn_box[I][2], tgt_morphology[i]);
    }

    I = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE;
    bttn_box[I][0] = bttn_box[T_FLUX_SECTION+MRP_SECT        ][0];
    bttn_box[I][1] = bttn_box[T_FLUX_SECTION+MRP_SECT+1      ][0] + 0.220;
    bttn_box[I][2] = bttn_box[T_FLUX_SECTION+FLX_SECT+tgt_morpho_num][2] - 0.04;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    if (src[0].morphology == SRC_MULTI_COMP) {
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              0.0, 1.0, ch_file[0].mcm);
      cpgsci(1);
      TV_menu_hatch(0.025, 0.40,
                    bttn_box[T_FLUX_SECTION+FLX_SECT][2],
                    bttn_box[T_FLUX_SECTION+FLX_SECT][3], 7, 4);
    }

    I = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 1;
    bttn_box[I][0] = bttn_box[T_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][0];
    bttn_box[I][1] = bttn_box[T_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][1];
    bttn_box[I][2] = bttn_box[T_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][2];
    bttn_box[I][3] = bttn_box[T_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][3];
    if (src[0].morphology == SRC_CC_COMP) {
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              0.0, 1.0, ch_file[0].cct);
      cpgsci(1);
      TV_menu_hatch(0.025, 0.40,
                    bttn_box[T_FLUX_SECTION+FLX_SECT][2],
                    bttn_box[T_FLUX_SECTION+FLX_SECT][3], 7, 4);
    }

    I = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 2;
    bttn_box[I][0] = bttn_box[T_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][0];
    bttn_box[I][1] = bttn_box[T_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][1];
    bttn_box[I][2] = bttn_box[T_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][2];
    bttn_box[I][3] = bttn_box[T_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][3];
    if (src[0].morphology == SRC_BHS_MOD) {
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              0.0, 1.0, ch_file[0].bhs);
      TV_menu_hatch(0.025, 0.40,
                    bttn_box[T_FLUX_SECTION+FLX_SECT][2],
                    bttn_box[T_FLUX_SECTION+FLX_SECT][3], 7, 4);
    }

/*
-----------------
*/

    cpgsci(1);
    cpgsfs(2);
    cpgsci(1);
    cpgrect(0.715, 1.400, y_pos-0.125, y_pos+0.031);
    cpgsfs(1);

    I = R_FLUX_SECTION;
    bttn_box[I][0] = 0.930;
    bttn_box[I][1] = 0.940;
    bttn_box[I][2] = y_pos + 0.3 * pitch;
    bttn_box[I][3] = bttn_box[I][2] + 0.01;
    if (src[1].positionID == 0) {
      _on_button (&Bdum, "", bttn_box[I]);
    } else if (src[1].positionID == 1) {
      _off_button(&Bdum, "", bttn_box[I]);
    }
    cpgsci(1);
    cpgtext(bttn_box[I][0]+0.02, bttn_box[I][2], "Source-1\0");
    I++;
    bttn_box[I][0] = 1.050;
    bttn_box[I][1] = 1.060;
    bttn_box[I][2] = y_pos + 0.3 * pitch;
    bttn_box[I][3] = bttn_box[I][2] + 0.01;
    if (src[1].positionID == 0) {
      _off_button(&Bdum, "", bttn_box[I]);
    } else if (src[1].positionID == 1) {
      _on_button (&Bdum, "", bttn_box[I]);
    }
    cpgsci(1);
    cpgtext(bttn_box[I][0]+0.02, bttn_box[I][2], "Source-2\0");
    cpgtext(0.735, bttn_box[I][2], "Reference Position\0");

    I = R_FLUX_SECTION + FLX_SECT;
    bttn_box[I][0] = 0.960;
    bttn_box[I][1] = 1.030;
    bttn_box[I][2] = bttn_box[R_FLUX_SECTION][2] - 1.3 * pitch;
    bttn_box[I][3] = bttn_box[I][2] + pitch;

    cpgtext(0.735, bttn_box[I][2] + 0.3 * pitch, "Reference Total Flux [Jy]\0");
    cpgsci(1);
    cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
    cpgsci(0);
    cpgptxt(bttn_box[I][1]-0.015, text_bottom(bttn_box[I][2], bttn_box[I][3]),
            0.0, 1.0, ch_rflux);
    cpgsci(1);

    for (i=0; i<ref_morpho_num; i++) {
      I = R_FLUX_SECTION + MRP_SECT + i;
      bttn_box[I][0] = 0.740 + 0.2 * (float)(i%3);
      bttn_box[I][1] = 0.750 + 0.2 * (float)(i%3);
      bttn_box[I][2] = bttn_box[R_FLUX_SECTION][2] - 2.0 * pitch
                       - pitch * (float)(i/3);
      bttn_box[I][3] = bttn_box[I][2] + 0.01;
      if (i == src[1].morphology) {
        _on_button( &ref_proc_flag[i], "", bttn_box[I]);
      } else {
        _off_button(&ref_proc_flag[i], "", bttn_box[I]);
      }
      cpgsci(1);
      cpgtext(bttn_box[I][0]+0.02, bttn_box[I][2], ref_morphology[i]);
    }

    I = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE;
    bttn_box[I][0] = bttn_box[R_FLUX_SECTION+MRP_SECT        ][0];
    bttn_box[I][1] = bttn_box[R_FLUX_SECTION+MRP_SECT+1      ][0] + 0.220;
    bttn_box[I][2] = bttn_box[R_FLUX_SECTION+FLX_SECT+tgt_morpho_num][2] - 0.04;
    bttn_box[I][3] = bttn_box[I][2] + pitch;
    if (src[1].morphology == SRC_MULTI_COMP) {
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              0.0, 1.0, ch_file[1].mcm);
      cpgsci(1);
      TV_menu_hatch(0.725, 0.90,
                    bttn_box[R_FLUX_SECTION+FLX_SECT][2],
                    bttn_box[R_FLUX_SECTION+FLX_SECT][3], 7, 4);
    }

    I = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 1;
    bttn_box[I][0] = bttn_box[R_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][0];
    bttn_box[I][1] = bttn_box[R_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][1];
    bttn_box[I][2] = bttn_box[R_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][2];
    bttn_box[I][3] = bttn_box[R_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][3];
    if (src[1].morphology == SRC_CC_COMP) {
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              0.0, 1.0, ch_file[1].cct);
      TV_menu_hatch(0.725, 0.90,
                    bttn_box[R_FLUX_SECTION+FLX_SECT][2],
                    bttn_box[R_FLUX_SECTION+FLX_SECT][3], 7, 4);
    }

    I = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 2;
    bttn_box[I][0] = bttn_box[R_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][0];
    bttn_box[I][1] = bttn_box[R_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][1];
    bttn_box[I][2] = bttn_box[R_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][2];
    bttn_box[I][3] = bttn_box[R_FLUX_SECTION+3+SRC_MORPH_CANDIDATE][3];
    if (src[1].morphology == SRC_BHS_MOD) {
      cpgsci(1);
      cpgrect(bttn_box[I][0], bttn_box[I][1], bttn_box[I][2], bttn_box[I][3]);
      cpgsci(0);
      cpgptxt(bttn_box[I][1]-0.015, text_bottom(bttn_box[I][2], bttn_box[I][3]),
              0.0, 1.0, ch_file[0].bhs);
      TV_menu_hatch(0.725, 0.90,
                    bttn_box[R_FLUX_SECTION+FLX_SECT][2],
                    bttn_box[R_FLUX_SECTION+FLX_SECT][3], 7, 4);
    }
    cpgebuf();

/*
------------------------------------
*/

    while (1) {
      cpgcurs(cursor_pos, cursor_pos+1, string);

      if (_button_chk(cursor_pos, ctrl_bttn_box[0]) == true) {
        _on_button(&Bdum, "GO (baseline_base)\0", ctrl_bttn_box[0]);
        *reference_phase_process_mode = 0;
        break;
      } else if (_button_chk(cursor_pos, ctrl_bttn_box[1]) == true) {
        _on_button(&Bdum, "GO (antenna_base)\0",  ctrl_bttn_box[1]);
        *reference_phase_process_mode = 1;
        break;
      } else if (_button_chk(cursor_pos, ctrl_bttn_box[2]) == true) {
        _on_button(&Bdum, "EXIT (without saving parameters)\0",
                  ctrl_bttn_box[2]);
        sprintf(string, "Exit ARIS. Good by!");
        comment_disp(cmnt, comment, string, true);
        return (-1);
      }

/*
---------------------------
*/

      for (I=0; I<blsel_num; I++) {
        if (_button_chk(cursor_pos, blsel_bttn_box[I]) == true) {
          BL_SELECT = I;
          for (i=0; i<blsel_num; i++) {
            if (i == BL_SELECT) {
              _on_button (&blsel_status[i].code, (blsel_status+i)->name,
                          blsel_bttn_box[I]);
            } else {
              J = i;
              _off_button(&blsel_status[i].code, (blsel_status+i)->name,
                          blsel_bttn_box[J]);
            }
          }
        }
      }

/*
---------------------------
*/

      if (ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true) {
        for (I=TROPOS_SECTION; I<TROPOS_SECTION+trp_num; I++) {
          if (_button_chk(cursor_pos, bttn_box[I]) == true) {
            *TRP_CONDITION = I - TROPOS_SECTION;
            for (i=0; i<trp_num; i++) {
              if (i == *TRP_CONDITION) {
               _on_button(&trp_status[i].code, (trp_status+i)->name,
                          bttn_box[I]);
              } else {
                J = TROPOS_SECTION + i;
                _off_button(&trp_status[i].code, (trp_status+i)->name,
                            bttn_box[J]);
              }
            }

            if (*TRP_CONDITION < trp_status_num) {
              break;
            } else if (*TRP_CONDITION == trp_tune) {

/** INDIVIDUAL TUNING ON **/

              if ((ant_trp = (struct ant_trp_tune *)
                calloc(GRT_NUM, sizeof(struct ant_trp_tune))) == NULL) {
                printf("ERROR: err_parameter: ");
                printf("tropospheric tuning parameter memry allocation.\n");
                exit (-1);
              } else {
                for (iant=0; iant<GRT_NUM; iant++) {
                  strcpy((ant_trp+iant)->IDC, (ant_prm+iant)->IDC);
                }
              }

              if ((fp = fopen( "aris_input/ant_trp_tuning.prm", "r")) != NULL) {
                while (1) {
                  if (fgets(string, sizeof(string), fp) == NULL) {
                    break;
                  } else {
                    sscanf(string, "%s %d", &ant_tmp, &idum);
                    for (iant=0; iant<GRT_NUM; iant++) {
                      if (strncmp((ant_trp+iant)->IDC, ant_tmp,
                                  strlen(ant_tmp)) == 0) {
                        ant_trp[iant].code = idum;
                        break;
                      }
                    }
                  }
                }
                fclose (fp);
              } else {
                for (iant=0; iant<GRT_NUM; iant++) {
                  ant_trp[iant].code = 0;
                }
              }

              cpgid2 = (int)cpgopen("/xs");
              if (cpgid2 < 0) {
                cpgask(-1);
              }
              cpgslct(cpgid2);

              cpgpap(2.00*pgpap_prm, 1.0/1.4);
              cpgsch(1.5*pgpap_prm/13.0);
              cpgsvp(0.0,  1.0, 0.0, 1.0);
              cpgswin(0.0, 1.4, 0.0, 1.0);

              for (i=0; i<trp_status_num; i++) {
                cpgptxt(0.335+(float)i*0.150, 0.950, 0.0, 0.5,
                        (trp_status+i)->name);
              }

              for (iant=0; iant<GRT_NUM; iant++) {
                ant_bttn_box[iant][0] = 0.010 + 0.097;
                ant_bttn_box[iant][1] = 0.110 + 0.097;
                ant_bttn_box[iant][2] = 0.90 - pitch * (float)iant;
                ant_bttn_box[iant][3] = ant_bttn_box[iant][2] + pitch;
                _off_button(&Bdum, (ant_trp+iant)->IDC,
                                    &ant_bttn_box[iant][0]);
                for (i=0; i<trp_status_num; i++) {
                  ftmp1 = 0.325 + (float)i*0.150;
                  ftmp2 = ant_bttn_box[iant][2] + 0.2 * pitch;
                  trp_bttn_box[iant][i][0] = ftmp1;
                  trp_bttn_box[iant][i][1] = ftmp1 + 0.02;
                  trp_bttn_box[iant][i][2] = ftmp2;
                  trp_bttn_box[iant][i][3] = ftmp2 + 0.02;

                  if (ant_trp[iant].code == i) {
                    _on_button(&Bdum, "", trp_bttn_box[iant][i]);
                  } else {
                    _off_button(&Bdum, "", trp_bttn_box[iant][i]);
                  }
                }
              }

              ant_bttn_box[GRT_NUM][0] = 0.600;
              ant_bttn_box[GRT_NUM][1] = 0.660;
              ant_bttn_box[GRT_NUM][2] = 0.050;
              ant_bttn_box[GRT_NUM][3] = ant_bttn_box[GRT_NUM][2] + pitch;
              _off_button(&Bdum, "QUIT", &ant_bttn_box[GRT_NUM][0]);

              QUIT_SWT = false;
              while (1) {
                cpgcurs(cursor_pos, cursor_pos+1, string);
                if (_button_chk(cursor_pos, ant_bttn_box[GRT_NUM]) == true) {
                  QUIT_SWT = true;
                  break;
                } else {
                  for (iant=0; iant<GRT_NUM; iant++) {
                    for (i=0; i<trp_status_num; i++) {
                      if (_button_chk(cursor_pos,
                            trp_bttn_box[iant][i]) == true) {
                        for (j=0; j<trp_status_num; j++) {
                          if (j == i) {
                            _on_button(&Bdum, "", trp_bttn_box[iant][j]);
                          } else {
                            _off_button(&Bdum, "", trp_bttn_box[iant][j]);
                          }
                        }
                        ant_trp[iant].code = i;
                        break;
                      }
                    }
                  }
                }
              }

              if ((fp = fopen(
                    "aris_input/ant_trp_tuning.prm", "w")) != NULL) {
                for (iant=0; iant<GRT_NUM; iant++) {
                  fprintf(fp, "%s %d\n",
                       (ant_trp+iant)->IDC, ant_trp[iant].code);
                }
                fclose (fp);
              }
              cpgclos();
              cpgslct(cpgid1);

              for (iant=0; iant<GRT_NUM; iant++) {
                ant_prm[iant].WVturb = ant_trp[iant].code;
              }
              free (ant_trp);
            }
          }
        }
      }

/*
---------------------------
*/

      if (ERROR_FLAG[IONTRB] == true) {
        for (I=IONOS_SECTION; I<IONOS_SECTION+ion_num; I++) {
          if (_button_chk(cursor_pos, bttn_box[I]) == true) {
            *ION_CONDITION = I - IONOS_SECTION;
            for (i=0; i<ion_num; i++) {
              if (i == *ION_CONDITION) {
                _on_button(&ion_status[i].code,  (ion_status+i)->name,
                           bttn_box[I]);
              } else {
                J = IONOS_SECTION + i;
                _off_button(&ion_status[i].code, (ion_status+i)->name,
                            bttn_box[J]);
              }
            }
          }
        }
      }

/*
---------------------------
*/

      if (SRT_NUM >= 1) {
        for (I=TRK_SECTION; I<TRK_SECTION+*TRK_NUM+TRK_GEN_CMD; I++) {
          if (_button_chk(cursor_pos, bttn_box[I]) == true) {

            if (I == TRK_SECTION) {
              for (itrk=0; itrk<*TRK_NUM; itrk++) {
                J = TRK_SECTION + TRK_GEN_CMD + itrk;
                trk_priority[itrk] = 1;
                sprintf(string, "%s(%d)", trk_name[itrk], trk_priority[itrk]);
                _on_button(&Bdum, string, bttn_box[J]);
              }
            } else if (I == TRK_SECTION + 1) {
              for (itrk=0; itrk<*TRK_NUM; itrk++) {
                J = TRK_SECTION + TRK_GEN_CMD + itrk;
                trk_priority[itrk] = 0;
                _off_button(&Bdum, trk_name[itrk], bttn_box[J]);
              }
            } else {
              itrk = I - TRK_GEN_CMD - TRK_SECTION;
              if (trk_priority[itrk] >= 1) {
                trk_priority[itrk] = 0;
              } else {
                trk_priority[itrk] = *TRK_NUM;
              }
              trk_priority_check(*TRK_NUM, trk_priority, trk_name);
              I = TRK_SECTION + TRK_GEN_CMD;
              for (itrk=0; itrk<*TRK_NUM; itrk++) {
                if (trk_priority[itrk] >= 1) {
                  sprintf(string, "%s(%d)", trk_name[itrk], trk_priority[itrk]);
                  _on_button(&Bdum, string, bttn_box[I]);
                } else {
                  _off_button(&Bdum, trk_name[itrk], bttn_box[I]);
                }
                I++;
              }

#ifdef __DEBUG__
              for (j=1; j<=*TRK_NUM; j++) {
                for (k=0; k<*TRK_NUM; k++) {
                  if (trk_priority[k] == j) {
                    printf("__DEBUG__ %d. %s\n", j, trk_pos[k].IDC);
                  }
                }
              }
              printf("\n");
#endif /* __DEBUG__ */
            }

            trk_num = 0;
            for (itrk=0; itrk<*TRK_NUM; itrk++) {
              if (trk_priority[itrk] >= 1) {
                trk_num++;
              }
            }

            srtatt_button_disp(SRT_NUM, srtatt_num, srtatt_code, trk_num,
                               srtatt_y_pos, pitch, srtatt_name,
                               bttn_box+SRTATT_SECTION);
          }
        }
      }

/*
---------------------------
*/

/********
      for (I=OBSTIM_SECTION; I<OBSTIM_SECTION+2; I++) {
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_timer[I-OBSTIM_SECTION], 0.0, 0.0);
          sscanf(ch_timer[I-OBSTIM_SECTION], "%d/%d:%d:%d", &id, &ih, &im, &is);
        }
      }
********/

/*
---------------------------
*/

      for (I=T_WAVE_SECTION; I<T_WAVE_SECTION+wave_num; I++) {
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          wave_id[0] = wave_code[I - T_WAVE_SECTION];
          wave_id[1] = wave_id[0];
          for (i=0; i<wave_num; i++) {
            if (wave_id[0] == wave_code[i]) {
              _on_button(&wave_swt[0][i], wave_name__on[i], bttn_box[I]);
              J = R_WAVE_SECTION + i;
              _on_button(&wave_swt[1][i], wave_name__on[i], bttn_box[J]);
            } else {
              J = T_WAVE_SECTION + i;
              _off_button(&wave_swt[0][i], wave_name_off[i], bttn_box[J]);
              J = R_WAVE_SECTION + i;
              _off_button(&wave_swt[1][i], wave_name_off[i], bttn_box[J]);
            }
          }
        }
      }

/*
---------------------------
*/

      for (I=R_WAVE_SECTION; I<R_WAVE_SECTION+wave_num; I++) {
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          wave_id[1] = wave_code[I - R_WAVE_SECTION];
          for (i=0; i<wave_num; i++) {
            if (wave_id[1] == wave_code[i]) {
              _on_button(&wave_swt[1][i], wave_name__on[i], bttn_box[I]);
            } else {
              J = R_WAVE_SECTION + i;
              _off_button(&wave_swt[1][i], wave_name_off[i], bttn_box[J]);
            }
          }
        }
      }

/*
---------------------------
*/

      for (I=BW_SECTION; I<BW_SECTION+bw_num; I++) {
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          bw_id = I - BW_SECTION;
          for (i=0; i<bw_num; i++) {
            if (i == bw_id) {
              _on_button(&bw_code[i], bw_name[i], bttn_box[I]);
            } else {
              J = BW_SECTION + i;
              _off_button(&bw_code[i], bw_name[i], bttn_box[J]);
            }
          }
        }
      }

/*
---------------------------
*/

      for (I=FRCHAN_SECTION; I<FRCHAN_SECTION+frchan_num; I++) {
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          frchan_id = I - FRCHAN_SECTION;
          for (i=0; i<frchan_num; i++) {
            if (i == frchan_id) {
              _on_button(&frchan_code[i], frchan_name[i], bttn_box[I]);
            } else {
              J = FRCHAN_SECTION + i;
              _off_button(&frchan_code[i], frchan_name[i], bttn_box[J]);
            }
          }
          *FRCHAN_NUM = channel_num_set(frchan_id);
        }
      }

/*
---------------------------
*/

      for (I=CFACT_SECTION; I<CFACT_SECTION+cfact_num; I++) {
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          cfact_id = I - CFACT_SECTION;
          for (i=0; i<cfact_num; i++) {
            if (i == cfact_id) {
              _on_button(&cfact_code[i], cfact_name[i], bttn_box[I]);
            } else {
              J = CFACT_SECTION + i;
              _off_button(&cfact_code[i], cfact_name[i], bttn_box[J]);
            }
          }
        }
      }

/*
---------------------------
*/

      if (ERROR_FLAG[TDSECZ] == true && GRT_NUM > 0
       && array->TYPE == _VLBI_ARRAY_) {
        I = TROPOS_SECTION + trp_num;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_tdscz, 0.0, 0.0);
          sscanf(ch_tdscz, "%lf", &tdscz);
        }
      }

/*
---------------------------
*/

      if (ERROR_FLAG[IDSECZ] == true && GRT_NUM > 0) {
        I = IONOS_SECTION + ion_num;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_idscz, 0.0, 0.0);
          sscanf(ch_idscz, "%lf", &idscz);
        }
      }

/*
---------------------------
*/

      if (ERROR_FLAG[APOSER] == true && SRT_NUM > 0) {
        I = ORBIT_SECTION;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_orbae, 0.0, 0.0);
          sscanf(ch_orbae, "%lf", &orbit_error_ap);
          orbit_error_ap *= 1.0e-2;
        }

        I = ORBIT_SECTION + 1;
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          tv_get_param("double", cursor_pos, bttn_box[I],
                       pitch, ch_orbpe, 0.0, 0.0);
          sscanf(ch_orbpe, "%lf", &orbit_error_pe);
          orbit_error_pe *= 1.0e-2;
        }
      }

/*
---------------------------
*/

      if (SRT_NUM >= 1 && trk_num >= 1) {
        for (iant=0; iant<SRT_NUM; iant++) {
          for (I=SRTATT_SECTION; I<SRTATT_SECTION+srtatt_num; I++) {
            if (_button_chk(cursor_pos, bttn_box[srtatt_num*iant+I]) == true) {
              srtatt_code[iant] = I - SRTATT_SECTION;
              for (i=0; i<srtatt_num; i++) {
                if (i == srtatt_code[iant]) {
                  _on_button (&Bdum, srtatt_name[i],
                                    bttn_box[srtatt_num*iant+I]);
                } else {
                  J = SRTATT_SECTION + i;
                  _off_button(&Bdum, srtatt_name[i],
                                    bttn_box[srtatt_num*iant+J]);
                }
              }
            }
          }
        }
      }

/*
---------------------------
*/

      if (_button_chk(cursor_pos, bttn_box[T_FLUX_SECTION  ]) == true) {
        src[0].positionID = 0;
        _on_button (&Bdum, "", bttn_box[T_FLUX_SECTION  ]);
        _off_button(&Bdum, "", bttn_box[T_FLUX_SECTION+1]);
      }

      if (_button_chk(cursor_pos, bttn_box[T_FLUX_SECTION+1]) == true) {
        src[0].positionID = 1;
        _off_button(&Bdum, "", bttn_box[T_FLUX_SECTION  ]);
        _on_button (&Bdum, "", bttn_box[T_FLUX_SECTION+1]);
      }

/*
----
*/

      if (_button_chk(cursor_pos, bttn_box[T_FLUX_SECTION+FLX_SECT]) == true &&
          src[0].morphology != SRC_CC_COMP &&
          src[0].morphology != SRC_BHS_MOD) {
        tv_get_param("double", cursor_pos, bttn_box[T_FLUX_SECTION+FLX_SECT],
                     pitch, ch_tflux, 0.0, 0.0);
        sscanf(ch_tflux, "%lf", &src[0].flux);
      }

/*
----
*/

      for (I=T_FLUX_SECTION+MRP_SECT;
           I<T_FLUX_SECTION+MRP_SECT+tgt_morpho_num; I++) {
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          src[0].morphology = I - T_FLUX_SECTION - MRP_SECT;
          for (i=0; i<tgt_morpho_num; i++) {
            J = i + T_FLUX_SECTION + MRP_SECT;
            if (i == src[0].morphology) {
              _on_button( &tgt_proc_flag[i], "", bttn_box[J]);
            } else {
              _off_button(&tgt_proc_flag[i], "", bttn_box[J]);
            }
          }
          if (src[0].morphology == SRC_MULTI_COMP) {
            J = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE;
            cpgsci(1);
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);
            cpgsci(0);
            cpgptxt(bttn_box[J][1]-0.015,
                    text_bottom(bttn_box[J][2], bttn_box[J][3]),
                    0.0, 1.0, ch_file[0].mcm);
            cpgsci(1);
            TV_menu_hatch(0.025, 0.40,
                          bttn_box[T_FLUX_SECTION+FLX_SECT][2],
                          bttn_box[T_FLUX_SECTION+FLX_SECT][3], 7, 4);
          } else if (src[0].morphology == SRC_CC_COMP) {
            J = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 1;
            cpgsci(1);
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);
            cpgsci(0);
            cpgptxt(bttn_box[J][1]-0.015,
                    text_bottom(bttn_box[J][2], bttn_box[J][3]),
                    0.0, 1.0, ch_file[0].cct);
            TV_menu_hatch(0.025, 0.40,
                          bttn_box[T_FLUX_SECTION+FLX_SECT][2],
                          bttn_box[T_FLUX_SECTION+FLX_SECT][3], 7, 4);
          } else if (src[0].morphology == SRC_BHS_MOD) {
            J = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 2;
            cpgsci(1);
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);
            cpgsci(0);
            cpgptxt(bttn_box[J][1]-0.015,
                    text_bottom(bttn_box[J][2], bttn_box[J][3]),
                    0.0, 1.0, ch_file[0].bhs);
            TV_menu_hatch(0.025, 0.40,
                          bttn_box[T_FLUX_SECTION+FLX_SECT][2],
                          bttn_box[T_FLUX_SECTION+FLX_SECT][3], 7, 4);
          }

          if (src[0].morphology == SRC_POINT         ||
              src[0].morphology == SRC_DISK_JETCJET  ||
              src[0].morphology == SRC_DISK_VSOP2JET) {
            J = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE;
            cpgsci(0);
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);

            J = T_FLUX_SECTION + FLX_SECT;
            TV_menu_hatch(0.025, 0.40,
                          bttn_box[J][2], bttn_box[J][3], 0, 1);
            cpgtext(0.035, bttn_box[J][2] + 0.3 * pitch,
                    "Target Total Flux [Jy]\0");
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);
            cpgsci(0);
            cpgptxt(bttn_box[J][1]-0.015,
                    text_bottom(bttn_box[J][2], bttn_box[J][3]),
                    0.0, 1.0, ch_tflux);
            cpgsci(1);
          }
          break;
        }
      }

/*
----
*/

      J = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE;
      if (src[0].morphology == SRC_MULTI_COMP &&
          _button_chk(cursor_pos, bttn_box[J]) == true) {
        tv_get_param("char", cursor_pos, bttn_box[J],
                     pitch, ch_file[0].mcm, 0.0, 0.0);
      }

      J = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 1;
      if (src[0].morphology == SRC_CC_COMP &&
          _button_chk(cursor_pos, bttn_box[J]) == true) {
        tv_get_param("char", cursor_pos, bttn_box[J],
                     pitch, ch_file[0].cct, 0.0, 0.0);
      }

      J = T_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 2;
      if (src[0].morphology == SRC_BHS_MOD &&
          _button_chk(cursor_pos, bttn_box[J]) == true) {
        tv_get_param("char", cursor_pos, bttn_box[J],
                     pitch, ch_file[0].bhs, 0.0, 0.0);
      }

/*
---------------------------
*/

      if (_button_chk(cursor_pos, bttn_box[R_FLUX_SECTION  ]) == true) {
        src[1].positionID = 0;
        _on_button (&Bdum, "", bttn_box[R_FLUX_SECTION  ]);
        _off_button(&Bdum, "", bttn_box[R_FLUX_SECTION+1]);
      }

      if (_button_chk(cursor_pos, bttn_box[R_FLUX_SECTION+1]) == true) {
        src[1].positionID = 1;
        _off_button(&Bdum, "", bttn_box[R_FLUX_SECTION  ]);
        _on_button (&Bdum, "", bttn_box[R_FLUX_SECTION+1]);
      }

/*
----
*/

      if (_button_chk(cursor_pos, bttn_box[R_FLUX_SECTION+FLX_SECT]) == true &&
          src[1].morphology != SRC_CC_COMP &&
          src[1].morphology != SRC_BHS_MOD) {
        tv_get_param("double", cursor_pos, bttn_box[R_FLUX_SECTION+FLX_SECT],
                     pitch, ch_rflux, 0.0, 0.0);
        sscanf(ch_rflux, "%lf", &src[1].flux);
      }

/*
----
*/

      for (I=R_FLUX_SECTION+MRP_SECT;
           I<R_FLUX_SECTION+MRP_SECT+ref_morpho_num; I++) {
        if (_button_chk(cursor_pos, bttn_box[I]) == true) {
          src[1].morphology = I - R_FLUX_SECTION - MRP_SECT;
          for (i=0; i<ref_morpho_num; i++) {
            J = i + R_FLUX_SECTION + MRP_SECT;
            if (i == src[1].morphology) {
              _on_button( &ref_proc_flag[i], "", bttn_box[J]);
            } else {
              _off_button(&ref_proc_flag[i], "", bttn_box[J]);
            }
          }
          if (src[1].morphology == SRC_MULTI_COMP) {
            J = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE;
            cpgsci(1);
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);
            cpgsci(0);
            cpgptxt(bttn_box[J][1]-0.015,
                    text_bottom(bttn_box[J][2], bttn_box[J][3]),
                    0.0, 1.0, ch_file[1].mcm);
            cpgsci(1);
            TV_menu_hatch(0.725, 1.10,
                          bttn_box[R_FLUX_SECTION+FLX_SECT][2],
                          bttn_box[R_FLUX_SECTION+FLX_SECT][3], 7, 4);
          } else if (src[1].morphology == SRC_CC_COMP) {
            J = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 1;
            cpgsci(1);
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);
            cpgsci(0);
            cpgptxt(bttn_box[J][1]-0.015,
                    text_bottom(bttn_box[J][2], bttn_box[J][3]),
                    0.0, 1.0, ch_file[1].cct);
            TV_menu_hatch(0.725, 1.10,
                          bttn_box[R_FLUX_SECTION+FLX_SECT][2],
                          bttn_box[R_FLUX_SECTION+FLX_SECT][3], 7, 4);
          } else if (src[1].morphology == SRC_BHS_MOD) {
            J = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 2;
            cpgsci(1);
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);
            cpgsci(0);
            cpgptxt(bttn_box[J][1]-0.015,
                    text_bottom(bttn_box[J][2], bttn_box[J][3]),
                    0.0, 1.0, ch_file[1].bhs);
            TV_menu_hatch(0.725, 1.10,
                          bttn_box[R_FLUX_SECTION+FLX_SECT][2],
                          bttn_box[R_FLUX_SECTION+FLX_SECT][3], 7, 4);
          }

          if (src[1].morphology == SRC_POINT         ||
              src[1].morphology == SRC_DISK_JETCJET  ||
              src[1].morphology == SRC_DISK_VSOP2JET) {
            cpgsci(0);
            J = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE;
            cpgsci(0);
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);

            J = R_FLUX_SECTION + FLX_SECT;
            TV_menu_hatch(0.725, 1.10,
                          bttn_box[J][2], bttn_box[J][3], 0, 1);
            cpgtext(0.735, bttn_box[J][2] + 0.3 * pitch,
                    "Reference Total Flux [Jy]\0");
            cpgrect(bttn_box[J][0], bttn_box[J][1],
                    bttn_box[J][2], bttn_box[J][3]);
            cpgsci(0);
            cpgptxt(bttn_box[J][1]-0.015,
                    text_bottom(bttn_box[J][2], bttn_box[J][3]),
                    0.0, 1.0, ch_rflux);
            cpgsci(1);
          }
          break;
        }
      }

/*
----
*/

      J = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE;
      if (src[1].morphology == SRC_MULTI_COMP &&
          _button_chk(cursor_pos, bttn_box[J]) == true) {
        tv_get_param("char", cursor_pos, bttn_box[J],
                     pitch, ch_file[1].mcm, 0.0, 0.0);
      }

      J = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 1;
      if (src[1].morphology == SRC_CC_COMP &&
          _button_chk(cursor_pos, bttn_box[J]) == true) {
        tv_get_param("char", cursor_pos, bttn_box[J],
                     pitch, ch_file[1].cct, 0.0, 0.0);
      }

      J = R_FLUX_SECTION + MRP_SECT + SRC_MORPH_CANDIDATE + 2;
      if (src[1].morphology == SRC_BHS_MOD &&
          _button_chk(cursor_pos, bttn_box[J]) == true) {
        tv_get_param("char", cursor_pos, bttn_box[J],
                     pitch, ch_file[1].bhs, 0.0, 0.0);
      }

/*
---------------------------
*/

      if (_button_chk(cursor_pos, bttn_box[SWTCYC_SECTION]) == true) {
        tv_get_param("double", cursor_pos, bttn_box[SWTCYC_SECTION],
                     pitch, ch_swtim, 0.0, 0.0);
        sscanf(ch_swtim, "%lf", &swt_cyc_time);

/*
--------
*/
        *nswt    = (int)lrint(swt_cyc_time);
        if (*nswt % 2 == 1) {
          (*nswt)++;
          sprintf(string,
                  "Swtching cycle time is changed to %d [sec].", *nswt);
          comment_disp(cmnt, comment, string, true);
        }
        swt_cyc_time = (double)(*nswt);

        for (i=0; i<sizeof(string); i++) {
          string[i] = 0;
        }
      }

/*
---------------------------
*/

      if (_button_chk(cursor_pos, bttn_box[SWTCYC_SECTION+1]) == true) {
        tv_get_param("double", cursor_pos, bttn_box[SWTCYC_SECTION+1],
                     pitch, ch_appon, 0.0, 0.0);
        sscanf(ch_appon, "%lf", apparent_tgt_on);
      }

/*
---------------------------
*/

    }
  }

/*
------------------------------------------------------
*/

  sprintf(string, "ARIS continues...");
  if (TV_SWT == false) {
    printf("%s\n", string);
  } else if (TV_SWT == true) {
    comment_disp(cmnt, comment, string, true);
  }

/*
------------------------------------------------------
*/

  if (BL_SELECT == ALLBL) {
    *BGN_ANT_I = 0;
    *END_ANT_I = ANT_NUM;
    *BGN_ANT_J = 0;
    *END_ANT_J = ANT_NUM;
  } else if (BL_SELECT == GG_BL) {
    *BGN_ANT_I = 0;
    *END_ANT_I = GRT_NUM;
    *BGN_ANT_J = 0;
    *END_ANT_J = GRT_NUM;
  } else if (BL_SELECT == GS_BL) {
    *BGN_ANT_I = 0;
    *END_ANT_I = GRT_NUM;
    *BGN_ANT_J = GRT_NUM;
    *END_ANT_J = ANT_NUM;
  } else if (BL_SELECT == SS_BL) {
    *BGN_ANT_I = GRT_NUM;
    *END_ANT_I = ANT_NUM;
    *BGN_ANT_J = GRT_NUM;
    *END_ANT_J = ANT_NUM;
  }

/*
------
*/

  for (iant=0; iant<GRT_NUM; iant++) {
    if (ant_prm[iant].WVturb == 0) {
      Cw[iant] = pow(sqrt(pow(wvc[0].v[0], 2.0)
                    + pow(wvc[0].v[1], 2.0)) * 120.0,
                                                   -0.5*wvc[0].o_expon)
              / sqrt(wvc[0].i_scale[0])   / 12.0;
      Cw[iant] *= 4.0;
    } else if (ant_prm[iant].WVturb == 1) {
      Cw[iant] = pow(sqrt(pow(wvc[0].v[0], 2.0) + pow(wvc[0].v[1], 2.0)) * 120.0,
                                                   -0.5*wvc[0].o_expon)
              / sqrt(wvc[0].i_scale[0])   /  9.0;
      Cw[iant] *= 4.0;
    } else if (ant_prm[iant].WVturb == 2) {
      Cw[iant] = pow(sqrt(pow(wvc[0].v[0], 2.0) + pow(wvc[0].v[1], 2.0)) * 120.0,
                                                   -0.5*wvc[0].o_expon)
              / sqrt(wvc[0].i_scale[0])   /  6.0;
      Cw[iant] *= 4.0;
    } else if (ant_prm[iant].WVturb == 3) {
      Cw[iant] = 0.5e-7 * sqrt(1.4 * 1.0e3);
      Cw[iant] = 0.189e-3 / pow(1.0e3, 0.5) / pow(10.0e3, 1.0/3.0) * 10.0;
    } else if (ant_prm[iant].WVturb == 4) {
      Cw[iant] = 1.0e-7 * sqrt(1.4 * 1.0e3);
      Cw[iant] = 2.835e-3 / pow(1.0e3, 0.5) / pow(10.0e3, 1.0/3.0) * 10.0;
    } else if (ant_prm[iant].WVturb == 5) {
      Cw[iant] = 2.0e-7 * sqrt(1.4 * 1.0e3);
    } else if (ant_prm[iant].WVturb == 6) {
      Cw[iant] = 4.0e-7 * sqrt(1.4 * 1.0e3);
    }
    Cw[iant] /= speed_of_light;
  }

  /*                                                        */
  /*                                                        */
  /*  1.67 X 10^{-16/3} [m^{1-5/6}] :  Asaki et al., 1996   */
  /*         ----> Cn = 2.1e-7 [m^{-1/3}]                   */
  /*  0.67 X 10^{-16/3} [m^{1-5/6}] :  Asaki et al., 1998   */
  /*         ----> Cn = 0.8e-7 [m^{-1/3}]                   */
  /*                                                        */
  /*  In Beasley and Conway, the coefficient is defined by  */
  /*  CW^2 = 1.4*Cn^2*L, where L is the inner scale. Then,  */
  /*  a factor sqrt(1.4 * 1000[m]) has to be multiplied     */
  /*  to change their defined coefficient to the exact      */
  /*  structure coefficient defined Drvskii and Finkelstein.*/
  /*                                                        */
  /*  VERY GOOD SITUATION:                                  */
  /*  The value is described in Asaki et al.                */
  /*  (ALMA Memo, No. 535, 2005). [50-% condition]          */
  /*                                                        */
  /*                                                        */

/*
------
*/

  if (*ION_CONDITION == 0) {
    *CI = 0.5;
  } else if (*ION_CONDITION == 1) {
    *CI = 1.0;
  } else if (*ION_CONDITION == 2) {
    *CI = 2.0;
  }

/*
------
*/

  if (SRT_NUM > 0) {
    if (ERROR_FLAG[APOSER] == true) {
      for (iant=0; iant<SRT_NUM; iant++) {
        srt[iant].ODDA = orbit_error_ap;
        srt[iant].ODDP = orbit_error_pe;
      }
    } else {
      for (iant=0; iant<SRT_NUM; iant++) {
        srt[iant].ODDA = 0.0;
        srt[iant].ODDP = 0.0;
      }
    }

    for (iant=0; iant<SRT_NUM; iant++) {
      if (srtatt_code[iant] == 0) {
        srt[iant].BODY_X_SUN = +1;
      } else if (srtatt_code[iant] == 1) {
        srt[iant].BODY_X_SUN = -1;
      } else if (srtatt_code[iant] == 2) {
        srt[iant].BODY_X_SUN =  0;
      }
    }
  }

/*
------
*/

  if (SRT_NUM >= 1) {
    for (itrk=0; itrk<*TRK_NUM; itrk++) {
      trk_pos[itrk].UFL      = false;
      trk_pos[itrk].priority = 0;
    }

    trk_num = 0;
    for (I=1; I<=*TRK_NUM; I++) {
      for (itrk=0; itrk<*TRK_NUM; itrk++) {
        if (trk_priority[itrk] == I) {
          array_config(0,      ALL_ANT,          wave_id,       0,     &idum,
                       trk_pos+trk_num,
                       trk_name[itrk],
                       "aris_input/tracking_network.prm",  false,   false);
          trk_pos[trk_num].UFL      = true;
          trk_pos[trk_num].priority = trk_priority[itrk];
#ifdef __DEBUG__
          printf(__DEBUG__ "%d. %s\n",
                 trk_pos[trk_num].priority,  trk_pos[trk_num].IDC);
#endif /* __DEBUG__ */
          trk_num++;
        }
      }
    }
    *TRK_NUM = trk_num;
  } else {
    for (itrk=0; itrk<trk_num; itrk++) {
      sprintf((trk_pos+itrk)->IDC, "%s", trk_name[itrk]);
      trk_pos[itrk].UFL      = true;
      trk_pos[itrk].priority = trk_priority[itrk];
    }
    *TRK_NUM = 0;
  }

/*
---------------
*/

  if (ERROR_FLAG[TDSECZ] == true) {
    if (       array->TYPE == __CONNECTED_) {
      for (iant=0; iant<GRT_NUM; iant++) {
        dz[iant].trp /= speed_of_light;
      }
    } else if (array->TYPE == _VLBI_ARRAY_) {
      for (iant=0; iant<GRT_NUM; iant++) {
        dz[iant].trp *= (tdscz * 1.0e-3 / speed_of_light);
      }
    }
  }

/*
---------------
*/

  if (ERROR_FLAG[IDSECZ] == true) {
    for (iant=0; iant<GRT_NUM; iant++) {
      dz[iant].tec *= (idscz * 1.0e+16);
    }
  }

/*
---------------
*/

  wave_select(wave_id[0], &wave_length[0], &nu[0]);
  wave_select(wave_id[1], &wave_length[1], &nu[1]);
  wave_id[2]     = wave_id[0];
  nu[2]          = nu[0];
  wave_length[2] = wave_length[0];

/*
------
*/

  *band_width = band_width_set(bw_id);
  *cfact      = cfact_set(cfact_id);

/*
------
*/

  for (iant=0; iant<GRT_NUM; iant++) {
    if (*TRP_CONDITION == 0 || *TRP_CONDITION == 1 || *TRP_CONDITION == 2) {
      Cw[iant] *= wave_length[0];
    }
  }

/*
------
*/

  number_char_cut(ch_tdscz);
  number_char_cut(ch_idscz);
  number_char_cut(ch_orbae);
  number_char_cut(ch_orbpe);
  number_char_cut(ch_swtim);
  number_char_cut(ch_tflux);
  number_char_cut(ch_rflux);
  number_char_cut(ch_appon);

  if ((fp = fopen("aris_input/sim.prm", "w")) != NULL) {
    fprintf(fp, "BASELINE SELECT     %d\n", BL_SELECT);
    fprintf(fp, "TROPOS CONDITION    %d\n", *TRP_CONDITION);
    fprintf(fp, "IONOS CONDITION     %d\n", *ION_CONDITION);
    fprintf(fp, "TROPOS ZENITH ERROR %s\n", ch_tdscz);
    fprintf(fp, "IONOS ZENITH ERROR  %s\n", ch_idscz);
    fprintf(fp, "ORBIT ERROR APOGEE  %s\n", ch_orbae);
    fprintf(fp, "ORBIT ERROR PERIGEE %s\n", ch_orbpe);
    for (iant=0; iant<SRT_NUM; iant++) {
      fprintf(fp, "X_AXIS_SUN_ANGLE    %d,%d\n", iant+1, srt[iant].BODY_X_SUN);
    }
    for (itrk=0; itrk<trk_num; itrk++) {
      fprintf(fp, "TRACKING NETWORK    %s,%d\n",
              trk_pos[itrk].IDC, trk_pos[itrk].priority);
    }
    fprintf(fp, "START TIME (UTC)    %s\n", ch_timer[0]);
    fprintf(fp, "STOP TIME  (UTC)    %s\n", ch_timer[1]);
    fprintf(fp, "TGT WAVE ID         %d\n", wave_id[0]);
    fprintf(fp, "REF WAVE ID         %d\n", wave_id[1]);
    fprintf(fp, "BAND WIDTH          %d\n", bw_id);
    fprintf(fp, "FREQUENCY CHANNEL   %d\n", frchan_id);
    fprintf(fp, "COHERENCE FACTOR    %d\n", cfact_id);
    fprintf(fp, "SWITCHING CYC TIME  %s\n", ch_swtim);
    fprintf(fp, "APPARENT TARGET ON  %s\n", ch_appon);
    fprintf(fp, "TGT FLUX DENSITY    %s\n", ch_tflux);
    fprintf(fp, "TGT POSITION        %d\n", src[0].positionID);
    fprintf(fp, "TGT MORPHOLOGY      %d\n", src[0].morphology);
    fprintf(fp, "TGT MULTI COMPONENT %s\n", ch_file[0].mcm);
    fprintf(fp, "TGT CC TABLE        %s\n", ch_file[0].cct);
    fprintf(fp, "TGT BHS MODEL       %s\n", ch_file[0].bhs);
    fprintf(fp, "REF FLUX DENSITY    %s\n", ch_rflux);
    fprintf(fp, "REF POSITION        %d\n", src[1].positionID);
    fprintf(fp, "REF MORPHOLOGY      %d\n", src[1].morphology);
    fprintf(fp, "REF MULTI COMPONENT %s\n", ch_file[1].mcm);
    fprintf(fp, "REF CC TABLE        %s\n", ch_file[1].cct);
    fprintf(fp, "REF BHS MODEL       %s\n", ch_file[1].bhs);
    fprintf(fp, "REF PHASE PROCESS   %d\n", *reference_phase_process_mode);
    fclose (fp);
  } else {
    printf("CAUTION: Sim parameters cannot be saved. ");
    printf("Make directory \"./aris_input/\".\n");
  }

  return (1);
}



double  band_width_set(int  bw_id)
{
  double band_width = 0.0;
  if (bw_id == 0) {
    band_width =   1.953125e3;
  } else if (bw_id == 1) {
    band_width =    4.0e6;
  } else if (bw_id == 2) {
    band_width =    8.0e6;
  } else if (bw_id == 3) {
    band_width =   16.0e6;
  } else if (bw_id == 4) {
    band_width =   32.0e6;
  } else if (bw_id == 5) {
    band_width =  128.0e6;
  } else if (bw_id == 6) {
    band_width =  256.0e6;
  } else if (bw_id == 7) {
    band_width =  512.0e6;
  } else if (bw_id == 8) {
    band_width = 1024.0e6;
  } else if (bw_id == 9) {
    band_width = 2048.0e6;
  }
  return band_width;
}


double  cfact_set(int  cfact_id)
{
  double cfact = 0.0;
  if (cfact_id == 0) {         /* Two-level sampling  + three-level FR */
    cfact = 0.637 * 0.960;
  } else if (cfact_id == 1) {  /* four-level sampling + three-level FR */
    cfact = 0.881 * 0.960;
  } else if (cfact_id == 2) {  /* eight-level sampling + three-level FR */
    cfact = 0.960 * 0.960;
  }
  return cfact;
}


int      channel_num_set(int  frchan_id)
{
  if (frchan_id == 0) {
    return (1);
  } else if (frchan_id == 1) {
    return (4);
  } else if (frchan_id == 2) {
    return (8);
  } else if (frchan_id == 3) {
    return (16);
  } else if (frchan_id == 4) {
    return (32);
  }
}
