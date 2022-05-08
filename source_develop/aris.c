#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <unistd.h>
#include <aris.h>

/****
#define __DEMO__
#define __ANT_DEBUG__
#define __UV_DEBUG__
#define __CF_DEBUG__
#define __RANDOM_SEED__
****/

#define __RANDOM_SEED__

#define __CMG_FULL_SPEC__
#ifndef __CMG_FULL_SPEC__
  #define __CMG_HALF_SPEC__
#endif



/****
#define __ALAM_ANTPOS_RANDOM__
****/
#ifndef __ALMA_ANTPOS_RANDOM__
  #define __ALMA_ANTPOS_SLOPE__
#else
/****
  #define __ALMA_ANTPOS_NOW__
****/
  #define __ALMA_ANTPOS_NOW__
  #ifndef __ALMA_ANTPOS_NOW__
    #define __ALMA_ANTPOS_FTR__
  #endif
#endif


 
/****
  src[0]: target
  src[1]: refernce calibrator
****/

int  main(int argc, char **argv)
{
  FILE   *fp;
  int    i, j, k, I, J, K, L;
  int    ns, SRCID[SRC_NUM];
  int    itmp, proc_mode;
  _Bool  ana__mode;
  int    wave_id[SRC_NUM_P1];
  int    ARRAY_TYPE, ARRAY_ID;
  int    iant, jant, ANT_NUM, GRT_NUM;
  int    ndata, nant[SRC_NUM];
  int    BGN_ANT_I, BGN_ANT_J, END_ANT_I, END_ANT_J;
  int    isite;
  int    NLOOP;
  int    IX[ANTMAX], IY[ANTMAX], TX, TY;

  _Bool  ERROR_FLAG[ERROR_NUM];
  int    iobs, nobs;
  int    iave, nswt=0, nslewt=0, noffset;
  double apparent_tgt_on;
  int    ON_TIME[SRC_NUM];
  int    ibase, nbase;
  int    iaxis;
  char   string[200];
  float  *mapr, *mapi;
  float  src_flag[SRC_NUM];
  double pix_uvl[SRC_NUM], pix_mas[SRC_NUM];
  double wave_length[SRC_NUM_P1], nu[SRC_NUM_P1];
  double uv_max_tmp[SRC_NUM], uv_max, uv_factor[SRC_NUM];
  double bl_max;
  double SEPANG;
  double azl=0.0, ell=0.0;
  double amperr;
  double lftmp;
  double pslope;
  double DS[SRC_NUM];
  double sefd[SRC_NUM];
  double soffx, soffy;
  int    nseries, NELEM, SITE_NUM=1; 
  double *D_S;

  char   comment[40][NCOMLEN];
  struct comment_param cmnt;
  char   c_dummy = 0;

  int    TRP_CONDITION, ION_CONDITION;
  struct phase_screen_parameter wvc[ANTMAX], ion[ANTMAX], dry[ANTMAX];
  struct phase_screen_parameter ALMA_ant;
  double Cw[ANTMAX], Cd[ANTMAX], CW, CD, CI;
  double Ci[2][86400];
  double alpha1, alpha2, wind_v;
  int    corner_position[2][2];

  double amp, SIGMA=0.0, sigma=0.0, f1, f2;
  double dSP1, dSP2;
  double Sv, Sh;
  struct EOP_data EOP, EOP_e;
  double dW, dWX, dWY, dUT1, dPN;
  double W[3][3], R1[3][3];
  double xyz_tmp[3];
  double d_Delta_psi, d_Delta_eps;
  double pnt_err_tau = 0.0;

  int    TimUTC[6], timUTC[6];
  int    DOY;
  double UT1_UTC,  ut1_utc;
  double obs_duration;
  double OBS_T, OBS_p;
  double grt_elevation_limit;
  struct srt_data_link  *srt_link=NULL;
  struct atmospheric_zenith_error *dz[2];
  struct TID  MSTID;

  struct fringe *frng[SRC_NUM_P1];
  float  *fringe_weight[SRC_NUM_P1];

  struct st_observable     *int_obs[SRC_NUM];
  double *wvc_ds[SRC_NUM], *ion_ds[SRC_NUM];
  double *dry_ds[SRC_NUM];
  double *fqs_ds=NULL;

  int    COUNT_NOD;
  double *seed_dist=NULL;
  double seed_end[2];

  struct baseline_uvw  *bluvw[SRC_NUM];

  struct source_parameter src[SRC_NUM], src_proc[SRC_NUM], sun;
  struct morphology_file_name ch_file[SRC_NUM];

  int    IFREF=FLDMAX/2, JFREF=FLDMAX/2;
  int    ifrq, nfrq;
  double ar, ai, fr, fi, F, tau_err, ion_err, phs_err, fai_err;
  double d_lo_phs;
  double band_width, d_band_width, cfact;
  double inttim=1.0;
  struct antenna_parameter ant_prm[ANTMAX];
  struct antenna_error_parameter ant_err[ANTMAX];
  char   antenna_code[ANTMAX][10];
  char   antenna_list_file[500];
  struct array_parameter array;
  struct data_number data_num;
  int    reference_phase_process_mode;

  int    SRT_NUM;
  double sep_angle_limit_from_earth_limb;
  struct srt_orbit_parameter srt[SRTMAX];
  double init_l[SRTMAX];

  int    TRK_NUM;
  struct antenna_parameter  trk_pos[TRKMAX];

  _Bool  AZEL_FIX = false;
  char   fits_fname[100], fname[100];

  _Bool  TV_SWT = false;
  int    pgid[10];
  float  cursor_pos[2];
  float  bttn_box[2][4];
  float  pitch = 0.03;

  double fov_st;

  double frq_coh_factor[SRC_NUM];

/**** ACA DEMO PARAMETER ****/

  double pwv_exp, pwv_rms;

/*
===================================================================
*/

  for (i=0; i<10; i++) {
    pgid[i] = -1;
  }

  for (i=0; i<SRC_NUM; i++) {
    sefd[i] = 0.0;
  }

/*
===================================================================
*/

  for (i=0; i<ERROR_NUM; i++) {
    ERROR_FLAG[i]      = false;
  }

  for (ns=0; ns<SRC_NUM; ns++) {
    src_flag[ns] = 100.0 * (float)(ns + 1);
    frq_coh_factor[ns] = 1.0;
  }

/*
===================================================================
*/

  for (iant=0; iant<SRTMAX; iant++) {
    srt[iant].BODY_X_SUN = 0;
  }

/*
===================================================================
*/

  OBS_T               =   290.0;
  OBS_p               =  1013.0;
  grt_elevation_limit =     7.0  / 180.0 * dpi;
  sep_angle_limit_from_earth_limb
                      =     0.25 / 180.0 * dpi;

/*
===================================================================
*/

#ifdef __RANDOM_SEED__
  seed_random(true);
#else
  seed_random(false);
#endif /* __RANDOM_SEED__ */

/*
===================================================================
*/

  if (argc >= 2 && strncmp(argv[1], "-g", 2) == 0) {
    TV_SWT = true;
  } else {
    TV_SWT = false;
  }

/*
===================================================================
*/

/***********************************************
  cpgslct(cpgid2);
  cpgbbuf();
  cpgsave();
  Gaussian_noise_check();
  cpgebuf();
  cpgunsa();
***********************************************/

/*
===================================================================
*/

  if (TV_SWT == true) {
    if (pgid[0] < 0) {
      pgid[0] = (int)cpgopen("/xs");
      if (pgid[0] < 0) {
        cpgask(-1);
        TV_SWT = false;
      }
    }
    cmnt.ncol  = 12;
    cmnt.xmin  = 0.020;
    cmnt.xmax  = 0.980;
    cmnt.ymin  = 0.005;
    cmnt.ymax  = 0.180;
    cmnt.pitch = 0.030;
  }

/*
----------------
*/

  for (iant=0; iant<ANTMAX; iant++) {
    for (i=0; i<10; i++) {
      ant_prm[iant].IDC[i]  = 0;
      ant_err[iant].IDC[i]  = 0;
      antenna_code[iant][i] = 0;
    }
  }

  NLOOP = 1;
  while (NLOOP) {
    if (TV_SWT == true) {
      cpgslct(pgid[0]);
    }
    SRT_NUM  = 0;
    array.ID = ALL_ANT;
    for (ns=0; ns<SRC_NUM; ns++) {
      wave_id[ns] = -1;
    }
    if (TV_SWT == true) {
      cpgscr(0, 0.8, 0.8, 0.8);
      cpgscr(1, 0.0, 0.0, 0.0);

      cpgbbuf();

      cpgslct(pgid[0]);
      cpgpap(2.00*pgpap_prm, 1.0/1.4);
      cpgsch(1.5*pgpap_prm/13.0);
      cpgsvp(0.0,  1.0, 0.0, 1.0);
      cpgswin(0.0, 1.4, 0.0, 1.0);
      comment_init(&cmnt, comment, true);
  
      sprintf(string,
       "ARIS: Observation Parameter Setup --GRAPHICAL USER INTERFACE MODE--");
      comment_disp(&cmnt, comment, string, true);
      sprintf(string, "Welcome!! Please set the observing parameters.");
      comment_disp(&cmnt, comment, string, true);
      if (NLOOP > 1) {
        sprintf(string, "Antenna list file reloaded.");
        comment_disp(&cmnt, comment, string, true);
      }

    } else if (TV_SWT == false) {
      printf("ARIS: Observation Parameter Setup --TERMINAL USER INTERFACE MODE--\n");
      printf("Welcome!! Please set the observing parameters.\n");
      if (NLOOP > 1) {
        printf("Antenna list file reloaded.\n");
      }
    }

    proc_mode = obs_param_input(ERROR_FLAG, &array, wave_id,
                    &ANT_NUM, &GRT_NUM, &SRT_NUM, srt, &grt_elevation_limit,
                    sep_angle_limit_from_earth_limb,
                    TimUTC, &UT1_UTC, &obs_duration, src, &sun,
                    antenna_list_file, ant_prm, antenna_code,
                    &cmnt, comment,
                    TV_SWT, cursor_pos, pgid);

    if (proc_mode == __NG__) {
      sprintf(string, "ERROR: ARIS: in OBS_PARAM_INPUT: exit.");
      printf("%s\n", string);
      if (TV_SWT == true) {
        comment_disp(&cmnt, comment, string, true);
      }
      return (0);
    } else if (proc_mode == _EXIT_) {
      sprintf(string, "See you again!!");
      if (TV_SWT == true) {
        comment_disp(&cmnt, comment, string, true);
      }
      printf("EXIT.\n");
      return (0);
    } else if (proc_mode == ANT_VIS) {
      proc_mode = antenna_visibility(
                    ANT_NUM, GRT_NUM, SRT_NUM, ant_prm, antenna_code, src, sun,
                    wvc, dry, ion, Cw, Cd, Ci, fqs_ds, src_flag, srt,
                    trk_pos, sep_angle_limit_from_earth_limb,
                    &cmnt, comment, TV_SWT, pgid);
      if (proc_mode == -1 || proc_mode == 0) {
        return (0);
      }
    } else if (proc_mode == PHS_SCR) {
      phase_screen_check(NOISE_MTRX, wvc[0], cursor_pos, TV_SWT, pgid);
    } else if (proc_mode == __GO__) {
      antenna_selection(&ANT_NUM, &GRT_NUM, &SRT_NUM, wave_id,
                        grt_elevation_limit, ant_prm, antenna_code,
                        antenna_list_file, true);
#ifdef __ANT_DEBUG__
      for (iant=0; iant<GRT_NUM; iant++) {
        printf("__ANT_DEBUG__   %d   %d    %s    %s  %lf  %lf  %lf\n",
              iant+1, strlen(antenna_code[iant]),
              antenna_code[iant],
              ant_prm[iant].IDC,
              ant_prm[iant].LLH[0],
              ant_prm[iant].LLH[1],
              ant_prm[iant].LLH[2]);
      }
#endif
      break;
    }
    NLOOP++;
  }

  bl_max = 0;
  for (iant=0; iant<GRT_NUM; iant++) {
    for (jant=iant+1; jant<GRT_NUM; jant++) {
      lftmp = pow(ant_prm[iant].XYZ[0] - ant_prm[jant].XYZ[0], 2.0)
            + pow(ant_prm[iant].XYZ[1] - ant_prm[jant].XYZ[1], 2.0)
            + pow(ant_prm[iant].XYZ[2] - ant_prm[jant].XYZ[2], 2.0);
      if (lftmp > bl_max) {
        bl_max = lftmp;
      }
    } 
  }
  bl_max = sqrt(bl_max);

/*
===================================================================
*/

  nobs = (int)lrint(obs_duration);
  data_num.sobs =    0;
  data_num.eobs = nobs;
  data_num.nobs = nobs;

/*
------------------------------
*/

  if (SRT_NUM >= 1) {
    if ((srt_link = (struct srt_data_link *)
                        calloc(1, sizeof(struct srt_data_link))) == NULL) {
      printf("ARIS: ERROR: memory allocation: srt_link.\n");
      return (0);
    }
    if (get_srt_link(srt_link, &fov_st) == __NG__) {
      printf("ARIS: ERROR: GET_SRT_LINK\n");
      return (0);
    }
  }

  for (ns=0; ns<SRC_NUM; ns++) {
    uv_factor[ns] = 4.0;
  }
  /*
     For wider imaging area, uv_factor would be a small value. 
     (for example, 1.0.) For small imaging area, it would be 
     3.0 - 4.0. 
  */

  if (       array.TYPE == __CONNECTED_) {
    SITE_NUM = 1;
    NELEM    = GRT_NUM;
    nseries  = SRC_NUM * GRT_NUM;
  } else if (array.TYPE == _VLBI_ARRAY_) {
    SITE_NUM = GRT_NUM;
    NELEM    = 1;
    nseries  = 2;
  }

/*
==============================================================
*/

  if (       array.TYPE == __CONNECTED_) {
    alpha1  = 2.0 * 0.60;   /* Matsushita et al. (2017), Holdaway (2004) */
    alpha2  = 2.0 * 0.20;   /* Matsushita et al. (2017) */
    wind_v  =  -6.0;
  } else if (array.TYPE == _VLBI_ARRAY_) {
    alpha1  = 5.0 / 3.0;    /* Kolmogorov turbulence */
    alpha2  = 2.0 / 3.0;    /* Kolmogorov turbulence */
    wind_v  = -10.0;
  }
  pwv_exp = 1.5;   /* [mm] */
  pwv_rms = pwv_exp * 0.02 / pow(100.0, 0.5*alpha1);

/*
===================================================================
*/

  if (ERROR_FLAG[TWVTRB] == true) {
    for (isite=0; isite<SITE_NUM; isite++) {
      wvc[isite].pixel      = 1.0;
      if (array.TYPE == __CONNECTED_) {
        if (bl_max > 800.0) {
          wvc[isite].pixel      = (bl_max / 800.0) * 2.0;
          wvc[isite].pixel      = bl_max / 8000.0;
        }
        if (wvc[isite].pixel < 1.0) {
          wvc[isite].pixel = 1.0;
        }
      }
      wvc[isite].H_d        = 1.0e3;
      wvc[isite].H_s        = 1.0e3;
      wvc[isite].i_scale[0] = 1024.0;
      fit_power_two_number(wvc[isite].i_scale[0],  wvc[isite].pixel,
                         &(wvc[isite].i_scale[0]), &itmp);
      wvc[isite].o_scale[0] = 8192.0;
      fit_power_two_number(wvc[isite].o_scale[0],  wvc[isite].pixel,
                         &(wvc[isite].o_scale[0]), &itmp);
      wvc[isite].i_expon    = alpha1;
      wvc[isite].o_expon    = alpha2;
      wvc[isite].v[0]       = wind_v;
      wvc[isite].v[1]       = 0.0;
      wvc[isite].i_coeffi   = 1.0;
      wvc[isite].o_coeffi   = 0.0;
      wvc[isite].c_coeffi   = 0.0;
      Cw[isite] = wvc[isite].i_coeffi;
    }
  }

/*
==============================================================
*/

  if (ERROR_FLAG[DRYTRB] == true && array.TYPE == __CONNECTED_) {
    for (isite=0; isite<SITE_NUM; isite++) {
      dry[isite].H_d        = 20.0e+3;
      dry[isite].H_s        = 20.0e+3;
      dry[isite].pixel      = 20.0;
      dry[isite].i_scale[0] = dry[isite].pixel *  1024.0;
      dry[isite].o_scale[0] = dry[isite].pixel *  8192.0;
      dry[isite].i_expon    = alpha1;
      dry[isite].o_expon    = 1.0e-10;
      dry[isite].i_coeffi   = 1.0;
      dry[isite].o_coeffi   = 0.0;
      dry[isite].c_coeffi   = 0.0;
      dry[isite].v[0]       = 0.0;
      if (ant_prm[isite].XYZ[2] >= 0.0) {
        dry[isite].v[1]     = -20.0;
      } else {
        dry[isite].v[1]     =  20.0;
      }
      Cd[isite] = dry[isite].i_coeffi;
      CD        = 6.0e-2 / 3.0e-3;
/**** 
6.0e-2 is an appropriate scale to fit to the ALMA LB test resutls. 
3.0e-3 (= 3 mm) is the wave lenght scaling factor with respect to 
Band3. Later, the observing wavelength is multiplied.
****/
    }
  }

/*
==============================================================
*/

  if (ERROR_FLAG[IONTRB] == true) {
    for (isite=0; isite<SITE_NUM; isite++) {
      ion[isite].H_d        = 300.0e+3;
      ion[isite].H_s        = 450.0e+3;
      ion[isite].pixel      = 150.0;
      ion[isite].pixel      = 300.0;
      ion[isite].i_scale[0] = ion[isite].pixel *   512.0;
      ion[isite].o_scale[0] = ion[isite].pixel *  8192.0;
      ion[isite].i_expon    = 5.0 / 3.0;
      ion[isite].o_expon    = -40.0;
      ion[isite].i_coeffi   = 1.0e12;
      ion[isite].o_coeffi   = 0.0;
      ion[isite].c_coeffi   = 0.0;
      ion[isite].v[0]       = 0.0;
      if (ant_prm[isite].XYZ[2] >= 0.0) {
        ion[isite].v[1]       = -100.0;
      } else {
        ion[isite].v[1]       =  100.0;
      }
    }

    if (TEC_fluctuation_amp(Ci, TimUTC) == __NG__) {
      printf("ERROR: ARIS: TEC_fluctuation_amp.\n");
      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      return (0);
    }

    MSTID.amp    = 1.0e16;
    MSTID.lambda = 200.0e3;
    MSTID.v[0]   =   0.0;
    MSTID.v[1]   = 100.0;
  }

/*
=====================================================
*/

  Sv   = 0.0;    /* Antenna Position: Vertical Comp.    */
  Sh   = 0.0;    /* Antenna Position: Horizontal Comp.  */
  dSP1 = 0.0;
  dSP2 = 0.0;
  dUT1 = 0.0;
  dW   = 0.0;
  dPN  = 0.0;

  if ((fp=fopen("aris_input/fixed_parameter.prm", "r")) == NULL) {
    printf("ERROR: MAIN PROGRAM: ./aris_input/fixed_parameter.prm.\n");
    if (SRT_NUM >= 1) {
      free (srt_link);
    }
    return -1;
  }
  while (1) {
    if (fgets(string, sizeof(string), fp) == NULL) {
      break;
    }
    if (string[0] != '#') {
      if (strncmp(string, "SOURCE-1 POSITION ERROR", 23) == 0) {
        i = 0;
        while (1) {
          if (string[i++] == ':') {
            break;
          }
        }
        sscanf(string+i, "%lf", &dSP1);
      } else if (strncmp(string, "SOURCE-2 POSITION ERROR", 23) == 0) {
        i = 0;
        while (1) {
          if (string[i++] == ':') {
            break;
          }
        }
        sscanf(string+i, "%lf", &dSP2);
      } else if (strncmp(string, "ANTENNA POSITION ERROR-H", 24) == 0) {
        i = 0;
        while (1) {
          if (string[i++] == ':') {
            break;
          }
        }
        sscanf(string+i, "%lf", &Sh);
      } else if (strncmp(string, "ANTENNA POSITION ERROR-V", 24) == 0) {
        i = 0;
        while (1) {
          if (string[i++] == ':') {
            break;
          }
        }
        sscanf(string+i, "%lf", &Sv);
      } else if (strncmp(string, "UT1 ERROR", 9) == 0) {
        i = 0;
        while (1) {
          if (string[i++] == ':') {
            break;
          }
        }
        sscanf(string+i, "%lf", &dUT1);
      } else if (strncmp(string, "POLAR MOTION ERROR", 18) == 0) {
        i = 0;
        while (1) {
          if (string[i++] == ':') {
            break;
          }
        }
        sscanf(string+i, "%lf", &dW);
      } else if (strncmp(string, "PRECESSION-NUTATION ERROR", 25) == 0) {
        i = 0;
        while (1) {
          if (string[i++] == ':') {
            break;
          }
        }
        sscanf(string+i, "%lf", &dPN);
      }
    }
  }
  fclose (fp);

/*
-------------------------------
*/

  EOP.WX   = 0.0;
  EOP.WY   = 0.0;
  EOP_e.WX = 0.0;
  EOP_e.WY = 0.0;

  ALMA_ant.pixel      = 4.0;
  ALMA_ant.H_d        = 1.0e3;
  ALMA_ant.H_s        = 1.0e3;
  ALMA_ant.i_scale[0] = 1024.0;
  ALMA_ant.o_scale[0] = 16384.0;
  ALMA_ant.i_expon    = alpha1;
  ALMA_ant.o_expon    = alpha2;
  ALMA_ant.v[0]       =   0.0;
  ALMA_ant.v[1]       =   0.0;
  ALMA_ant.i_coeffi   = 0.0;
  ALMA_ant.o_coeffi   = 0.0;
  ALMA_ant.c_coeffi   = 1.0;

/*
---- Antenna Position Error ----
*/

  if (ERROR_FLAG[APOSER] == true) {
    if (array.TYPE == __CONNECTED_) {

/****
#### __ALMA_ANTENNA_POSITION_ERROR__ ####

According to Hunter et al. (2016, Proc. SPIE 9906, 99062M), the ALMA antenna 
position has a antenna-distance dependency with respect to the array center 
position. The following code is an emprical one based on Hunter et al. (2016).
((x, y, z) = (0.071, 0.054, 0.198) mm/km)
****/

      lftmp = 0.5 * dpi;
      TX = 0;
      TY = 0;
      for (iant=0; iant<GRT_NUM; iant++) {
        position_on_screen(0.0, ant_prm[iant].OFS, ALMA_ant, 0.0,
                           lftmp, &soffx, &soffy);
        IX[iant] = soffx / ALMA_ant.pixel;
        IY[iant] = soffy / ALMA_ant.pixel;
        if (IX[iant] < TX) {
          TX = IX[iant];
        }
        if (IY[iant] < TY) {
          TY = IY[iant];
        }
      }
      for (iant=0; iant<GRT_NUM; iant++) {
        IX[iant] -= TX;
        IY[iant] -= TY;
      }

      corner_position[0][0] = 0.0;
      corner_position[0][1] = 0.0;
      corner_position[1][0] = 0.0;
      corner_position[0][1] = 0.0;
      D_S       = (double *)calloc(NOISE_MTRX+1, sizeof(double));
      seed_dist = (double *)calloc(NOISE_MTRX+1, sizeof(double));

      for (iaxis=0; iaxis<3; iaxis++) {
        if (iaxis == 0) {
          sprintf(string, "ALMA antenna position error (E-W) calculating.");
        } else if (iaxis == 1) {
          sprintf(string, "ALMA antenna position error (N-S) calculating.");
        } else if (iaxis == 2) {
          sprintf(string, "ALMA antenna position error (Zenith) calculating.");
        }
        printf("%s\n", string);
        if (TV_SWT == true) {
          comment_disp(&cmnt, comment, string, true);
        }


#ifdef __ALAM_ANTPOS_RANDOM__
  #ifdef __ALMA_ANTPOS_NOW__
        if (iaxis == 0 || iaxis == 1) {
          ALMA_ant.c_coeffi   = 1.0e-3;
        } else if (iaxis == 2) {
          ALMA_ant.c_coeffi   = 2.0e-3;
        }
  #elif defined __ALMA_ANTPOS_FTR__
        if (iaxis == 0 || iaxis == 1) {
          ALMA_ant.c_coeffi   = 0.7e-3;
        } else if (iaxis == 2) {
          ALMA_ant.c_coeffi   = 1.0e-3;
        }
  #endif
        if ((COUNT_NOD = turbulent_phase_screen
                   (NOISE_MTRX, 0, seed_dist, &ALMA_ant,
                   false, GRT_NUM, IX, IY, D_S, corner_position, 0)) == -1) {
            printf("ERROR: ARIS: ");
            printf("ALMA antenna position errors could not be calculated ");
            printf("due to errors in ATMOSPHERIC_FLUCTUATION.\n");
            exit (-1);
        }
        for (iant=0; iant<GRT_NUM; iant++) {
          ant_err[iant].ERR[iaxis] = D_S[iant];
        }
#elif defined __ALMA_ANTPOS_SLOPE__
        if (iaxis == 0) {
          pslope = 0.2e-3 / 8.0e3;
        } else if (iaxis == 1) {
          pslope = 0.2e-3 / 8.0e3;
        } else if (iaxis == 2) {
          pslope = 4.0e-3 / 8.0e3;
        }
        for (iant=0; iant<GRT_NUM; iant++) {
          ant_err[iant].ERR[iaxis] = pslope * ant_prm[iant].OFS[0];
        }
#endif


      }
      free (D_S);

      for (iant=0; iant<GRT_NUM; iant++) {
        xyz_tmp[0] = ant_err[iant].ERR[0];
        xyz_tmp[1] = ant_err[iant].ERR[1];
        xyz_tmp[2] = ant_err[iant].ERR[2];
        drotate(xyz_tmp, -ant_prm[iant].LLH[1] + 0.5*dpi, "x");
        drotate(xyz_tmp,  ant_prm[iant].LLH[0] + 0.5*dpi, "z");
        ant_err[iant].ERR[0] = xyz_tmp[0];
        ant_err[iant].ERR[1] = xyz_tmp[1];
        ant_err[iant].ERR[2] = xyz_tmp[2];
/****
        printf("%11.6e, %11.6e, %11.6e, \n",
                xyz_tmp[0], xyz_tmp[1], xyz_tmp[2]);
****/
      }

/****
#### __ALMA_ANTENNA_POSITION_ERROR__ ####
****/

/*
---- VLBI ----
*/

    } else if (array.TYPE == _VLBI_ARRAY_) {
      for (iant=0; iant<GRT_NUM; iant++) {
        if (ant_prm[iant].AAE[0] == 0.0 &
            ant_prm[iant].AAE[1] == 0.0 &
            ant_prm[iant].AAE[2] == 0.0) {

/****
#### __VLBI_ANTENNA_POSITION_ERROR__ ####

The VLBI antenna position errors are set randomly with accordance 
with Asaki et al. (2007).
****/

          f1 = 0.5 * dpi * random_val1();
          f2 =       dpi * random_val1();
          amp = gauss_dev() / sqrt(pow(cos(f1)/Sh, 2.0) + pow(sin(f1)/Sv , 2.0));
          xyz_tmp[0] = amp * cos(f1) * cos(f2);
          xyz_tmp[1] = amp * cos(f1) * sin(f2);
          xyz_tmp[2] = amp * sin(f1);
          drotate(xyz_tmp, -dpi/2.0, "z");
          drotate(xyz_tmp,
            dpi/2.0-atan2(ant_prm[iant].XYZ[2], vlen2(ant_prm[iant].XYZ)), "y");
          drotate(xyz_tmp,
                  atan2(ant_prm[iant].XYZ[1], ant_prm[iant].XYZ[0]), "z");
          ant_err[iant].ERR[0] = xyz_tmp[0];
          ant_err[iant].ERR[1] = xyz_tmp[1];
          ant_err[iant].ERR[2] = xyz_tmp[2];
        } else {
          ant_err[iant].ERR[0] = ant_prm[iant].AAE[0];
          ant_err[iant].ERR[1] = ant_prm[iant].AAE[1];
          ant_err[iant].ERR[2] = ant_prm[iant].AAE[2];
        }
      }
    }
  } else if (ERROR_FLAG[APOSER] == false) {
    for (iant=0; iant<GRT_NUM; iant++) {
      ant_err[iant].ERR[0] = 0.0;
      ant_err[iant].ERR[1] = 0.0;
      ant_err[iant].ERR[2] = 0.0;
    }
  }

/*
---- Target Position Error ----
*/

  if (ERROR_FLAG[TPOSER] == true) {
    DS[0] = dSP1 * dpi / 180.0 / 3600.0 * gauss_dev();
    ref_pos_shift(src,   src+1, DS[0], &cmnt, comment, TV_SWT);
  }

/*
---- Reference Calibrator Position Error ----
*/

  if (ERROR_FLAG[RPOSER] == true) {
    DS[1] = dSP2 * dpi / 180.0 / 3600.0 * gauss_dev();
    ref_pos_shift(src+1, src,   DS[1], &cmnt, comment, TV_SWT);
  }

/*
---- EOP ----
*/

  transformation_matrices(TimUTC, (&EOP)->Qt,   (&EOP)->W,
                          0.0,     0.0,   EOP.WX,  EOP.WY);
  transformation_matrices(TimUTC, (&EOP_e)->Qt, (&EOP_e)->W,
                          0.0,     0.0,   EOP.WX,  EOP.WY);
  for (iant=0; iant<GRT_NUM; iant++) {
    dvector_calc(EOP.W, ant_prm[iant].XYZ);
  }

/*
---- EOP Error ----
*/

  if (ERROR_FLAG[EOPERR] == true) {
    dUT1        *=    gauss_dev();
                      /* UT1 Error                           : 0.02 ms */
    dWX         =     dW  / sqrt(2.0) / 3600.0 / 180.0 * dpi * gauss_dev();
                      /* Polar Motion Error [x]              : 0.3 mas */
    dWY         =     dW  / sqrt(2.0) / 3600.0 / 180.0 * dpi * gauss_dev();
                      /* Polar Motion Error [y]              : 0.3 mas */
    d_Delta_psi =     dPN / sqrt(2.0) / 3600.0 / 180.0 * dpi * gauss_dev();
    d_Delta_eps =     dPN / sqrt(2.0) / 3600.0 / 180.0 * dpi * gauss_dev();
                      /* Precession and Nutation             : 0.3 mas */

    transformation_matrices(TimUTC, (&EOP_e)->Qt, (&EOP_e)->W,
                            d_Delta_psi,  d_Delta_eps,
                            EOP.WX+dWX,  EOP.WY+dWY);

    for (iant=0; iant<GRT_NUM; iant++) {
      dvector_calc(EOP_e.W, ant_prm[iant].ERR);
    }
  } else {
    for (iant=0; iant<GRT_NUM; iant++) {
      dvector_calc(EOP.W,   ant_prm[iant].ERR);
    }
  }

/*
---------------
*/

  for (ns=0; ns<SRC_NUM; ns++) {
    dvector_calc(EOP.Qt, src[ns].s);
    src[ns].RA     = atan2(src[ns].s[1], src[ns].s[0]);
    src[ns].DC     = atan2(src[ns].s[2], vlen2(src[ns].s));
  }

  for (ns=0; ns<SRC_NUM; ns++) {
    if (ERROR_FLAG[EOPERR] == true) {
      dvector_calc(EOP_e.Qt, src[ns].s_e);
    } else {
      dvector_calc(EOP.Qt,   src[ns].s_e);
    }
    src[ns].RA_e   = atan2(src[ns].s_e[1], src[ns].s_e[0]);
    src[ns].DC_e   = atan2(src[ns].s_e[2], vlen2(src[ns].s_e));
  }

/*
---------------
*/

  if (ERROR_FLAG[EOPERR] == true || ERROR_FLAG[APOSER] == true ||
      ERROR_FLAG[RPOSER] == true) {
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        R1[j][i] = EOP.W[i][j];
      }
    }
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
        W[i][j] = 0.0;
        for (k=0; k<3; k++) {
          W[i][j] += EOP_e.W[i][k] * R1[k][j];
        }
      }
    }
  }

/*
---------------
*/

  for (i=0; i<2; i++) {
    if ((dz[i] = (struct atmospheric_zenith_error *)
      calloc(GRT_NUM, sizeof(struct atmospheric_zenith_error))) == NULL) {
      printf("ARIS: ERROR: memory allocation: dz.\n");
      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, i);
      return (0);
    }
  }

  if (ERROR_FLAG[TDSECZ] == true) {
    if (array.TYPE == __CONNECTED_) {

/****
#### __ALMA_ANTENNA_POSITION_ERROR__ ####

According to Hunter et al. (2016, Proc. SPIE 9906, 99062M), the ALMA antenna 
position has a antenna-distance dependency with respect to the array center 
position. The following code is an emprical one based on Hunter et al. (2016).
((x, y, z) = (0.071, 0.054, 0.198) mm/km)
****/

      lftmp = 0.5 * dpi;
      TX = 0;
      TY = 0;
      for (iant=0; iant<GRT_NUM; iant++) {
        position_on_screen(0.0, ant_prm[iant].OFS, ALMA_ant, 0.0,
                           lftmp, &soffx, &soffy);
        IX[iant] = soffx / ALMA_ant.pixel;
        IY[iant] = soffy / ALMA_ant.pixel;
        if (IX[iant] < TX) {
          TX = IX[iant];
        }
        if (IY[iant] < TY) {
          TY = IY[iant];
        }
      }
      for (iant=0; iant<GRT_NUM; iant++) {
        IX[iant] -= TX;
        IY[iant] -= TY;
      }

      corner_position[0][0] = 0.0;
      corner_position[0][1] = 0.0;
      corner_position[1][0] = 0.0;
      corner_position[0][1] = 0.0;
      D_S       = (double *)calloc(NOISE_MTRX+1, sizeof(double));
      seed_dist = (double *)calloc(NOISE_MTRX+1, sizeof(double));

      sprintf(string, "Atacama zenith path delay calculating.");
      printf("%s\n", string);
      if (TV_SWT == true) {
        comment_disp(&cmnt, comment, string, true);
      }

#ifdef __ALMA_ANTPOS_NOW__
      ALMA_ant.c_coeffi   = 2.0e-3;
#elif defined __ALMA_ANTPOS_FTR__
      ALMA_ant.c_coeffi   = 1.0e-3;
#endif
      ALMA_ant.c_coeffi   = 0.5e-3;
      if ((COUNT_NOD = turbulent_phase_screen
                 (NOISE_MTRX, 0, seed_dist, &ALMA_ant,
                 false, GRT_NUM, IX, IY, D_S, corner_position, 0)) == -1) {
        printf("ERROR: ARIS: ");
        printf("ALMA antenna position errors could not be calculated ");
        printf("due to errors in ATMOSPHERIC_FLUCTUATION.\n");
        exit (-1);
      }
      for (iant=0; iant<GRT_NUM; iant++) {
        dz[0][iant].trp = D_S[iant];
      }
      free (D_S);

/*
---- VLBI ----
*/

    } else if (array.TYPE == _VLBI_ARRAY_) {
      for (iant=0; iant<GRT_NUM; iant++) {
        dz[0][iant].trp = gauss_dev();
      }
    }
  }

  if (ERROR_FLAG[IDSECZ] == true) {
    if (array.TYPE == _VLBI_ARRAY_) {
      for (iant=0; iant<GRT_NUM; iant++) {
        dz[0][iant].tec = gauss_dev();
      }
    }
  }

/*
---- Local Oscilator Phase Error ----
*/

  if (ERROR_FLAG[LOPOFS] == false && ERROR_FLAG[LOPJMP] == false) {
    for (iant=0; iant<ANT_NUM; iant++) {
      for (i=0; i<N_WAVE; i++) {
        ant_err[iant].LOPHS[i] = 0.0;
      }
    }
  } else if (ERROR_FLAG[LOPOFS] == true) {
    for (iant=0; iant<ANT_NUM; iant++) {
      for (i=0; i<N_WAVE; i++) {
        ant_err[iant].LOPHS[i] = dpi * random_val1();
      }
    }
  }

/*
---- Antenna Gain Loss Error ----
*/

  if (ERROR_FLAG[AMPERR] == true) {
    for (iant=0; iant<ANT_NUM; iant++) {
      ant_err[iant].d_gain = 0.05 * gauss_dev();
    }
  } else if (ERROR_FLAG[AMPERR] == false) {
    for (iant=0; iant<ANT_NUM; iant++) {
      ant_err[iant].d_gain = 0.00;
    }
  }

/*
==============================================================
*/

  I = data_num.nobs * ANT_NUM;
  for (ns=0; ns<SRC_NUM; ns++) {
    if ((int_obs[ns] = (struct st_observable *)
                   calloc(I, sizeof(struct st_observable))) == NULL) {
      printf("ARIS: ERROR: memory allocation: obs[%d]\n", ns);

      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, ns);
      return (0);
    }
  }
  for (ns=0; ns<SRC_NUM; ns++) {
    for (i=0; i<I; i++) {
      int_obs[ns][i].wt        = src_flag[ns];
      int_obs[ns][i].amp_error = 1.0;
    }
  }

  if (ERROR_FLAG[TWVTRB] == true) {
    for (ns=0; ns<SRC_NUM; ns++) {
      if ((wvc_ds[ns] = (double *)calloc(I, sizeof(double))) == NULL) {
        printf("ARIS: ERROR: calloc for wvc_ds[%d].\n", ns);

        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        free_memory_block((void *)wvc_ds,  ns);
        return (0);
      }
    }
  }

  if (ERROR_FLAG[DRYTRB] == true) {
    for (ns=0; ns<SRC_NUM; ns++) {
      if ((dry_ds[ns] = (double *)calloc(I, sizeof(double))) == NULL) {
        printf("ARIS: ERROR: calloc for dry_ds[%d].\n", ns);

        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        free_memory_block((void *)wvc_ds,  ns);
        free_memory_block((void *)dry_ds,  ns);
        return (0);
      }
    }
  }

  if (ERROR_FLAG[IONTRB] == true) {
    for (ns=0; ns<SRC_NUM; ns++) {
      if ((ion_ds[ns] = (double *)calloc(I, sizeof(double))) == NULL) {;
        printf("AIRS: ERROR: calloc for ion_ds[%d].\n", ns);

        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        if (ERROR_FLAG[TWVTRB] == true) {
          free_memory_block((void *)wvc_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[DRYTRB] == true) {
          free_memory_block((void *)dry_ds,  SRC_NUM);
        }
        free_memory_block((void *)ion_ds,  ns);
        return (0);
      }
    }
  }
  if (ERROR_FLAG[FQSERR] == true) {
    if ((fqs_ds = (double *)calloc(I, sizeof(double))) == NULL) {;
      printf("AIRS: ERROR: calloc for fqs_ds.\n");

      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, SRC_NUM);
      if (ERROR_FLAG[TWVTRB] == true) {
        free_memory_block((void *)wvc_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[DRYTRB] == true) {
        free_memory_block((void *)dry_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[IONTRB] == true) {
        free_memory_block((void *)ion_ds,  SRC_NUM);
      }
      return (0);
    }
  }

/*
==============================================================
---------------------- Niell Mapping Function
*/

  DOY = (int)(MJD(TimUTC[0], TimUTC[1], TimUTC[2],
                  TimUTC[3], TimUTC[4], TimUTC[5],     0.0)
            - MJD(TimUTC[0],         1,         0,
                          0,         0,         0,     0.0));
  for (iant=0; iant<GRT_NUM; iant++) {
    nmf20_coeffi(ant_prm[iant].LLH[2], DOY,
                 &ant_prm[iant].nmfh[0], &ant_prm[iant].nmfh_hcr[0],
                 &ant_prm[iant].nmfw[0]);
  }

/*
==============================================================
*/

  sprintf(string, "**** UVW Calculation ****");
  if (TV_SWT == false) {
    printf("%s\n", string);
  } else if (TV_SWT == true) {
    comment_disp(&cmnt, comment, string, true);
  }

  for (ns=0; ns<SRC_NUM; ns++) {
    nant[ns] = 0;
    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;
    for (iant=0; iant<SRT_NUM; iant++) {
      init_l[iant] = dpi;
    }

    for (iobs=0; iobs<data_num.nobs; iobs++) {
      timUTC[5] = TimUTC[5] + iobs;
      i = uvw_calc(100.0 * (float)(ns + 1),
              ANT_NUM, GRT_NUM,  data_num.nobs, iobs,
              timUTC, UT1_UTC, ant_prm, src[ns], sun,
              wvc, dry, ion, Cw, Cd, Ci, CI,
              wvc_ds[ns], dry_ds[ns], ion_ds[ns], fqs_ds,
              false, ERROR_FLAG, W, dUT1, int_obs[ns],
              OBS_T, OBS_p, dz[1], false,
              srt, init_l, TRK_NUM, trk_pos, srt_link,
              sep_angle_limit_from_earth_limb, MSTID);
/****
      usleep(1);
Although the above ``usleep'' is not needed, the performance could 
be stable with this usleep in order to disturb SEGMENTATION FAULT. 
****/

      nant[ns] += i;
    }
    az_rotation(GRT_NUM, data_num, int_obs[ns], ant_prm);
  }
  az_adjustment(GRT_NUM, data_num, int_obs, ant_prm);

  for (ns=0; ns<SRC_NUM; ns++) {
    if (nant[ns] == 0) {
      printf("ERROR: ARIS: NO valid data because the elevation ");
      printf("angle is below the limit.\n");

      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, SRC_NUM);
      if (ERROR_FLAG[TWVTRB] == true) {
        free_memory_block((void *)wvc_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[DRYTRB] == true) {
        free_memory_block((void *)dry_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[IONTRB] == true) {
        free_memory_block((void *)ion_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[FQSERR] == true) {
        free(fqs_ds);
      }
      return (0);
    }
  }

  for (ns=0; ns<SRC_NUM; ns++) {
    uv_max_tmp[ns] = 0.0;
  }
  if (uv_display(data_num.nobs, uv_max_tmp, ANT_NUM, 0, ANT_NUM, 0, ANT_NUM,
                 TimUTC, UT1_UTC,
                 ant_prm, src, src_flag,
                 int_obs, wave_length, nswt, false, 1.25, &c_dummy, NO_STRUCTURE) == -1) {
    if (SRT_NUM >= 1) {
      free (srt_link);
    }
    free_memory_block((void *)dz, 2);
    free_memory_block((void *)int_obs, SRC_NUM);
    if (ERROR_FLAG[TWVTRB] == true) {
      free_memory_block((void *)wvc_ds,  SRC_NUM);
    }
    if (ERROR_FLAG[DRYTRB] == true) {
      free_memory_block((void *)dry_ds,  SRC_NUM);
    }
    if (ERROR_FLAG[IONTRB] == true) {
      free_memory_block((void *)ion_ds,  SRC_NUM);
    }
    if (ERROR_FLAG[FQSERR] == true) {
      free(fqs_ds);
    }
    return (0);
  }
/****
  if (uv_max_tmp[0] > uv_max_tmp[1]) {
    uv_max = uv_max_tmp[0];
  } else {
    uv_max = uv_max_tmp[1];
  }
****/
  uv_max = uv_max_tmp[0];

/*
---- ACA_DEMO ----
*/

/******** Asaki et al. ALMA MEMO 535 (2005)
  if (array.ID == ACA) {
    for (ns=0; ns<SRC_NUM; ns++) {
      for (iant=0; iant<GRT_NUM; iant++) {
        i = data_num.nobs * iant;
        for (iobs=0; iobs<data_num.nobs; iobs++) {
          int_obs[ns][i].az = fabs(src[ns].RA2k);
          int_obs[ns][i].el = fabs(src[ns].DC2k);
          i++;
        }
      }
    }
  }
********/

/*
==============================================================
*/

  if (ERROR_FLAG[TWVTRB] == true) {
    sprintf(string, "**** Tropospheric Phase Screen ****");
    if (TV_SWT == false) {
      printf("%s\n", string);
    } else if (TV_SWT == true) {
      comment_disp(&cmnt, comment, string, true);
    }

    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;
    if (atmospheric_fluctuation(
           SRC_NUM,
           NOISE_MTRX, GRT_NUM, data_num.nobs,  wvc_ds, timUTC, ut1_utc,
           ant_prm, wvc, src, OBS_T, OBS_p, false, nseries, SITE_NUM,
           NELEM, AZEL_FIX,
           &cmnt, comment, TV_SWT) == -1) {
      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, SRC_NUM);
      free_memory_block((void *)wvc_ds,  SRC_NUM);
      if (ERROR_FLAG[DRYTRB] == true) {
        free_memory_block((void *)dry_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[IONTRB] == true) {
        free_memory_block((void *)ion_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[FQSERR] == true) {
        free(fqs_ds);
      }
      return 0;
    }
  }

/*
==============================================================
*/

  if (ERROR_FLAG[DRYTRB] == true) {
    sigma = 4.0e-6 * sqrt(10.0) / speed_of_light;
    if (atmospheric_fluctuation(
           SRC_NUM,
           NOISE_MTRX, GRT_NUM, data_num.nobs,  dry_ds, timUTC, ut1_utc,
           ant_prm, dry, src, OBS_T, OBS_p, false, nseries, SITE_NUM,
           NELEM, AZEL_FIX,
           &cmnt, comment, TV_SWT) == -1) {
      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, SRC_NUM);
      free_memory_block((void *)wvc_ds,  SRC_NUM);
      free_memory_block((void *)dry_ds,  SRC_NUM);
      if (ERROR_FLAG[IONTRB] == true) {
        free_memory_block((void *)ion_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[FQSERR] == true) {
        free(fqs_ds);
      }
      return 0;
    }

    for (ns=0; ns<SRC_NUM; ns++) {
      for (i=0; i<data_num.nobs * ANT_NUM; i++) {
        dry_ds[ns][i] *= 2.0e-14;
      }
    }
  }

/*
==============================================================
*/

  if (ERROR_FLAG[IONTRB] == true) {
    sprintf(string, "**** Ionospheric Phase Screen ****");
    if (TV_SWT == false) {
      printf("%s\n", string);
    } else if (TV_SWT == true) {
      comment_disp(&cmnt, comment, string, true);
    }

    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;
    if (atmospheric_fluctuation(
           SRC_NUM,
           NOISE_MTRX, GRT_NUM, data_num.nobs,  ion_ds, timUTC, ut1_utc,
           ant_prm,
           ion, src, OBS_T, OBS_p, false,
           nseries, SITE_NUM, NELEM, AZEL_FIX,
           &cmnt, comment, TV_SWT) == -1) {
      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, SRC_NUM);
      if (ERROR_FLAG[TWVTRB] == true) {
        free_memory_block((void *)wvc_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[DRYTRB] == true) {
        free_memory_block((void *)dry_ds,  SRC_NUM);
      }
      free_memory_block((void *)ion_ds,  SRC_NUM);
      if (ERROR_FLAG[FQSERR] == true) {
        free(fqs_ds);
      }
      return 0;
    }

/*
--------------------------------------------------------------
*/

    sprintf(string, "**** TEC MS-TID Simulation ****");
    if (TV_SWT == false) {
      printf("%s\n", string);
    } else if (TV_SWT == true) {
      comment_disp(&cmnt, comment, string, true);
    }

    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;
    GRT_TEC_TID(GRT_NUM, data_num.nobs, ion_ds, ant_prm, ion, int_obs, MSTID);
  }

/*
==============================================================
*/

  if (ERROR_FLAG[FQSERR] == true) {
    sprintf(string, "**** Frequency Standard Error ****");
    if (TV_SWT == false) {
      printf("%s\n", string);
    } else if (TV_SWT == true) {
      comment_disp(&cmnt, comment, string, true);
    }

    for (iant=0; iant<ANT_NUM; iant++) {
      i = data_num.nobs * iant;
      for (iobs=0; iobs<data_num.nobs; iobs++) {
        fqs_ds[i] = 0.0;
        i++;
      }

      if (ant_prm[iant].FRQSTD != 0) {
        if (freq_std_instability(data_num.nobs, fqs_ds+data_num.nobs*iant,
                               ant_prm[iant].FRQSTD) == -1) {
          if (SRT_NUM >= 1) {
            free (srt_link);
          }
          free_memory_block((void *)dz, 2);
          free_memory_block((void *)int_obs, SRC_NUM);
          if (ERROR_FLAG[TWVTRB] == true) {
            free_memory_block((void *)wvc_ds,  SRC_NUM);
          }
          if (ERROR_FLAG[DRYTRB] == true) {
            free_memory_block((void *)dry_ds,  SRC_NUM);
          }
          free_memory_block((void *)ion_ds,  SRC_NUM);
          if (ERROR_FLAG[FQSERR] == true) {
            free(fqs_ds);
          }
          return 0;
        }
      }
    }
  }

/*
==============================================================
*/

  while (1) {

/*
----------------------------------------------------
*/

    if (ERROR_FLAG[TDSECZ] == true) {
      for (iant=0; iant<GRT_NUM; iant++) {
        dz[1][iant].trp = dz[0][iant].trp;
      }
    }

    if (ERROR_FLAG[IDSECZ] == true) {
      for (iant=0; iant<GRT_NUM; iant++) {
        dz[1][iant].tec = dz[0][iant].tec;
      }
    }

/*
----------------------------------------------------
*/

    for (i=0; i<6; i++) {
      timUTC[i] = TimUTC[i];
    }
    ut1_utc = UT1_UTC;
    data_num.sobs =    0;
    data_num.eobs = nobs;
    data_num.nobs = nobs;
    if (TV_SWT == true) {
      cmnt.ncol  = 11;
      cmnt.xmin  = 0.02;
      cmnt.xmax  = 0.98;
      cmnt.ymin  = 0.01;
      cmnt.ymax  = 0.16;
      cmnt.pitch = 0.03;
    }
    if (err_parameter_set(ANT_NUM, GRT_NUM, SRT_NUM,
              &BGN_ANT_I, &END_ANT_I, &BGN_ANT_J, &END_ANT_J,
              &array, &TRP_CONDITION, &ION_CONDITION,
              &CW, &CI, wvc,
              dz[1], src_proc, ch_file,
              wave_id, wave_length, nu,
              &band_width, &nfrq, &cfact,
              &reference_phase_process_mode,
              ERROR_FLAG, srt, &data_num,
              timUTC, &ut1_utc,
              &nswt,  &apparent_tgt_on, ant_prm, &TRK_NUM, trk_pos,
              &cmnt, comment, TV_SWT, cursor_pos, pgid[0]) == -1) {
      return 1;
    }

    if (TV_SWT == true) {
      cpgslct(pgid[0]);
      cpgpap(1.5*pgpap_prm, 1.00);
      cpgsch(1.5*pgpap_prm/12.0);
      cpgsvp(0.0, 1.0, 0.0, 1.0);
      cpgswin(0.0, 1.0, 0.0, 1.0);
      cpgsci(11);
      cpgrect(0.07, 0.47, 0.50, 0.98);
      cpgrect(0.57, 0.97, 0.50, 0.98);

      cmnt.ncol  = 14;
      cmnt.xmin  = 0.02;
      cmnt.xmax  = 0.98;
      cmnt.ymin  = 0.05;
      cmnt.ymax  = 0.48;
      cmnt.pitch = 0.03;
      comment_init(&cmnt, comment, true);
    }

    for (ns=0; ns<SRC_NUM; ns++) {
      SRCID[ns]           = src_proc[ns].positionID;
      src_flag[ns]        = 100.0 * (float)(SRCID[ns] + 1);

      src_proc[ns].RA2k   = src[SRCID[ns]].RA2k;
      src_proc[ns].DC2k   = src[SRCID[ns]].DC2k;
      src_proc[ns].s2k[0] = src[SRCID[ns]].s2k[0];
      src_proc[ns].s2k[1] = src[SRCID[ns]].s2k[1];
      src_proc[ns].s2k[2] = src[SRCID[ns]].s2k[2];

      src_proc[ns].RA     = src[SRCID[ns]].RA;
      src_proc[ns].DC     = src[SRCID[ns]].DC;
      src_proc[ns].s[0]   = src[SRCID[ns]].s[0];
      src_proc[ns].s[1]   = src[SRCID[ns]].s[1];
      src_proc[ns].s[2]   = src[SRCID[ns]].s[2];

      src_proc[ns].RA_e   = src[SRCID[ns]].RA_e;
      src_proc[ns].DC_e   = src[SRCID[ns]].DC_e;
      src_proc[ns].s_e[0] = src[SRCID[ns]].s_e[0];
      src_proc[ns].s_e[1] = src[SRCID[ns]].s_e[1];
      src_proc[ns].s_e[2] = src[SRCID[ns]].s_e[2];
    }

/*
----------------------------------------------------
*/

    if (ERROR_FLAG[TWVTRB] == true || ERROR_FLAG[DRYTRB] == true) {
      for (iant=0; iant<GRT_NUM; iant++) {
        Cw[iant] = CW;
      }
    }

/*
----------------------------------------------------
*/

    if (ERROR_FLAG[DRYTRB] == true) {
      for (iant=0; iant<GRT_NUM; iant++) {
        Cd[iant] = CD * 3.0e-3 * (wave_length[0]/3.0e-3);
        if (TRP_CONDITION == 0) {
          Cd[iant] *=  2.0;
        } else if (TRP_CONDITION == 1) {
          Cd[iant] *=  3.0;
        } else if (TRP_CONDITION == 2) {
          Cd[iant] *=  4.0;
        }
      }
    }

/*
----------------------------------------------------
*/

    antenna_selection(&ANT_NUM, &GRT_NUM, &SRT_NUM, wave_id,
                      grt_elevation_limit, ant_prm, antenna_code,
                      antenna_list_file, false);
    for (iant=0; iant<ANT_NUM; iant++) {
      ant_prm[iant].ERR[0] = ant_prm[iant].XYZ[0] + ant_err[iant].ERR[0];
      ant_prm[iant].ERR[1] = ant_prm[iant].XYZ[1] + ant_err[iant].ERR[1];
      ant_prm[iant].ERR[2] = ant_prm[iant].XYZ[2] + ant_err[iant].ERR[2];
      ant_prm[iant].d_gain = ant_err[iant].d_gain;
      for (i=0; i<N_WAVE; i++) {
        ant_prm[iant].LOPHS[i] = ant_err[iant].LOPHS[i];
      }
    }

#ifdef __ANT_DEBUG__
    for (iant=0; iant<ANT_NUM; iant++) {
      printf("__ANT_DEBUG__  %d  %s  %d  %d\n",
        iant, ant_prm[iant].IDC, ant_prm[iant].ARRAY, ant_prm[iant].UFL);
      printf("__ANT_DEBUG__  %lf  %lf  %lf\n",
        ant_prm[iant].XYZ[0], ant_prm[iant].XYZ[1], ant_prm[iant].XYZ[2]);
      printf("__ANT_DEBUG__  %lf  %lf  %lf\n",
        ant_prm[iant].LLH[0], ant_prm[iant].LLH[1], ant_prm[iant].LLH[2]);
      printf("__ANT_DEBUG__  %lf  %lf  %lf\n",
        ant_prm[iant].ERR[0], ant_prm[iant].ERR[1], ant_prm[iant].ERR[2]);
      printf("__ANT_DEBUG__  %lf  %lf  %lf\n",
        ant_prm[iant].OFS[0], ant_prm[iant].OFS[1], ant_prm[iant].OFS[2]);
      printf("__ANT_DEBUG__  %lf  %lf  %lf  %lf  %lf\n",
            ant_prm[iant].AZSV, ant_prm[iant].AZSA,
            ant_prm[iant].ELSV, ant_prm[iant].ELSA, ant_prm[iant].FQSST);
      for (ns=0; ns<SRC_NUM; ns++) {
        printf("__ANT_DEBUG__  %lf  %lf  %lf  %lf\n",
          ant_prm[iant].Dm[ns],  ant_prm[iant].Ae[ns],
          ant_prm[iant].Trx[ns], ant_prm[iant].Tsky[ns]);
      }
      printf("__ANT_DEBUG__\n");
    }
#endif

    for (iant=0; iant<ANT_NUM; iant++) {
      ant_prm[iant].NOSTA = iant + 1;
    }
    if (ERROR_FLAG[THRMNS] == true) {
      for (iant=0; iant<ANT_NUM; iant++) {
        for (ns=0; ns<SRC_NUM; ns++) {
          if (ant_prm[iant].WID[ns] >= 0) {
            sefd[ns] = SEFD(ant_prm[iant].Trx[ns], ant_prm[iant].Tsky[ns],
                            ant_prm[iant].Dm[ns],  ant_prm[iant].Ae[ns]);
          }
        }

        if (ant_prm[iant].WID[0] >= 0 && ant_prm[iant].WID[1] >= 0) {
          sprintf(string, "SEFD (%s) : %8.2lf Jy    %8.2lf Jy",
                  ant_prm[iant].IDC, sefd[0] * 1.0e26, sefd[1] * 1.0e26);
        } else if (ant_prm[iant].WID[0] >= 0 && ant_prm[iant].WID[1] < 0) {
          sprintf(string, "SEFD (%s) : %8.2lf Jy         --  Jy",
                  ant_prm[iant].IDC, sefd[0] * 1.0e26);
        } else if (ant_prm[iant].WID[0] < 0 && ant_prm[iant].WID[1] >= 0) {
          sprintf(string, "SEFD (%s) :      --  Jy    %8.2lf Jy",
                  ant_prm[iant].IDC, sefd[1] * 1.0e26);
        }
        if (ant_prm[iant].WID[0] >= 0 || ant_prm[iant].WID[1] >= 0) {
          if (TV_SWT == true) {
            comment_disp(&cmnt, comment, string, true);
          } else {
            printf("%s\n", string);
          }
        }
      }
    }

/*
---- SRT Slew Time Setting ----  2007.10.15
*/

    SEPANG = sepang(src_proc[0].s2k, src_proc[1].s2k) * 180.0 / dpi;
    for (iant=GRT_NUM; iant<ANT_NUM; iant++) {
#ifdef __CMG_FULL_SPEC__
      if (SEPANG <= 1.0) {
        ant_prm[iant].slewt =  5.0;
      } else {
        ant_prm[iant].slewt =  5.0 + (SEPANG - 1.0) *  5.0;
      }
#elif defined __CMG_HALF_SPEC__
      if (SEPANG <= 1.0) {
        ant_prm[iant].slewt = 10.0;
      } else {
        ant_prm[iant].slewt = 10.0 + (SEPANG - 1.0) * 10.0;
      }
#endif
    }

/*
----------------------------------------------------
*/

    nbase =  (ANT_NUM * (ANT_NUM - 1)) / 2;
    I = data_num.nobs * ANT_NUM;
    for (ns=0; ns<SRC_NUM; ns++) {
      for (i=0; i<I; i++) {
        int_obs[ns][i].u         = (double)0.0;
        int_obs[ns][i].v         = (double)0.0;
        int_obs[ns][i].w         = (double)0.0;
        int_obs[ns][i].grp_delay = (double)0.0;
        int_obs[ns][i].ion_delay = (double)0.0;
        int_obs[ns][i].local_phs = (double)0.0;
        int_obs[ns][i].az        = (double)0.0;
        int_obs[ns][i].el        = (double)0.0;
        int_obs[ns][i].wt        = (float)src_flag[ns];
        int_obs[ns][i].amp_error = 1.0;
      }
    }

/*
---- Switching Flagging -----------------------------
*/

    for (ns=0; ns<SRC_NUM; ns++) {
      for (i=0; i<6; i++) {
        timUTC[i] = TimUTC[i];
      }
      ut1_utc = UT1_UTC;
      for (iant=0; iant<SRT_NUM; iant++) {
        init_l[iant] = dpi;
      }

      for (iobs=0; iobs<data_num.nobs; iobs++) {
        timUTC[5] = TimUTC[5] + iobs;
        for (j=0; j<ERROR_NUM; j++) {
          if (ERROR_FLAG[j] == true) {
            ERROR_FLAG[j] = true;
          } else {
            ERROR_FLAG[j] = false;
          }
        }
        uvw_calc(src_flag[ns],
                ANT_NUM, GRT_NUM,  data_num.nobs, iobs, timUTC, ut1_utc,
                ant_prm, src_proc[ns], sun,
                wvc, dry, ion, Cw, Cd, Ci, CI,
                wvc_ds[SRCID[ns]], dry_ds[SRCID[ns]], ion_ds[SRCID[ns]], fqs_ds,
                false,  ERROR_FLAG, W, dUT1, int_obs[ns],
                OBS_T, OBS_p, dz[1], false,
                srt, init_l, TRK_NUM, trk_pos, srt_link,
                sep_angle_limit_from_earth_limb, MSTID);
      }
      az_rotation(GRT_NUM, data_num, int_obs[ns], ant_prm);
    }
    az_adjustment(GRT_NUM, data_num, int_obs, ant_prm);

#ifndef __DEMO__
    if (apparent_tgt_on == 0.0) {
      ON_TIME[0] = nswt / 2;
      ON_TIME[1] = nswt / 2;
    } else {
      ON_TIME[0] = (int)lrint(apparent_tgt_on * 1.0e-2 * (double)nswt);
      ON_TIME[1] = nswt - ON_TIME[0];
    }
    if (nswt != 0) {
      for (ns=0; ns<SRC_NUM; ns++) {
        noffset = (1 - ns) * ON_TIME[1 - ns];
        for (iant=0; iant<ANT_NUM; iant++) {
          if (iant < GRT_NUM) {
            azl = ant_prm[iant].AZSV * ant_prm[iant].AZSV / ant_prm[iant].AZSA;
            ell = ant_prm[iant].ELSV * ant_prm[iant].ELSV / ant_prm[iant].ELSA;
          } else if (iant >= GRT_NUM) {
            nslewt = (int)lrint(ant_prm[iant].slewt);
/*XXXXXXXX*/
/**
            nslewt = 0;
**/
/*XXXXXXXX*/
          }

          if (nslewt >= ON_TIME[ns]) {
            printf("wARNING: Switching Cycle Time (%d s) is too short\n", nswt);
            printf("WARNING: to assign ON-SOURCE (%d s) ", ON_TIME[ns]);
            printf("to Source-%d for %s.\n", ns+1, (ant_prm+iant)->IDC);
            printf("wARNING: It is suggested to lengthen Switching Cycle Time ");
            printf("and/or adjust TARGET ON ratio.\n");
          }

          iobs = 0;
          while (iobs < data_num.nobs) {
            iave = (iobs - noffset) % nswt;
            i = iant * data_num.nobs + iobs;

            if (iave == 0) {
              if (iant < GRT_NUM) {
                nslewt = (int)lrint(
                          slew_time(int_obs[0][i].az, int_obs[1][i].az, azl,
                                    ant_prm[iant].AZSV, ant_prm[iant].AZSA,
                                    int_obs[0][i].el, int_obs[1][i].el, ell,
                                    ant_prm[iant].ELSV, ant_prm[iant].ELSA));
              } else if (iant >= GRT_NUM && ERROR_FLAG[SRTAER] == true) {
                pnt_err_tau = 0.684e-12 * gauss_dev();

/*
----
               [2009/07/31]
                 1-sigma Attitude Error : 1.9/1000 deg
             --> 1-sigma Delay Error due to the Attitude Error : 0.684 ps
               [2007/01/01]
                 5/1000-deg Attitude Error : 1.8-ps Delay Error
             --> 1/1000-deg Attitude Error : 0.36-ps Delay Error
----
*/

/*XXXXXXXX*/
/**
                nslewt = 0;
**/
/*XXXXXXXX*/
              }
            }

            if (nu[0] != nu[1] && ant_prm[iant].FQSST > nslewt) {
              nslewt = ant_prm[iant].FQSST;
            }

            if        (iave <  nslewt || iave >= ON_TIME[ns]) {
              int_obs[ns][i].wt = src_flag[ns] + (float)OFF_SOURCE;
            } else if (iave >= nslewt && iave <  ON_TIME[ns]) {
              int_obs[ns][i].wt = src_flag[ns];
              if (iant >= GRT_NUM && ERROR_FLAG[SRTAER] == true) {
                int_obs[ns][i].local_phs += pnt_err_tau;
              }
            }
            iobs++;
          }
        }
      }
    }

/*
------------------------------
*/

    for (ns=0; ns<SRC_NUM; ns++) {
      pixel_calc(&pix_uvl[ns], &pix_mas[ns], wave_length[ns],
                 uv_max, &uv_factor[ns], FLDMAX, 1);
      lftmp = log10(pix_mas[ns]);
      if (lftmp < 0.0) {
        lftmp = pow(10.0, fabs(floor(lftmp)));
        pix_mas[ns] = lrint(pix_mas[ns] * lftmp) / lftmp;
        pixel_calc(&pix_uvl[ns], &pix_mas[ns], wave_length[ns],
                   uv_max, &uv_factor[ns], FLDMAX, -1);
      }
    }

/*
----------------------------------------------------
*/

/**** ALMA B2B Phase Jump Simulation ****/

    if (ERROR_FLAG[LOPJMP] == true) {
      for (ns=0; ns<1; ns++) {
        for (iant=0; iant<ANT_NUM; iant++) {
          if (ant_prm[iant].lo_phs_jmp_val != 0.0) {
            I = 0;
            J = 0;
            if (ant_prm[iant].lo_phs_jmp_tim >= 0.0) {
              I = (int)rint(ant_prm[iant].lo_phs_jmp_tim * (float)data_num.nobs);
              J = data_num.nobs;
            } else if (ant_prm[iant].lo_phs_jmp_val < 0.0) {
              I = 0;
              J = (int)rint((1.0 + ant_prm[iant].lo_phs_jmp_tim) * (float)data_num.nobs);
            }
            i = data_num.nobs * iant;
            for (iobs=I; iobs<J; iobs++) {
              int_obs[ns][i].local_phs += ant_prm[iant].lo_phs_jmp_val / 180.0 * dpi;
              i++;
            }
          }
        }
      }
    }

/*
----------------------------------------------------
*/

    for (ns=0; ns<SRC_NUM; ns++) {

      src[SRCID[ns]].flux        = src_proc[ns].flux;
      src[SRCID[ns]].morphology  = src_proc[ns].morphology;

#ifdef __UV_DEBUG__
      sprintf(string, "uvtxt-%d.txt", ns+1);
      fp = fopen(string, "w");
#endif /*__UV_DEBUG__*/
      if (ns == 0) {
        sprintf(string, "**** Target data processed.      ****");
      } else {
        sprintf(string, "**** Reference data processed.   ****");
      }
      if (TV_SWT == false) {
        printf("%s\n", string);
      } else if (TV_SWT == true) {
        comment_disp(&cmnt, comment, string, true);
      }

      I = data_num.nobs * nbase;
      if ((bluvw[ns] = (struct baseline_uvw *)
           calloc(I, sizeof(struct baseline_uvw))) == NULL) {
        printf("ARIS: ERROR: memory allocation: bluvw\n");

        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        if (ERROR_FLAG[TWVTRB] == true) {
          free_memory_block((void *)wvc_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[DRYTRB] == true) {
          free_memory_block((void *)dry_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[IONTRB] == true) {
          free_memory_block((void *)ion_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[FQSERR] == true) {
          free(fqs_ds);
        }
        return (0);
      }

      if ((fringe_weight[ns] = (float *)calloc(I, sizeof(float))) == NULL) {
        printf("ARIS: ERROR: memory allocation: fringe_weight.\n");

        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        if (ERROR_FLAG[TWVTRB] == true) {
          free_memory_block((void *)wvc_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[DRYTRB] == true) {
          free_memory_block((void *)dry_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[IONTRB] == true) {
          free_memory_block((void *)ion_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[FQSERR] == true) {
          free(fqs_ds);
        }
        free_memory_block((void *)bluvw,         ns);
        return (0);
      }

      for (i=0; i<I; i++) {
        bluvw[ns][i].u = (double)0.0;
        bluvw[ns][i].v = (double)0.0;
        bluvw[ns][i].w = (double)0.0;
        fringe_weight[ns][i] = 0.0;
      }

      I = data_num.nobs * nbase * nfrq;
      if ((frng[ns]  = (struct fringe *)
           calloc(I, sizeof(struct fringe))) == NULL) {
        printf("ARIS: ERROR: memory allocation: frng.\n");

        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        if (ERROR_FLAG[TWVTRB] == true) {
          free_memory_block((void *)wvc_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[DRYTRB] == true) {
          free_memory_block((void *)dry_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[IONTRB] == true) {
          free_memory_block((void *)ion_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[FQSERR] == true) {
          free(fqs_ds);
        }
        free_memory_block((void *)bluvw,         ns);
        free_memory_block((void *)fringe_weight, ns);
        return (0);
      }
      for (i=0; i<I; i++) {
        frng[ns][i].rl = (double)0.0;
        frng[ns][i].im = (double)0.0;
        frng[ns][i].wt = (double)0.0;
      }

/*
----------------------------------------------------
*/

      for (i=0; i<6; i++) {
        timUTC[i] = TimUTC[i];
      }
      ut1_utc = UT1_UTC;
      for (iant=0; iant<SRT_NUM; iant++) {
        init_l[iant] = dpi;
      }

      for (iobs=0; iobs<data_num.nobs; iobs++) {
        timUTC[5] = TimUTC[5] + iobs;
        uvw_calc(src_flag[ns],
                ANT_NUM, GRT_NUM,  data_num.nobs, iobs, timUTC, ut1_utc,
                ant_prm, src_proc[ns], sun,
                wvc, dry, ion, Cw, Cd, Ci, CI,
                wvc_ds[SRCID[ns]], dry_ds[SRCID[ns]], ion_ds[SRCID[ns]], fqs_ds,
                true,  ERROR_FLAG, W, dUT1, int_obs[ns],
                OBS_T, OBS_p, dz[1], true,
                srt, init_l, TRK_NUM, trk_pos, srt_link,
                sep_angle_limit_from_earth_limb, MSTID);
      }
      az_rotation(GRT_NUM, data_num, int_obs[ns], ant_prm);
      if (ns == 1) {
        az_adjustment(GRT_NUM, data_num, int_obs, ant_prm);
      }

/*
-------------------------------------------------------------
*/

      if ((mapr = (float *)calloc(FLDMAX*FLDMAX, sizeof(float))) == NULL) {
        printf("ERROR: ARIS: memory allocation: mapr.\n");
        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        if (ERROR_FLAG[TWVTRB] == true) {
          free_memory_block((void *)wvc_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[DRYTRB] == true) {
          free_memory_block((void *)dry_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[IONTRB] == true) {
          free_memory_block((void *)ion_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[FQSERR] == true) {
          free(fqs_ds);
        }
        free_memory_block((void *)bluvw,         SRC_NUM);
        free_memory_block((void *)fringe_weight, SRC_NUM);
        free_memory_block((void *)frng,          SRC_NUM);
        return (0);
      }
      if ((mapi = (float *)calloc(FLDMAX*FLDMAX, sizeof(float))) == NULL) {
        printf("ERROR: ARIS: memory allocation: mapi.\n");
        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        if (ERROR_FLAG[TWVTRB] == true) {
          free_memory_block((void *)wvc_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[DRYTRB] == true) {
          free_memory_block((void *)dry_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[IONTRB] == true) {
          free_memory_block((void *)ion_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[FQSERR] == true) {
          free(fqs_ds);
        }
        free_memory_block((void *)bluvw,         SRC_NUM);
        free_memory_block((void *)fringe_weight, SRC_NUM);
        free_memory_block((void *)frng,          SRC_NUM);
        free (mapr);
        return (0);
      }

      if (visibility_calc(mapr, mapi, FLDMAX, &pix_uvl[ns], &pix_mas[ns],
                          wave_length[ns], uv_max, &uv_factor[ns],
                          src_proc[ns], ch_file+ns,
                          0.05+(float)ns*0.50,
                          0.45+(float)ns*0.50, 0.60, 0.60+0.40,
                          true, ns, pgid[0]) == __NG__) {
        printf("ERROE: aris: visibility_calc returns ERROR.\n");
        printf("ERROR: aris: abnormaly ended.\n");
        if (SRT_NUM >= 1) {
          free (srt_link);
        }
        free_memory_block((void *)dz, 2);
        free_memory_block((void *)int_obs, SRC_NUM);
        if (ERROR_FLAG[TWVTRB] == true) {
          free_memory_block((void *)wvc_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[DRYTRB] == true) {
          free_memory_block((void *)dry_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[IONTRB] == true) {
          free_memory_block((void *)ion_ds,  SRC_NUM);
        }
        if (ERROR_FLAG[FQSERR] == true) {
          free(fqs_ds);
        }
        free_memory_block((void *)bluvw,         SRC_NUM);
        free_memory_block((void *)fringe_weight, SRC_NUM);
        free_memory_block((void *)frng,          SRC_NUM);
        free (mapr);
        free (mapi);
        return (0);
      }

      if (TV_SWT == true) {
        cpgsvp(0.0, 1.0, 0.0, 1.0);
        cpgswin(0.0, 1.0, 0.0, 1.0);
      }

/*
----------------------------------------------------------------
*/

      d_band_width = band_width / (double)nfrq;

      for (iant=0; iant<ANT_NUM; iant++) {
        for (jant=iant+1; jant<ANT_NUM; jant++) {
          if (baseline_check(iant, jant, ns,
                             BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J,
                             ant_prm) == true) {

/********XXXXXXXXXX********/
            frq_coh_factor[ns] 
               = coherence_factor_calc(nu[ns], 1.0,
                                    ant_prm[iant].FRQSTD, 
                                    ant_prm[jant].FRQSTD);
#ifdef __CF_DEBUG__
              printf("## __CF_DEBUG__   %lf\n", frq_coh_factor[ns]);
#endif /* __CF_DEBUG__ */
/********XXXXXXXXXX********/

            ibase = baseline_number(ANT_NUM, iant, jant);
            if (ERROR_FLAG[THRMNS] == true) {
              sigma = thermal_noise_cal(1.0, cfact,
                             ant_prm[iant].Dm[ns], ant_prm[jant].Dm[ns],
                             ant_prm[iant].Ae[ns], ant_prm[jant].Ae[ns],
                             1.0,  1.0,   d_band_width,  inttim);
            }

            d_lo_phs = 0.0;
            if (ERROR_FLAG[LOPOFS] == false || ERROR_FLAG[LOPJMP] == false) {
              d_lo_phs = diff(ant_prm[iant].LOPHS[wave_id[ns]],
                              ant_prm[jant].LOPHS[wave_id[ns]]);
            }

            for (iobs=0; iobs<data_num.nobs; iobs++) {
              I = iant  * data_num.nobs + iobs;
              J = jant  * data_num.nobs + iobs;
              K = ibase * data_num.nobs + iobs;

              if (int_obs[ns][I].wt >= src_flag[ns] &&
                  int_obs[ns][J].wt >= src_flag[ns]) {
                bluvw[ns][K].u = diff(int_obs[ns][I].u, int_obs[ns][J].u);
                bluvw[ns][K].v = diff(int_obs[ns][I].v, int_obs[ns][J].v);
                bluvw[ns][K].w = diff(int_obs[ns][I].w, int_obs[ns][J].w);

#ifdef __UV_DEBUG__
                printf("%2d %2d %5d %15.4lf\n",
                         iant+1, jant+1, iobs, bluvw[ns][K].w / speed_of_light);
                fprintf(fp, "%2d %2d %5d %15.4lf %15.4lf\n",
                         iant+1, jant+1, iobs, bluvw[ns][K].u, bluvw[ns][K].v);
#endif /* __UV_DEBUG__ */

                vis_smoothing(&fr, &fi, mapr, mapi,
                              bluvw[ns][K].u, bluvw[ns][K].v,
                              FLDMAX, IFREF, JFREF, pix_uvl[ns]);
                tau_err = diff(int_obs[ns][I].grp_delay,
                               int_obs[ns][J].grp_delay);
                ion_err = diff(int_obs[ns][I].ion_delay,
                               int_obs[ns][J].ion_delay);
                phs_err = diff(int_obs[ns][I].local_phs,
                               int_obs[ns][J].local_phs);

                if (ERROR_FLAG[THRMNS] == true) {
                  SIGMA = sqrt(
                          (ant_prm[iant].Trx[ns]
                         + ant_prm[iant].Tsky[ns] / sin(int_obs[ns][I].el))
                        * (ant_prm[jant].Trx[ns]
                         + ant_prm[jant].Tsky[ns] / sin(int_obs[ns][J].el))
                                );
                }

                for (ifrq=0; ifrq<nfrq; ifrq++) {
                  L = (ibase * data_num.nobs + iobs) * nfrq + ifrq;
                  F = nu[ns] + (double)ifrq  * d_band_width;
                  fai_err = 2.0 * dpi * (F * tau_err
                                       + F * phs_err
                                       + ion_err / F) + d_lo_phs;

                  L = (ibase * data_num.nobs + iobs) * nfrq + ifrq;
                  ar = cos(fai_err);
                  ai = sin(fai_err);
                  frng[ns][L].rl = fr * ar - fi * ai;
                  frng[ns][L].im = fr * ai + fi * ar;
                  frng[ns][L].wt = 1.0;

                  if (ERROR_FLAG[THRMNS] == true) {
                    frng[ns][L].rl += SIGMA * sigma * gauss_dev();
                    frng[ns][L].im += SIGMA * sigma * gauss_dev();
                  }
                  if (ERROR_FLAG[AMPERR] == true) {
                    amperr = frq_coh_factor[ns] * sqrt(
                      int_obs[ns][I].amp_error * int_obs[ns][J].amp_error);
                    frng[ns][L].rl *= amperr;
                    frng[ns][L].im *= amperr;
                  }
                }
                fringe_weight[ns][K] = src_flag[ns] + 1.0;
              } else {
#ifdef __DEBUG__
                printf("__DEBUG__   %f  %f  %f  %f  %d  %d\n",
                        int_obs[ns][I].wt, int_obs[ns][J].wt, 
                        src_flag[ns], src_flag[ns],
                        SRCID[ns], SRCID[ns]);
#endif /* __DEBUG__ */
                fringe_weight[ns][K] = 0.0;
              }
            }
          }
        }
      }
      free (mapr);
      free (mapi);
#ifdef __UV_DEBUG__
      fclose (fp);
#endif /*__UV_DEBUG__*/
    }

/*
----------------------------------------
*/

    for (ns=0; ns<SRC_NUM; ns++) {
      ndata = 0;
      for (iant=0; iant<ANT_NUM; iant++) {
        for (jant=iant+1; jant<ANT_NUM; jant++) {
          if (baseline_check(iant, jant, ns,
                             BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J,
                             ant_prm) == true) {
            ibase = baseline_number(ANT_NUM, iant, jant);
            for (iobs=0; iobs<data_num.nobs; iobs++) {
              I = ibase * data_num.nobs + iobs;
              if (fringe_weight[ns][I] > 0.0) {
                ndata++;
              }
            }
          }
        }
      }
      if (ndata == 0 && ANT_NUM >= 2) {
        sprintf(string, "WARNING: ARIS: NO valid data for the source.");
        printf("%s\n", string);
        if (TV_SWT == true) {
          comment_disp(&cmnt, comment, string, true);
        }
        sprintf(string, "WARNING: Please check the observation situation.");
        printf("%s\n", string);
        if (TV_SWT == true) {
          comment_disp(&cmnt, comment, string, true);
        }
        sprintf(string, "WARNING: The problem may be solved, for example, ");
        printf("%s\n", string);
        if (TV_SWT == true) {
          comment_disp(&cmnt, comment, string, true);
        }
        sprintf(string, "WARNING: if the observation time becomes longer.");
        printf("%s\n", string);
        if (TV_SWT == true) {
          comment_disp(&cmnt, comment, string, true);
        }
      }
    }

/*
----------------------------------------------------
*/

    I = data_num.nobs * nbase;
    if ((fringe_weight[SRC_NUM_P1-1]  = (float *)calloc(I, sizeof(float)))
        == NULL) {
      printf("ARIS: ERROR: memory allocation: fringe_weight[%d]\n",
             SRC_NUM_P1 - 1);

      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, SRC_NUM);
      if (ERROR_FLAG[TWVTRB] == true) {
        free_memory_block((void *)wvc_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[DRYTRB] == true) {
        free_memory_block((void *)dry_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[IONTRB] == true) {
        free_memory_block((void *)ion_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[FQSERR] == true) {
        free(fqs_ds);
      }
      free_memory_block((void *)bluvw,         SRC_NUM);
      free_memory_block((void *)fringe_weight, SRC_NUM);
      free_memory_block((void *)frng,          SRC_NUM);
      return (0);
    }
    for (i=0; i<I; i++) {
      fringe_weight[SRC_NUM_P1-1][i] = (double)0.0;
    }

    I = data_num.nobs * nbase * nfrq;
    if ((frng[SRC_NUM_P1-1]  = (struct fringe *)
         calloc(I, sizeof(struct fringe))) == NULL) {
      printf("memory allocation error: frng[%d]\n", SRC_NUM_P1 - 1);

      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, SRC_NUM);
      if (ERROR_FLAG[TWVTRB] == true) {
        free_memory_block((void *)wvc_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[DRYTRB] == true) {
        free_memory_block((void *)dry_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[IONTRB] == true) {
        free_memory_block((void *)ion_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[FQSERR] == true) {
        free(fqs_ds);
      }
      free_memory_block((void *)bluvw,         SRC_NUM);
      free_memory_block((void *)fringe_weight, SRC_NUM_P1);
      free_memory_block((void *)frng,          SRC_NUM);
      return (0);
    }
    for (i=0; i<I; i++) {
      frng[SRC_NUM_P1-1][i].rl = (double)0.0;
      frng[SRC_NUM_P1-1][i].im = (double)0.0;
    }


    if (reference_phase_process_mode == 1) {
      sprintf(string, "Antenna-base phase solution for the reference source.");
      printf("%s\n", string);
      if (TV_SWT == true) {
        comment_disp(&cmnt, comment, string, true);
      }
      antenna_base_phase(ANT_NUM, nfrq,
                    BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J,
                    ant_prm, data_num,
                    data_num.nobs, nswt, nu, frng, fringe_weight);
    } else {
      sprintf(string,
             "Baseline-base phase solution for the reference source.");
      if (TV_SWT == true) {
        comment_disp(&cmnt, comment, string, true);
      } else {
        printf("%s\n", string);
      }
    }

    if (phase_reference(ANT_NUM, nfrq,
                    BGN_ANT_I, END_ANT_I, BGN_ANT_J, END_ANT_J,
                    data_num,
                    data_num.nobs, nswt, nu, frng, fringe_weight) == -1) {
      printf("ARIS: PHASE_REFERNCE: ERROR\n");
      free_memory_block((void *)bluvw,         SRC_NUM);
      free_memory_block((void *)fringe_weight, SRC_NUM_P1);
      free_memory_block((void *)frng,          SRC_NUM_P1);

      if (SRT_NUM >= 1) {
        free (srt_link);
      }
      free_memory_block((void *)dz, 2);
      free_memory_block((void *)int_obs, SRC_NUM);
      if (ERROR_FLAG[TWVTRB] == true) {
        free_memory_block((void *)wvc_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[DRYTRB] == true) {
        free_memory_block((void *)dry_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[IONTRB] == true) {
        free_memory_block((void *)ion_ds,  SRC_NUM);
      }
      if (ERROR_FLAG[FQSERR] == true) {
        free(fqs_ds);
      }

      exit (-1);
    }

/*
------------------------------------------------------------
*/

#endif /* __DEMO__ */

    ana__mode = true;
    if (TV_SWT == true) {
      cpgsvp(0.0, 1.0, 0.0, 1.0);
      cpgswin(0.0, 1.0, 0.0, 1.1);

      ana__mode = false;

      bttn_box[0][0] = 0.020;
      bttn_box[0][1] = 0.450;
      bttn_box[0][2] = 0.020;
      bttn_box[0][3] = bttn_box[0][2] + pitch;
      cursor_pos[0] = 0.5 * (bttn_box[0][0] + bttn_box[0][1]);
      cursor_pos[1] = 0.5 * (bttn_box[0][2] + bttn_box[0][3]);
      off_button(&i, "Continue", bttn_box[0]);

      bttn_box[1][0] = 0.550;
      bttn_box[1][1] = 0.980;
      bttn_box[1][2] = bttn_box[0][2];
      bttn_box[1][3] = bttn_box[0][3];
      off_button(&i, "Return",   bttn_box[1]);

      while (1) {
        cpgcurs(cursor_pos, cursor_pos+1, string);
        if (_button_chk(cursor_pos, bttn_box[0]) == true) {
          ana__mode = true;
          break;
        } else if (_button_chk(cursor_pos, bttn_box[1]) == true) {
          ana__mode = false;
          break;
        }
      }
    }

    if (ana__mode == true) {
      for (j=0; j<ERROR_NUM; j++) {
        if (ERROR_FLAG[j] == true) {
          ERROR_FLAG[j] = true;
        } else {
          ERROR_FLAG[j] = false;
        }
      }
      proc_mode = menu_config(ANT_NUM, GRT_NUM,  SRT_NUM,  TRK_NUM,
                              BGN_ANT_I, END_ANT_I,
                              BGN_ANT_J, END_ANT_J,
                              data_num.nobs, TimUTC, UT1_UTC, ERROR_FLAG, nswt,
                              itmp,
                              ant_prm, ant_err, 
                              dz[1], int_obs, grt_elevation_limit,
                              wave_length, nu, EOP, nbase, src_flag,
                              frng, fringe_weight, inttim, nfrq, band_width,
                              bluvw, src_proc, sun,
                              pix_uvl, pix_mas, srt, trk_pos, srt_link,
                              sep_angle_limit_from_earth_limb, OBS_T, OBS_p,
                              wvc_ds, ion_ds,
                              TV_SWT, cursor_pos, pgid[0], pgid[1]);
    }

/*
-------------------------------
*/

    free_memory_block((void *)bluvw,         SRC_NUM);
    free_memory_block((void *)fringe_weight, SRC_NUM_P1);
    free_memory_block((void *)frng,          SRC_NUM_P1);

    if (proc_mode == _EXIT_ || proc_mode == __NG__) {
      break;
    }
  }

/*
===================================================================
*/

  cpgend();

/*
===================================================================
*/

  if (SRT_NUM >= 1) {
    free (srt_link);
  }
  free_memory_block((void *)dz, 2);
  free_memory_block((void *)int_obs, SRC_NUM);
  if (ERROR_FLAG[TWVTRB] == true) {
    free_memory_block((void *)wvc_ds,  SRC_NUM);
  }
  if (ERROR_FLAG[DRYTRB] == true) {
    free_memory_block((void *)dry_ds,  SRC_NUM);
  }
  if (ERROR_FLAG[IONTRB] == true) {
    free_memory_block((void *)ion_ds,  SRC_NUM);
  }
  if (ERROR_FLAG[FQSERR] == true) {
    free(fqs_ds);
  }

/*
===================================================================
*/
 
  return 1;
}
