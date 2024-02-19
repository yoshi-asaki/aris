#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <mathtools.h>
#include <astrotools.h>

#define __NG__              -10
#define __GO__               10
#define _EXIT_               20
#define PHS_SCR              30
#define ANT_VIS              40
#define RE_LOAD              50

#define __READ__              0
#define _WRITE__              1

#define EPOCH              2000

#define SRC_NUM               2
#define SRC_NUM_P1            3
#define ANTMAX              120
#define SRTMAX                2
#define TRKMAX               20
/****
#define NOISE_MTRX         8192
#define NOISE_MTRX        16384
****/
#define NOISE_MTRX        32768

#define SRC_MORPH_CANDIDATE   6
#define SRC_POINT             0
#define SRC_MULTI_COMP        1
#define SRC_DISK_JETCJET      2
#define SRC_DISK_VSOP2JET     3
#define SRC_CC_COMP           4
#define SRC_BHS_MOD           5

#define LFDMAX              512
#define FLDMAX             2048

#define DAY                   1
#define NIGHT                 0

#define SSF_PHASE             0
#define SSF_DELAY             1

#define CTRL_C                3
#define CTRL_D                4
#define BACKSPACE_KEY         8
#define RETURN_KEY           13

#define NO_DEF_ARRAY          0
#define _VLBI_ARRAY_          1
#define __CONNECTED_          2

#define ANT_RESET            -1
#define ALL_ANT               0
#define VLBA                  1
#define EVN                   2
#define HSA                   3
#define VERA                  4
#define JVN                   5
#define KVN                   6
#define LBA                   7
#define KAVA                  8
#define EALMA                 9
#define ALMA                 10
#define ACA                  11
#define NG_VLA               12
#define STAND_ALONE          20
#define TRACKING_NETWORK     30
#define ORBITING             40
#define SLR_NETWORK          50

#define ONBOARD               0
#define TRACKING              1

#define N_WAVE               20

/**** HEMT ****/
#define ALLBAND           0

#define L_BAND            1
#define S_BAND            2
#define C_BAND            3
#define X_BAND            4
#define KU_BAND           5
#define K_BAND            6
#define Q_BAND            7
#define W_BAND            8

/**** SIS  ****/
#define BAND01           11
#define BAND02           12
#define BAND03           13  /*  84 - 119 GHz */
#define BAND04           14  /* 125 - 163 GHz */
#define BAND05           15  /* 163 - 211 GHz */
#define BAND06           16  /* 211 - 275 GHz */
#define BAND07           17  /* 275 - 370 GHz */
#define BAND08           18  /* 385 - 500 GHz */
#define BAND09           19  /* 602 - 720 GHz */
#define BAND10           20  /* 787 - 950 GHz */

#define TGT               0
#define REF               1

#define ERROR_NUM     15
#define ERROR_NUM_P2  17
#define APOSER      0   /* Antenna POSition ERror              */
#define TPOSER      1   /* Target POSition ERror               */
#define RPOSER      2   /* Reference POSition ERror            */
#define EOPERR      3   /* Earth Orientation Parameter ERRor   */
#define TDSECZ      4   /* Tropospheric Delta SECZ             */
#define IDSECZ      5   /* Ionospheric Delta SECZ              */
#define TWVTRB      6   /* Tropospheric Water Vapor TuRBulence */
#define DRYTRB      7   /* DRY air TuRBulence                  */
#define IONTRB      8   /* IONospheric TuRBulence              */
#define THRMNS      9   /* THeRMal NoiSe                       */
#define FQSERR     10   /* FreQuency Standard ERRor            */
#define LOPOFS     11   /* Local Oscillator Phase OFfSet       */
#define LOPJMP     12   /* Local Oscillator Phase JuMP         */
#define AMPERR     13   /* AMPlitude ERRor                     */
#define SRTAER     14   /* SRT Attitude ERror                  */

#define OFF_SOURCE                    -1.0
#define TELESCOPE_SLEW_SPEED_LIMIT    -2.0
#define GRT_ELEVATION_LIMIT           -4.0
#define SRT_TRACKING_CONDITION_LIMIT  -8.0
#define SRT_EARTH_ECLIPS_SRC         -16.0
#define SRT_EARTH_ECLIPS_SUN         -32.0

#define ALAZ        0   /* ALT AZIMUTH  */
#define EQUA        1   /* EQUATORIAL   */
#define ORBI        2   /* ORBITING     */

#define OFF_TRK   100
#define ON_TRK    200

/*#define TIME_STEP 300*/  /* TIME STEP [sec] */
#define TIME_STEP 30  /* TIME STEP [sec] */

#define NCOMLEN  100

/*
---- TROPOSPHERIC CONDITION ----
*/

/*
VERY_GODD_WVR  : ALMA site condition + WVR phase correction
VERY_GOOD      : ALMA site condition
GOOD           : Good tropospheric condition (Beasely anc Conway)
TYPICAL        : Typical tropospheric condition (Beasely anc Conway)
POOR           : Poor tropospheric condition (Beasely anc Conway)
*/

/******
#define TRP_NUM          7
#define FANTASTIC_ALMA   0
#define VERY_GOOD_ALMA   1
#define GOOD_ALMA        2
#define VERY_GOOD        3
#define GOOD             4
#define TYPICAL          5
#define POOR             6
********/

/*
---- IONOSPHERIC CONDITION ----
*/

#define HALF             0
#define NOMINAL          1
#define DOUBLE           2

/*
---- FREQUENCY STANDARD ----
*/

/*
H_M     : H maser
CSO_10  : CryoCooled Sapphire Oscillator (10-MHz Reference)
CCSO    : CryoCooled Sapphire Oscillator (100-MHz Reference)
TRP1    : Very good tropospheric condition
TRP2    : Very good tropospheric condition
TRP3    : Good tropospheric condition
TRP4    : Typical tropospheric condition
TRP5    : Bad tropospheric condition
*/

#define NONE        0
#define H_M         2
#define CSO_10      8
#define CSO_100    32
#define TRP1      128
#define TRP2      512
#define TRP3     2048
#define TRP4     8196
#define TRP5    32784

/*
---------------------------------------------------------
      STRUCTURES
---------------------------------------------------------
*/

struct array_parameter
        {
          int    ID;
          int    TYPE;
          int    SITE_NUM;
        };

struct antenna_position
        {
          char   IDC[10];  /** IDentificatio Code of antenna **/
          int    ARRAY;    /** Array Code **/
          int    MNTSTA;   /** Mount Type **/
          int    NOSTA;    /** Station ID Number in Process **/
          _Bool  UFL;      /** Use FLag **/
          int    FRQSTD;   /** FReQuency STanDard **/
          double XYZ[3];   /** Position in Cartesian coordinate **/
          double LLH[3];   /** Latitude, Longitude, and Height **/
          double OFS[3];   /** position offset for some reasons **/
          double AAE[3];   /** A priori Antenna position Error **/
          double AZSV;     /** Azimuth Slew Velocity       [rad/s] */
          double ELSV;     /** Elevation Slew Velocity     [rad/s] */
          double AZSA;     /** Azimuth Slew Accerelation   [rad/s/s] */
          double ELSA;     /** Elevation Slew Accerelation [rad/s/s] */
          double ELLIM;    /** Elevation Limit [rad] */

          double ERR[3];   /** ERRor in position               */
        };
struct antenna_parameter
        {
          char   IDC[10];  /** IDentificatio Code of antenna **/
          int    ARRAY;    /** Array Code **/
          int    MNTSTA;   /** Mount Type **/
          int    NOSTA;    /** Station ID Number in Process **/
          _Bool  UFL;      /** Use FLag **/
          int    FRQSTD;   /** FReQuency STanDard **/
          double XYZ[3];   /** Position in Cartesian coordinate **/
          double LLH[3];   /** Latitude, Longitude, and Height **/
          double OFS[3];   /** position offset for some reasons **/
          double AAE[3];   /** A priori Antenna position Error **/
          double AZSV;     /** Azimuth Slew Velocity       [rad/s] */
          double ELSV;     /** Elevation Slew Velocity     [rad/s] */
          double AZSA;     /** Azimuth Slew Accerelation   [rad/s/s] */
          double ELSA;     /** Elevation Slew Accerelation [rad/s/s] */
          double ELLIM;    /** Elevation Limit [rad] */
          double FQSST;    /** Frequency Switch Settle Time [s]  */
          double STAXOF;   /** STAXOF **/

          double ERR[3];   /** ERRor in position               */
          double LOPHS[N_WAVE]; /** Local Oscilator Phase      */
          double d_gain;   /** Delta GAIN Error                */

          int    WVturb;   /** Water Vapor turnulence          */
          int    DAturb;   /** Dry Air turbulence              */
          int    IOturb;   /** Ionosphere turbulence           */

          double Cw;       /** Water vapor structure coefficient */

          double lo_phs_jmp_val; /** LO phase jump value (deg)      */
          double lo_phs_jmp_tim; /** LO phase jump timing (0-1)     */

/****
          double WAVE[SRC_NUM];
****/
          int    WID[SRC_NUM];   /** Rx Use Flag (Wave ID)          */
          int    IF_FLG[16];     /** IF Channel use/nouse Flag      */
          int    priority;       /** tracking priority              */
          double Dm[SRC_NUM];    /** Diameter                       */
          double Ae[SRC_NUM];    /** Aperture Efficiency            */
          double Trx[SRC_NUM];   /** RX Noise Temperature           */
          double Tsky[SRC_NUM];  /** Sky Noise Temperature          */
          double slewt;          /** Frequency Switch Time Interval */
          double nmfh[3];
          double nmfh_hcr[3];
          double nmfw[3];
        };
struct antenna_error_parameter
        {
          char   IDC[10];
          double ERR[3];
          double LOPHS[N_WAVE];
          double d_gain;
        };
struct atmospheric_zenith_error
        {
          double trp;
          double tec;
        };
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
struct char_obs_time
        {
          char  start_t[6][10];
          char  obsd[20];
        };
struct char_srt_info
        {
          char  apo[20];
          char  per[20];
          char  inc[20];
          char  OMG[20];
          char  omg[20];
          char  t_0[20];
          char  d_OMG[20];
          char  d_omg[20];
        };
struct char_src_info
        {
          char  tgt_ra[40];
          char  tgt_dc[40];
          char  ref_ra[40];
          char  ref_dc[40];
          char  mid_ra[40];
          char  mid_dc[40];
          char  dlt_ra[40];
          char  dlt_dc[40];
          char  sepang[40];
          char  posang[40];
        };
struct pair_src_info
        {
          double mid_ra, mid_dc;
          double dlt_ra, dlt_dc;
          double sepang, posang;
        };
struct TID
        {
          double amp;
          double lambda;
          double v[2];
        };
struct srt_data_link
        {
          double FOV_rot_axis[3];
          double FOV_rot_angle;
          int    nzenith;
          _Bool  mask[54000];
        };
struct srt_orbit_parameter
        {
          double a;
          double apogee;
          double perigee;
          double inclination;
          double e;
          double t0;
          double n;
          double Omega;
          double omega;
          double d_Omega;
          double d_omega;
          double ODDA;    /* OD Displacement at Apogee  */
          double ODDP;    /* OD Displacement at Perigee */
          double initial_phase;
          double ef;
          double uplink_freq;
          double downlink_freq;
          int    BODY_X_SUN;
        };
struct source_parameter
        {
          char   name[20];

          double RA2k;
          double DC2k;
          double s2k[3];

          double RA;
          double DC;
          double epoch;
          double s[3];

          double RA_e;
          double DC_e;
          double epoch_e;
          double s_e[3];

          double flux;
          double weight;
          int    morphology;
          int    positionID;
        };
struct st_observable
        {
          double u;
          double v;
          double w;
          double grp_delay;
          double ion_delay;
          double local_phs;
          double amp_error;
          double az;
          double el;
          float  wt;
        };
/****
wt : TARGET          : 100 + weight  / 100 + FLAG
     REFERENCE       : 200 + weight  / 100 + FLAG
     (weight = 0 - 1)
     NO FLAG                      :   0
     OFF SOURCE                   : - 1
     TELESCOPE SLEW SPEED LIMIT   : - 2
     GRT ELEVATION LIMIT          : - 4
     SRT TRACKING CONDITION LIMIT : - 8
     SRT EARTH ECLIPS (SOURCE)    : -16
     SRT EARTH ECLIPS (SUN)       : -32
****/
struct fringe
        {
          double rl;
          double im;
          double wt;
        };
struct baseline_uvw
        {
          double u;
          double v;
          double w;
        };
struct EOP_data
        { 
          double WX, WY;
          double Qt[3][3];
          double W[3][3];
        };
struct earth_shape
        {
          int    n_day, n_ngt, n_shd;
          float  s_day[2][3000];
          float  s_ngt[2][3000];
          float  s_shd[2][3000];
        };
struct comment_param
        {
          int    ncol;
          float  xmin, xmax, ymin, ymax;
          float  pitch;
        };
struct data_number
        {
          int    sobs;
          int    eobs;
          int    nobs;
        };
struct morphology_file_name
        {
          char   mcm[100];
          char   cct[100];
          char   bhs[100];
          char   asc[100];
        };

/*
---------------------------------------------------------
      FUNCTIONS
---------------------------------------------------------
*/

void    free_memory_block
        (
          void   *[],    int  
        );
void    free_memory_block_char
        (
          char    **,    int  
        );
void    free_memory_block_float
        (
          float   **,    int  
        );
void    free_memory_block_double
        (
          double  **,    int  
        );
int     obs_param_set
        (
          _Bool     *,    int        ,
          struct char_srt_info      *,
          struct srt_orbit_parameter                *,
          char      *,    double    *,
          double     ,
          struct char_obs_time      *,    int       *,    double    *,
          double    *,
          int        ,
          struct char_src_info      *,
          struct pair_src_info      *,
          struct source_parameter   *,
          struct source_parameter   *
        );
int     obs_param_input
        (
          _Bool     *,    struct array_parameter    *,
          int       *,
          int       *,    int       *,    int       *,
          struct srt_orbit_parameter                *,
          double    *,
          double     ,
          int       *,    double    *,    double    *,
          struct source_parameter   *,
          struct source_parameter   *,
          char      *,
          struct antenna_parameter  *,    char [][10],
          struct comment_param      *,    char   [][NCOMLEN],
          _Bool      ,    int        ,    float     *,    int       *
        );
int     obs_param_file_io(
          _Bool     *,    char      *,
          struct array_parameter    *,
          int       *,    int       *,    int       *,
          struct srt_orbit_parameter                *,
          double    *,
          double     ,
          int       *,    double    *,    double    *,  int       *,
          struct source_parameter   *,
          struct source_parameter   *,
          char [][10],
          char      *,
          struct  pair_src_info     *,
          struct  char_src_info     *,
          struct  char_srt_info     *,
          struct  char_obs_time     *,
          int
        );
_Bool   baseline_check
        (
          int        ,    int        ,   int        ,   int        ,
          int        ,    int        ,   int        ,
          struct antenna_parameter  *
        );
void    in__src_proc
        (
          int        ,    int       *,   int       *
        );
int     out_src_proc
        (
          int        ,    int
        );
int     ch_time_set
        (
          int       *,    struct char_obs_time  *
        );
int     antenna_selection
        (
          int       *,    int       *,   int    *,
          int       *,    double     ,
          struct antenna_parameter  *,
          char [][10],    char      *,   _Bool
        );
int     number_char_cut
        (
          char      *
        );
void    srt_info_disp
        (
          int        ,    float      ,    float [][4],    float      ,
          char      *,    _Bool     *,
          struct char_srt_info      *
        );
int     err_parameter_set
        (
          int        ,    int        ,    int        ,
          int       *,    int       *,
          int       *,    int       *,
          struct array_parameter    *,
          int       *,    int       *,
          double    *,    double    *,
          struct phase_screen_parameter             *,
          struct atmospheric_zenith_error           *,
          struct source_parameter   *,
          struct morphology_file_name               *,
          int       *,    double    *,    double    *,
          double    *,    int       *,    double    *,
          int       *,
          _Bool     *,
          struct srt_orbit_parameter                *,
          struct data_number        *,
          int       *,    double    *,
          int       *,    double    *,
          struct antenna_parameter                  *,
          int       *,    struct antenna_parameter                  *,
          struct comment_param      *,
          char   [][NCOMLEN],
          _Bool      ,    float     *,    int        ,    int        
        );
int     trk_priority_check
        (
          int        ,    int       *,    char [][10]
        );
int     tracking_init
        (
          int        ,    int       *,    int       *,
          char [][10],    struct antenna_parameter  *
        );
int     tracking_button_disp
        (
          int        ,    float      ,    float      , float [][4],
          int        ,    int       *,    char [][10],
          struct antenna_parameter  *
        );
void    srtatt_button_disp
        (
          int        ,     int       ,    int       *, int        ,
          float      ,     float     ,    char [][10], float [][4]
        );
double  sun_angle_check
        (
          int       *,    double     ,
          struct source_parameter   *,    struct source_parameter   *
        );
void    source_info_disp
        (
          int        ,    int        ,
          float      ,    float [][4],    float      ,
          struct char_src_info      *
        );
int     menu_config
        (
          int        ,    int        ,    int        ,    int        ,
          int        ,
          int        ,    int        ,    int        ,    int        ,
          int       *,    double     ,    _Bool     *,
          int        ,    int        ,
          struct antenna_parameter                  *,
          struct antenna_error_parameter            *,
          struct atmospheric_zenith_error           *,
          struct st_observable    *[],
          double     ,    double    *,    double    *,
          struct EOP_data            ,
          int        ,    float     *,
          struct fringe            **,    float    **,
          double     ,
          int        ,    double     ,
          struct baseline_uvw     *[],
          struct source_parameter   *,
          struct source_parameter    ,
          double    *,    double    *,
          struct srt_orbit_parameter                *,
          struct antenna_parameter                  *,
          struct srt_data_link      *,
          double     ,
          double     ,    double     ,
          double   **,    double   **,
          _Bool      ,    float     *,    int        ,    int
        );
int     on_source_disp
        ( 
          int        ,    int        ,    int        ,
          int       *,    double     ,    int        ,
          struct antenna_parameter  *,
          struct st_observable    *[],    double 
        );
int     TV_menu_hatch
        (
          float      ,    float      ,    float      ,    float      ,
          int        ,    int        
        );
int     transformation_matrices
        (
          int       *,
          double [][3],
          double [][3],
          double     ,    double     ,
          double     ,    double
        );
int     fits_data_select
        (
          char [][40],    _Bool     *,
          float      ,    float     *,    _Bool
        );
float   text_bottom
        (
          float      ,    float
        );
void    comment_init
        (
          struct comment_param      *,    char   [][NCOMLEN],
          _Bool      
        );
void    comment_disp
        (
          struct comment_param      *,    char   [][NCOMLEN],
          char      *,    _Bool      
        );
int     source_position
        (
          struct source_parameter   *,
          struct pair_src_info      *,    struct char_src_info      *,
          int        ,    int
        );
void    input_star_position
        (
          char      *,    char      *,    double    *,    double    *
        );
_Bool   output_star_position
        (
          char      *,    char      *,    double    *
        );
void    on_button
        (
          int       *,    char      *,    float     *
        );
void    off_button
        (
          int       *,    char      *,    float     *
        );
void    toggle_button
        (
          int       *,    char      *,    float     *
        );
int     button_chk
        (
          float     *,    float     *
        );
void    _on_button
        (
          _Bool     *,    char      *,    float     *
        );
void    _off_button
        (
          _Bool     *,    char      *,    float     *
        );
void    _toggle_button
        (
          _Bool     *,    char      *,    float     *
        );
_Bool   _button_chk
        (
          float     *,    float     *
        );
int     tv_get_param
        (
          char      *,    float     *,    float     *,
          float      ,    char      *,    double     ,    double
        );
void    tv_button_disp
        (
          float     *,    char      *
        );
void    str_init
        (
          char      *,    int
        );
void    char_copy
        (
          char      *,    char      *
        );
void    char_ncopy
        (
          char      *,    char      *,    int
        );
int     minor_shift_refpos
        (
          double    *,    double    *,    double    *,    double    *,
          double    *,    double    *
        );
int     ref_pos_shift
        (
          struct source_parameter   *,
          struct source_parameter   *,
          double     ,
          struct comment_param      *,
          char   [][NCOMLEN],
          _Bool
        );
int     array_config
        (
          int        ,    int        ,    int       *,    int        ,
          int       *,
          struct antenna_parameter  *,
          char      *,
          char      *,    _Bool      ,    _Bool
        );
int     wave_select
        (
          int        ,    double    *,    double    *
        );

void    Gaussian_Noise1D
        (
          int        ,    double     ,    double    *
        );
double  gauss_dev
        (
          void
        );
int     atmospheric_fluctuation
        (
          int        ,
          int        ,    int        ,    int        ,    double  *[],
          int       *,    double     ,
          struct antenna_parameter  *,
          struct phase_screen_parameter             *,
          struct source_parameter                   *,
          double     ,    double     ,    _Bool      ,
          int        ,    int        ,    int        ,    _Bool      ,
          struct comment_param      *,
          char            [][NCOMLEN],    _Bool
        );
int     turbulent_phase_screen
        (
          int        ,    int        ,    double    *,
          struct phase_screen_parameter             *,
          _Bool      ,    int        ,
          int       *,    int       *,    double    *,
          int   [][2],    int
        );
int     phase_screen_check
        (
          int        ,    struct phase_screen_parameter              ,
          float     *,    _Bool      ,    int       *
        );
int     antenna_visibility
        (
          char      *,    int        ,    int        ,    int        ,
          struct antenna_parameter  *,    char [][10],
          struct source_parameter   *,    struct source_parameter    ,
          struct phase_screen_parameter             *,
          struct phase_screen_parameter             *,
          struct phase_screen_parameter             *,
          double    *,    double    *,
          double            [][86400],
          double    *,
          float     *,
          struct srt_orbit_parameter                *,
          struct antenna_parameter  *,
          double     ,
          struct comment_param      *,    char   [][NCOMLEN],
          _Bool      ,    int       *
        );
int     GRT_TEC_TID
        (
          int        ,    int        ,    double  *[],
          struct antenna_parameter  *,
          struct phase_screen_parameter             *,
          struct st_observable    *[],
          struct TID                 
        );
int     TEC_fluctuation_amp
        (
          double [][86400],    int       *
        );
void    Gaussian_noise_check
        (
          void 
        );
void    source_model
        (
          float     *,    int        ,    float      ,    int        ,
          double    *,    float     *,    float     *,    float     *,
          char      *
        );
void    source_model3
        (
          float     *,    int        ,    float      ,    int        ,
          float     *,    float     *,
          float     *,    float     *,
          float     *,    float     *,
          char      *
        );
void    jet_cjet
        (
          float     *,    int        ,    double     ,    double     ,
          double
        );
void    PA_rotate
        (
          int        ,    float     *,    float     *,    double
        );
void    ADAF_disk
        (
          float     *,    int        ,    double     ,    double     ,
          double
        );
void    VSOP2_logo
        (
          float     *,    int        ,    double     ,    double     ,
          double
        );
void    set_color
        (
          float      ,    float     *,    float     *,    float     *,
          _Bool      ,    char      *
        );
double  random_val0
        (
          void
        );
double  random_val1
        (
          void
        );
void    fft2d
        (
          float     *,    float     *,    int        ,    int
        );
void    wfft2d
        (
          float     *,    float     *,    float     *,    int        ,
          int
        );
double  snrcal
        (
          double     ,    double     ,    double     ,    double     ,
          double     ,    double     ,    double     ,    double     ,
          double     ,    double
        );
double  thermal_noise_cal
        (
          double     ,    double     ,    double     ,    double     ,
          double     ,    double     ,    double     ,    double     ,
          double     ,    double
        );
void    parabo
        (
          float     *,    float     *,    float      ,    float      ,
          float     *,    float     *,    float     *,    float     *
        );



int     brightness_disp
        (
          int        ,    int        ,    int        ,    int        ,
          float      ,    float      ,    float      ,    float      ,
          float      ,    float      ,
          _Bool      ,
          float     *,    float     *,    float     *,    float     *,
          float   *[],    char      *,    char    *[],
          _Bool      ,    _Bool      ,    _Bool      ,    _Bool      ,
          _Bool      ,    int        ,
          char      *,
          float     *,    float     *,    float     *,
          float     *,    float     *,    float     *,    float     *
        );


void    get_vis
        (
          float     *,    float     *,    int        ,    float     *,
          int
        );
int     visibility_calc
        (
          float     *,    float     *,    int        ,
          double    *,    double    *,    double     ,
          double     ,    double    *,
          struct source_parameter    ,
          struct morphology_file_name               *,
          float      ,    float      ,    float      ,    float      ,
          _Bool      ,    int        ,    int       
        );
int     bhs_model
        (
          char      *,    float     *,    int        ,    int        ,
          double    *,
          double    *,    double     ,    double     ,    double    *
        );
int     ccm_read
        (
          char      *,    float     *,    int        ,    int        ,
          double    *,
          double    *,    double     ,    double     ,    double    *
        );
int     mcm_read
        (
          char      *,    float     *,    int        ,    double      
        );
void    vis_smoothing
        (
          double    *,    double    *,    float     *,    float     *,
          double     ,    double     ,    int        ,    int        ,
          int        ,    double
        );
void    seed_random
        (
          _Bool
        );
int     qlook_imager
        (
          int        ,    int        ,
          int        ,    int        ,    int        ,    int        ,
          int        ,    int        ,    struct baseline_uvw     *[],
          struct fringe           *[],    float   *[],
          double    *
        );
int     imager
        (
          int        ,    struct baseline_uvw       *,    double    *,
          int        ,    struct fringe             *,    float     *,
          int
        );
int     pixel_calc
        (
          double    *,    double    *,    double     ,    double     ,
          double    *,    int        ,    int
        );
int     phase_reference
        (
          int        ,    int        ,
          int        ,    int        ,    int        ,    int        ,
          struct data_number         ,
          int        ,    int        ,    double    *,
          struct fringe           *[],    float   *[]
        );
int     antenna_base_phase
        (
          int        ,    int        ,
          int        ,    int        ,    int        ,    int        ,
          struct antenna_parameter  *,
          struct data_number         ,
          int        ,    int        ,    double    *,
          struct fringe           *[],    float   *[]
        );
void    pg_color_bar
        (
          float      ,    float      ,    float      ,    float      ,
          float      ,    float      ,    _Bool      ,    char      *
        );
void    pg_color_map
        (
          int        ,    int        ,    int        ,
          float      ,    float      ,    float      ,    float      ,
          float      ,    float      ,    float      ,    float      ,
          float     *,    float      ,
          float      ,    float      ,    float      ,
          char      *,    char      *,    char      *,
          _Bool      ,    _Bool      ,    char      *
        );
void    palett
        (
          int        ,    float      ,    float
        );
void    pgcont_map
        (
          int        ,    int        ,    int        ,
          float      ,    float      ,    float      ,    float      ,
          float      ,    float      ,    float      ,    float      ,
          float      ,
          float     *,    float      ,
          float      ,    float      ,    float      ,
          char      *,    char      *,    char      *,
          _Bool      ,    _Bool      
        );
void    peak_normalize
        (
          int        ,
          float     *,    float     *,    float     *,
          float     *,    float     *,
          float     *,    float     *,
          float      ,
          float      ,
          int        ,    int        ,    int        ,    int        ,
          float     *
        );
int     fitsidi_save
        (
          char      *,    int        ,
          int        ,    int        ,    int        ,
          int        ,    int        ,    int        ,    int        ,
          struct antenna_parameter  *,
          int       *,    double     ,    double     ,
          int        ,    double     ,    double     ,
          struct EOP_data            ,
          double    *,
          struct baseline_uvw       *,
          struct fringe             *,    float     *,
          int        ,
          int        ,    int        ,    int        ,    int        ,
          struct antenna_parameter  *,
          struct source_parameter    ,
          struct srt_orbit_parameter                *
        );
int     uvw_calc
        (
          float      ,
          int        ,
          int        ,    int        ,    int        ,
          int       *,    double     ,
          struct antenna_parameter                  *,
          struct source_parameter                    ,
          struct source_parameter                    ,
          struct phase_screen_parameter             *,
          struct phase_screen_parameter             *,
          struct phase_screen_parameter             *,
          double    *,    double    *,
          double [][86400],               double     ,
          double    *,
          double    *,
          double    *,
          double    *,
          _Bool      ,    _Bool     *,    double [][3],
          double     ,
          struct st_observable      *,
          double     ,    double     ,
          struct atmospheric_zenith_error           *,
          _Bool      ,
          struct srt_orbit_parameter                *,
          double    *,
          int        ,
          struct antenna_parameter                  *,
          struct srt_data_link      *,
          double     ,                    struct TID 
        );
int     az_rotation
        (
          int        ,
          struct data_number         ,
          struct st_observable      *,
          struct antenna_parameter  *
        );
int     az_adjustment
        (
          int        ,
          struct data_number         ,
          struct st_observable    *[],
          struct antenna_parameter  *
        );
void    weight_flag_add
        (
          float      ,    float     *,    float
        );
_Bool   weight_flag_chk
        (
          float      ,    float     *,    float
        );
void    attitude_Q
        (
          int        ,
          struct source_parameter    ,
          struct source_parameter    ,
          double    *
        );
int     tracking_condition
        (
          int       *,    double     ,
          struct antenna_parameter   ,
          double     ,    double     ,    _Bool      ,
          double    *,    double    *,
          double    *,    double    *,
          double    *,    double    *,
          struct srt_data_link      *,
          struct source_parameter    ,
          struct source_parameter    ,
          double    *,    double    *,    double    *,    double    *,
          double    *,    double    *,
          int
        );
_Bool   tracking_status
        (
          _Bool      ,    float     *,    int       *,
          int       *,    double     ,
          struct srt_orbit_parameter ,
          struct source_parameter    ,
          struct source_parameter    ,
          double     ,    double     ,
          int        ,
          double    *,    double    *,    double    *,
          struct antenna_parameter  *,
          struct srt_data_link      *,
          double    *,    double    *,    double    *,    double    *,
          double    *,    double    *,
          double     ,
          double    *,    double    *,    double    *
        );
int     tracking_station_name_read(
          int        ,    char      *,    char   [][10],  int       *
        );

int     get_srt_link
        (
          struct srt_data_link      *,    double    *
        );
int     allocate_uvcalc_param
        (
          int        ,    double  *[]
        );
int     uv_display
        (
          int        ,    double    *,    int        ,
          int        ,    int        ,    int        ,    int        ,
          int       *,    double     ,
          struct antenna_parameter  *,
          struct source_parameter   *,
          float     *,
          struct st_observable    *[],
          double    *,    int        ,    _Bool      ,    float      ,
          char      *,    int
        );
int     ant_pos_err_disp
        (
          int        ,
          struct antenna_parameter  *,
          struct antenna_error_parameter            *,
          struct atmospheric_zenith_error           *
        );
double  airmass
        (
          double
        );
int     nmf20_coeffi
        (
          double     ,    int        ,
          double    *,    double    *,    double    *
        );
int     nmf20
        (
          double    *,    double    *,
          double    *,    double    *,    double    *,
          double     ,    double     
        );
void    spherical_geometry
        (
          double     ,    double     ,    double    *,
          double    *,    double    *,    double    *,    double    *
        );
int     station_select
        (
          int        ,
          int        ,    int        ,    int        ,    int        ,
          int       *,    int       *,
          struct antenna_parameter  *,    float      ,
          float     *,    _Bool      ,    int        
        );
void    EPL_disp
        (
          int        ,    int        ,    int        ,
          struct antenna_parameter  *,    double    *,
          int       *,    struct st_observable    *[]
        );
int     spectl_disp
        (
          int        ,    int        ,    int        ,
          struct antenna_parameter  *,    double    *,
          int       *,    struct st_observable    *[]
        );
int     allanv_disp
        (
          int        ,    int        ,    int        ,
          struct antenna_parameter  *,    double    *,
          int       *,    struct st_observable    *[]
        );
int     fringe_disp
        (
          int        ,    int        ,
          int        ,    int        ,
          struct antenna_parameter  *,
          int        ,    int        ,
          int       *,    double     ,
          struct fringe           *[],    float   *[],
          float     *,    float     *
        );
int     coherence_disp
        (
          int        ,    int        ,
          int        ,    int        ,
          struct antenna_parameter  *,
          int        ,    int        ,
          int       *,    double     ,
          struct fringe           *[],    float   *[],
          float     *
        );
void    position2uvw
        (
          double    *,    double    *,    double    *,    double    *,
          double    *,    double    *,    double    *,
          double     ,    double     ,    double    *
        );
void    spacecraft_position
        (
          struct srt_orbit_parameter ,
          int       *,    double     ,    double     ,
          double    *,    double    *,
          double    *,    double    *,    double    *
        );
void    orbit_propagation
        (
          struct srt_orbit_parameter ,
          int       *,    double     ,    double     ,
          double    *,    double    *,
          double    *,    double    *,    double    *,    int        ,
          double    *,    double    *,    double    *
        );
void    RK04_init
        (
          double    *,    double    *,    double    *
        );
void    RK07_init
        (
          double    *,    double    *,    double    *
        );
void    RK_end
        (
          double    *,    double    *,    double    *
        );
double  f1
        (
          double     ,    double     ,    double     ,    double
        );
double  f2
        (
          double     ,    double     ,    double     ,    double
        );
double  Force_Model
        (
          double    *,    int       *,    double
        );
int     slr_config
        (
          struct antenna_parameter  *
        );
void    intput_star_position
        (
          double    *,    double    *
        );
void    azel_disp
        (
          int        ,    int        ,
          struct antenna_parameter  *,
          struct st_observable    *[],
          int       *,    double     ,    char      *
        );
int     SSF_disp
        (
          int       ,    int         ,
          int       ,    int         ,    int       ,    int        ,
          int       ,    int         ,
          int       *,   double      ,    double   *,
          struct antenna_parameter  *,
          struct fringe           *[],
          float   *[],
          struct st_observable    *[],
          char       *,  int
        );
void    orbit_disp
        (
          int        ,    int        ,
          struct srt_orbit_parameter                *,
          int       *,    double     ,    double     ,
          struct srt_data_link      *,
          double     ,    double     ,
          int        ,    int        ,
          struct st_observable    *[],
          struct antenna_parameter  *,
          struct source_parameter    ,
          struct source_parameter    ,

          char      *
        );
void    draw_earth
        (
          double     ,    double     ,    double     ,
          struct source_parameter    ,
          struct source_parameter    ,
          struct earth_shape        *,
          int       *,
          _Bool      ,    _Bool      ,    _Bool
        );
void    data_link_schedule
        (
          int        ,    int        ,
          struct srt_orbit_parameter                *,
          struct antenna_parameter                  *,
          struct source_parameter    ,
          struct source_parameter    ,
          struct srt_data_link      *,
          int       *,    double     ,    double     ,
          double     ,    double     ,
          int        ,    int        
        );
void    link_TRK_sim
        (
          int        ,    int        ,
          struct srt_orbit_parameter                *,
          struct antenna_parameter                  *,
          struct source_parameter    ,
          struct source_parameter    ,
          struct srt_data_link      *,
          int       *,    double     ,    double     ,
          double     ,    double     ,
          int        ,    int        
        );
void    link_SLR_sim
        (
          int        ,    int        ,
          struct srt_orbit_parameter                *,
          struct antenna_parameter                  *,
          struct source_parameter    ,
          struct source_parameter    ,
          struct srt_data_link      *,
          int       *,    double     ,    double     ,
          double     ,    double     ,
          int        ,    int        
        );
void    orbit_info_print
        (
          FILE      *,    int        ,
          struct srt_orbit_parameter *,
          int       *,    double
        );
void    delay_rate_disp
        (
          int        ,    int        ,
          struct antenna_parameter  *,
          int        ,
          int       *,
          struct st_observable    *[],    double  *[]
        );
float   ordinate_axis
        (
          int        ,    float     *,    float     *,
          float     *,    char      *
        );
void    position_on_screen
        (
          double     ,    double    *,
          struct phase_screen_parameter              ,
          double     ,    double     ,    double    *,    double    *
        );
int     earth_eclipse
        (
          double    *,    double    *,    double    *,    double    *,
          double
        );
double  diff
        (
          double     ,    double
        );
double  slew_time
        (
          double     ,    double     ,    double     ,
          double     ,    double     ,
          double     ,    double     ,    double     ,
          double     ,    double     
        );
double  SEFD
        (
          double     ,    double     ,    double     ,    double    
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
int     freq_std_instability
        (
          int        ,    double    *,    int
        );
double  asd_model
        (
          double     ,    int        
        );
double  modified_asd_model
        (
          double     ,    int        
        );
double  coherence_factor_calc
        (
          double     ,    double     ,    int        ,    int 
        );
int     fit_power_two_number
        (
          double     ,    double     ,    double    *,    int       *
        );


/****
#define __RANDOM_SEED__
****/



#define NSRCPROC    10
#define SRCPROC1     0
#define SRCPROC2     1
#define SRCPROC3     2
#define SRCPROC4     3
#define SRCPROC5     4
#define SRCPROC6     5
#define SRC__RA__DEC 0   /* RA         -          Dec */
#define SRC_dRA_dDEC 1   /* Delta RA   -    Delta Dec */
#define SRC_SEP_POSA 2   /* Separation - Position ang */
#define SRC_POS1     0   /* Reference : Source-1      */
#define SRC_POS2     1   /* Reference : Mid-Point     */

#define NO_STRUCTURE  0   /* UV_display No matter               */
#define F__STRUCTURE  1   /* UV_display Fine Structure Option   */
#define C__STRUCTURE  2   /* UV_display Coarse Structure Option */

/****
static float  pgpap_prm      = 5.75;
static float  pgpap_prm      = 3.00;
static float  pgpap_prm      = 4.00;
static float  pgpap_prm      = 5.00;
static float  pgpap_prm      = 5.50;
static float  pgpap_prm      = 5.75;
static float  pgpap_prm      = 6.25;
static float  pgpap_prm      = 5.75;
****/
static float  pgpap_prm      = 5.75;

static double dpi            = 3.141592653589793238462643;
static double speed_of_light = 2.99792458e8;
static double OMEGA          = 7.29211585791599e-5;
/**  OMEGA = 2PI / (23*3600 + 56*60 + 4.0905)  **/
static double BOLTZ          = 1.3806503e-23;
static double TEC_CONST      = 40.308;
static double earth_radius   = 6.378137e6;
static double GM             = 398600.4415e9;
/********
#define   G   6.67259e-11
#define   M   5.974e24
********/

