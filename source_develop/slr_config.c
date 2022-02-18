#include <stdio.h>
#include <math.h>
#include <aris.h>

/****
#define __DEBUG__
****/

int  slr_config(struct antenna_parameter *slr_pos)
{
  int    SLR_NUM;
  int    islr;
  double DPI;

  DPI = dpi / 180.0;

/****
#  station_name station_id epoch          x             y             z     
 0 Tanegashima     7358 1997/01/01 -3607651.5503  4147873.9561  3223722.8677
 1 Koganei         7308 1997/01/01 -3942020.1228  3368097.5570  3702191.0797
 2 Yarragadee      7090 1997/01/01 -2389006.9238  5043329.3395 -3078524.8896
 3 Mount_Stromlo_2 7825 1997/01/01 -4467064.3092  2683034.8862 -3667007.9248
 4 Zimmerwald      7810 1997/01/01  4331283.6759   567549.7428  4633140.2669
 5 Graz            7839 1997/01/01  4194426.5167  1162694.0327  4647246.6497
 6 Wettzell        8834 1997/01/01  4075576.8504   931785.4556  4801583.5588
 7 San_Juan        7406 1997/01/01  1984104.3220 -5068867.1983 -3314482.5237
 8 Riyadh          7832 1997/01/01  3992100.9957  4192172.5015  2670410.6889
 9 Monument_Peak   7110 1997/01/01 -2386278.2108 -4802354.1447  3444881.5980
10 Herstmonceux    7840 1997/01/01  4033463.7122    23662.4784  4924305.1712
11 Changchun       7237 1997/01/01 -2674386.7386  3757189.3075  4391508.4012
12 Matera_MLRO     7941 1997/01/01  4641978.8587  1393067.4785  4133249.4277
13 Hartebeesthoek  7501 1997/01/01  5085401.1179  2668330.0746 -2768688.8727
14 Potsdam_3       7841 1997/01/01  3800432.2932   881691.9670  5029030.0329
15 Simosato        7838 1997/01/01 -3822388.3468  3699363.5675  3507573.1060
16 Greenbelt       7105 1997/01/01  1130719.6324 -4831350.5774  3994106.5389
17 San_Fernando    7824 1997/01/01  5105473.9094  -555110.7233  3769892.7442
18 McDonald        7080 1997/01/01 -1330021.0667 -5328401.8557  3236480.7817
19 Beijing         7249 1997/01/01 -2148760.3719  4426759.6707  4044509.7266
20 Shanghai_2      7821 1997/01/01 -2830744.1715  4676580.3937  3275072.9259
21 Maidanak_1      1864 1997/01/01  1953285.7944  4588973.8346  3966767.6833
****/

  for (islr=0; islr<ANTMAX; islr++) {
    slr_pos[islr].UFL = false;
  }

  sprintf(slr_pos[ 0].IDC, "TANEGASH");
  slr_pos[ 0].ARRAY  = SLR_NETWORK;
  slr_pos[ 0].MNTSTA = ALAZ;
  slr_pos[ 0].UFL    = true;
  slr_pos[ 0].XYZ[0] = -3607651.5503;
  slr_pos[ 0].XYZ[1] =  4147873.9561;
  slr_pos[ 0].XYZ[2] =  3223722.8677;

  sprintf(slr_pos[ 1].IDC, "KOGANEI ");
  slr_pos[ 1].ARRAY  = SLR_NETWORK;
  slr_pos[ 1].MNTSTA = ALAZ;
  slr_pos[ 1].UFL    = true;
  slr_pos[ 1].XYZ[0] = -3942020.1228;
  slr_pos[ 1].XYZ[1] =  3368097.5570;
  slr_pos[ 1].XYZ[2] =  3702191.0797;

  sprintf(slr_pos[ 2].IDC, "YARRAGAD");
  slr_pos[ 2].ARRAY  = SLR_NETWORK;
  slr_pos[ 2].MNTSTA = ALAZ;
  slr_pos[ 2].UFL    = true;
  slr_pos[ 2].XYZ[0] = -2389006.9238;
  slr_pos[ 2].XYZ[1] =  5043329.3395;
  slr_pos[ 2].XYZ[2] = -3078524.8896;

  sprintf(slr_pos[ 3].IDC, "MOUNT_ST");
  slr_pos[ 3].ARRAY  = SLR_NETWORK;
  slr_pos[ 3].MNTSTA = ALAZ;
  slr_pos[ 3].UFL    = true;
  slr_pos[ 3].XYZ[0] = -4467064.3092;
  slr_pos[ 3].XYZ[1] =  2683034.8862;
  slr_pos[ 3].XYZ[2] = -3667007.9248;

  sprintf(slr_pos[ 4].IDC, "ZIMMERWA");
  slr_pos[ 4].ARRAY  = SLR_NETWORK;
  slr_pos[ 4].MNTSTA = ALAZ;
  slr_pos[ 4].UFL    = true;
  slr_pos[ 4].XYZ[0] =  4331283.6759;
  slr_pos[ 4].XYZ[1] =   567549.7428;
  slr_pos[ 4].XYZ[2] =  4633140.2669;

  sprintf(slr_pos[ 5].IDC, "GRAZ    ");
  slr_pos[ 5].ARRAY  = SLR_NETWORK;
  slr_pos[ 5].MNTSTA = ALAZ;
  slr_pos[ 5].UFL    = true;
  slr_pos[ 5].XYZ[0] =  4194426.5167;
  slr_pos[ 5].XYZ[1] =  1162694.0327;
  slr_pos[ 5].XYZ[2] =  4647246.6497;

  sprintf(slr_pos[ 6].IDC, "WETTZELL");
  slr_pos[ 6].ARRAY  = SLR_NETWORK;
  slr_pos[ 6].MNTSTA = ALAZ;
  slr_pos[ 6].UFL    = true;
  slr_pos[ 6].XYZ[0] =  4075576.8504;
  slr_pos[ 6].XYZ[1] =   931785.4556;
  slr_pos[ 6].XYZ[2] =  4801583.5588;

  sprintf(slr_pos[ 7].IDC, "SAN_JUAN");
  slr_pos[ 7].ARRAY  = SLR_NETWORK;
  slr_pos[ 7].MNTSTA = ALAZ;
  slr_pos[ 7].UFL    = true;
  slr_pos[ 7].XYZ[0] =  1984104.3220;
  slr_pos[ 7].XYZ[1] = -5068867.1983;
  slr_pos[ 7].XYZ[2] = -3314482.5237;

  sprintf(slr_pos[ 8].IDC, "RIYADH  ");
  slr_pos[ 8].ARRAY  = SLR_NETWORK;
  slr_pos[ 8].MNTSTA = ALAZ;
  slr_pos[ 8].UFL    = true;
  slr_pos[ 8].XYZ[0] =  3992100.9957;
  slr_pos[ 8].XYZ[1] =  4192172.5015;
  slr_pos[ 8].XYZ[2] =  2670410.6889;

  sprintf(slr_pos[ 9].IDC, "MT_PEAK ");
  slr_pos[ 9].ARRAY  = SLR_NETWORK;
  slr_pos[ 9].MNTSTA = ALAZ;
  slr_pos[ 9].UFL    = true;
  slr_pos[ 9].XYZ[0] = -2386278.2108;
  slr_pos[ 9].XYZ[1] = -4802354.1447;
  slr_pos[ 9].XYZ[2] =  3444881.5980;

  sprintf(slr_pos[10].IDC, "HERSTMON");
  slr_pos[10].ARRAY  = SLR_NETWORK;
  slr_pos[10].MNTSTA = ALAZ;
  slr_pos[10].UFL    = true;
  slr_pos[10].XYZ[0] =  4033463.7122;
  slr_pos[10].XYZ[1] =    23662.4784;
  slr_pos[10].XYZ[2] =  4924305.1712;

  sprintf(slr_pos[11].IDC, "CHANGCHU");
  slr_pos[11].ARRAY  = SLR_NETWORK;
  slr_pos[11].MNTSTA = ALAZ;
  slr_pos[11].UFL    = true;
  slr_pos[11].XYZ[0] = -2674386.7386;
  slr_pos[11].XYZ[1] =  3757189.3075;
  slr_pos[11].XYZ[2] =  4391508.4012;

  sprintf(slr_pos[12].IDC, "MATERA_M");
  slr_pos[12].ARRAY  = SLR_NETWORK;
  slr_pos[12].MNTSTA = ALAZ;
  slr_pos[12].UFL    = true;
  slr_pos[12].XYZ[0] =  4641978.8587;
  slr_pos[12].XYZ[1] =  1393067.4785;
  slr_pos[12].XYZ[2] =  4133249.4277;

  sprintf(slr_pos[13].IDC, "HARTEBEE");
  slr_pos[13].ARRAY  = SLR_NETWORK;
  slr_pos[13].MNTSTA = ALAZ;
  slr_pos[13].UFL    = true;
  slr_pos[13].XYZ[0] =  5085401.1179;
  slr_pos[13].XYZ[1] =  2668330.0746;
  slr_pos[13].XYZ[2] = -2768688.8727;

  sprintf(slr_pos[14].IDC, "POTSDAM ");
  slr_pos[14].ARRAY  = SLR_NETWORK;
  slr_pos[14].MNTSTA = ALAZ;
  slr_pos[14].UFL    = true;
  slr_pos[14].XYZ[0] =  3800432.2932;
  slr_pos[14].XYZ[1] =   881691.9670;
  slr_pos[14].XYZ[2] =  5029030.0329;

  sprintf(slr_pos[15].IDC, "SIMOSATO");
  slr_pos[15].ARRAY  = SLR_NETWORK;
  slr_pos[15].MNTSTA = ALAZ;
  slr_pos[15].UFL    = true;
  slr_pos[15].XYZ[0] = -3822388.3468;
  slr_pos[15].XYZ[1] =  3699363.5675;
  slr_pos[15].XYZ[2] =  3507573.1060;

  sprintf(slr_pos[16].IDC, "GREENBLT");
  slr_pos[16].ARRAY  = SLR_NETWORK;
  slr_pos[16].MNTSTA = ALAZ;
  slr_pos[16].UFL    = true;
  slr_pos[16].XYZ[0] =  1130719.6324;
  slr_pos[16].XYZ[1] = -4831350.5774;
  slr_pos[16].XYZ[2] =  3994106.5389;

  sprintf(slr_pos[17].IDC, "SAN_FERN");
  slr_pos[17].ARRAY  = SLR_NETWORK;
  slr_pos[17].MNTSTA = ALAZ;
  slr_pos[17].UFL    = true;
  slr_pos[17].XYZ[0] =  5105473.9094;
  slr_pos[17].XYZ[1] =  -555110.7233;
  slr_pos[17].XYZ[2] =  3769892.7442;

  sprintf(slr_pos[18].IDC, "MCDONALD");
  slr_pos[18].ARRAY  = SLR_NETWORK;
  slr_pos[18].MNTSTA = ALAZ;
  slr_pos[18].UFL    = true;
  slr_pos[18].XYZ[0] = -1330021.0667;
  slr_pos[18].XYZ[1] = -5328401.8557;
  slr_pos[18].XYZ[2] =  3236480.7817;

  sprintf(slr_pos[19].IDC, "BEIJING ");
  slr_pos[19].ARRAY  = SLR_NETWORK;
  slr_pos[19].MNTSTA = ALAZ;
  slr_pos[19].UFL    = true;
  slr_pos[19].XYZ[0] = -2148760.3719;
  slr_pos[19].XYZ[1] =  4426759.6707;
  slr_pos[19].XYZ[2] =  4044509.7266;

  sprintf(slr_pos[20].IDC, "SHANGHAI");
  slr_pos[20].ARRAY  = SLR_NETWORK;
  slr_pos[20].MNTSTA = ALAZ;
  slr_pos[20].UFL    = true;
  slr_pos[20].XYZ[0] = -2830744.1715;
  slr_pos[20].XYZ[1] =  4676580.3937;
  slr_pos[20].XYZ[2] =  3275072.9259;

  sprintf(slr_pos[21].IDC, "MAIDANAK");
  slr_pos[21].ARRAY  = SLR_NETWORK;
  slr_pos[21].MNTSTA = ALAZ;
  slr_pos[21].UFL    = true;
  slr_pos[21].XYZ[0] =  1953285.7944;
  slr_pos[21].XYZ[1] =  4588973.8346;
  slr_pos[21].XYZ[2] =  3966767.6833;

  SLR_NUM = 0;
  for (islr=0; islr<ANTMAX; islr++) {
    if (slr_pos[islr].UFL == true) {
      slr_pos[SLR_NUM] = slr_pos[islr];
      SLR_NUM++;
    }
  }

  for (islr=0; islr<SLR_NUM; islr++) {
    J_system_geocentric_equatorial_rectangular_coordinate2llh(
        &slr_pos[islr].LLH[0],
        &slr_pos[islr].LLH[1],
        &slr_pos[islr].LLH[2], slr_pos[islr].XYZ);
#ifdef __DEBUG__
    printf("__DEBUG__  %lf  %lf  %lf\n",
        slr_pos[islr].LLH[0],
        slr_pos[islr].LLH[1],
        slr_pos[islr].LLH[2]);
#endif /* __DEBUG__ */
  }

  return SLR_NUM;
}
