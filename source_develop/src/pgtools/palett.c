#include <stdio.h>
#include <cpgplot.h>
#include "./pgtools.h"

void  palett(int MODE_SWT, float contra, float bright)
{
  float  BL[2] = { 0.0,  1.0};
  float  BR[2] = { 1.0,  0.0};
  float  BG[2] = { 1.0,  0.0};
  float  BB[2] = { 1.0,  0.0};

  float  GL[2] = { 0.0,  1.0};
  float  GR[2] = { 0.0,  1.0};
  float  GG[2] = { 0.0,  1.0};
  float  GB[2] = { 0.0,  1.0};

  float  RL[9] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
  float  RR[9] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
  float  RG[9] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
  float  RB[9] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};

  float  OHL[5] = { 0.0, 0.2, 0.4, 0.6, 1.0};
  float  OHR[5] = { 0.0, 0.5, 1.0, 1.0, 1.0};
  float  OHG[5] = { 0.0, 0.0, 0.5, 1.0, 1.0};
  float  OHB[5] = { 0.0, 0.0, 0.0, 0.3, 1.0};

  float  RHL[5] = { 0.0, 0.2, 0.4, 0.6, 1.0};
  float  RHR[5] = { 1.0, 1.0, 1.0, 0.5, 0.0};
  float  RHG[5] = { 1.0, 1.0, 0.5, 0.0, 0.0};
  float  RHB[5] = { 1.0, 0.3, 0.0, 0.0, 0.0};

  float  WL[10] = { 0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
  float  WR[10] = { 0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
  float  WG[10] = { 0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
  float  WB[10] = { 0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};

  float  OAL[20] = { 0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
                     0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
  float  OAR[20] = { 0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  float  OAG[20] = { 0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
                     0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
  float  OAB[20] = { 0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  float  RAL[20] = { 0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
                     0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
  float  RAR[20] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.3, 0.3, 0.0, 0.0};
  float  RAG[20] = { 0.0, 0.0, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 0.6, 0.6,
                     0.8, 0.8, 0.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.0, 0.0};
  float  RAB[20] = { 0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
                     0.9, 0.9, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.0, 0.0};

  int    I;
  if        (MODE_SWT == 1) {
/* -- gray scale (B->W) */
    I =  2;
    cpgctab(GL, GR, GG, GB, I, contra, bright);

  } else if (MODE_SWT == 2) {
/* -- rainbow */
    I =  9;
    cpgctab(RL, RR, RG, RB, I, contra, bright);

  } else if (MODE_SWT == 3) {
/* -- heat */
    I =  5;
    cpgctab(OHL, OHR, OHG, OHB, I, contra, bright);

  } else if (MODE_SWT == -3) {
/* -- heat */
    I =  5;
    cpgctab(RHL, RHR, RHG, RHB, I, contra, bright);

  } else if (MODE_SWT == 4) {
/* -- weired IRAF */
    I = 10;
    cpgctab(WL, WR, WG, WB, I, contra, bright);

  } else if (MODE_SWT == 5) {
/* -- AIPS */
    I = 20;
    cpgctab(OAL, OAR, OAG, OAB, I, contra, bright);

  } else if (MODE_SWT == -5) {
/* -- AIPS */
    I = 20;
    cpgctab(OAL, RAR, RAG, RAB, I, contra, bright);

  } else {
/* -- gray scale (W->B) */
    I =  2;
    cpgctab(BL, BR, BG, BB, I, contra, bright);

  }

  return;
}
