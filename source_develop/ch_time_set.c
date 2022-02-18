#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


int  ch_time_set(int  *TimUTC, struct char_obs_time *ch_obs_t)
{
  int    i;

/*
-------------------------------------------
*/

  sprintf(ch_obs_t->start_t[0], "%4d", TimUTC[0]);
  for (i=0; i<4; i++) {
    if (ch_obs_t->start_t[0][i] == ' ') {
      ch_obs_t->start_t[0][i] = '0';
    }
  }
  sprintf(ch_obs_t->start_t[1], "%2d", TimUTC[1]);
  for (i=0; i<2; i++) {
    if (ch_obs_t->start_t[1][i] == ' ') {
      ch_obs_t->start_t[1][i] = '0';
    }
  }
  sprintf(ch_obs_t->start_t[2], "%2d", TimUTC[2]);
  for (i=0; i<2; i++) {
    if (ch_obs_t->start_t[2][i] == ' ') {
      ch_obs_t->start_t[2][i] = '0';
    }
  }
  sprintf(ch_obs_t->start_t[3], "%2d", TimUTC[3]);
  for (i=0; i<2; i++) {
    if (ch_obs_t->start_t[3][i] == ' ') {
      ch_obs_t->start_t[3][i] = '0';
    }
  }
  sprintf(ch_obs_t->start_t[4], "%2d", TimUTC[4]);
  for (i=0; i<2; i++) {
    if (ch_obs_t->start_t[4][i] == ' ') {
      ch_obs_t->start_t[4][i] = '0';
    }
  }
  sprintf(ch_obs_t->start_t[5], "%2d", TimUTC[5]);
  for (i=0; i<2; i++) {
    if (ch_obs_t->start_t[5][i] == ' ') {
      ch_obs_t->start_t[5][i] = '0';
    }
  }

  return 1;
}
