#include "astrotools.h"

int   baseline2antenna_number(int  baseline_num,  int  ant_num,
                              int  *ant_X,        int  *ant_Y)
{
  int    ibase;

/*
  This tool calculate the antenna number starting from 0.
  
  baseline_num  : the number of the order of a baseline
  ant_num       : the total number of antennas
  ant_X         : (RETURN) the number of the antenna X
  ant_Y         : (RETURN) the number of the antenna Y

  In the case of ant_num=6 and baseline_num=9,
    *ant_X = 2 (C)
    *ant_Y = 3 (D)
*/

  *ant_X = 0;
  *ant_Y = 0;

  ant_num--;
  while (1) {
    if ((ibase = baseline_num - ant_num) < 0) {
      break;
    } else {
      baseline_num = ibase;
      ant_num--;
      (*ant_X)++;
    }
  }
  *ant_Y = *ant_X + 1 + baseline_num;

  return 1;
}
