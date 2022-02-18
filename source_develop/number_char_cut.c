#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <aris.h>


int  number_char_cut(char *ch_num)
{
  int     i, j, k;
  char    *c_tmp;

/*
-----------------------------------------------------
*/

  c_tmp = (char *)calloc(strlen(ch_num) + 1, 1);
  for (i=0; i<strlen(ch_num); i++) {
    if (ch_num[i] != ' ') {
      break;
    }
  }
  if (i != 0) {
    char_copy(c_tmp, ch_num+i);
    char_copy(ch_num, c_tmp);
  }
  free (c_tmp);

  for (i=0; i<strlen(ch_num); i++) {
    if (ch_num[i] == ' ') {
      ch_num[i] = '\0';
      break;
    }
  }

  for (i=0; i<strlen(ch_num); i++) {
    if (ch_num[i] == '.') {
      break;
    }
  }
  k = i;
  j = strlen(ch_num) - 1;
  for (i=j; i>=k+2; i--) {
    if (ch_num[i] != '0') {
      break;
    }
    ch_num[i] = '\0';
  }
  if (ch_num[k+2] == '\0' && ch_num[k+1] == '0') {
    ch_num[k] = '\0';
  }

  return 1;
}
