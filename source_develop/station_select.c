#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <mathtools.h>
#include <astrotools.h>
#include <aris.h>


int  station_select(int ANT_NUM,
                    int BGN_ANT_I, int END_ANT_I,
                    int BGN_ANT_J, int END_ANT_J,
                    int *NX, int *NY,
                    struct antenna_parameter *ant_prm,
                    float sel_pos, float *cursor_pos,
                    _Bool TV_SWT,  int ST_NUM)
{
  int    i, j, I, idum, END_ANT;
  int    PR_SWT, CNTL_BTN_NUM;
  char   string[20];
  float  pitch = 0.03;
  float  **bttn_box;
  char   **bttn_lab;

/*
---- ANT_NUM+1  : `+1' is for QUIT button. ----------
*/

  if (TV_SWT == false) {
    if (ST_NUM == 1) {
      while (1) {
        *NX = 0;
        printf("Input station number.\n");
        printf("(Just numbers. CR -> quit) : ");
        if (fgets(string, sizeof(string), stdin) == NULL) {
          printf("ERROR: STATION_SELECT: invalid input.\n");
          return (__NG__);
        }
        if (string[0] == '\n') {
          return (_EXIT_);
        }

        sscanf(string, "%d", NX);
        (*NX)--;

        if (*NX >= BGN_ANT_I && *NX < END_ANT_I &&
            *NX >= BGN_ANT_J && *NX < END_ANT_J) {
          return (__GO__);
        } else {
          printf("Selected number is wrong.");
        }
      }
    } else if (ST_NUM >= 2) {
      while (1) {
        *NX = 0; *NY = 0;
        printf("Input X- and Y-station number.\n");
        printf("(Just numbers. CR -> quit) : ");
        if (fgets(string, sizeof(string), stdin) == NULL) {
          printf("ERROR: STATION_SELECT: invalid input.\n");
          return (__NG__);
        }
        if (string[0] == '\n') {
          return (_EXIT_);
        }

        for (i=0; i<strlen(string); i++)
          if (string[i] == ',')
            string[i] = ' ';
        sscanf(string, "%d %d", NX, NY);
        (*NX)--;
        (*NY)--;

        if (*NX > *NY) {
          i = *NX;
          *NX = *NY;
          *NY = i;
        }

        if (*NX >= BGN_ANT_I && *NX < END_ANT_I &&
            *NY >= BGN_ANT_J && *NY < END_ANT_J) {
          if (ST_NUM == 3) {
            printf("1. No P-R   2. With P-R   : ");
            if (fgets(string, sizeof(string), stdin) == NULL) {
              printf("ERROR: STATION_SELECT: invalid input.\n");
              return (__NG__);
            }
            if (string[0] == '1') {
              return ( 2);
            } else if (string[0] == '2') {
              return ( 3);
            }
          } else {
            return (__GO__);
          }
        } else {
          printf("Selected numbers are wrong combination.");
        }
      }
    }

/*
-----------------------------------------------------
*/

  } else if (TV_SWT == true) {
    if (ST_NUM == 3) {
      CNTL_BTN_NUM = 3;
    } else {
      CNTL_BTN_NUM = 1;
    }

    if ((bttn_box = (float **)malloc((ANT_NUM+CNTL_BTN_NUM) * sizeof(float *)))
                                                              == NULL) {
      printf("ERROR: BASELINE_SELECT: memory alloc for **bttn_box.\n");
      return (__NG__);
    }
    for (i=0; i<ANT_NUM+CNTL_BTN_NUM; i++) {
      if ((bttn_box[i] = (float *)calloc(4, sizeof(float))) == NULL) {
        printf("ERROR: BASELINE_SELECT: memory alloc for *bttn_box.\n");
        free_memory_block_float(bttn_box, i);
        return (__NG__);
      }
    }

    if ((bttn_lab = (char **)malloc((ANT_NUM+CNTL_BTN_NUM) * sizeof(char *)))
                                                              == NULL) {
      printf("ERROR: BASELINE_SELECT: memory alloc for **bttn_lab.\n");
      free_memory_block_float(bttn_box, ANT_NUM+CNTL_BTN_NUM);
      return (__NG__);
    }
    for (i=0; i<ANT_NUM+CNTL_BTN_NUM; i++) {
      if ((bttn_lab[i] = (char *)calloc(10, sizeof(char))) == NULL) {
        printf("ERROR: BASELINE_SELECT: memory alloc for *bttn_lab.\n");
        free_memory_block_char(bttn_lab, i);
        free_memory_block_float(bttn_box, ANT_NUM+CNTL_BTN_NUM);
        return (__NG__);
      }
    }

    if (END_ANT_I < END_ANT_J) {
      END_ANT = END_ANT_J;
    } else {
      END_ANT = END_ANT_I;
    }

    i = 0;
    bttn_box[i][0] = 0.00 + 0.10 * (float)(i % 10);
    bttn_box[i][1] = 0.10 + 0.10 * (float)(i % 10);
    bttn_box[i][2] = sel_pos - pitch * (float)(i / 10);
    bttn_box[i][3] = bttn_box[i][2] + pitch;
    sprintf(bttn_lab[i], "QUIT");
    off_button(&idum, bttn_lab[i], bttn_box[i]);
    for (i=1; i<=END_ANT; i++) {
      bttn_box[i][0] = 0.00 + 0.10 * (float)(i % 10);
      bttn_box[i][1] = 0.10 + 0.10 * (float)(i % 10);
      bttn_box[i][2] = sel_pos - pitch * (float)(i / 10);
      bttn_box[i][3] = bttn_box[i][2] + pitch;
      sprintf(bttn_lab[i], "%s", ant_prm[i-1].IDC);
    }
    for (i=1; i<=END_ANT; i++) {
      off_button(&idum, bttn_lab[i], bttn_box[i]);
    }
    if (ST_NUM == 3) {
      bttn_box[i][0] = 0.00 + 0.10 * (float)(i % 10);
      bttn_box[i][1] = 0.10 + 0.10 * (float)(i % 10);
      bttn_box[i][2] = sel_pos - pitch * (float)(i / 10);
      bttn_box[i][3] = bttn_box[i][2] + pitch;
      sprintf(bttn_lab[i], "No P-R");
      off_button(&idum, bttn_lab[i], bttn_box[i]);
      i++;
      bttn_box[i][0] = 0.00 + 0.10 * (float)(i % 10);
      bttn_box[i][1] = 0.10 + 0.10 * (float)(i % 10);
      bttn_box[i][2] = sel_pos - pitch * (float)(i / 10);
      bttn_box[i][3] = bttn_box[i][2] + pitch;
      sprintf(bttn_lab[i], "With P-R");
      off_button(&idum, bttn_lab[i], bttn_box[i]);
    }

    I = 0;
    *NX = -1;
    *NY = -1;
    PR_SWT = 1;
    while (1) {
      cpgcurs(cursor_pos, cursor_pos+1, string);
      if (_button_chk(cursor_pos, bttn_box[0]) == true) {
        free_memory_block_float(bttn_box, ANT_NUM+CNTL_BTN_NUM);
        free_memory_block_char(bttn_lab,  ANT_NUM+CNTL_BTN_NUM);
        return (_EXIT_);
      } else {
        for (i=1; i<END_ANT+CNTL_BTN_NUM; i++) {
          if (_button_chk(cursor_pos, bttn_box[i]) == true) {
            if (ST_NUM == 1) {
              *NX = i - 1;
              free_memory_block_float(bttn_box, ANT_NUM+CNTL_BTN_NUM);
              free_memory_block_char(bttn_lab,  ANT_NUM+CNTL_BTN_NUM);
              return (__GO__);
            } else if (ST_NUM >= 2) {
              if (ST_NUM == 3) {
                if (i == END_ANT + 1) {
                  on_button(&idum, bttn_lab[i], bttn_box[i]);
                  PR_SWT = 2;
                  I++;
                } else if (i == END_ANT + 2) {
                  on_button(&idum, bttn_lab[i], bttn_box[i]);
                  PR_SWT = 3;
                  I++;
                }
              }

              if (*NX == -1 && i < ANT_NUM + 1) {
                on_button(&idum, bttn_lab[i], bttn_box[i]);
                *NX = i - 1;
                I++;
              }

              if (*NX != -1 && *NY == -1 && i - 1 != *NX && i < ANT_NUM + 1) {
                on_button(&idum, bttn_lab[i], bttn_box[i]);
                *NY = i - 1;
                if (*NX > *NY) {
                  j = *NX;
                  *NX = *NY;
                  *NY = j;
                }
                I++;
              }

              if ((ST_NUM == 2 && I == 2) || (ST_NUM == 3 && I == 3)) {
                if (*NX >= BGN_ANT_I && *NX < END_ANT_I &&
                    *NY >= BGN_ANT_J && *NY < END_ANT_J) {
                  off_button(&idum, bttn_lab[*NX+1], bttn_box[*NX+1]);
                  off_button(&idum, bttn_lab[*NY+1], bttn_box[*NY+1]);
                  free_memory_block_float(bttn_box, ANT_NUM+CNTL_BTN_NUM);
                  free_memory_block_char(bttn_lab,  ANT_NUM+CNTL_BTN_NUM);
                  return (PR_SWT);
                } else {
                  I = 0;
                  *NX = -1;
                  *NY = -1;
                  PR_SWT = 1;
                  for (i=0; i<=END_ANT; i++) {
                    off_button(&idum, bttn_lab[i], bttn_box[i]);
                  }
                  if (ST_NUM == 3) {
                    off_button(&idum, bttn_lab[i], bttn_box[i]);
                    i++;
                    off_button(&idum, bttn_lab[i], bttn_box[i]);
                  }
                }
              }
            }
            break;
          }
        }
      }
    }
  }
}
