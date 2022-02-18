#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <aris.h>

#define CTRL_C           3
#define CTRL_D           4
#define BACKSPACE_KEY    8
#define RETURN_KEY      13


void red_curs(float , float , float , float );
void show_character(char *, float , float , float );

int   tv_get_param(char   *data_type,
                   float  *cursor_pos,
                   float  *bttn_box,
                   float  pitch,  char   *string,
                   double lim_min,  double  lim_max)
{
  int    i;
  char   c;
  float  cxmin, cxmax, cymin, cymax;
  double SET_PARAM;

/*
--------------------
*/

  cxmin = bttn_box[1] - 0.40 * pitch;
  cxmax = bttn_box[1] - 0.20 * pitch;
  cymin = 0.6*bttn_box[2] + 0.4*bttn_box[3] - 0.05 * pitch;
  cymax = 0.6*bttn_box[2] + 0.4*bttn_box[3] + 0.50 * pitch;

/*
--------------------
*/

  cpgsci(1);
  cpgrect(bttn_box[0], bttn_box[1], bttn_box[2], bttn_box[3]);
  cpgsci(2);
  cpgsfs(2);
  cpgrect(bttn_box[0], bttn_box[1], bttn_box[2], bttn_box[3]);
  cpgsfs(1);
  cpgsci(0);
  cpgptxt(bttn_box[1]-0.015, 0.6*bttn_box[2]+0.4*bttn_box[3], 0.0, 1.0, string);

  red_curs(cxmin, cxmax, cymin, cymax);

  i = strlen(string);
  while (1) {
    cpgband(0, 1, cursor_pos[0], cursor_pos[1], cursor_pos, cursor_pos+1, &c);
    if (c == CTRL_C || c == CTRL_D) {
      exit (-1);
    } else if (c == RETURN_KEY) {
      break;
    } else if (c == BACKSPACE_KEY) {
      string[i] = '\0';
      if (i != 0) {
        i--;
        string[i] = '\0';
      }

      cpgsci(1);
      cpgrect(bttn_box[0], bttn_box[1], bttn_box[2], bttn_box[3]);
      cpgsci(2);
      cpgsfs(2);
      cpgrect(bttn_box[0], bttn_box[1], bttn_box[2], bttn_box[3]);
      cpgsfs(1);
      cpgsci(0);
      show_character(string, bttn_box[1], bttn_box[2], bttn_box[3]);

    } else {
      string[i] = c;
      i++;
      string[i] = '\0';

      cpgsci(1);
      cpgrect(bttn_box[0], bttn_box[1], bttn_box[2], bttn_box[3]);
      cpgsci(2);
      cpgsfs(2);
      cpgrect(bttn_box[0], bttn_box[1], bttn_box[2], bttn_box[3]);
      cpgsfs(1);
      cpgsci(0);
      show_character(string, bttn_box[1], bttn_box[2], bttn_box[3]);
      cpgsci(1);
    }
    red_curs(cxmin, cxmax, cymin, cymax);
  }

  cpgsci(1);
  cpgrect(bttn_box[0], bttn_box[1], bttn_box[2], bttn_box[3]);

  if (strncmp(data_type, "char", 4) == 0) {
    cpgptxt(bttn_box[1]-0.015, 0.6*bttn_box[2]+0.4*bttn_box[3],
            0.0, 1.0, string);
  } else {
    sscanf(string, "%lf", &SET_PARAM);
    if (strlen(string) == 0) {
      sprintf(string, "0");
    }

    if (lim_min != lim_max) {
      if (SET_PARAM < lim_min || SET_PARAM > lim_max) {
        sprintf(string, "Invalid!");
        SET_PARAM = 0;
      }
    }
  }
  cpgsci(0);
  cpgptxt(bttn_box[1]-0.015, 0.6*bttn_box[2]+0.4*bttn_box[3], 0.0, 1.0, string);
  cpgsci(1);

  return 1;
}



void red_curs(float pgxmin, float pgxmax, float pgymin, float pgymax)
{
  cpgsci(2);
  cpgrect(pgxmin, pgxmax, pgymin, pgymax);
  cpgsci(1);
}



void show_character(char *string, float pgxmax, float pgymin, float pgymax)
{
  int   i;

  i = 1;
  while (string[strlen(string)-i] == ' ') {
    i++;
  }
  cpgptxt(pgxmax-(float)i*0.012, 0.6*pgymin+0.4*pgymax,
          0.0, 1.0, string);
}
