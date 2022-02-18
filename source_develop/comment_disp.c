#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mathtools.h>
#include <cpgplot.h>
#include <astrotools.h>
#include <aris.h>


void comment_disp(struct comment_param *cmnt, char comment[][NCOMLEN],
                  char *add_comment, _Bool NEW_ADD)
{
  int    i, N;
  float  pgx, pgyr, pgyt;

  N = (int)((cmnt->ymax - cmnt->ymin) / cmnt->pitch);

  pgx = cmnt->xmin + 0.02;

  if (NEW_ADD == true) {
    if (strlen(comment[N-1]) == 0) {
      for (i=0; i<N; i++) {
        if (strlen(comment[i]) == 0) {
          char_copy(comment[i], add_comment);
          pgyr = cmnt->ymax - (i + 1) * cmnt->pitch + 0.4 * cmnt->pitch;
          cpgsci(7);
          cpgtext(pgx, pgyr, comment[i]);
          break;
        }
      }
    } else if (strlen(comment[N-1]) != 0) {
      for (i=0; i<N-1; i++) {
        pgyt = cmnt->ymax - (float)(i + 1) * cmnt->pitch + 0.4 * cmnt->pitch;
        cpgsci(4);
        pgyr = cmnt->ymax - (float)(i + 1) * cmnt->pitch;
        cpgrect(cmnt->xmin, cmnt->xmax, pgyr, pgyr + cmnt->pitch);
        if (strlen(comment[i+1]) != 0) {
          char_copy(comment[i], comment[i+1]);
          cpgsci(7);
          cpgtext(pgx, pgyt, comment[i]);
        } else {
          break;
        }
      }

      char_copy(comment[i], add_comment);
      pgyt = cmnt->ymax - (float)(i + 1) * cmnt->pitch + 0.4 * cmnt->pitch;
      cpgsci(4);
      pgyr = cmnt->ymax - (float)(i + 1) * cmnt->pitch;
      cpgrect(cmnt->xmin, cmnt->xmax, pgyr, pgyr + cmnt->pitch);
      cpgsci(7);
      cpgtext(pgx, pgyt, comment[i]);
      cpgsci(1);
    }
  } else if (NEW_ADD == false) {
    if (strlen(comment[N-1]) == 0) {
      for (i=0; i<N; i++) {
        if (strlen(comment[i]) == 0) {
          strcat(comment[i-1], add_comment);
          pgyt = cmnt->ymax - (float)i * cmnt->pitch + 0.4 * cmnt->pitch;
          cpgsci(7);
          cpgtext(pgx, pgyt, comment[i-1]);
          cpgsci(1);
          break;
        }
      }
    } else if (strlen(comment[N-1]) != 0) {
      strcat(comment[N-1], add_comment);
      pgyt = cmnt->ymax - (float)N * cmnt->pitch + 0.4 * cmnt->pitch;
      cpgsci(7);
      cpgtext(pgx, pgyt, comment[N-1]);
      cpgsci(1);
    }
  }

  return;
}
