#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <aris.h>


int      get_srt_link(struct srt_data_link *srt_link, double *fov_st)
{
  int    i, j, n;
  int    dtmp1, dtmp2, dtmp3;
  char   string[60];
  FILE   *srt_link_fp;

/*
---------------------------
*/

  for (i=0; i<360*150; i++) {
    *(srt_link->mask + i) = false;
  }

  if ((srt_link_fp = fopen("aris_input/srt_mask.dat", "r")) != NULL) {
    if (fgets(string, sizeof(string), srt_link_fp) == NULL) {
      printf("ERROR: GET_SRT_LINK: Invalid input.\n");
      fclose (srt_link_fp);
      return (__NG__);
    }
    sscanf(string, "%lf %lf %lf %lf",
           srt_link->FOV_rot_axis,
           srt_link->FOV_rot_axis + 1,
           srt_link->FOV_rot_axis + 2,
           &srt_link->FOV_rot_angle);
    if (fgets(string, sizeof(string), srt_link_fp) == NULL) {
      printf("ERROR: GET_SRT_LINK: Invalid input.\n");
      fclose (srt_link_fp);
      return (__NG__);
    }
    sscanf(string, "%d", &srt_link->nzenith);

    for (i=0; i<360; i++) {
      for (j=0; j<srt_link->nzenith; j++) {
        if (fgets(string, sizeof(string), srt_link_fp) == NULL) {
          printf("ERROR: GET_SRT_LINK: Invalid input.\n");
          fclose (srt_link_fp);
          return (__NG__);
        }
        sscanf(string, "%d %d %d", &dtmp1, &dtmp2, &dtmp3);
        if (dtmp3 != 0) {
          *(srt_link->mask + srt_link->nzenith*dtmp1 + dtmp2) = true;
        }
      }
    }
    fclose (srt_link_fp);
  }

/*
--------------- FOV Steradian ----------------
*/

  *fov_st = 0.0;
  for (j=0; j<srt_link->nzenith; j++) {
    n = 0;
    for (i=0; i<360; i++) {
      if (*(srt_link->mask + srt_link->nzenith*i + j) == false) {
        n++;
      }
    }
    *fov_st += (double)n / 180.0 * dpi *
               (cos((double)j/180.0*dpi) - cos((double)(j+1)/180.0*dpi));
  }

  return (__GO__);
}
