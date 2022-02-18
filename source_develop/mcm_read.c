#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <aris.h>


int  mcm_read(char *mcm_file, float *dist, int dmax, double pix_mas)
{
  int    n;
  char   string[500], comp;
  float  semi_maj, semi_min, pa;
  float  xoff,  yoff;
  float  tflux;
  FILE   *fp;

/*
-----------------------------------
*/

  if ((fp = fopen(mcm_file, "r")) == NULL) {
    printf("ERROR: mcm_read: file cannot be found : %s\n", mcm_file);
    return (__NG__);
  }
  while (1) {
    if (fgets(string, sizeof(string), fp) == NULL) {
      break;
    }
    if (string[0] != '#') {
      sscanf(string, "%f, %f, %f, %f, %f, %f, %d",
             &tflux, &semi_maj, &semi_min, &pa, &xoff, &yoff, &n);
      semi_maj *= 0.5;
      semi_min *= 0.5;
      pa *= (dpi / 180.0);
      if (n == 0) {
        comp = 'P';
      } else if (n == 1) {
        comp = 'S';
      } else if (n == 2) {
        comp = 'G';
      }
      source_model3(dist, dmax, pix_mas, 1, &tflux,
                    &semi_maj, &semi_min, &pa, &xoff, &yoff, &comp);
    }
  }
  fclose (fp);

  return (__GO__);
}
