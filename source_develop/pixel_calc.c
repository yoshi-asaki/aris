#include <stdio.h>
#include <aris.h>

int  pixel_calc(double *pix_uvl, double *pix_mas,
                double wave_length, double uv_max,
                double *uv_factor, int FIELD_SIZE,
                int    CALC_ORDER_SWT)
{
  if (CALC_ORDER_SWT == 1) {
    *pix_uvl = *uv_factor * uv_max / (float)(FIELD_SIZE/2);
    *pix_mas = 1.0 / (2.0 * *uv_factor * uv_max / wave_length)
                          / dpi * 180.0 * 3600.0e3;
    return 1;
  } else if (CALC_ORDER_SWT == -1) {
    if (*pix_mas == 0.0) {
      printf("ERROR: PIXEL_CALC: pix_mas is not set.\n");
      return -1;
    } else {
      *pix_uvl = 1.0 / (2.0 * *pix_mas / wave_length)
                            / dpi * 180.0 * 3600.0e3 / (float)(FIELD_SIZE/2);
      *uv_factor = *pix_uvl / uv_max * (float)(FIELD_SIZE/2);
      return 1;
    }
  } else {
    printf("ERROR: PIXEL_CALC: SWITCH is invalid: %d.\n", CALC_ORDER_SWT);
    return -1;
  }
}
