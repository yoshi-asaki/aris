#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <aris.h>

int    button_chk(float *cursor_pos, float *button_box)
{
  if (cursor_pos[0] >= button_box[0] && cursor_pos[0] <= button_box[1] &&
      cursor_pos[1] >= button_box[2] && cursor_pos[1] <= button_box[3]) {
    return (1);
  } else {
    return (0);
  }
}
