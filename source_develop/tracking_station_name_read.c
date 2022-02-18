#include <stdio.h>
#include <stdlib.h>
#include <aris.h>


int  tracking_station_name_read(int trk_num,
                                char *station_name,
                                char trk_name[][10],
                                int *priority)
{
  int    i, j, idum;

  idum = 0;
  i = 0;
  while (station_name[i] != '\0') {
    if (station_name[i] == ',' || station_name[i] == ' ') {
      station_name[i] = '\0';

      j = 0;
      while (station_name[i+1+j] != '\0') {
        if (station_name[i+1+j] == '\n') {
          idum = 0;
          break;
        } else if (station_name[i+1+j] != ' ') {
          sscanf(station_name + i + 1 + j, "%d", &idum);
          break;
        }
        j++;
      }
      break;
    } else if (station_name[i] == '\n') {
      station_name[i] = '\0';
      idum = 0;
      break;
    }
    i++;
  }
  sprintf(trk_name[trk_num], "%s", station_name);
  priority[trk_num] = idum;
  trk_num++;

  return trk_num;
}
