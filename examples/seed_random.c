#include <stdio.h>
#include <phase_screen.h>

#include <sys/time.h>
#include <unistd.h>
struct timeval  TV;
struct timezone TZ;


void seed_random(int RAND_SWT)
{
  if (RAND_SWT == ON) {
    gettimeofday(&TV, &TZ);
    srand(TV.tv_usec);
  } else {
    srand(17);
  }
}
