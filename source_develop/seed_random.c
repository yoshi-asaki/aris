#include <stdio.h>
#include <stdbool.h>
#include <aris.h>
#include <unistd.h>

/****
#define __CLOCK_GETTIME__
****/
#ifndef __CLOCK_GETTIME__
  #define __GETTIMEOFDAY__
#endif


#ifdef __CLOCK_GETTIME__
  #include <time.h>
  struct timespec ts;
#elif defined __GETTIMEOFDAY__
  #include <sys/time.h>
  struct timeval  TV;
  /**** 2010.08.31  timezone structure is too old.
  struct timezone TZ;
  ****/
#endif


void seed_random(_Bool RAND_SWT)
{
  if (RAND_SWT == true) {
#ifdef __CLOCK_GETTIME__
/**
    clock_gettime(CLOCK_REALTIME, &ts);
    srand(ts.tv_nsec);
**/
    srand((unsigned)time(NULL));

#elif defined __GETTIMEOFDAY__
/**** 2010.08.31  timezone structure is too old.
    gettimeofday(&TV, &TZ);
****/
    gettimeofday(&TV, NULL);
    srand(TV.tv_usec);
#endif
  } else {
    srand(17);
  }
}
