#include <stdio.h>
#include <stdlib.h>
#include <cpgplot.h>


int   pgarea(
         float        , float        , float       *, float       *);
int   pgarea_set(
         float       *, float       *, float       *, float       *);
int   pgxregion(
         float        , float        , float       *, float       *);
void  pg_skyplot_frame(                                            );
void  pg_skyplot_trajectory(
         int          , float       *, float       *, int          ,
         int                                                       );
void  palett(
         int          , float        , float                       );
void  pg_panel_tile(
         int          , int          , float       *,
         int         *, int         *, int         *, int         *,
         float       *, float       *, float       *, float       *,
         float       *, float       *, float       *, float       *);
