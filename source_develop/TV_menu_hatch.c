#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cpgplot.h>
#include <aris.h>


int  TV_menu_hatch(float pgxmin, float pgxmax, float pgymin, float pgymax,
                   int color,    int fill_style)
{
  cpgsci(color);
  cpgsfs(fill_style);

  cpgrect(pgxmin, pgxmax, pgymin, pgymax);

  cpgsci(1);
  cpgsfs(1);

  return 1;
}
