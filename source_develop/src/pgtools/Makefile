.c.o:
		$(CC) $(CFLAG) -c -I$(INCTO) $*.c

SH		= /bin/sh

AR		= ar rcv

CC		= gcc

CFLAG           = -O2

LIBS		=

LIBDIR		= ../../lib

INCDIR		= ../../include

INCTO		= /usr/local/pgplot

HEADER		= pgtools.h

RANLIB		= ranlib

OBJS            = pgarea.o \
                  palett.o \
                  pgxregion.o \
                  pgarea_set.o \
                  pg_panel_tile.o \
                  pg_skyplot_frame.o \
                  pg_skyplot_trajectory.o

PROGRAM		= libpgtools.a

all		: $(PROGRAM)

$(PROGRAM)	: $(OBJS)
		  @echo "Archive...."
		  $(AR) $(PROGRAM) $(OBJS) $(LIBS)
		  cp $(PROGRAM) $(LIBDIR)
		  @echo "Ranlib...."
		  $(RANLIB) $(LIBDIR)/$(PROGRAM)
		  cp $(HEADER) $(INCDIR)
		  @echo "done"

clean           :; rm -f $(OBJS) $(PROGRAM) core

