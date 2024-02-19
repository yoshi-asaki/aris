.c.o:
		$(CC) $(CFLAG) -c $*.c

SH		= /bin/sh

AR		= ar rcv

CC		= gcc

CFLAG           = -O2

LIBS		=

LIBDIR		= ../../lib

INCDIR		= ../../include

HEADER		= astrotools.h

RANLIB		= ranlib

OBJS		= MJD.o \
		  MJD2date.o \
		  UTC_minus_TAI.o \
		  GC_CC.o \
		  B1950toJ2000.o \
		  vlen2.o \
		  vlen3.o \
		  ET.o \
		  EPSIRON.o \
		  sepang.o \
		  Q_tools.o \
		  sun_position.o \
		  heirocentric_equatorial_rectangular_coordinates_earth_position.o \
		  luna_position.o \
		  azel_position.o \
		  nutation.o \
		  nutation_calc.o \
		  mean_obliquity_of_the_ecliptic.o \
		  precession.o \
		  precession_calc.o \
		  LST.o \
		  GST.o \
		  earth_rotation.o \
		  topocentric_equatorial_rectangular_coordinate.o \
		  annual_parallax.o \
		  solar_motion_J2000.o \
		  baseline2antenna_number.o \
		  baseline_number.o \
		  baseline_base2antenna_base_solution.o \
		  minor_shift_refpos.o

PROGRAM		= libastrotools.a

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

