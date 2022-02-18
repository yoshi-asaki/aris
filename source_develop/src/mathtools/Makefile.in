.c.o:
		$(CC) $(CFLAG) -c $*.c

SH		= /bin/sh

AR		= ar rcv

CC		= gcc

CFLAG           = -O2

LIBS		=

LIBDIR		= ../../lib

INCDIR		= ../../include

HEADER		= mathtools.h

RANLIB		= ranlib

OBJS		= allanv.o \
		  drotate.o \
		  drotation_matrix_set.o \
		  fftrfm.o \
		  inner_product3.o \
		  cross_product3.o \
		  pol_fit.o \
		  pol_fit_d.o \
		  frotate.o \
		  lstsqr.o \
		  lstsqr2.o \
		  fmaxmin.o \
		  dmaxmin.o \
		  dvec2scalar.o \
		  dvector_calc.o \
		  vec2scalar.o

PROGRAM		= libmathtools.a

all		: $(PROGRAM)

$(PROGRAM)	: $(OBJS)
		  @echo "Archive...."
		  $(AR) $(PROGRAM) $(OBJS) $(LIBS)
		  cp $(PROGRAM) $(LIBDIR)
		  @echo "Ranlib...."
		  $(RANLIB) $(LIBDIR)/$(PROGRAM)
		  cp $(HEADER) $(INCDIR)
		  @echo "done"

clean           :; rm -f $(OBJS) $(PROGRAM) core; cd $(LIBDIR); rm -f $(PROGRAM)

