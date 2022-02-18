#
# Makefile
#

CC = gcc
FC = gfortran
DEFS = -DHAVE_CONFIG_H
LIBS = 

SRCS = sample.c
OBJS = $(SRCS:.c=.o)
PROG = sample

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LIBS)

.c.o:
	$(CC) $(DEFS) -c $< -o $@

clean:
	rm -f $(OBJS)
