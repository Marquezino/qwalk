# QWalk
# Version: 1.4
# Date: Nov 2012
# Author: F.L. Marquezino

SHELL = /bin/sh

prefix = ..
exec_prefix = $(prefix)
includedir = $(exec_prefix)/include
bindir = $(exec_prefix)/bin
libdir = $(exec_prefix)/lib
docdir = $(exec_prefix)/doc

CC = gcc
AR = ar
RANLIB = ranlib

LINK = -L$(libdir) -lqwalk -lm
CFLAGS = -I$(includedir) -O2 -std=c99

_MEM_OBJS = qwmem_int.o qwmem_real.o qwmem_complex.o
_QW_OBJS = qwcoin.o qwstate.o qwprob.o qwstatistics.o qwlinks.o \
	qwscreen.o qwmeasure.o
_QWIO_OBJS = qwcoin_io.o qwstate_io.o qwprob_io.o qwstatistics_io.o \
	qwoptions_io.o qwoptions_io_read.o qwextra_io.o 

MEM_OBJS = $(patsubst %,$(libdir)/%,$(_MEM_OBJS))
QW_OBJS = $(patsubst %,$(libdir)/%,$(_QW_OBJS))
QWIO_OBJS = $(patsubst %,$(libdir)/%,$(_QWIO_OBJS))

LIB_FILE = $(libdir)/libqwalk.a
LIB_OBJS = $(MEM_OBJS) $(QW_OBJS) $(QWIO_OBJS)

PROG_QW1D = $(bindir)/qw1d
QW1D_OBJS = $(libdir)/qw1d.o $(libdir)/qw1d_sub.o
#I could have used patsubst here

PROG_QW2D = $(bindir)/qw2d
QW2D_OBJS = $(libdir)/qw2d.o $(libdir)/qw2d_sub.o

PROG_AKR = $(bindir)/akr
AKR_OBJS = $(libdir)/akr.o $(libdir)/qw2d_sub.o

PROG_AKRR = $(bindir)/akrr
AKRR_OBJS = $(libdir)/akrr.o $(libdir)/qw2d_sub.o

PROG_QWAMPL = $(bindir)/qwamplify
QWAMPL_OBJS = $(libdir)/qwamplify.o

PROGS = $(PROG_QW1D) $(PROG_QW2D) $(PROG_QWAMPL) $(PROG_AKR) $(PROG_AKRR)

all:	$(LIB_FILE) $(PROGS)

$(PROG_QW2D):	$(QW2D_OBJS) $(LIB_FILE)
	$(CC) $(QW2D_OBJS) $(LINK) $(CFLAGS) -o $(PROG_QW2D)

$(PROG_QW1D):	$(QW1D_OBJS) $(LIB_FILE)
	$(CC) $(QW1D_OBJS) $(LINK) $(CFLAGS) -o $(PROG_QW1D)

$(PROG_AKR):	$(AKR_OBJS) $(LIB_FILE)
	$(CC) $(AKR_OBJS) $(LINK) $(CFLAGS) -o $(PROG_AKR)

$(PROG_AKRR):	$(AKRR_OBJS) $(LIB_FILE)
	$(CC) $(AKRR_OBJS) $(LINK) $(CFLAGS) -o $(PROG_AKRR)


$(PROG_QWAMPL):	$(QWAMPL_OBJS)
	$(CC) $(QWAMPL_OBJS) $(LINK) $(CFLAGS) -o $(PROG_QWAMPL)

$(LIB_FILE):	$(LIB_OBJS)
	$(AR) rc $(LIB_FILE) $(LIB_OBJS)
	$(RANLIB) $(LIB_FILE)

$(libdir)/%.o:	%.c $(includedir)/%.h $(includedir)/qwconsts.h Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

# some files may not have a corresponding header (e.g., qw1d.c and qw2d.c)
$(libdir)/%.o:	%.c $(includedir)/qwconsts.h Makefile 
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY:	clean doc

clean:
	rm $(libdir)/*.o

doc:
	a2ps *.c $(includedir)/*.h -b"QWalk - Source Code" --toc -4 -o $(docdir)/QWalk-Source.ps
