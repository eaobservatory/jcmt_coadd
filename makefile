# Makefile for jcmt_coadd
# Quick bodge

# Author: Tim J (JACH)

# Last written by         Time-stamp: <10 Jan 96 0941 timj>


PKG_NAME = jcmt_coadd.exe

# Source files

FFILES = jcmt_coadd.f
CFILES = 

# Header files

HEADFILES = gsd_pars.inc


# Paths

INC     = -I/star/include -I.
LIB     = -L/star/lib -lgsd -lcnf

#LIB     = -lnag -L/jcmt_sw/specx/test/gsd -lgsd -L/star/lib -lcnf
NLIBS   = `./pgplotlink`

# Define default compilations

.f.o:
	$(FC) -c -e -w $(INC) $(@:.o=.f)

.c.o:
	$(CC) -c -w $(CFLAGS) $(INC) $(@:.o=.c)


# make options

default: build

build:  $(FFILES:.f=.o) $(CFILES:.c=.o) $(HEADFILES)
	echo "Make sure the SYSTEM variable is set"
	$(FC) $(FFILES:.f=.o) $(CFILES:.c=.o) -o $(PKG_NAME) \
        $(LIB) $(NLIBS)

clean:
	rm -f *.o
	rm -f $(PKG_NAME)
