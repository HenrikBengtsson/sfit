INCLUDES = ${HOME}/include
LIBS = ${HOME}/lib
PKG_CFLAGS = -Wall -O3

all: cfit
	cp cfit dummy.so
	mkdir -p ../inst/bin/
	mv cfit ../inst/bin/

cfit: cfit.o spa.o cli.c spa.h
	$(CC) $(CFLAGS) $(PKG_CFLAGS) cli.c cfit.o spa.o -o cfit -lm 
  
test: cfit dye2.dat
	cfit 5 dye2.dat -a .01 -d 1 --outM stdout -q 10 -Q 90

clean:
	rm -f ../inst/bin/cfit
	rm -f dummy.so
	rm -f *.o

######################################################################
# HISTORY:
# 2011-05-15
# o Now ../inst/bin/ is created, in case it's missing.
# 2006-12-05
# o Update to make use of CC, CFLAGS and PKG_CFLAGS.
# 2006-05-20
# o Added a clean rule.
# o Added a rule to move the build binary to ../inst/bin/ when done.
# 2006-05-16
# o Removed -liberty from gcc in cfit target. /HB
######################################################################

