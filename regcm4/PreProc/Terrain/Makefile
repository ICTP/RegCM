#Makefile for RegCM@  run on PC_LINUX (Intel 64bit CPU)
#  by using Intel FORTRAN compiler version 8 or 9

include ../../Makefile.inc

SRC = terrain.f

OBJ = terrain.o

all: terrain

.f.o:
	$(FC) -c $(FFLAGS) $*.f

terrain: $(OBJ)
	$(LD) -o $@ $(OBJ) $(NETCDFLIB)

clean:
	rm -f *.o terrain
# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

terrain.o: domain.param terrain.f 
