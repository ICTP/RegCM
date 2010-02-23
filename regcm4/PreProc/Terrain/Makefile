# Makefile for RegCM4
# for Terrain code.
  
include ../../Makefile.inc

SRC = terrain.f

OBJ = terrain.o

all: terrain

.f.o:
	$(FC) -c $(FFLAGS) $*.f

terrain: $(OBJ)
	$(LD) -o $@ $(OBJ)

clean:
	rm -f *.o terrain
# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

terrain.o: domain.param terrain.f 
