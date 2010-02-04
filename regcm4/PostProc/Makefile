#LINUX
include ../Makefile.inc

SRC = postproc.f netcdf.f

OBJ = postproc.o netcdf.o

all: postproc

postproc: $(OBJ)
	$(LD) -o $@ $(OBJ) $(NETCDFINC) $(NETCDFLIB)

clean:
	rm -f *.o postproc
# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

postproc.o: postproc.param postproc.f
