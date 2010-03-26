#LINUX
include ../Makefile.inc
SRC = postproc.f netcdf.f
OBJ = postproc.o netcdf.o
all: postproc
doSRF_average:
	f2py --verbose  -m  average -c  --f77exec=ifort  averages.f --f77flags=" -assume byterecl  -convert big_endian" --f90flags=" -assume byterecl -convert big_endian"

postproc: $(OBJ)
	$(LD) -o $@ $(OBJ) $(NETCDFINC) $(NETCDFLIB)
clean:
	rm -f *.o postproc
# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
postproc.o: postproc.param postproc.f
