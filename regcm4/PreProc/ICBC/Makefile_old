#Makefile for RegCM@  run on PC_LINUX (Intel 64bit CPU)
#  by using Intel FORTRAN compiler version 8 or 9

include ../../Makefile.inc

SRC1a = SST_1DEG.f
OBJ1a = SST_1DEG.o

SRC1b = SST_FVGCM.f
OBJ1b = SST_FVGCM.o

SRC1c = SST_EH5OM.f
OBJ1c = SST_EH5OM.o

SRC1d = SST_ERSST.f
OBJ1d = SST_ERSST.o

SRC2 = ICBC.f
OBJ2 = ICBC.o

all: SST SST_FVGCM SST_EH5OM SST_ERSST ICBC

.f.o:
	$(FC) -c $(FFLAGS) $(NETCDFINC) $*.f

SST: $(OBJ1a)
	$(LD) -o $@ $(OBJ1a) $(NETCDFLIB)

SST_FVGCM: $(OBJ1b)
	$(LD) -o $@ $(OBJ1b)

SST_EH5OM: $(OBJ1c)
	$(LD) -o $@ $(OBJ1c)

SST_ERSST: $(OBJ1d)
	$(LD) -o $@ $(OBJ1d) $(NETCDFLIB)

ICBC: $(OBJ2)
	$(LD) -o $@ $(OBJ2) $(NETCDFLIB)

clean:
	rm -f *.o SST SST_FVGCM SST_EH5OM SST_ERSST ICBC
# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

SST_1DEG.o: icbc.param SST_1DEG.f
SST_FVGCM.o: icbc.param SST_FVGCM.f
SST_EH5OM.o: icbc.param SST_EH5OM.f
SST_ERSST.o: icbc.param SST_ERSST.f
ICBC.o: icbc.param ICBC.f
