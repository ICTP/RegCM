# Makefile for RegCM4
# for CLM interface code.
  
include ../../Makefile.inc

SVNDEF := -D'SVN_REV="$(shell svnversion -n .)"'

all: clm2rcm

MODOBJ = mod_param_clm.o

CLM2RCMOBJ = clm2rcm.o bilinx4d.o clm3grid1.o clm3grid2.o comp.o \
             fexist.o julian.o makectl.o maskme.o param.o readcdfr4_360.o \
             readcdfr4.o readcdfr4_iso.o xlon360.o nclib.o

SHRMOD = mod_param.o

clm2rcm: $(SHRMOD) $(MODOBJ) $(CLM2RCMOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(CLM2RCMOBJ) $(LIB) \
               $(MODOBJ) $(SHRMOD) $(NETCDFLIB)

mod_param.o: ../Terrain/mod_param.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o clm2rcm

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

bilinx4d.o: bilinx4d.f90
clm2rcm.o: clm2rcm.f90 mod_param_clm.o mod_param.o
clm3grid1.o: clm3grid1.f90
clm3grid2.o: clm3grid2.f90
comp.o: comp.f90
fexist.o: fexist.f90
julian.o: julian.f90
makectl.o: makectl.f90
maskme.o: maskme.f90
mod_param_clm.o: mod_param_clm.f90
nclib.o: nclib.f90
param.o: param.f90
readcdfr4_360.o: readcdfr4_360.f90
readcdfr4.o: readcdfr4.f90
readcdfr4_iso.o: readcdfr4_iso.f90
xlon360.o: xlon360.f90
