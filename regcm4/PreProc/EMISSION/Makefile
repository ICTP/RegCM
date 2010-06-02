# Makefile for RegCM4
# for CLM interface code.
  
include ../../Makefile.inc

SVNDEF := -D'SVN_REV="$(shell svnversion -n .)"'

all: emission report

MODOBJ = mod_emission.o

EMISSOBJ = emission.o bilinxo.o griddef.o julian2.o

REPORTOBJ = report.o bilinxo.o julian2.o

SHRMOD = mod_dynparam.o

report: $(SHRMOD) $(MODOBJ) $(REPORTOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(REPORTOBJ) $(LIB) \
               $(MODOBJ) $(SHRMOD) $(NETCDFLIB)

emission: $(SHRMOD) $(MODOBJ) $(EMISSOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(EMISSOBJ) $(LIB) \
               $(MODOBJ) $(SHRMOD) $(NETCDFLIB)

mod_dynparam.o: ../../Config/mod_dynparam.F90
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o emission report

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_emission.o: mod_emission.f90
bilinxo.o: bilinxo.f90
julian2.o: julian2.f90
griddef.o: griddef.f90 mod_emission.o
emission.o: emission.f90 mod_dynparam.o mod_emission.o
report.o: mod_emission.o mod_dynparam.o
