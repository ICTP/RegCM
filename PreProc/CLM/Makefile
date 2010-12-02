# Makefile for RegCM4
# for CLM interface code.
  
include ../../Makefile.inc

all: clm2rcm

MODOBJ = mod_param_clm.o mod_read_domain.o mod_nclib.o mod_clm3grid.o

CLM2RCMOBJ = clm2rcm.o

SHRMOD = mod_dynparam.o mod_constants.o mod_date.o

clm2rcm: $(SHRMOD) $(MODOBJ) $(CLM2RCMOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(CLM2RCMOBJ) \
               $(MODOBJ) $(SHRMOD) $(NETCDFLIB)
	@cp $@ $(BINDIR)

mod_constants.o: ../../Config/mod_constants.F90
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

mod_dynparam.o: ../../Config/mod_dynparam.F90 mod_constants.o
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

mod_date.o: ../ICBC/mod_date.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c -o $@ $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

%.o: %.F90
	$(F90) $(CPPFLAGS) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o clm2rcm

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

clm2rcm.o: clm2rcm.F90 mod_dynparam.o mod_read_domain.o mod_date.o \
           mod_clm3grid.o
mod_clm3grid.o: mod_clm3grid.f90
mod_param_clm.o: mod_param_clm.f90
mod_read_domain.o: mod_read_domain.f90
mod_nclib.o: mod_nclib.f90
