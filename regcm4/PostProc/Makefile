# Makefile for RegCM4
# for PostProc code.
  
include ../Makefile.inc

all: average regrid sigma2p
	@cp average regrid sigma2p $(BINDIR)
	@echo "Programs all installed in $(BINDIR)"

mod_constants.o: ../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_derived.o: mod_derived.f90 mod_constants.o
	$(F90) -c $(LDFLAGS) $(F90FLAGS) -o $@ $<

average: average.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $<

regrid: regrid.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $<

sigma2p: sigma2p.F90 mod_constants.o
	$(F90) $(CPPFLAGS) $(NETCDFINC) $(LDFLAGS) $(F90FLAGS) \
               -o $@ $< $(NETCDFLIB)

clean:
	rm -f average regrid sigma2p *.o *.mod *__genmod*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

