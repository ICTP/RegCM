# Makefile for RegCM4
# for PostProc code.
  
include ../Makefile.inc

all: sigma2p
	@cp sigma2p $(BINDIR)
	@cp regrid.sh $(BINDIR)/regrid
	@cp average.sh $(BINDIR)/average
	@echo "Programs all installed in $(BINDIR)"

mod_constants.o: ../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_derived.o: mod_derived.f90 mod_constants.o
	$(F90) -c $(LDFLAGS) $(F90FLAGS) -o $@ $<

sigma2p: sigma2p.F90 mod_constants.o
	$(F90) $(CPPFLAGS) $(NETCDFINC) $(LDFLAGS) $(F90FLAGS) \
               -o $@ $< $(NETCDFLIB)

clean:
	rm -f sigma2p *.o *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

