# Makefile for RegCM4
# for PostProc code.
  
include ../Makefile.inc

all: average regrid sigma2p
	@cp average regrid sigma2p $(BINDIR)
	@echo "Programs all installed in $(BINDIR)"

average: average.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $<

regrid: regrid.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $<

sigma2p: sigma2p.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $<

clean:
	rm -f average regrid sigma2p *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

