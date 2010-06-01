# Makefile for RegCM4
# for PostProc code.
  
include ../Makefile.inc

all: average regrid sigma2p
	@cp average regrid sigma2p $(BINDIR)
	@echo "Programs all installed in $(BINDIR)"

average: average.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $< $(LIB)

regrid: regrid.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $< $(LIB)

sigma2p: sigma2p.f90
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $< $(LIB)

clean:
	rm -f average regrid sigma2p

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

