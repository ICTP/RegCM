
include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:.

bin_PROGRAMS = ReadRegcmio
ReadRegcmio_SOURCES = ReadRegcmio.F90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod

distclean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

ReadRegcmio.o: ReadRegcmio.F90
