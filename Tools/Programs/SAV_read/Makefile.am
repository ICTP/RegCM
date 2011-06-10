
include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = reSAV
reSAV_SOURCES = reSAV.F90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod

distclean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

reSAV.o: reSAV.F90
