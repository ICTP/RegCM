
include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = CheckSun
CheckSun_SOURCES = CheckSun.F90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

CheckSun.o: CheckSun.F90
