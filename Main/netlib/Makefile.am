
noinst_LIBRARIES = libnetlib.a
libnetlib_a_SOURCES = linpack.f90 eispack.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

linpack.o: linpack.f90
eispack.o: eispack.f90
