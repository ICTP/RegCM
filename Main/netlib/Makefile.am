
noinst_LIBRARIES = libnetlib.a
libnetlib_a_SOURCES = linpack.f90 eispack.f90

%.o: %.f90
	$(FC) -c $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) -c $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean:
		@echo -n Removing objects, mods and $(noinst_LIBRARIES) ... 
		@rm -f *.o *.mod $(noinst_LIBRARIES)
		@echo Done

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

linpack.o: linpack.f90
eispack.o: eispack.f90
