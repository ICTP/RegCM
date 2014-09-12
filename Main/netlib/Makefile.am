
include $(REGCM_ROOT)/makeinc

noinst_LIBRARIES = libnetlib.a
libnetlib_a_SOURCES = linpack.F90 eispack.F90 lapack_dgbsv.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

linpack.o : linpack.F90
eispack.o : eispack.F90
lapack_dgbsv.o : lapack_dgbsv.F90
