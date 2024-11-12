include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = emcre_grid
emcre_grid_SOURCES = emcre_grid.F90 emcre_netcdf.F90 emcre_tools.F90 \
                     spline.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

emcre_grid.o: emcre_grid.F90 emcre_tools.o emcre_netcdf.o
emcre_netcdf.o : emcre_netcdf.F90 emcre_tools.o
emcre_tools.o : emcre_tools.F90
