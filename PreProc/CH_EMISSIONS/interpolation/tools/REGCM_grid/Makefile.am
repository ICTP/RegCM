include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = emcre_grid
emcre_grid_SOURCES = emcre_grid.f90 emcre_netcdf.f90 emcre_tools.f90 \
                     mo_f2kcli.F90 spline.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.f90

distclean-local:
	rm -f *.mod *__genmod.f90

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

emcre_grid.o: emcre_grid.f90 mo_f2kcli.o emcre_tools.o emcre_netcdf.o
emcre_netcdf.o emcre_netcdf.mod: emcre_netcdf.f90 emcre_tools.o
emcre_tools.o emcre_tools.mod: emcre_tools.f90
