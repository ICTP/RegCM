include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = emcre_grid
emcre_grid_SOURCES = emcre_grid.F90 emcre_netcdf.F90 emcre_tools.F90 \
                     mo_f2kcli.F90 spline.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.F90

distclean-local:
	rm -f *.mod *__genmod.F90

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

emcre_grid.o: emcre_grid.F90 mo_f2kcli.o emcre_tools.o emcre_netcdf.o
emcre_netcdf.o : emcre_netcdf.F90 emcre_tools.o
emcre_tools.o : emcre_tools.F90
mo_f2kcli.o : mo_f2kcli.F90
