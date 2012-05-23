
include $(REGCM_ROOT)/makeinc

vpath %.mod $(REGCM_ROOT)/Share:.

bin_PROGRAMS = clm2rcm
clm2rcm_SOURCES = mod_param_clm.F90 mod_grid.F90 mod_nclib.f90  \
                  mod_clm3grid.f90 clm2rcm.F90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

%.mod: %.f90 %.o
	@true

%.mod: %.F90 %.o
	@true

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
clm2rcm.o clm2rcm.mod: clm2rcm.F90 mod_nclib.mod mod_dynparam.mod mod_message.mod mod_grid.mod mod_param_clm.mod mod_date.mod mod_clm3grid.mod mod_memutil.mod mod_realkinds.mod mod_stdio.mod mod_domain.mod mod_nchelper.mod
mod_clm3grid.o mod_clm3grid.mod: mod_clm3grid.f90 mod_realkinds.mod
mod_grid.o mod_grid.mod: mod_grid.f90 mod_realkinds.mod mod_dynparam.mod mod_memutil.mod
mod_nclib.o mod_nclib.mod: mod_nclib.f90 mod_memutil.mod mod_message.mod mod_stdio.mod mod_realkinds.mod
mod_param_clm.o mod_param_clm.mod: mod_param_clm.f90 mod_realkinds.mod mod_message.mod mod_stdio.mod
