
include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = clm2rcm
clm2rcm_SOURCES = mod_param_clm.F90 mod_grid.F90 mod_nclib.f90  \
                  mod_clm3grid.f90 clm2rcm.F90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

clm2rcm.o: clm2rcm.F90 mod_nclib.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_param_clm.mod mod_grid.mod \
 $(REGCM_ROOT)/Share/mod_date.mod $(REGCM_ROOT)/Share/mod_nchelper.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/Share/mod_message.mod \
 $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_stdio.mod \
 $(REGCM_ROOT)/Share/mod_domain.mod mod_clm3grid.mod
mod_clm3grid.o mod_clm3grid.mod: mod_clm3grid.f90 \
 $(REGCM_ROOT)/Share/mod_realkinds.mod
mod_nclib.o mod_nclib.mod: mod_nclib.f90 $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/Share/mod_stdio.mod $(REGCM_ROOT)/Share/mod_message.mod \
 $(REGCM_ROOT)/Share/mod_realkinds.mod
mod_param_clm.o mod_param_clm.mod: mod_param_clm.f90 \
 $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_stdio.mod \
 $(REGCM_ROOT)/Share/mod_message.mod
mod_grid.o mod_grid.mod: mod_grid.f90 \
 $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod
