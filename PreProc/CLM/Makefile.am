
include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = clm2rcm
clm2rcm_SOURCES = mod_param_clm.F90 mod_read_domain.f90 mod_nclib.f90  \
                  mod_clm3grid.f90 clm2rcm.F90

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

clm2rcm.o: clm2rcm.F90 mod_nclib.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod mod_read_domain.mod mod_param_clm.mod \
 $(REGCM_ROOT)/Share/mod_date.mod mod_clm3grid.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_die.mod
mod_clm3grid.o mod_clm3grid.mod: mod_clm3grid.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod
mod_nclib.o mod_nclib.mod: mod_nclib.f90 $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_realkinds.mod
mod_param_clm.o mod_param_clm.mod: mod_param_clm.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/external/m_stdio.mod \
 $(REGCM_ROOT)/external/m_die.mod
mod_read_domain.o mod_read_domain.mod: mod_read_domain.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_die.mod
