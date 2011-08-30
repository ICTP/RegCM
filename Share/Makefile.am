
noinst_LIBRARIES = librcmlib.a
librcmlib_a_SOURCES = mod_constants.F90 mod_date.f90 mod_space.f90 \
                      mod_message.f90 mod_memutil.F90 mod_projections.f90 \
                      mod_dynparam.F90 mod_service.F90 mod_ensemble.f90

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

mod_constants.o mod_constants.mod: mod_constants.F90
mod_date.o mod_date.mod: mod_date.f90 $(REGCM_ROOT)/external/m_realkinds.mod \
 $(REGCM_ROOT)/external/m_stdio.mod $(REGCM_ROOT)/external/m_die.mod
mod_dynparam.o mod_dynparam.mod: mod_dynparam.F90 mod_constants.mod \
 mod_date.mod
mod_memutil.o mod_memutil.mod: mod_memutil.F90 mod_space.mod \
 mod_message.mod $(REGCM_ROOT)/external/m_realkinds.mod mod_constants.mod \
 mod_date.mod
mod_message.o mod_message.mod: mod_message.f90 \
 $(REGCM_ROOT)/external/m_die.mod $(REGCM_ROOT)/external/m_stdio.mod
mod_projections.o mod_projections.mod: mod_projections.f90 \
 $(REGCM_ROOT)/external/m_realkinds.mod mod_constants.mod
mod_space.o mod_space.mod: mod_space.f90 \
  $(REGCM_ROOT)/external/m_realkinds.mod mod_date.mod
testdate.o: testdate.f90 mod_date.mod
mod_service.o mod_service.mod : mod_service.F90 \
  $(REGCM_ROOT)/external/m_realkinds.mod mod_dynparam.mod
mod_ensemble.o mod_ensemble.mod: mod_ensemble.f90 \
  mod_memutil.mod $(REGCM_ROOT)/external/m_realkinds.mod mod_constants.mod
