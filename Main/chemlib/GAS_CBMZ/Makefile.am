
noinst_LIBRARIES = libcbmz.a
libcbmz_a_SOURCES = mod_cbmz_boxvars.F90 mod_cbmz_chemmech.F90 \
  mod_cbmz_chemvars.F90 mod_cbmz_chemlocal.F90 mod_cbmz_molwg.f90

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

mod_cbmz_boxvars.o mod_cbmz_boxvars.mod: mod_cbmz_boxvars.F90 \
  $(REGCM_ROOT)/external/m_realkinds.mod
mod_cbmz_chemmech.o mod_cbmz_chemmech.mod: mod_cbmz_chemmech.F90 \
  $(REGCM_ROOT)/external/m_realkinds.mod
mod_cbmz_chemlocal.o mod_cbmz_chemlocal.mod: mod_cbmz_chemlocal.F90 \
  $(REGCM_ROOT)/external/m_realkinds.mod mod_cbmz_chemmech.o
mod_cbmz_chemvars.o mod_cbmz_chemvars.mod: mod_cbmz_chemvars.F90 \
  $(REGCM_ROOT)/external/m_realkinds.mod mod_cbmz_chemmech.o
mod_cbmz_molwg.o mod_cbmz_molwg.mod: mod_cbmz_molwg.f90 \
  $(REGCM_ROOT)/external/m_realkinds.mod
