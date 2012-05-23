
noinst_LIBRARIES = librcmlib.a
librcmlib_a_SOURCES = mod_constants.F90 mod_date.f90 mod_space.f90 \
                      mod_message.f90 mod_memutil.F90 mod_projections.f90 \
                      mod_dynparam.F90 mod_ensemble.f90 mod_stdio.f90 \
                      mod_realkinds.f90 mod_interp.f90 mod_nchelper.F90 \
                      mod_vertint.F90 mod_spline.f90 mod_domain.F90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

%.mod: %.f90 %.o
	@true

%.mod: %.F90 %.o
	@true

distclean-local:
distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_constants.o mod_constants.mod: mod_constants.F90
mod_date.o mod_date.mod: mod_date.f90 mod_realkinds.mod mod_stdio.mod mod_message.mod
mod_domain.o mod_domain.mod: mod_domain.F90 mod_realkinds.mod mod_dynparam.mod mod_memutil.mod mod_nchelper.mod mod_stdio.mod mod_message.mod
mod_dynparam.o mod_dynparam.mod: mod_dynparam.F90 mod_constants.mod mod_date.mod
mod_ensemble.o mod_ensemble.mod: mod_ensemble.f90 mod_realkinds.mod mod_memutil.mod mod_constants.mod
mod_interp.o mod_interp.mod: mod_interp.f90 mod_memutil.mod mod_realkinds.mod mod_stdio.mod mod_message.mod mod_constants.mod
mod_memutil.o mod_memutil.mod: mod_memutil.F90 mod_space.mod mod_message.mod mod_realkinds.mod mod_constants.mod mod_date.mod
mod_message.o mod_message.mod: mod_message.f90 mod_stdio.mod
mod_nchelper.o mod_nchelper.mod: mod_nchelper.F90 mod_realkinds.mod mod_stdio.mod mod_constants.mod mod_memutil.mod mod_dynparam.mod mod_message.mod
mod_projections.o mod_projections.mod: mod_projections.f90 mod_realkinds.mod mod_constants.mod
mod_realkinds.o mod_realkinds.mod: mod_realkinds.f90
mod_space.o mod_space.mod: mod_space.f90 mod_realkinds.mod mod_date.mod
mod_spline.o mod_spline.mod: mod_spline.f90 mod_realkinds.mod mod_constants.mod
mod_stdio.o mod_stdio.mod: mod_stdio.f90
mod_vertint.o mod_vertint.mod: mod_vertint.F90 mod_realkinds.mod mod_stdio.mod mod_constants.mod mod_message.mod
