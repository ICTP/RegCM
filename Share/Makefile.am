
include $(REGCM_ROOT)/makeinc

noinst_LIBRARIES = librcmlib.a
librcmlib_a_SOURCES = mod_constants.F90 mod_date.f90 mod_space.f90 \
                      mod_message.f90 mod_memutil.F90 mod_projections.f90 \
                      mod_dynparam.F90 mod_ensemble.f90 mod_stdio.f90 \
                      mod_realkinds.f90 mod_interp.f90 mod_nchelper.F90 \
                      mod_vertint.F90 mod_spline.f90 mod_domain.F90 \
                      mod_intkinds.f90

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

mod_constants.o: mod_constants.F90 mod_realkinds.o
mod_date.o: mod_date.f90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_message.o
mod_domain.o: mod_domain.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_nchelper.o mod_stdio.o mod_message.o
mod_dynparam.o: mod_dynparam.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_date.o
mod_ensemble.o: mod_ensemble.f90 mod_intkinds.o mod_realkinds.o mod_memutil.o mod_constants.o
mod_interp.o: mod_interp.f90 mod_intkinds.o mod_realkinds.o mod_memutil.o mod_stdio.o mod_message.o mod_constants.o
mod_memutil.o: mod_memutil.F90 mod_intkinds.o mod_realkinds.o mod_space.o mod_message.o mod_constants.o mod_date.o
mod_message.o: mod_message.f90 mod_intkinds.o mod_stdio.o
mod_nchelper.o: mod_nchelper.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_constants.o mod_memutil.o mod_dynparam.o mod_message.o
mod_projections.o: mod_projections.f90 mod_intkinds.o mod_realkinds.o mod_constants.o
mod_realkinds.o: mod_realkinds.f90
mod_intkinds.o: mod_intkinds.f90
mod_space.o: mod_space.f90 mod_intkinds.o mod_realkinds.o mod_date.o
mod_spline.o: mod_spline.f90 mod_intkinds.o mod_realkinds.o mod_constants.o
mod_stdio.o: mod_stdio.f90 mod_intkinds.o
mod_vertint.o: mod_vertint.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_constants.o mod_message.o
