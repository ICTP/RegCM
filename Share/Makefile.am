
noinst_LIBRARIES = librcmlib.a
librcmlib_a_SOURCES = mod_constants.F90 mod_date.f90 mod_space.f90 \
                      mod_message.f90 mod_memutil.F90 mod_projections.f90 \
                      mod_dynparam.F90 mod_service.F90 mod_ensemble.f90 \
                      mod_stdio.f90 mod_realkinds.f90 mod_mppparam.F90 \
                      mod_mpmessage.f90

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
mod_date.o mod_date.mod: mod_date.f90 mod_realkinds.o mod_stdio.o mod_message.o
mod_dynparam.o mod_dynparam.mod: mod_dynparam.F90 mod_constants.o mod_date.o
mod_memutil.o mod_memutil.mod: mod_memutil.F90 mod_space.o \
  mod_message.o mod_realkinds.o mod_constants.o mod_date.o
mod_message.o mod_message.mod: mod_message.f90 mod_stdio.o
mod_projections.o mod_projections.mod: mod_projections.f90 \
  mod_realkinds.o mod_constants.o
mod_space.o mod_space.mod: mod_space.f90 mod_realkinds.o mod_date.o
mod_service.o mod_service.mod : mod_service.F90 mod_realkinds.o mod_dynparam.o
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
mod_ensemble.o mod_ensemble.mod: mod_ensemble.f90 \
  mod_memutil.o mod_realkinds.o mod_constants.o
mod_stdio.o mod_stdio.mod: mod_stdio.f90
mod_realkinds.o mod_realkinds.mod: mod_realkinds.f90
mod_mpmessage.o mod_mpmessage.mod: mod_mpmessage.f90 mod_stdio.o mod_dynparam.o
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
mod_mppparam.o mod_mppparam.mod: mod_mppparam.F90 mod_dynparam.o mod_date.o \
  mod_message.o mod_realkinds.o
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
