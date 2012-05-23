
vpath %.mod $(REGCM_ROOT)/Share:.

noinst_LIBRARIES = libmpp.a
libmpp_a_SOURCES = mod_mppparam.F90 mod_mpmessage.f90 mod_service.F90 \
                   mod_runparams.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

%.mod: %.f90 %.o
	@true

%.mod: %.F90 %.o
	@true

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_mpmessage.o mod_mpmessage.mod: mod_mpmessage.f90 mod_message.mod
mod_mppparam.o mod_mppparam.mod: mod_mppparam.F90 mod_realkinds.mod mod_dynparam.mod mod_mpmessage.mod mod_memutil.mod mod_date.mod mod_stdio.mod
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_runparams.o mod_runparams.mod: mod_runparams.f90 mod_realkinds.mod mod_date.mod mod_dynparam.mod mod_memutil.mod
mod_service.o mod_service.mod: mod_service.F90 mod_realkinds.mod mod_dynparam.mod
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
