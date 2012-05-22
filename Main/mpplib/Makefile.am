
noinst_LIBRARIES = libmpp.a
libmpp_a_SOURCES = mod_mppparam.F90 mod_mpmessage.f90 mod_service.F90

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

mod_service.o mod_service.mod : mod_service.F90 \
 $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_dynparam.mod
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
mod_mpmessage.o mod_mpmessage.mod: mod_mpmessage.f90 \
 $(REGCM_ROOT)/Share/mod_message.mod $(REGCM_ROOT)/Share/mod_dynparam.mod
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
mod_mppparam.o mod_mppparam.mod: mod_mppparam.F90 \
 $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_date.mod \
 $(REGCM_ROOT)/Share/mod_message.mod $(REGCM_ROOT)/Share/mod_realkinds.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod mod_mpmessage.mod
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
