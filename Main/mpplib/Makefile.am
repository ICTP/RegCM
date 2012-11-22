
vpath %.o $(REGCM_ROOT)/Share:.

noinst_LIBRARIES = libmpp.a
libmpp_a_SOURCES = mod_mppparam.F90 mod_mpmessage.f90 mod_service.F90 \
                   mod_runparams.F90 mod_outvars.f90 mod_ncout.F90

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
mod_mpmessage.o : mod_mpmessage.f90 mod_message.o
mod_mppparam.o : mod_mppparam.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_mpmessage.o mod_memutil.o mod_date.o mod_stdio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_outvars.o: mod_outvars.f90 mod_realkinds.o
mod_runparams.o : mod_runparams.F90 mod_intkinds.o mod_realkinds.o mod_date.o mod_dynparam.o mod_memutil.o
mod_service.o : mod_service.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o mod_mpmessage.o mod_mppparam.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_ncout.o: mod_ncout.F90 mod_dynparam.o mod_runparams.o mod_ncstream_types.o mod_ncstream.o mod_outvars.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
