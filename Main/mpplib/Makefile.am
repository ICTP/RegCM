
include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:.

noinst_LIBRARIES = libmpp.a
libmpp_a_SOURCES = mod_mppparam.F90 mod_mpmessage.F90 mod_service.F90 \
                   mod_runparams.F90 mod_outvars.F90 mod_ncout.F90 \
                   mod_regcm_types.F90 mod_ipcc_scenario.F90 \
                   mod_timer.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_mpmessage.o : mod_mpmessage.F90 mod_intkinds.o mod_message.o mod_realkinds.o mod_stdio.o
mod_ipcc_scenario.o : mod_ipcc_scenario.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_mpmessage.o mod_runparams.o mod_memutil.o mod_stdio.o mod_mppparam.o
mod_mppparam.o : mod_mppparam.F90 mod_date.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_stdio.o mod_constants.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_ncout.o : mod_ncout.F90 mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mppparam.o mod_ncstream.o mod_ncstream_types.o mod_outvars.o mod_realkinds.o mod_runparams.o mod_service.o mod_stdio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_outvars.o : mod_outvars.F90 mod_realkinds.o
mod_runparams.o : mod_runparams.F90 mod_date.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_timer.o mod_spline.o mod_mpmessage.o mod_constants.o
mod_service.o : mod_service.F90 mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_stdio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_regcm_types.o : mod_regcm_types.F90 mod_intkinds.o mod_realkinds.o
mod_timer.o : mod_timer.F90 mod_realkinds.o mod_intkinds.o mod_constants.o mod_date.o
