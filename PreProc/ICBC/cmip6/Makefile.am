
vpath %.o $(REGCM_ROOT)/Share:.

include $(REGCM_ROOT)/makeinc

noinst_LIBRARIES = libcmip6.a
libcmip6_a_SOURCES = mod_cmip6_helper.F90 mod_cmip6_cnrm.F90 \
                     mod_cmip6_cesm.F90 mod_cmip6_miroc6.F90 \
                     mod_cmip6_hadmm.F90 mod_cmip6_normm.F90 \
                     mod_cmip6_mpihr.F90 mod_cmip6_ecea.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_cmip6_helper.o : mod_cmip6_helper.F90 mod_realkinds.o mod_intkinds.o mod_constants.o mod_date.o mod_message.o mod_dynparam.o mod_kdinterp.o mod_stdio.o
mod_cmip6_cesm.o : mod_cmip6_cesm.F90 mod_cmip6_helper.o mod_realkinds.o mod_intkinds.o mod_message.o mod_date.o mod_stdio.o mod_dynparam.o mod_memutil.o
mod_cmip6_cnrm.o : mod_cmip6_cnrm.F90 mod_cmip6_helper.o mod_realkinds.o mod_intkinds.o mod_message.o mod_date.o mod_stdio.o mod_dynparam.o mod_memutil.o
mod_cmip6_ecea.o : mod_cmip6_ecea.F90 mod_cmip6_helper.o mod_realkinds.o mod_intkinds.o mod_message.o mod_date.o mod_stdio.o mod_dynparam.o mod_memutil.o
mod_cmip6_hadmm.o : mod_cmip6_hadmm.F90 mod_cmip6_helper.o mod_realkinds.o mod_intkinds.o mod_message.o mod_date.o mod_stdio.o mod_dynparam.o mod_memutil.o
mod_cmip6_miroc6.o : mod_cmip6_miroc6.F90 mod_cmip6_helper.o mod_realkinds.o mod_intkinds.o mod_message.o mod_date.o mod_stdio.o mod_dynparam.o mod_memutil.o
mod_cmip6_mpihr.o : mod_cmip6_mpihr.F90 mod_cmip6_helper.o mod_realkinds.o mod_intkinds.o mod_message.o mod_date.o mod_stdio.o mod_dynparam.o mod_memutil.o
mod_cmip6_normm.o : mod_cmip6_normm.F90 mod_cmip6_helper.o mod_realkinds.o mod_intkinds.o mod_message.o mod_date.o mod_stdio.o mod_dynparam.o mod_memutil.o
