
vpath %.inc .

include $(REGCM_ROOT)/makeinc

noinst_LIBRARIES = librcmlib.a
librcmlib_a_SOURCES = mod_constants.F90 mod_date.F90 mod_space.F90 mod_hgt.F90 \
                      mod_message.F90 mod_memutil.F90 mod_projections.F90 \
                      mod_dynparam.F90 mod_ensemble.F90 mod_stdio.F90 \
                      mod_humid.F90 mod_realkinds.F90 mod_interp.F90 \
                      mod_nchelper.F90 mod_vertint.F90 mod_spline.F90 \
                      mod_domain.F90 mod_intkinds.F90 mod_ncstream_types.F90 \
                      mod_ncstream.F90 mod_sunorbit.F90 mod_nhinterp.F90 \
                      mod_sigma.F90 mod_capecin.F90 mod_stdatm.F90 \
                      mod_rdldtr.F90 mod_intldtr.F90 mod_kdtree2_priority.F90 \
                      mod_kdtree2.F90 mod_earth.F90 mod_kdinterp.F90 \
                      mod_spbarcoord.F90 mod_posix.F90 mod_zita.F90 \
                      mod_isaatm.F90 mod_hash.F90 mod_sort.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_capecin.o : mod_capecin.F90 mod_realkinds.o mod_intkinds.o mod_constants.o
mod_constants.o: mod_constants.F90 mod_realkinds.o
mod_date.o: mod_date.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_message.o
mod_domain.o: mod_domain.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_nchelper.o mod_stdio.o mod_message.o
mod_dynparam.o: mod_dynparam.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_date.o
mod_ensemble.o: mod_ensemble.F90 mod_intkinds.o mod_realkinds.o mod_memutil.o mod_constants.o
mod_hgt.o: mod_hgt.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_message.o
mod_humid.o : mod_humid.F90 mod_intkinds.o mod_realkinds.o mod_constants.o pfesat.inc pfwsat.inc sig2p.inc
mod_interp.o: mod_interp.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o mod_memutil.o mod_constants.o mod_dynparam.o mod_earth.o
mod_intldtr.o : mod_intldtr.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o mod_earth.o
mod_memutil.o: mod_memutil.F90 mod_intkinds.o mod_realkinds.o mod_space.o mod_message.o mod_constants.o mod_date.o
mod_message.o: mod_message.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o
mod_nchelper.o: mod_nchelper.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_constants.o mod_memutil.o mod_dynparam.o mod_message.o $(REGCM_ROOT)/version
mod_nhinterp.o : mod_nhinterp.F90 mod_constants.o mod_intkinds.o mod_realkinds.o mod_stdio.o mod_stdatm.o
mod_projections.o: mod_projections.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_memutil.o
mod_realkinds.o: mod_realkinds.F90
mod_intkinds.o: mod_intkinds.F90
mod_space.o: mod_space.F90 mod_intkinds.o mod_realkinds.o mod_date.o
mod_spline.o: mod_spline.F90 mod_intkinds.o mod_realkinds.o mod_constants.o
mod_stdio.o: mod_stdio.F90 mod_intkinds.o
mod_vertint.o: mod_vertint.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_constants.o mod_message.o mod_interp.o mod_stdatm.o
mod_ncstream.o: mod_ncstream.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_constants.o mod_dynparam.o mod_message.o mod_date.o mod_ncstream_types.o mod_dynparam.o mod_zita.o $(REGCM_ROOT)/version
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
mod_ncstream_types.o: mod_ncstream_types.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_date.o
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
mod_posix.o : mod_posix.F90 mod_stdio.o
mod_rdldtr.o : mod_rdldtr.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_constants.o mod_memutil.o mod_nchelper.o mod_earth.o mod_message.o mod_kdinterp.o
mod_sigma.o : mod_sigma.F90 mod_constants.o mod_intkinds.o mod_memutil.o mod_message.o mod_realkinds.o mod_stdio.o
mod_sunorbit.o : mod_sunorbit.F90 mod_constants.o mod_message.o mod_intkinds.o mod_realkinds.o mod_stdio.o mod_dynparam.o
mod_stdatm.o : mod_stdatm.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o
mod_kdtree2_priority.o : mod_kdtree2_priority.F90 mod_realkinds.o mod_message.o
mod_kdtree2.o : mod_kdtree2.F90 mod_realkinds.o mod_message.o mod_stdio.o mod_kdtree2_priority.o
mod_earth.o : mod_earth.F90 mod_realkinds.o mod_intkinds.o mod_constants.o
mod_spbarcoord.o : mod_spbarcoord.F90 mod_realkinds.o mod_intkinds.o mod_constants.o
mod_kdinterp.o : mod_kdinterp.F90 mod_realkinds.o mod_intkinds.o mod_constants.o mod_stdio.o mod_message.o mod_earth.o mod_kdtree2.o mod_spbarcoord.o mod_dynparam.o
	$(FC) $(CPPFLAGS) $(GITVER) $(OMPFLAGS) $(FCFLAGS) -c $<
mod_zita.o : mod_zita.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o
mod_hash.o : mod_hash.F90 mod_intkinds.o
mod_isaatm.o : mod_isaatm.F90 mod_realkinds.o mod_constants.o
mod_sort.o : mod_sort.F90 mod_realkinds.o mod_constants.o
