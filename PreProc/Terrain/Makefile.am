
include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:.

bin_PROGRAMS = terrain
terrain_SOURCES = mod_header.F90 mod_block.F90 mod_fudge.F90 mod_maps.F90 \
                mod_write.F90 mod_smooth.F90 mod_maputils.F90 mod_intldtr.F90 \
                mod_rdldtr.F90 mod_snow.F90 terrain.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_block.o : mod_block.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_message.o mod_constants.o
mod_fudge.o : mod_fudge.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_message.o
mod_header.o : mod_header.F90 mod_intkinds.o mod_stdio.o
mod_intldtr.o : mod_intldtr.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o
mod_maps.o : mod_maps.F90 mod_intkinds.o mod_realkinds.o mod_memutil.o
mod_maputils.o : mod_maputils.F90 mod_constants.o mod_projections.o mod_projections.o mod_projections.o mod_projections.o
mod_rdldtr.o : mod_rdldtr.F90 mod_block.o mod_stdio.o mod_constants.o mod_memutil.o mod_nchelper.o
mod_smooth.o : mod_smooth.F90 mod_intkinds.o mod_realkinds.o mod_constants.o
mod_write.o : mod_write.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_ncstream_types.o mod_ncstream.o
mod_snow.o: mod_snow.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_interp.o mod_stdio.o
terrain.o : terrain.F90 mod_realkinds.o mod_stdio.o mod_constants.o mod_message.o mod_date.o mod_dynparam.o mod_space.o mod_memutil.o mod_maps.o mod_block.o mod_smooth.o mod_projections.o mod_maputils.o mod_intldtr.o mod_fudge.o mod_nchelper.o mod_rdldtr.o mod_ncstream_types.o mod_ncstream.o mod_write.o mod_header.o mod_interp.o mod_snow.o
