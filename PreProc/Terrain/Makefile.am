
include $(REGCM_ROOT)/makeinc

vpath %.mod $(REGCM_ROOT)/Share:.

bin_PROGRAMS = terrain
terrain_SOURCES = mod_header.F90 mod_block.f90 mod_fudge.f90 mod_maps.f90 \
                mod_write.F90 mod_smooth.f90 mod_maputils.f90 mod_intldtr.f90 \
                mod_rdldtr.F90 terrain.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

%.mod: %.f90 %.o
	@true

%.mod: %.F90 %.o
	@true

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_block.o mod_block.mod: mod_block.f90 mod_stdio.mod mod_realkinds.mod mod_message.mod mod_constants.mod
mod_fudge.o mod_fudge.mod: mod_fudge.f90 mod_stdio.mod mod_realkinds.mod mod_message.mod
mod_header.o mod_header.mod: mod_header.F90 mod_stdio.mod
mod_intldtr.o mod_intldtr.mod: mod_intldtr.f90 mod_realkinds.mod mod_constants.mod mod_stdio.mod mod_message.mod
mod_maps.o mod_maps.mod: mod_maps.f90 mod_realkinds.mod mod_memutil.mod
mod_maputils.o mod_maputils.mod: mod_maputils.f90 mod_constants.mod mod_projections.mod mod_projections.mod mod_projections.mod mod_projections.mod
mod_rdldtr.o mod_rdldtr.mod: mod_rdldtr.F90 mod_block.mod mod_stdio.mod mod_constants.mod mod_memutil.mod mod_nchelper.mod
mod_smooth.o mod_smooth.mod: mod_smooth.f90 mod_realkinds.mod mod_stdio.mod
mod_write.o mod_write.mod: mod_write.F90 mod_realkinds.mod mod_dynparam.mod mod_nchelper.mod
terrain.o terrain.mod: terrain.f90 mod_dynparam.mod mod_constants.mod mod_maps.mod mod_block.mod mod_smooth.mod mod_maputils.mod mod_intldtr.mod mod_fudge.mod mod_rdldtr.mod mod_write.mod mod_header.mod mod_stdio.mod mod_message.mod mod_memutil.mod
