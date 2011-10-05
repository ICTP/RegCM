
include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = terrain
terrain_SOURCES = mod_header.F90 mod_block.f90 mod_fudge.f90 mod_maps.f90 \
                  mod_write.F90 mod_smooth.f90 mod_maputils.f90 mod_interp.f90 \
                  mod_rdldtr.F90 terrain.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod

distclean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_block.o mod_block.mod: mod_block.f90 $(REGCM_ROOT)/Share/mod_stdio.mod \
 $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_message.mod \
 $(REGCM_ROOT)/Share/mod_constants.mod
mod_fudge.o mod_fudge.mod: mod_fudge.f90 $(REGCM_ROOT)/Share/mod_stdio.mod \
 $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_message.mod
mod_header.o mod_header.mod: mod_header.F90 $(REGCM_ROOT)/Share/mod_stdio.mod
mod_interp.o mod_interp.mod: mod_interp.f90 \
 $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_constants.mod \
 $(REGCM_ROOT)/Share/mod_stdio.mod $(REGCM_ROOT)/Share/mod_message.mod
mod_maps.o mod_maps.mod: mod_maps.f90 $(REGCM_ROOT)/Share/mod_realkinds.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod
mod_maputils.o mod_maputils.mod: mod_maputils.f90 \
 $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_projections.mod \
 $(REGCM_ROOT)/Share/mod_projections.mod $(REGCM_ROOT)/Share/mod_projections.mod \
 $(REGCM_ROOT)/Share/mod_projections.mod
mod_rdldtr.o mod_rdldtr.mod: mod_rdldtr.F90 mod_block.mod \
 $(REGCM_ROOT)/Share/mod_stdio.mod $(REGCM_ROOT)/Share/mod_constants.mod \
 $(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/Share/mod_message.mod
mod_smooth.o mod_smooth.mod: mod_smooth.f90 \
 $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_stdio.mod
mod_write.o mod_write.mod: mod_write.F90 $(REGCM_ROOT)/Share/mod_stdio.mod \
 $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_memutil.mod \
 $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_message.mod \
 mod_block.mod mod_maps.mod
terrain.o: terrain.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
 $(REGCM_ROOT)/Share/mod_constants.mod mod_maps.mod mod_block.mod \
 mod_smooth.mod mod_maputils.mod mod_interp.mod mod_fudge.mod mod_rdldtr.mod \
 mod_write.mod mod_header.mod $(REGCM_ROOT)/Share/mod_stdio.mod \
 $(REGCM_ROOT)/Share/mod_message.mod $(REGCM_ROOT)/Share/mod_memutil.mod
