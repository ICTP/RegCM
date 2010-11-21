
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/Config:.
vpath %.f90 $(REGCM_ROOT)/Config:.

bin_PROGRAMS = terrain
terrain_SOURCES = mod_constants.F90 mod_dynparam.F90 mod_header.F90      \
                  mod_block.f90 mod_fudge.f90 mod_maps.f90 mod_write.F90 \
                  mod_smooth.f90 mod_projections.f90 mod_maputils.f90    \
                  mod_interp.f90 mod_rdldtr.f90 terrain.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean:
		@echo -n Removing objects, mods and $(bin_PROGRAMS) ... 
		@rm -f *.o *.mod $(bin_PROGRAMS)
		@echo Done

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_constants.o: mod_constants.F90
mod_dynparam.o: mod_dynparam.F90 mod_constants.o
mod_block.o: mod_block.f90
mod_fudge.o: mod_fudge.f90
mod_maps.o: mod_maps.f90
mod_write.o: mod_write.F90 mod_block.o mod_maps.o mod_dynparam.o
mod_smooth.o: mod_smooth.f90
mod_projections.o: mod_projections.f90 mod_constants.o
mod_maputils.o: mod_maputils.f90 mod_constants.o mod_block.o mod_maps.o \
                mod_projections.o
mod_interp.o: mod_interp.f90 mod_block.o mod_maps.o
mod_rdldtr.o: mod_rdldtr.f90 mod_block.o
terrain.o: terrain.f90 mod_dynparam.o mod_projections.o \
           mod_block.o mod_maps.o mod_interp.o mod_rdldtr.o \
           mod_write.o mod_header.o
