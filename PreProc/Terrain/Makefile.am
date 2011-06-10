
include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = terrain
terrain_SOURCES = mod_header.F90 mod_block.f90 mod_fudge.f90 mod_maps.f90 \
									mod_write.F90 mod_smooth.f90 mod_maputils.f90 mod_interp.f90 \
									mod_rdldtr.f90 terrain.f90

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

mod_block.o: mod_block.f90
mod_fudge.o: mod_fudge.f90
mod_maps.o: mod_maps.f90
mod_write.o: mod_write.F90 mod_block.o mod_maps.o
mod_smooth.o: mod_smooth.f90
mod_maputils.o: mod_maputils.f90 mod_block.o mod_maps.o
mod_interp.o: mod_interp.f90 mod_block.o mod_maps.o
mod_rdldtr.o: mod_rdldtr.f90 mod_block.o
terrain.o: terrain.f90 mod_block.o mod_maps.o mod_interp.o \
	         mod_rdldtr.o mod_write.o mod_header.o
