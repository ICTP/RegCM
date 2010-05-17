# Makefile for RegCM4
# for Terrain code.
  
include ../../Makefile.inc


all: terrain install

SHRMODS =  mod_dynparam.o mod_constants.o

MODOBJS =  mod_block.o     mod_interp.o    mod_maps.o        mod_rdldtr.o \
           mod_smooth.o    mod_fudge.o     mod_interfaces.o  mod_projections.o 

OBJS =     output.o        surf.o          terrain.o  \
           xyobsll.o       setup.o         header_terrain.o

mod_dynparam.o: ../../Config/mod_dynparam.F90
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_constants.o: ../../Config/mod_constants.F90
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

terrain: $(SHRMODS) $(MODOBJS) $(OBJS)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SHRMODS) $(MODOBJS) $(OBJS) \
               $(LIB) $(NETCDFLIB)

install: terrain 
	@cp terrain $(BINDIR)
	@chmod 755  $(BINDIR)/terrain
	@echo "Executable terrain is installed now in $(BINDIR)"

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

header_terrain.o: header_terrain.F90 mod_interfaces.o
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
  
clean:
	rm -f *.mod *.o terrain

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_block.o: mod_block.f90
mod_fudge.o: mod_fudge.f90
mod_interfaces.o: mod_interfaces.f90
mod_maps.o: mod_maps.f90
mod_smooth.o: mod_smooth.f90
mod_projections.o: mod_projections.f90 mod_constants.o
mod_interp.o: mod_interp.f90 mod_block.o mod_maps.o
output.o: output.f90 mod_interfaces.o
mod_rdldtr.o: mod_rdldtr.f90 mod_block.o
setup.o: setup.f90 mod_block.o mod_interfaces.o
surf.o: surf.f90 mod_block.o mod_interp.o mod_interfaces.o
xyobsll.o: xyobsll.f90 mod_block.o mod_constants.o mod_projections.o \
           mod_interfaces.o
terrain.o: terrain.f90 mod_dynparam.o mod_projections.o \
           mod_block.o mod_maps.o mod_interp.o mod_rdldtr.o mod_interfaces.o
