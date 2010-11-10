# Makefile for RegCM4
# for Terrain code.
  
include ../../Makefile.inc

all: terrain install

SHRMODS =  mod_dynparam.o mod_constants.o

MODOBJS =  mod_block.o     mod_interp.o    mod_maps.o        mod_rdldtr.o \
           mod_smooth.o    mod_fudge.o     mod_projections.o mod_write.o  \
           mod_header.o    mod_maputils.o

OBJS =     terrain.o

mod_constants.o: ../../Config/mod_constants.F90
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_dynparam.o: ../../Config/mod_dynparam.F90 mod_constants.o
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

terrain: $(SHRMODS) $(MODOBJS) $(OBJS)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SHRMODS) $(MODOBJS) $(OBJS) \
               $(NETCDFLIB)

install: terrain 
	@cp terrain $(BINDIR)
	@chmod 755  $(BINDIR)/terrain
	@echo "Executable terrain is installed now in $(BINDIR)"

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

%.o: %.F90
	$(F90) $(CPPFLAGS) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

mod_header.o: mod_header.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<
  
clean:
	rm -f *.mod *.o terrain

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

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
