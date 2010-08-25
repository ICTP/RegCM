# Makefile for RegCM4
  
include ../../../Makefile.inc

# CPPFLAGS += -DNETCDF4_HDF5

all: GrADSNcPlot

OBJS = mod_projections.o mod_constants.o

GrADSNcPlot: GrADSNcPlot.f90 $(OBJS)
	$(F90) $(LDFLAGS) $(NETCDFINC) $(F90FLAGS) -o $@ $< $(OBJS) \
               $(NETCDFLIB)

mod_constants.o: ../../../Config/mod_constants.F90
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_projections.o: ../../../PreProc/Terrain/mod_projections.f90 mod_constants.o
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c -o $@ $<

install: GrADSNcPlot
	@cp GrADSNcPlot $(BINDIR)
	@chmod 755  $(BINDIR)/GrADSNcPlot
	@echo "Executable GrADSNcPlot is installed now in $(BINDIR)"

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

%.o: %.F90
	$(F90) $(CPPFLAGS) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o GrADSNcPlot

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
