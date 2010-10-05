# Makefile for RegCM4
  
include ../../../Makefile.inc

# CPPFLAGS += -DNETCDF4_HDF5

all: GrADSNcPlot GrADSNcPrepare

OBJS = mod_projections.o mod_constants.o

GrADSNcPlot: GrADSNcPlot.F90 $(OBJS)
	$(F90) $(LDFLAGS) $(NETCDFINC) $(CPPFLAGS) $(F90FLAGS) -o $@ $< $(OBJS) \
               $(NETCDFLIB)

GrADSNcPrepare: GrADSNcPrepare.F90 $(OBJS)
	$(F90) $(LDFLAGS) $(NETCDFINC) $(CPPFLAGS) $(F90FLAGS) -o $@ $< $(OBJS) \
               $(NETCDFLIB)

mod_constants.o: ../../../Config/mod_constants.F90
	$(F90) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_projections.o: ../../../PreProc/Terrain/mod_projections.f90 mod_constants.o
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c -o $@ $<

install: GrADSNcPlot GrADSNcPrepare
	@cp GrADSNcPlot GrADSNcPrepare $(BINDIR)
	@chmod 755  $(BINDIR)/GrADSNcPlot GrADSNcPrepare
	@echo "Executables are installed now in $(BINDIR)"

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

%.o: %.F90
	$(F90) $(CPPFLAGS) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o GrADSNcPlot GrADSNcPrepare

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
