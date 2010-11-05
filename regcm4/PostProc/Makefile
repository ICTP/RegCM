# Makefile for RegCM4 for PostProc code.
  
include ../Makefile.inc

all: sigma2p GrADSNcPlot GrADSNcPrepare
	@cp sigma2p GrADSNcPlot GrADSNcPrepare $(BINDIR)
	@cp regrid.sh $(BINDIR)/regrid
	@cp average.sh $(BINDIR)/average
	@echo "Programs all installed in $(BINDIR)"

OBJS = mod_projections.o mod_constants.o

mod_constants.o: ../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
               $(USCORING) -c -o $@ $<

mod_projections.o: ../PreProc/Terrain/mod_projections.f90 mod_constants.o
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c -o $@ $<

GrADSNcPlot: GrADSNcPlot.F90 $(OBJS)
	$(F90) $(LDFLAGS) $(NETCDFINC) $(CPPFLAGS) $(F90FLAGS) -o $@ $< $(OBJS) \
			   $(NETCDFLIB)

GrADSNcPrepare: GrADSNcPrepare.F90 $(OBJS)
	$(F90) $(LDFLAGS) $(NETCDFINC) $(CPPFLAGS) $(F90FLAGS) -o $@ $< $(OBJS) \
			   $(NETCDFLIB)

mod_derived.o: mod_derived.f90 mod_constants.o
	$(F90) -c $(LDFLAGS) $(F90FLAGS) -o $@ $<

sigma2p: sigma2p.F90 mod_constants.o
	$(F90) $(CPPFLAGS) $(NETCDFINC) $(LDFLAGS) $(F90FLAGS) \
               -o $@ $< $(NETCDFLIB)

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

%.o: %.F90
	$(F90) $(CPPFLAGS) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f sigma2p GrADSNcPlot GrADSNcPrepare *.o *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

