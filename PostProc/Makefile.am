
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/Config:$(REGCM_ROOT)/PreProc/Terrain:.
vpath %.f90 $(REGCM_ROOT)/Config:$(REGCM_ROOT)/PreProc/Terrain:$(REGCM_ROOT)/PreProc/ICBC:.

bin_PROGRAMS = sigma2p GrADSNcPlot GrADSNcPrepare average regrid
sigma2p_SOURCES = mod_constants.F90 sigma2p.F90
GrADSNcPlot_SOURCES = mod_constants.F90 mod_projections.f90 \
											mod_date.f90 GrADSNcPlot.F90
GrADSNcPrepare_SOURCES = mod_constants.F90 mod_projections.f90 \
												 mod_date.f90 GrADSNcPrepare.F90

average: average.sh
	cp $< $@

regrid: regrid.sh
	cp $< $@

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

mod_derived.o: mod_derived.f90 mod_constants.o
sigma2p.o: mod_constants.o sigma2p.F90
GrADSNcPlot.o: mod_constants.o mod_projections.o mod_date.o GrADSNcPlot.F90
GrADSNcPrepare.o: mod_constants.o mod_projections.o mod_date.o \
	                GrADSNcPrepare.F90
