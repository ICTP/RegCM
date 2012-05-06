
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/PreProc/Terrain:.
vpath %.f90 $(REGCM_ROOT)/PreProc/Terrain:.

bin_PROGRAMS = sigma2p GrADSNcPlot GrADSNcPrepare average regrid
sigma2p_SOURCES = sigma2p.F90 mod_derived.f90
GrADSNcPlot_SOURCES = GrADSNcPlot.F90
GrADSNcPrepare_SOURCES = GrADSNcPrepare.F90

average: average.sh
	cp $< $@

regrid: regrid.sh
	cp $< $@

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
		rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_derived.o: mod_derived.f90
sigma2p.o: sigma2p.F90 mod_derived.o
GrADSNcPlot.o: GrADSNcPlot.F90
GrADSNcPrepare.o: GrADSNcPrepare.F90
