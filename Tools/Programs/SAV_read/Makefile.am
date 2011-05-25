
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/Config:$(REGCM_ROOT)/PreProc/ICBC:.
vpath %.f90 $(REGCM_ROOT)/Config:$(REGCM_ROOT)/PreProc/ICBC:.

bin_PROGRAMS = reSAV
reSAV_SOURCES = mod_constants.F90 mod_dynparam.F90 reSAV.F90

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

reSAV.o: mod_constants.o mod_dynparam.o
