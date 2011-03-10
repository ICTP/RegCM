
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/Config:$(REGCM_ROOT)/PreProc/ICBC:.
vpath %.f90 $(REGCM_ROOT)/Config:$(REGCM_ROOT)/PreProc/ICBC:.

bin_PROGRAMS = ReadRegcmio
ReadRegcmio_SOURCES = mod_date.f90 ReadRegcmio.F90

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

ReadRegcmio.o: mod_date.o
