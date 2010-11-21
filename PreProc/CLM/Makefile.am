
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/Config:$(REGCM_ROOT)/PreProc/ICBC:.
vpath %.f90 $(REGCM_ROOT)/Config:$(REGCM_ROOT)/PreProc/ICBC:.

bin_PROGRAMS = clm2rcm
clm2rcm_SOURCES = mod_constants.F90 mod_dynparam.F90 mod_param_clm.F90 \
                  mod_date.f90 mod_read_domain.f90 mod_nclib.f90       \
                  mod_clm3grid.f90 clm2rcm.F90

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

mod_constants.o: mod_constants.F90
mod_dynparam.o: mod_constants.o mod_constants.F90
clm2rcm.o: clm2rcm.F90 mod_dynparam.o mod_read_domain.o mod_date.o \
           mod_clm3grid.o
mod_clm3grid.o: mod_clm3grid.f90
mod_param_clm.o: mod_param_clm.f90
mod_read_domain.o: mod_read_domain.f90
mod_nclib.o: mod_nclib.f90
