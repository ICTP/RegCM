
include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:.

CPPFLAGS += -DMAXPATCH_PFT=4

if OSNAME_AIX
  cpre = $(null)-WF,-D$(null)
  FPPFLAGS := $(patsubst -D%,$(cpre)%,$(CPPFLAGS))
else
  FPPFLAGS = $(CPPFLAGS)
endif

bin_PROGRAMS = clm2rcm
clm2rcm_SOURCES = mod_param_clm.F90 mod_grid.F90 mod_nclib.f90  \
                  mod_clm3grid.f90 clm2rcm.F90

%.o: %.f90
	$(FC) $(FPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
clm2rcm.o : clm2rcm.F90 mod_intkinds.o mod_realkinds.o mod_nclib.o mod_dynparam.o mod_message.o mod_grid.o mod_param_clm.o mod_date.o mod_clm3grid.o mod_memutil.o mod_stdio.o mod_domain.o mod_nchelper.o
mod_clm3grid.o : mod_clm3grid.f90 mod_intkinds.o mod_realkinds.o mod_stdio.o
mod_grid.o : mod_grid.f90 mod_realkinds.o mod_dynparam.o mod_memutil.o
mod_nclib.o : mod_nclib.f90 mod_intkinds.o mod_realkinds.o mod_memutil.o mod_message.o mod_stdio.o
mod_param_clm.o : mod_param_clm.f90 mod_intkinds.o mod_realkinds.o mod_message.o mod_stdio.o
