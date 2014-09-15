
include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:.

FPPFLAGS = $(CPPFLAGS) $(DEFINE)MAXPATCH_PFT=4

if DO_COMPILE_LIBCLM
bin_PROGRAMS = clm2rcm
clm2rcm_SOURCES = mod_param_clm.F90 mod_grid.F90 mod_nclib.F90  \
                  mod_clm3grid.F90 clm2rcm.F90

%.o: %.F90
	$(FC) $(FPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<
endif

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
clm2rcm.o : clm2rcm.F90 mod_intkinds.o mod_realkinds.o mod_nclib.o mod_dynparam.o mod_message.o mod_grid.o mod_param_clm.o mod_date.o mod_clm3grid.o mod_memutil.o mod_stdio.o mod_domain.o mod_nchelper.o
mod_clm3grid.o : mod_clm3grid.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o
mod_grid.o : mod_grid.F90 mod_realkinds.o mod_dynparam.o mod_memutil.o
mod_nclib.o : mod_nclib.F90 mod_intkinds.o mod_realkinds.o mod_memutil.o mod_message.o mod_stdio.o mod_dynparam.o
mod_param_clm.o : mod_param_clm.F90 mod_intkinds.o mod_realkinds.o mod_message.o mod_stdio.o
