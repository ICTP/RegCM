
include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:.

if DO_COMPILE_LIBCLM45
bin_PROGRAMS = mksurfdata interpinic

mksurfdata_SOURCES = mksurfdata.F90 mod_mkabm.F90 mod_mklaisai.F90 \
  mod_mksoitex.F90 mod_mkfmax.F90 mod_mkpeatf.F90 mod_mkurban.F90 \
  mod_mkgdp.F90 mod_mkpft.F90 mod_mkwetland.F90 mod_grid.F90 \
  mod_mkglacier.F90 mod_mksoilcol.F90 mod_mkvocef.F90 mod_mkorganic.F90 \
  mod_mklake.F90 mod_mklightning.F90 mod_mkpopd.F90 mod_mkharvest.F90 \
  mod_mkdynpft.F90 mod_mkvic.F90 mod_mkch4topm.F90 mod_mksoilph.F90 \
  mod_mkq10soil.F90 mod_mkndep.F90

interpinic_SOURCES = interpinic.F90

LDADD = $(OMPFLAGS)

%.o: %.F90
	$(FC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

endif

clean-local:
	rm -f *.mod *__genmod.*

distclean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
interpinic.o : interpinic.F90 mod_constants.o mod_intkinds.o mod_message.o mod_memutil.o mod_realkinds.o mod_kdinterp.o mod_stdio.o
mksurfdata.o : mksurfdata.F90 mod_constants.o mod_date.o mod_domain.o mod_dynparam.o mod_grid.o mod_intkinds.o mod_memutil.o mod_message.o mod_zita.o mod_mkabm.o mod_mkfmax.o mod_mkgdp.o mod_mkglacier.o mod_mklaisai.o mod_mklake.o mod_mkorganic.o mod_mkpeatf.o mod_mkpft.o mod_mksoilcol.o mod_mksoitex.o mod_mkurban.o mod_mkvocef.o mod_mkwetland.o mod_nchelper.o mod_realkinds.o mod_stdio.o mod_mklightning.o mod_mkpopd.o mod_mkharvest.o mod_mkdynpft.o mod_mkvic.o mod_mkch4topm.o mod_mksoilph.o mod_mkq10soil.o mod_mkndep.o $(REGCM_ROOT)/makeinc $(REGCM_ROOT)/version
mod_grid.o : mod_grid.F90 mod_dynparam.o mod_intkinds.o mod_memutil.o mod_message.o mod_realkinds.o mod_stdio.o
mod_mkabm.o : mod_mkabm.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mkfmax.o : mod_mkfmax.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mkgdp.o : mod_mkgdp.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mkglacier.o : mod_mkglacier.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mklaisai.o : mod_mklaisai.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mklake.o : mod_mklake.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mkorganic.o : mod_mkorganic.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mkpeatf.o : mod_mkpeatf.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mkpft.o : mod_mkpft.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mksoilcol.o : mod_mksoilcol.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mksoitex.o : mod_mksoitex.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o
mod_mkurban.o : mod_mkurban.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_message.o mod_rdldtr.o mod_constants.o
mod_mkwetland.o : mod_mkwetland.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_message.o mod_rdldtr.o mod_constants.o
mod_mkvocef.o : mod_mkvocef.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_intldtr.o mod_constants.o
mod_mklightning.o : mod_mklightning.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o $(REGCM_ROOT)/makeinc
mod_mkpopd.o : mod_mkpopd.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o $(REGCM_ROOT)/makeinc
mod_mkharvest.o : mod_mkharvest.F90 mod_dynparam.o mod_grid.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o $(REGCM_ROOT)/makeinc
mod_mkdynpft.o : mod_mkdynpft.F90 mod_dynparam.o mod_grid.o mod_message.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_mkvic.o : mod_mkvic.F90 mod_dynparam.o mod_grid.o mod_message.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o $(REGCM_ROOT)/makeinc
mod_mkch4topm.o : mod_mkch4topm.F90 mod_dynparam.o mod_grid.o mod_message.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o $(REGCM_ROOT)/makeinc
mod_mksoilph.o : mod_mksoilph.F90 mod_dynparam.o mod_grid.o mod_message.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o $(REGCM_ROOT)/makeinc
mod_mkq10soil.o : mod_mkq10soil.F90 mod_dynparam.o mod_grid.o mod_message.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o $(REGCM_ROOT)/makeinc
mod_mkndep.o : mod_mkndep.F90 mod_dynparam.o mod_grid.o mod_message.o mod_intkinds.o mod_realkinds.o mod_rdldtr.o mod_constants.o $(REGCM_ROOT)/makeinc
