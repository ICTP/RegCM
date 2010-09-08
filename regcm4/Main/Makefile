# Makefile for Main regcm Model

###########################################################################
# CONFIGURABLE PATAMETERS
###########################################################################

SCENARIO = A1B

###########################################################################
# END CONFIGURABLE PARAMETERS
###########################################################################

include ../Makefile.inc

destdir = $(BINDIR)
prgsuffix = Serial

# check if we need the MPI wrapper 
ifneq ($(origin PARALLEL), undefined)
  F90 = $(MPIF90)
  prgsuffix = MPI
endif  

ifeq ($(USECLM),1)
CLMPTH = ./clmlib
CLMLIB = $(CLMPTH)/libclm.a
all:: $(CLMLIB) program
prgsuffix = _clM
else
CLMLIB =
all:: program
endif

SHROBJS = mod_constants.o mod_dynparam.o mod_bats_param.o

MODOBJS = \
mod_mppio.o      mod_crdcae.o     mod_mainchem.o   mod_vecbats.o     \
mod_aerosol.o    mod_main.o       mod_rad.o        mod_bats.o        \
mod_vmodes.o     mod_message.o    mod_slice.o      mod_bdycod.o      \
mod_cvaria.o     mod_o3blk.o      mod_split.o      mod_holtbl.o      \
mod_date.o       mod_outrad.o     mod_cu_em.o      mod_cu_bm.o       \
mod_diagnosis.o  mod_runparams.o  mod_bndry.o      mod_precip.o      \
mod_dust.o       mod_scenarios.o  mod_colmod3.o    mod_tracer.o      \
mod_comozp.o     mod_pmoist.o     mod_trachem.o    mod_sun.o         \
mod_output.o     mod_pbldim.o     mod_lake.o       mod_drag.o        \
mod_zengocn.o    mod_ncio.o       mod_cu_grell.o   mod_radiation.o   \
mod_cu_kuo.o     mod_savefile.o   mod_leaftemp.o   mod_aero_sett_ddep.o

CLMOBJS = \
initializeMod.o  shr_orb_mod.o    shr_kind_mod.o  clm_varpar.o       \
clm_varsur.o     atmdrvMod.o      spmdMod.o       program_offMod.o   \
clm_comp.o       clmtype.o        perf_mod.o      clm_time_manager.o \
restFileMod.o

OBJS =       tractend2.o  tend.o       param.o      diffu.o      init.o       \
diffut.o     cumtran.o    bdyin.o      bdyuv.o      bdyval.o     htdiff.o     \
radout.o     tstep.o      vadv.o       spinit.o     conmas.o     nudge.o      \
regcm.o      inirad.o     condtq.o     sponge.o     tracbud.o    header.o     \
chsrfem.o    getdat.o     slice1D.o    tracdiag.o   slice3D.o    conadv.o     \
hadv.o       cldfrac.o

NETLIBOBJS = linpack.o eispack.o

EXTRAOBJS = mod_culookup.o mod_co2.o

ifeq ($(USECLM),1)
#  CPPFLAGS += -DCLM
  MODOBJS += mod_clm.o mod_mtrxclm.o
  LCLMLIB = -L$(CLMPTH) -lclm
else
ifeq ($(USESEAICE),1)
#  CPPFLAGS += -DSEAICE
  seaicesuffix = _si
endif
endif
ifeq ($(USEDCSST),1)
#  CPPFLAGS += -DDCSST
  dcsstsuffix = _ds
endif
ifeq ($(USEBAND),1)
#  CPPFLAGS += -DDCSST
  bandsuffix = _band
endif

program: regcm
	@cp regcm \
          $(destdir)/regcm$(prgsuffix)$(dcsstsuffix)$(seaicesuffix)$(bandsuffix)
	@chmod 755  $(destdir)/regcm$(prgsuffix)$(dcsstsuffix)$(seaicesuffix)$(bandsuffix) 
	@echo "Executable is installed now in $(destdir)"

$(CLMPTH)/libclm.a: mod_dynparam.o mod_clm.o
ifeq ($(origin PARALLEL), undefined)
	@echo "ERROR: CLM wants PARALLEL activated"
	@exit 1
endif
	$(MAKE) -C $(CLMPTH)

regcm: $(SHROBJS) $(CLMLIB) $(MODOBJS) $(OBJS) $(NETLIBOBJS)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(MODOBJS) $(OBJS) $(SHROBJS) \
               $(NETLIBOBJS) $(LCLMLIB) $(NETCDFLIB)

mod_constants.o: ../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS)  -c -o $@ $<

mod_dynparam.o: ../Config/mod_dynparam.F90 mod_constants.o
	$(F90) $(CPPFLAGS) $(F90FLAGS)  -c -o $@ $<

mod_bats_param.o: ../Config/mod_bats_param.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS)  -c -o $@ $<

%.o: %.F90
	$(F90) -c $(NETCDFINC) $(CPPFLAGS) $(F90FLAGS)  -c $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS)  -c $<

ifeq ($(USECLM),1)
clean:
	rm -f *.mod *.o regcm
	cd clmlib && $(MAKE) clean
else
clean:
	rm -f *.mod *.o regcm
endif

mod_aero_sett_ddep.o: mod_aero_sett_ddep.F90 mod_dynparam.o mod_message.o \
                      mod_dust.o
mod_aerosol.o: mod_aerosol.F90 mod_constants.o mod_main.o mod_mainchem.o \
               mod_message.o mod_runparams.o mod_dynparam.o mod_trachem.o \
               mod_slice.o
mod_bats.o: mod_bats.F90 mod_dynparam.o mod_bats_param.o mod_constants.o \
            mod_date.o mod_main.o mod_runparams.o mod_pbldim.o mod_slice.o
mod_bdycod.o: mod_bdycod.F90 mod_dynparam.o
mod_bndry.o: mod_bndry.F90 mod_constants.o mod_bats.o mod_leaftemp.o \
             mod_runparams.o mod_dynparam.o mod_drag.o
mod_clm.o: mod_clm.F90 mod_dynparam.o
mod_colmod3.o: mod_colmod3.F90 mod_constants.o mod_bats.o \
               mod_date.o mod_dynparam.o mod_radiation.o
mod_comozp.o: mod_comozp.f90
mod_crdcae.o: mod_crdcae.f90
mod_culookup.o: mod_culookup.f90 mod_constants.o
mod_cu_bm.o: mod_cu_bm.F90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
             mod_runparams.o mod_pmoist.o mod_rad.o mod_dynparam.o mod_trachem.o
mod_cu_em.o: mod_cu_em.f90 mod_bats.o mod_cvaria.o mod_date.o mod_main.o \
             mod_runparams.o mod_pmoist.o mod_rad.o mod_dynparam.o mod_slice.o
mod_cu_grell.o: mod_cu_grell.f90 mod_pmoist.o mod_rad.o mod_dynparam.o \
                mod_trachem.o mod_bats.o mod_constants.o mod_date.o mod_main.o \
                mod_runparams.o
mod_cu_kuo.o: mod_cu_kuo.f90 mod_bats.o mod_constants.o mod_cvaria.o \
              mod_date.o mod_main.o mod_runparams.o mod_pmoist.o mod_rad.o \
              mod_dynparam.o mod_trachem.o
mod_cvaria.o: mod_cvaria.F90 mod_dynparam.o
mod_date.o: mod_date.F90 mod_message.o mod_runparams.o mod_dynparam.o
mod_diagnosis.o: mod_diagnosis.f90 mod_dynparam.o
mod_drag.o: mod_drag.F90 mod_constants.o mod_bats.o mod_dynparam.o
mod_dust.o: mod_dust.F90 mod_trachem.o mod_dynparam.o mod_constants.o \
            mod_trachem.o mod_ncio.o
mod_holtbl.o: mod_holtbl.F90 mod_bats.o mod_constants.o mod_cvaria.o \
          mod_diagnosis.o mod_main.o mod_mainchem.o mod_runparams.o \
          mod_pbldim.o mod_pmoist.o mod_dynparam.o mod_slice.o \
          mod_trachem.o mod_mppio.o
mod_leaftemp.o: mod_leaftemp.f90 mod_constants.o mod_bats.o \
                mod_runparams.o mod_dynparam.o
mod_lake.o: mod_lake.F90 mod_bats.o mod_date.o mod_runparams.o mod_dynparam.o
mod_main.o: mod_main.F90 mod_dynparam.o
mod_mainchem.o: mod_mainchem.F90 mod_dynparam.o
mod_message.o: mod_message.F90 mod_dynparam.o
mod_mppio.o: mod_mppio.F90 mod_dynparam.o mod_dust.o
mod_ncio.o: mod_ncio.F90 mod_dynparam.o mod_message.o mod_constants.o \
            mod_date.o mod_runparams.o
	$(F90) $(NETCDFINC) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS)   -c $<
mod_o3blk.o: mod_o3blk.F90 mod_dynparam.o mod_main.o mod_runparams.o mod_rad.o
mod_outrad.o: mod_outrad.F90 mod_dynparam.o
ifeq ($(USECLM),1)
mod_output.o: mod_output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_runparams.o mod_pmoist.o mod_rad.o mod_diagnosis.o \
          mod_radiation.o mod_dynparam.o mod_split.o \
          mod_trachem.o mod_mppio.o mod_ncio.o mod_lake.o mod_savefile.o \
          $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_output.o: mod_output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_runparams.o mod_pmoist.o mod_rad.o mod_diagnosis.o \
          mod_radiation.o mod_dynparam.o mod_split.o \
          mod_trachem.o mod_mppio.o mod_lake.o mod_ncio.o mod_savefile.o
endif
mod_pbldim.o: mod_pbldim.F90 mod_dynparam.o
mod_pmoist.o: mod_pmoist.F90 mod_dynparam.o
mod_precip.o: mod_precip.f90 mod_bats.o mod_constants.o mod_cvaria.o \
              mod_date.o mod_main.o mod_runparams.o mod_pmoist.o \
              mod_dynparam.o mod_slice.o mod_trachem.o
mod_rad.o: mod_rad.F90 mod_dynparam.o
ifeq ($(USECLM),1)
mod_radiation.o: mod_radiation.F90 mod_comozp.o mod_crdcae.o mod_constants.o \
                 mod_date.o mod_message.o mod_dynparam.o mod_tracer.o \
                 mod_scenarios.o mod_aerosol.o mod_bats.o \
                 mod_runparams.o mod_clm.o
else
mod_radiation.o: mod_radiation.F90 mod_comozp.o mod_crdcae.o mod_constants.o \
                 mod_date.o mod_message.o mod_dynparam.o mod_tracer.o \
                 mod_scenarios.o mod_aerosol.o mod_bats.o \
                 mod_runparams.o
endif
mod_runparams.o: mod_runparams.f90 mod_dynparam.o
ifeq ($(USECLM),1)
mod_savefile.o: mod_savefile.F90 mod_dynparam.o mod_message.o mod_runparams.o \
                mod_bats.o mod_pmoist.o mod_main.o mod_mainchem.o mod_bdycod.o \
                mod_rad.o mod_trachem.o mod_date.o mod_radiation.o \
                mod_diagnosis.o mod_mppio.o mod_clm.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_savefile.o: mod_savefile.F90 mod_dynparam.o mod_message.o mod_runparams.o \
                mod_bats.o mod_pmoist.o mod_main.o mod_mainchem.o mod_bdycod.o \
                mod_rad.o mod_trachem.o mod_date.o mod_radiation.o \
                mod_diagnosis.o mod_mppio.o
endif
mod_slice.o: mod_slice.F90 mod_dynparam.o
mod_split.o: mod_split.F90 mod_dynparam.o mod_main.o mod_runparams.o \
             mod_constants.o
ifeq ($(USECLM),1)
mod_sun.o: mod_sun.F90 $(CLMLIB) mod_clm.o mod_message.o mod_date.o \
             mod_constants.o mod_dynparam.o mod_main.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_sun.o: mod_sun.F90 mod_constants.o mod_date.o mod_message.o \
           mod_dynparam.o mod_main.o
endif
mod_scenarios.o: mod_scenarios.F90
	$(F90) -c $(CPPFLAGS) -D$(SCENARIO) $(F90FLAGS)  -c $<
mod_tracer.o: mod_tracer.f90
mod_trachem.o: mod_trachem.F90 mod_dynparam.o mod_message.o
mod_vecbats.o: mod_vecbats.F90 mod_date.o mod_runparams.o mod_bats.o \
               mod_lake.o mod_main.o mod_bndry.o mod_drag.o
mod_vmodes.o: mod_vmodes.F90 mod_runparams.o mod_dynparam.o mod_split.o \
              mod_constants.o mod_message.o
ifeq ($(USECLM),1)
mod_zengocn.o: mod_zengocn.F90 mod_bats.o mod_constants.o mod_date.o \
               mod_main.o mod_runparams.o mod_pbldim.o \
               mod_slice.o mod_bdycod.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_zengocn.o: mod_zengocn.F90 mod_bats.o mod_constants.o mod_date.o \
               mod_main.o mod_runparams.o mod_pbldim.o \
               mod_slice.o mod_dynparam.o mod_bdycod.o
endif
#
bdyin.o: bdyin.F90 mod_bats.o mod_bdycod.o mod_date.o mod_main.o \
         mod_runparams.o mod_dynparam.o mod_mppio.o \
         mod_ncio.o
bdyuv.o: bdyuv.F90 mod_bdycod.o mod_cvaria.o mod_main.o mod_dynparam.o
bdyval.o: bdyval.F90 mod_bdycod.o mod_date.o mod_main.o mod_mainchem.o \
          mod_runparams.o mod_dynparam.o
ifeq ($(USECLM),1)
chsrfem.o: chsrfem.F90 mod_aerosol.o mod_dust.o mod_mainchem.o \
           mod_message.o mod_dynparam.o mod_trachem.o $(CLMLIB)
	$(F90) $(CPPFLAGS) $(NETCDFINC) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
chsrfem.o: chsrfem.F90 mod_aerosol.o mod_dust.o mod_mainchem.o \
           mod_message.o mod_dynparam.o mod_trachem.o
endif
cldfrac.o: cldfrac.f90 mod_pmoist.o mod_rad.o mod_dynparam.o mod_slice.o
co2.o: co2.f90 mod_bats.o mod_leaftemp.o mod_runparams.o mod_dynparam.o
conadv.o: conadv.F90 mod_constants.o mod_diagnosis.o mod_main.o \
          mod_runparams.o mod_dynparam.o
condtq.o: condtq.f90 mod_constants.o mod_cvaria.o mod_main.o mod_runparams.o \
          mod_pmoist.o mod_dynparam.o mod_slice.o
conmas.o: conmas.F90 mod_constants.o mod_date.o mod_diagnosis.o mod_main.o \
          mod_runparams.o mod_dynparam.o
cumtran.o: cumtran.F90 mod_mainchem.o mod_runparams.o mod_dynparam.o \
           mod_trachem.o
diffu.o: diffu.F90 mod_main.o mod_dynparam.o mod_slice.o
diffut.o: diffut.F90 mod_main.o mod_dynparam.o mod_slice.o
getdat.o: getdat.f90 mod_constants.o mod_date.o \
          mod_main.o mod_runparams.o mod_rad.o mod_dynparam.o
hadv.o: hadv.F90 mod_cvaria.o mod_main.o mod_dynparam.o
header.o: header.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS)   -c $<
htdiff.o: htdiff.F90 mod_pmoist.o mod_dynparam.o
ifport.o: ifport.f90
inirad.o: inirad.F90 mod_date.o mod_rad.o mod_dynparam.o mod_o3blk.o
ifeq ($(USECLM),1)
init.o: init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_main.o mod_mainchem.o mod_message.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o \
        mod_radiation.o mod_dynparam.o mod_mppio.o mod_trachem.o mod_lake.o \
        mod_ncio.o mod_sun.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
init.o: init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_main.o mod_mainchem.o mod_message.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o \
        mod_radiation.o mod_dynparam.o mod_mppio.o mod_trachem.o mod_ncio.o \
        mod_lake.o mod_sun.o
endif
nudge.o: nudge.F90 mod_bdycod.o mod_main.o mod_runparams.o mod_dynparam.o
param.o: param.F90 mod_constants.o mod_bats.o mod_cu_em.o mod_date.o \
         mod_main.o mod_message.o mod_runparams.o mod_o3blk.o \
         mod_pmoist.o mod_dynparam.o mod_trachem.o mod_mppio.o \
         mod_cu_bm.o mod_split.o mod_outrad.o mod_holtbl.o \
         mod_aero_sett_ddep.o mod_ncio.o mod_leaftemp.o \
         mod_aerosol.o mod_radiation.o mod_bdycod.o
radout.o: radout.F90 mod_bats.o mod_date.o mod_main.o mod_outrad.o \
          mod_runparams.o mod_rad.o mod_dynparam.o
outrad.o: outrad.F90 mod_date.o mod_outrad.o mod_dynparam.o mod_runparams.o \
          mod_main.o mod_mppio.o mod_ncio.o
ifeq ($(USECLM),1)
regcm.o: regcm.F90 mod_date.o mod_message.o mod_runparams.o \
         mod_dynparam.o mod_ncio.o mod_output.o mod_split.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
regcm.o: regcm.F90 mod_date.o mod_message.o mod_runparams.o \
         mod_dynparam.o mod_ncio.o mod_output.o
endif
slice1D.o: slice1D.f90 mod_bats.o mod_constants.o mod_main.o mod_pbldim.o \
           mod_dynparam.o
slice3D.o: slice3D.F90 mod_main.o mod_constants.o mod_mainchem.o \
           mod_runparams.o mod_pbldim.o mod_pmoist.o mod_dynparam.o mod_slice.o
spinit.o: spinit.F90 mod_bdycod.o mod_date.o mod_main.o \
          mod_runparams.o mod_dynparam.o mod_split.o mod_mppio.o \
          mod_constants.o mod_vmodes.o
sponge.o: sponge.F90 mod_bdycod.o mod_dynparam.o
ifeq ($(USECLM),1)
tend.o: tend.F90 mod_aerosol.o mod_bats.o mod_bdycod.o mod_constants.o \
        mod_cvaria.o mod_date.o mod_main.o mod_mainchem.o mod_message.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o mod_holtbl.o \
        mod_dynparam.o mod_slice.o mod_trachem.o mod_clm.o mod_cu_grell.o \
        mod_cu_kuo.o mod_cu_bm.o mod_cu_em.o mod_zengocn.o mod_sun.o \
        mod_precip.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
tend.o: tend.F90 mod_aerosol.o mod_bats.o mod_bdycod.o mod_constants.o \
        mod_cvaria.o mod_date.o mod_main.o mod_mainchem.o mod_message.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o mod_holtbl.o \
        mod_dynparam.o mod_slice.o mod_trachem.o mod_cu_grell.o mod_sun.o \
        mod_cu_kuo.o mod_cu_bm.o mod_cu_em.o mod_zengocn.o mod_precip.o
endif
tracbud.o: tracbud.F90 mod_constants.o mod_date.o mod_diagnosis.o mod_main.o \
           mod_mainchem.o mod_runparams.o mod_dynparam.o mod_trachem.o
tracdiag.o: tracdiag.F90 mod_constants.o mod_diagnosis.o mod_main.o \
            mod_mainchem.o mod_runparams.o mod_dynparam.o
tractend2.o: tractend2.f90 mod_constants.o mod_bats.o mod_cvaria.o mod_date.o \
             mod_dust.o mod_main.o mod_mainchem.o mod_message.o \
             mod_runparams.o mod_pbldim.o mod_pmoist.o mod_rad.o \
             mod_dynparam.o mod_trachem.o
tstep.o: tstep.f90 mod_date.o mod_runparams.o mod_dynparam.o
vadv.o: vadv.F90 mod_cvaria.o mod_main.o mod_runparams.o mod_dynparam.o

mod_mtrxclm.o: mod_mtrxclm.F90 mod_bats.o mod_dynparam.o mod_drag.o \
               mod_constants.o mod_runparams.o mod_mppio.o mod_date.o \
               mod_clm.o mod_main.o mod_pbldim.o mod_slice.o mod_date.o \
	$(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<

linpack.o: linpack.f90
eispack.o: eispack.f90

mod_co2.o: mod_co2.f90 mod_dynparam.o mod_runparams.o mod_bats.o mod_leaftemp.o
