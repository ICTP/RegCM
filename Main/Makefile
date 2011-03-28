# Makefile for Main regcm Model

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

SHROBJS = mod_constants.o mod_dynparam.o mod_bats_param.o mod_service.o

MODOBJS = \
mod_mppio.o      mod_param.o      mod_mainchem.o   mod_vecbats.o     \
mod_aerosol.o    mod_main.o       mod_rad.o        mod_bats.o        \
mod_vmodes.o     mod_message.o    mod_slice.o      mod_bdycod.o      \
mod_cvaria.o     mod_o3blk.o      mod_split.o      mod_holtbl.o      \
mod_date.o       mod_outrad.o     mod_cu_em.o      mod_cu_bm.o       \
mod_diagnosis.o  mod_runparams.o  mod_bndry.o      mod_precip.o      \
mod_dust.o       mod_scenarios.o  mod_colmod3.o    mod_tracer.o      \
mod_tstep.o      mod_pmoist.o     mod_trachem.o    mod_sun.o         \
mod_output.o     mod_pbldim.o     mod_lake.o       mod_drag.o        \
mod_zengocn.o    mod_ncio.o       mod_cu_grell.o   mod_radiation.o   \
mod_cu_kuo.o     mod_savefile.o   mod_leaftemp.o   mod_che_semdde.o  \
mod_cldfrac.o    mod_init.o       mod_cumtran.o    mod_header.o      \
mod_condtq.o     mod_diffusion.o  mod_advection.o  mod_nudge.o       \
mod_tendency.o   mod_che_tend.o   mod_chem.o

CLMOBJS = \
initializeMod.o  shr_orb_mod.o    shr_kind_mod.o  clm_varpar.o       \
clm_varsur.o     atmdrvMod.o      spmdMod.o       program_offMod.o   \
clm_comp.o       clmtype.o        perf_mod.o      clm_time_manager.o \
restFileMod.o

OBJS =  regcm.o

NETLIB = netlib/libnetlib.a

EXTRAOBJS = mod_co2.o

ifeq ($(USECLM),1)
#  CPPFLAGS += -DCLM
  MODOBJS += mod_clm.o mod_mtrxclm.o
  LCLMLIB = -L$(CLMPTH) -lclm
else
endif
ifeq ($(USEBAND),1)
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

$(NETLIB):
	$(MAKE) -C netlib

regcm: $(SHROBJS) $(CLMLIB) $(MODOBJS) $(OBJS) $(NETLIB)
	$(F90) $(F90FLAGS) $(LDFLAGS) -o $@ $(MODOBJS) $(OBJS) $(SHROBJS) \
               $(NETLIB) $(LCLMLIB) $(NETCDFLIB)

mod_constants.o: ../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS)  -c -o $@ $<

mod_dynparam.o: ../Config/mod_dynparam.F90 mod_constants.o
	$(F90) $(CPPFLAGS) $(F90FLAGS)  -c -o $@ $<

mod_bats_param.o: ../Config/mod_bats_param.f90
	$(F90) $(CPPFLAGS) $(F90FLAGS)  -c -o $@ $<

%.o: %.F90
	$(F90) -c $(NETCDFINC) $(CPPFLAGS) $(F90FLAGS)  -c $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS)  -c $<

ifeq ($(USECLM),1)
clean:
	rm -f *.mod *.o regcm
	$(MAKE) -C $(CLMPTH) clean
	$(MAKE) -C netlib clean
else
clean:
	rm -f *.mod *.o regcm
	$(MAKE) -C netlib clean
endif

#
# Modules
#
mod_advection.o: mod_advection.F90 mod_main.o mod_runparams.o mod_dynparam.o
mod_aerosol.o: mod_aerosol.F90 mod_constants.o mod_main.o mod_mainchem.o \
               mod_message.o mod_runparams.o mod_dynparam.o mod_trachem.o \
               mod_slice.o
mod_bats.o: mod_bats.F90 mod_dynparam.o mod_bats_param.o mod_runparams.o
mod_bdycod.o: mod_bdycod.F90 mod_dynparam.o mod_runparams.o mod_main.o \
              mod_bats.o mod_message.o mod_ncio.o mod_date.o mod_cvaria.o \
              mod_mppio.o mod_mainchem.o
mod_bndry.o: mod_bndry.F90 mod_constants.o mod_bats.o mod_leaftemp.o \
             mod_runparams.o mod_dynparam.o mod_drag.o
mod_chem.o: mod_chem.F90 mod_dynparam.o mod_date.o mod_message.o mod_mppio.o
ifeq ($(USECLM),1)
mod_che_semdde.o: mod_che_semdde.F90 mod_aerosol.o mod_dust.o mod_mainchem.o \
                  mod_message.o mod_dynparam.o mod_trachem.o mod_constants.o \
                  mod_ncio.o mod_mppio.o mod_bats_param.o $(CLMLIB)
	$(F90) $(CPPFLAGS) $(NETCDFINC) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_che_semdde.o: mod_che_semdde.F90 mod_aerosol.o mod_dust.o mod_mainchem.o \
                  mod_message.o mod_dynparam.o mod_trachem.o mod_constants.o \
                  mod_ncio.o mod_mppio.o mod_bats_param.o
endif
mod_che_tend.o: mod_che_tend.F90 mod_constants.o mod_bats.o mod_cvaria.o \
                mod_date.o mod_dust.o mod_main.o mod_mainchem.o mod_message.o \
                mod_runparams.o mod_pbldim.o mod_pmoist.o mod_rad.o \
                mod_dynparam.o mod_trachem.o mod_che_semdde.o mod_diffusion.o \
                mod_advection.o mod_diagnosis.o mod_mppio.o mod_slice.o
mod_cldfrac.o: mod_cldfrac.f90 mod_pmoist.o mod_rad.o mod_dynparam.o \
	             mod_slice.o mod_radiation.o
mod_clm.o: mod_clm.F90 mod_dynparam.o
mod_colmod3.o: mod_colmod3.F90 mod_constants.o mod_bats.o mod_runparams.o \
               mod_date.o mod_dynparam.o mod_radiation.o mod_main.o mod_rad.o \
               mod_outrad.o
mod_condtq.o: mod_condtq.f90 mod_constants.o mod_cvaria.o mod_main.o \
              mod_runparams.o mod_pmoist.o mod_dynparam.o mod_slice.o
mod_cumtran.o: mod_cumtran.F90 mod_mainchem.o mod_runparams.o mod_dynparam.o \
               mod_trachem.o
mod_cu_bm.o: mod_cu_bm.F90 mod_bats.o mod_constants.o mod_date.o mod_main.o \
             mod_runparams.o mod_pmoist.o mod_rad.o mod_dynparam.o mod_trachem.o
mod_cu_em.o: mod_cu_em.f90 mod_bats.o mod_cvaria.o mod_date.o mod_main.o \
             mod_runparams.o mod_pmoist.o mod_rad.o mod_dynparam.o mod_slice.o
mod_cu_grell.o: mod_cu_grell.F90 mod_pmoist.o mod_rad.o mod_dynparam.o \
                mod_trachem.o mod_bats.o mod_constants.o mod_date.o mod_main.o \
                mod_cvaria.o mod_runparams.o
mod_cu_kuo.o: mod_cu_kuo.F90 mod_bats.o mod_constants.o mod_cvaria.o \
              mod_date.o mod_main.o mod_runparams.o mod_pmoist.o mod_rad.o \
              mod_dynparam.o mod_trachem.o mod_advection.o
mod_cvaria.o: mod_cvaria.f90 mod_dynparam.o mod_main.o mod_runparams.o
mod_date.o: mod_date.F90 mod_message.o mod_runparams.o mod_dynparam.o
mod_diagnosis.o: mod_diagnosis.F90 mod_dynparam.o mod_constants.o \
                 mod_runparams.o mod_main.o mod_mppio.o mod_date.o \
                 mod_mainchem.o
mod_diffusion.o: mod_diffusion.F90 mod_dynparam.o mod_runparams.o
mod_drag.o: mod_drag.F90 mod_constants.o mod_bats.o mod_dynparam.o \
            mod_runparams.o
mod_dust.o: mod_dust.F90 mod_trachem.o mod_dynparam.o mod_constants.o \
            mod_trachem.o mod_ncio.o mod_mppio.o
mod_header.o: mod_header.F90 mod_date.o mod_constants.o
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS)   -c $<
mod_holtbl.o: mod_holtbl.F90 mod_bats.o mod_constants.o mod_cvaria.o \
          mod_diagnosis.o mod_main.o mod_mainchem.o mod_runparams.o \
          mod_pbldim.o mod_pmoist.o mod_dynparam.o mod_slice.o \
          mod_trachem.o mod_mppio.o
ifeq ($(USECLM),1)
mod_init.o: mod_init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_main.o mod_mainchem.o mod_message.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o mod_o3blk.o \
        mod_radiation.o mod_dynparam.o mod_mppio.o mod_trachem.o mod_lake.o \
        mod_ncio.o mod_sun.o mod_savefile.o mod_cu_bm.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_init.o: mod_init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_main.o mod_mainchem.o mod_message.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o mod_o3blk.o \
        mod_radiation.o mod_dynparam.o mod_mppio.o mod_trachem.o mod_ncio.o \
        mod_lake.o mod_sun.o mod_savefile.o
endif
mod_lake.o: mod_lake.F90 mod_bats.o mod_runparams.o mod_dynparam.o mod_mppio.o
mod_leaftemp.o: mod_leaftemp.f90 mod_constants.o mod_bats.o \
                mod_runparams.o mod_dynparam.o
mod_main.o: mod_main.f90 mod_dynparam.o mod_runparams.o
mod_mainchem.o: mod_mainchem.f90 mod_dynparam.o mod_runparams.o
mod_message.o: mod_message.F90 mod_dynparam.o
mod_mppio.o: mod_mppio.F90 mod_dynparam.o mod_message.o mod_runparams.o mod_main.o
mod_ncio.o: mod_ncio.F90 mod_dynparam.o mod_message.o mod_constants.o \
            mod_date.o mod_runparams.o mod_trachem.o
	$(F90) $(NETCDFINC) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS)   -c $<
mod_nudge.o: mod_nudge.F90 mod_bdycod.o mod_main.o mod_runparams.o \
             mod_dynparam.o
mod_o3blk.o: mod_o3blk.F90 mod_dynparam.o mod_main.o mod_runparams.o mod_rad.o
ifeq ($(USECLM),1)
mod_output.o: mod_output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_runparams.o mod_pmoist.o mod_rad.o mod_diagnosis.o \
          mod_radiation.o mod_dynparam.o mod_split.o mod_clm.o \
          mod_trachem.o mod_mppio.o mod_ncio.o mod_lake.o mod_savefile.o \
          mod_cu_bm.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_output.o: mod_output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_runparams.o mod_pmoist.o mod_rad.o mod_diagnosis.o \
          mod_radiation.o mod_dynparam.o mod_split.o mod_cu_bm.o \
          mod_trachem.o mod_mppio.o mod_lake.o mod_ncio.o mod_savefile.o
endif
ifeq ($(USECLM),1)
mod_outrad.o: mod_outrad.F90 mod_dynparam.o mod_runparams.o mod_bats.o \
              mod_rad.o mod_date.o mod_clm.o
else
mod_outrad.o: mod_outrad.F90 mod_dynparam.o mod_runparams.o mod_bats.o \
              mod_rad.o mod_date.o
endif
ifeq ($(USECLM),1)
mod_param.o: mod_param.F90 mod_constants.o mod_bats.o mod_cu_em.o mod_date.o \
             mod_main.o mod_message.o mod_runparams.o mod_o3blk.o \
             mod_pmoist.o mod_dynparam.o mod_trachem.o mod_mppio.o \
             mod_cu_bm.o mod_split.o mod_outrad.o mod_holtbl.o mod_chem.o \
             mod_che_semdde.o mod_ncio.o mod_leaftemp.o mod_clm.o \
             mod_aerosol.o mod_radiation.o mod_bdycod.o mod_scenarios.o \
						 mod_cu_grell.o
else
mod_param.o: mod_param.F90 mod_constants.o mod_bats.o mod_cu_em.o mod_date.o \
             mod_main.o mod_message.o mod_runparams.o mod_o3blk.o \
             mod_pmoist.o mod_dynparam.o mod_trachem.o mod_mppio.o \
             mod_cu_bm.o mod_split.o mod_outrad.o mod_holtbl.o \
             mod_che_semdde.o mod_ncio.o mod_leaftemp.o mod_chem.o \
             mod_aerosol.o mod_radiation.o mod_bdycod.o mod_scenarios.o \
						 mod_cu_grell.o
endif
mod_pbldim.o: mod_pbldim.f90 mod_dynparam.o
mod_pmoist.o: mod_pmoist.f90 mod_dynparam.o
mod_precip.o: mod_precip.f90 mod_bats.o mod_constants.o mod_cvaria.o \
              mod_date.o mod_main.o mod_runparams.o mod_pmoist.o \
              mod_dynparam.o mod_slice.o mod_trachem.o
mod_rad.o: mod_rad.f90 mod_dynparam.o
ifeq ($(USECLM),1)
mod_radiation.o: mod_radiation.F90 mod_constants.o mod_date.o mod_message.o \
                 mod_dynparam.o mod_tracer.o mod_scenarios.o mod_aerosol.o \
                 mod_bats.o mod_runparams.o mod_clm.o
else
mod_radiation.o: mod_radiation.F90 mod_constants.o mod_date.o mod_message.o \
                 mod_dynparam.o mod_tracer.o mod_scenarios.o mod_aerosol.o \
                 mod_bats.o mod_runparams.o
endif
mod_runparams.o: mod_runparams.F90 mod_message.o mod_dynparam.o mod_service.o
ifeq ($(USECLM),1)
mod_savefile.o: mod_savefile.F90 mod_dynparam.o mod_message.o mod_runparams.o \
                mod_bats.o mod_pmoist.o mod_main.o mod_mainchem.o mod_bdycod.o \
                mod_rad.o mod_trachem.o mod_date.o mod_radiation.o \
                mod_diagnosis.o mod_mppio.o mod_lake.o mod_cu_bm.o \
                mod_clm.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_savefile.o: mod_savefile.F90 mod_dynparam.o mod_message.o mod_runparams.o \
                mod_bats.o mod_pmoist.o mod_main.o mod_mainchem.o mod_bdycod.o \
                mod_rad.o mod_trachem.o mod_date.o mod_radiation.o \
                mod_diagnosis.o mod_mppio.o mod_lake.o mod_cu_bm.o
endif
mod_scenarios.o: mod_scenarios.F90 mod_message.o
mod_slice.o: mod_slice.F90 mod_dynparam.o mod_constants.o mod_main.o \
             mod_runparams.o mod_mainchem.o mod_pbldim.o mod_pmoist.o
mod_split.o: mod_split.F90 mod_dynparam.o mod_main.o mod_runparams.o \
             mod_constants.o mod_vmodes.o mod_date.o mod_savefile.o mod_service.o
ifeq ($(USECLM),1)
mod_sun.o: mod_sun.F90 $(CLMLIB) mod_clm.o mod_message.o mod_date.o \
             mod_constants.o mod_dynparam.o mod_main.o
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_sun.o: mod_sun.F90 mod_constants.o mod_date.o mod_message.o \
           mod_dynparam.o mod_main.o
endif
ifeq ($(USECLM),1)
mod_tendency.o: mod_tendency.F90 mod_aerosol.o mod_bats.o mod_bdycod.o \
                mod_constants.o mod_cvaria.o mod_date.o mod_main.o \
                mod_mainchem.o mod_message.o mod_runparams.o mod_pmoist.o \
                mod_rad.o mod_vecbats.o mod_holtbl.o mod_dynparam.o \
                mod_slice.o mod_trachem.o mod_cu_grell.o mod_cu_kuo.o \
                mod_cu_bm.o mod_cu_em.o mod_zengocn.o mod_sun.o \
                mod_precip.o mod_colmod3.o mod_diagnosis.o mod_condtq.o \
                mod_diffusion.o mod_advection.o mod_nudge.o mod_che_tend.o \
                mod_clm.o mod_mtrxclm.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
mod_tendency.o: mod_tendency.F90 mod_aerosol.o mod_bats.o mod_bdycod.o \
                mod_constants.o mod_cvaria.o mod_date.o mod_main.o \
                mod_mainchem.o mod_message.o mod_runparams.o mod_pmoist.o \
                mod_rad.o mod_vecbats.o mod_holtbl.o mod_dynparam.o \
                mod_slice.o mod_trachem.o mod_cu_grell.o mod_cu_kuo.o \
                mod_cu_bm.o mod_cu_em.o mod_zengocn.o mod_sun.o \
                mod_precip.o mod_colmod3.o mod_diagnosis.o mod_condtq.o \
                mod_diffusion.o mod_advection.o mod_nudge.o mod_che_tend.o
endif
mod_tracer.o: mod_tracer.f90 mod_constants.o mod_dynparam.o
mod_trachem.o: mod_trachem.F90 mod_dynparam.o mod_message.o mod_runparams.o
mod_tstep.o: mod_tstep.f90 mod_date.o mod_runparams.o mod_dynparam.o
mod_vecbats.o: mod_vecbats.F90 mod_constants.o mod_dynparam.o mod_date.o \
               mod_runparams.o mod_bats.o mod_slice.o mod_pbldim.o \
               mod_lake.o mod_main.o mod_bndry.o mod_drag.o mod_mppio.o \
               mod_message.o
mod_vmodes.o: mod_vmodes.F90 mod_runparams.o mod_dynparam.o \
              mod_constants.o mod_message.o $(NETLIB)
	$(F90) $(CPPFLAGS) -Inetlib $(F90FLAGS) -c $<
ifeq ($(USECLM),1)
mod_zengocn.o: mod_zengocn.F90 mod_bats.o mod_constants.o mod_date.o \
               mod_main.o mod_runparams.o mod_pbldim.o mod_dynparam.o \
               mod_slice.o mod_bdycod.o $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS) -c $<
else
mod_zengocn.o: mod_zengocn.F90 mod_bats.o mod_constants.o mod_date.o \
               mod_main.o mod_runparams.o mod_pbldim.o mod_dynparam.o \
               mod_slice.o mod_bdycod.o
endif
#
# RegCM Main Program
#
ifeq ($(USECLM),1)
regcm.o: regcm.F90 mod_date.o mod_message.o mod_runparams.o mod_init.o \
         mod_dynparam.o mod_ncio.o mod_output.o mod_split.o mod_header.o \
         mod_param.o mod_tendency.o mod_tstep.o mod_chem.o mod_mppio.o \
         $(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
else
regcm.o: regcm.F90 mod_date.o mod_message.o mod_runparams.o mod_init.o \
         mod_dynparam.o mod_ncio.o mod_output.o mod_split.o mod_header.o \
         mod_param.o mod_tendency.o mod_tstep.o mod_chem.o mod_mppio.o
endif
#
# CLM subroutines
#
mod_mtrxclm.o: mod_mtrxclm.F90 mod_bats.o mod_dynparam.o mod_drag.o \
               mod_constants.o mod_runparams.o mod_mppio.o mod_date.o \
               mod_clm.o mod_main.o mod_pbldim.o mod_slice.o mod_date.o \
	$(CLMLIB)
	$(F90) $(CPPFLAGS) -I$(CLMPTH) $(F90FLAGS)   -c $<
#
# Now unused, my be activated (?)
#
mod_co2.o: mod_co2.f90 mod_dynparam.o mod_runparams.o mod_bats.o mod_leaftemp.o
#
