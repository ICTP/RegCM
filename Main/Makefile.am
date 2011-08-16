
SUBDIRS = netlib clmlib chem

include $(REGCM_ROOT)/makeinc

CLMLIBPTH = ./clmlib
NETLIBPTH = ./netlib
CHELIBPTH = ./chem
KPPLIBPTH = chem/GAS_KPP_GEOS

CLMLIB = $(CLMLIBPTH)/libclm.a
NETLIB = $(NETLIBPTH)/libnetlib.a
CHELIB = $(CHELIBPTH)/libchem.a
KPPLIB = $(KPPLIBPTH)/libkpp.a

if REGCM_PARALLEL
rcsuffix = MPI
else
rcsuffix = Serial
endif  

if DO_COMPILE_LIBCLM
rcaddsuffix = _clm
endif

if DO_COMPILE_BAND
rcaddsuffix = _band
endif

program_transform_name = s&$$&$(rcsuffix)$(rcaddsuffix)&

LCLMLIB = -L$(CLMLIBPTH) -lclm
LNETLIB = -L$(NETLIBPTH) -lnetlib
LCHELIB = -L$(CHELIBPTH) -lchem
LKPPLIB = -L$(KPPLIBPTH) -lkpp

if COMPILER_SUN
  CPPFLAGS += -M$(CLMLIBPTH) -M$(NETLIBPTH) -M$(KPPLIBPTH) -M$(CHELIBPTH)
else
  CPPFLAGS += -I$(CLMLIBPTH) -I$(NETLIBPTH) -I$(KPPLIBPTH) -I$(CHELIBPTH)
endif

LIBS += $(LCLMLIB) $(LNETLIB) $(LCHELIB) $(LKPPLIB)

bin_PROGRAMS = regcm

regcm_SOURCES = mod_bats_param.f90 mod_diagnosis.F90   \
                mod_service.F90 mod_header.F90 mod_advection.F90      \
                mod_bats.F90 mod_bdycod.F90           \
                mod_bndry.F90 \
                mod_che_tend.F90 mod_cldfrac.f90 mod_colmod3.F90      \
                mod_condtq.f90 mod_cu_bm.F90 mod_cu_em.f90            \
                mod_cu_grell.F90 mod_cu_kuo.F90 mod_cumtran.F90       \
                mod_diffusion.F90 mod_drag.F90 \
                mod_holtbl.F90 mod_init.F90 mod_lake.F90              \
                mod_leaftemp.f90 mod_atm_interface.f90        \
                mod_mppio.F90 mod_ncio.F90 mod_cu_tables.f90          \
                mod_nudge.F90 mod_o3blk.F90 mod_output.F90            \
                mod_outrad.F90 mod_params.F90 mod_pbldim.f90          \
                mod_pmoist.f90 mod_precip.f90 mod_rad.f90             \
                mod_radiation.F90 mod_runparams.F90 mod_savefile.F90  \
                mod_scenarios.F90 mod_slice.F90 mod_split.F90         \
                mod_sun.F90 mod_tendency.F90 mod_tracer.f90           \
                mod_tstep.f90 mod_vecbats.F90         \
                mod_vmodes.F90 mod_zengocn.F90 mod_cu_tiedtke.F90     \
                mod_interface.F90 mod_ensemble.f90 mod_gridfuncs.f90  \
                mod_thetal.f90 mod_tcm_interface.F90 mod_uwtcm.F90    \
                regcm.F90

if DO_COMPILE_LIBCLM
regcm_SOURCES += mod_clm.F90 mod_mtrxclm.F90
endif

EXTRAOBJS = mod_co2.o

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

#
# Modules
#
mod_advection.o: mod_advection.F90 mod_atm_interface.o mod_runparams.o
mod_bats.o: mod_bats.F90 mod_bats_param.o mod_runparams.o
mod_bats_param.o: mod_bats_param.f90
mod_bdycod.o: mod_bdycod.F90 mod_runparams.o mod_atm_interface.o \
              mod_bats.o mod_ncio.o mod_mppio.o
mod_bndry.o: mod_bndry.F90 mod_bats.o mod_leaftemp.o \
             mod_runparams.o mod_drag.o
mod_che_tend.o: mod_che_tend.F90 mod_bats.o mod_atm_interface.o \
                mod_runparams.o mod_pbldim.o mod_pmoist.o mod_rad.o \
                mod_diffusion.o \
                mod_advection.o mod_diagnosis.o mod_mppio.o mod_slice.o
mod_cldfrac.o: mod_cldfrac.f90 mod_pmoist.o mod_rad.o \
               mod_slice.o
mod_clm.o: mod_clm.F90
mod_colmod3.o: mod_colmod3.F90 mod_bats.o mod_runparams.o \
               mod_radiation.o mod_atm_interface.o mod_rad.o \
               mod_outrad.o
mod_condtq.o: mod_condtq.f90 mod_atm_interface.o \
              mod_runparams.o mod_pmoist.o mod_slice.o
mod_cumtran.o: mod_cumtran.F90 mod_runparams.o
mod_cu_bm.o: mod_cu_bm.F90 mod_bats.o mod_atm_interface.o \
             mod_runparams.o mod_pmoist.o mod_rad.o
mod_cu_em.o: mod_cu_em.f90 mod_bats.o mod_atm_interface.o \
             mod_runparams.o mod_pmoist.o mod_rad.o mod_slice.o
mod_cu_grell.o: mod_cu_grell.F90 mod_pmoist.o mod_rad.o \
                mod_bats.o mod_atm_interface.o mod_runparams.o 
mod_cu_kuo.o: mod_cu_kuo.F90 mod_bats.o \
              mod_atm_interface.o mod_runparams.o mod_pmoist.o mod_rad.o \
              mod_advection.o
mod_cu_tables.o: mod_cu_tables.f90
mod_cu_tiedke.o: mod_cu_tiedtke.F90 mod_runparams.o \
                 mod_cu_tables.o
mod_diagnosis.o: mod_diagnosis.F90 mod_runparams.o mod_atm_interface.o mod_mppio.o
mod_diffusion.o: mod_diffusion.F90 mod_runparams.o
mod_drag.o: mod_drag.F90 mod_bats.o \
            mod_runparams.o
mod_ensemble.o: mod_ensemble.f90
mod_gridfuncs.o: mod_gridfuncs.f90 mod_atm_interface.o mod_runparams.o
mod_header.o: mod_header.F90
mod_holtbl.o: mod_holtbl.F90 mod_bats.o \
          mod_diagnosis.o mod_atm_interface.o mod_runparams.o \
          mod_pbldim.o mod_pmoist.o mod_slice.o \
          mod_mppio.o
if DO_COMPILE_LIBCLM
mod_init.o: mod_init.F90 mod_bats.o mod_bdycod.o mod_tcm_interface.o \
        mod_diagnosis.o mod_atm_interface.o mod_uwtcm.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o mod_o3blk.o \
        mod_radiation.o mod_mppio.o mod_lake.o \
        mod_ncio.o mod_sun.o mod_savefile.o mod_cu_bm.o $(CLMLIB)
else
mod_init.o: mod_init.F90 mod_bats.o mod_bdycod.o mod_tcm_interface.o \
        mod_diagnosis.o mod_atm_interface.o mod_uwtcm.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o mod_o3blk.o \
        mod_radiation.o mod_mppio.o mod_ncio.o \
        mod_lake.o mod_sun.o mod_savefile.o
endif
mod_lake.o: mod_lake.F90 mod_bats.o mod_runparams.o mod_mppio.o
mod_leaftemp.o: mod_leaftemp.f90 mod_bats.o \
                mod_runparams.o
mod_atm_interface.o: mod_atm_interface.f90 mod_runparams.o
mod_mppio.o: mod_mppio.F90 mod_runparams.o mod_tcm_interface.o \
             mod_atm_interface.o mod_outrad.o
mod_ncio.o: mod_ncio.F90 mod_runparams.o
mod_nudge.o: mod_nudge.F90 mod_bdycod.o mod_atm_interface.o mod_runparams.o
mod_o3blk.o: mod_o3blk.F90 mod_atm_interface.o mod_runparams.o mod_rad.o
if DO_COMPILE_LIBCLM
mod_output.o: mod_output.F90 mod_bats.o mod_bdycod.o \
          mod_atm_interface.o mod_outrad.o \
          mod_runparams.o mod_pmoist.o mod_rad.o mod_diagnosis.o \
          mod_radiation.o mod_split.o mod_clm.o \
          mod_mppio.o mod_ncio.o mod_lake.o mod_savefile.o \
          mod_cu_bm.o $(CLMLIB)
else
mod_output.o: mod_output.F90 mod_bats.o mod_bdycod.o \
          mod_atm_interface.o mod_outrad.o \
          mod_runparams.o mod_pmoist.o mod_rad.o mod_diagnosis.o \
          mod_radiation.o mod_split.o mod_cu_bm.o \
          mod_mppio.o mod_lake.o mod_ncio.o mod_savefile.o
endif
if DO_COMPILE_LIBCLM
mod_outrad.o: mod_outrad.F90 mod_runparams.o mod_bats.o \
              mod_rad.o mod_clm.o
else
mod_outrad.o: mod_outrad.F90 mod_runparams.o mod_bats.o \
              mod_rad.o
endif
if DO_COMPILE_LIBCLM
mod_params.o: mod_params.F90 mod_bats.o mod_cu_em.o \
              mod_atm_interface.o mod_runparams.o mod_o3blk.o \
              mod_pmoist.o mod_mppio.o \
              mod_cu_bm.o mod_split.o mod_outrad.o mod_holtbl.o \
              mod_ncio.o mod_leaftemp.o mod_clm.o \
              mod_radiation.o mod_bdycod.o mod_scenarios.o \
              mod_cu_tables.o mod_cu_tiedtke.o mod_cu_grell.o mod_tendency.o \
              mod_ncio.o mod_cldfrac.o
else
mod_params.o: mod_params.F90 mod_bats.o mod_cu_em.o \
              mod_atm_interface.o mod_runparams.o mod_o3blk.o \
              mod_pmoist.o mod_mppio.o \
              mod_cu_bm.o mod_split.o mod_outrad.o mod_holtbl.o \
              mod_ncio.o mod_leaftemp.o \
              mod_radiation.o mod_bdycod.o mod_scenarios.o \
              mod_cu_tables.o mod_cu_tiedtke.o mod_cu_grell.o mod_tendency.o \
              mod_ncio.o mod_cldfrac.o
endif
mod_pbldim.o: mod_pbldim.f90
mod_pmoist.o: mod_pmoist.f90
mod_precip.o: mod_precip.f90 mod_bats.o \
              mod_atm_interface.o mod_runparams.o mod_pmoist.o \
              mod_slice.o
mod_rad.o: mod_rad.f90
if DO_COMPILE_LIBCLM
mod_radiation.o: mod_radiation.F90 mod_tracer.o mod_scenarios.o \
                 mod_bats.o mod_runparams.o mod_clm.o
else
mod_radiation.o: mod_radiation.F90 mod_tracer.o mod_scenarios.o \
                 mod_bats.o mod_runparams.o
endif
mod_runparams.o: mod_runparams.F90 mod_service.o
if DO_COMPILE_LIBCLM
mod_savefile.o: mod_savefile.F90 mod_runparams.o \
                mod_bats.o mod_pmoist.o mod_atm_interface.o mod_bdycod.o \
                mod_rad.o mod_radiation.o \
                mod_diagnosis.o mod_mppio.o mod_lake.o mod_cu_bm.o \
                mod_clm.o $(CLMLIB)
else
mod_savefile.o: mod_savefile.F90 mod_runparams.o \
                mod_bats.o mod_pmoist.o mod_atm_interface.o mod_bdycod.o \
                mod_rad.o mod_radiation.o \
                mod_diagnosis.o mod_mppio.o mod_lake.o mod_cu_bm.o
endif
mod_scenarios.o: mod_scenarios.F90
mod_service.o: mod_service.F90
mod_slice.o: mod_slice.F90 mod_atm_interface.o \
             mod_runparams.o mod_pbldim.o mod_pmoist.o
mod_split.o: mod_split.F90 mod_atm_interface.o mod_runparams.o \
             mod_vmodes.o mod_savefile.o mod_service.o
if DO_COMPILE_LIBCLM
mod_sun.o: mod_sun.F90 $(CLMLIB) mod_clm.o \
             mod_atm_interface.o
else
mod_sun.o: mod_sun.F90 mod_atm_interface.o
endif
mod_tcm_interface.o: mod_tcm_interface.F90 mod_atm_interface.o mod_runparams.o \
                     mod_slice.o mod_rad.o mod_gridfuncs.o \
                     mod_service.o mod_pbldim.o
if DO_COMPILE_LIBCLM
mod_tendency.o: mod_tendency.F90 mod_bats.o mod_bdycod.o \
                mod_atm_interface.o \
                mod_runparams.o mod_pmoist.o \
                mod_rad.o mod_vecbats.o mod_holtbl.o \
                mod_slice.o mod_cu_grell.o mod_cu_kuo.o \
                mod_cu_bm.o mod_cu_em.o mod_sun.o mod_cumtran.o \
                mod_precip.o mod_colmod3.o mod_diagnosis.o mod_condtq.o \
                mod_diffusion.o mod_advection.o mod_nudge.o mod_che_tend.o \
                mod_clm.o mod_mtrxclm.o mod_cldfrac.o $(CLMLIB)
else
mod_tendency.o: mod_tendency.F90 mod_bats.o mod_bdycod.o \
                mod_atm_interface.o mod_runparams.o mod_pmoist.o \
                mod_rad.o mod_vecbats.o mod_holtbl.o \
                mod_slice.o mod_cu_grell.o mod_cu_kuo.o \
                mod_cu_bm.o mod_cu_em.o mod_sun.o mod_cumtran.o \
                mod_precip.o mod_colmod3.o mod_diagnosis.o mod_condtq.o \
                mod_diffusion.o mod_advection.o mod_nudge.o mod_che_tend.o \
                mod_cldfrac.o
endif
mod_thetal.o: mod_thetal.f90
mod_tracer.o: mod_tracer.f90
mod_tstep.o: mod_tstep.f90 mod_runparams.o
mod_uwtcm.o: mod_uwtcm.F90 mod_tcm_interface.o mod_atm_interface.o \
	           mod_runparams.o mod_thetal.o
mod_vecbats.o: mod_vecbats.F90 \
               mod_runparams.o mod_bats.o mod_slice.o mod_pbldim.o \
               mod_lake.o mod_atm_interface.o mod_bndry.o mod_drag.o mod_mppio.o \
               mod_zengocn.o
mod_vmodes.o: mod_vmodes.F90 mod_runparams.o
if DO_COMPILE_LIBCLM
mod_zengocn.o: mod_zengocn.F90 mod_bats.o \
               mod_atm_interface.o mod_runparams.o mod_pbldim.o \
               mod_slice.o mod_bdycod.o $(CLMLIB)
else
mod_zengocn.o: mod_zengocn.F90 mod_bats.o \
               mod_atm_interface.o mod_runparams.o mod_pbldim.o \
               mod_slice.o mod_bdycod.o
endif
#
# RegCM Main Program
#
if DO_COMPILE_LIBCLM
mod_interface.o: mod_interface.F90 mod_runparams.o mod_init.o mod_ncio.o \
                 mod_output.o mod_split.o mod_header.o mod_params.o mod_tendency.o \
                 mod_tstep.o mod_mppio.o $(CLMLIB)
else
mod_interface.o: mod_interface.F90 mod_runparams.o mod_init.o mod_ncio.o \
                 mod_output.o mod_split.o mod_header.o mod_params.o mod_tendency.o \
                 mod_tstep.o mod_mppio.o
endif
regcm.o: regcm.F90 regcm_nonesmf.F90 regcm_esmf.F90
#
# CLM subroutines
#
mod_mtrxclm.o: mod_mtrxclm.F90 mod_bats.o mod_drag.o \
               mod_runparams.o mod_mppio.o \
               mod_clm.o mod_atm_interface.o mod_pbldim.o mod_slice.o \
               $(CLMLIB)
#
# Now unused, my be activated (?)
#
mod_co2.o: mod_co2.f90 mod_runparams.o mod_bats.o mod_leaftemp.o
#
