
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/Config:.
vpath %.f90 $(REGCM_ROOT)/Config:.

CLMLIBPTH = ./clmlib
NETLIBPTH = ./netlib

CLMLIB = $(CLMLIBPTH)/libclm.a
NETLIB = $(NETLIBPTH)/libnetlib.a

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

if COMPILER_SUN
  CPPFLAGS += -M$(CLMLIBPTH) -M$(NETLIBPTH)
else
	CPPFLAGS += -I$(CLMLIBPTH) -I$(NETLIBPTH)
endif

LDFLAGS += $(LCLMLIB) $(LNETLIB)

bin_PROGRAMS = regcm

regcm_SOURCES = mod_constants.F90 mod_dynparam.F90 mod_bats_param.f90 \
								mod_service.F90 mod_header.F90 mod_advection.F90      \
                mod_aerosol.F90 mod_bats.F90 mod_bdycod.F90           \
                mod_bndry.F90 mod_chem.F90 mod_che_semdde.F90         \
                mod_che_tend.F90 mod_cldfrac.f90 mod_colmod3.F90      \
								mod_condtq.f90 mod_cu_bm.F90 mod_cu_em.f90            \
								mod_cu_grell.f90 mod_cu_kuo.F90 mod_cumtran.F90       \
								mod_cvaria.f90 mod_date.F90 mod_diagnosis.F90         \
								mod_diffusion.F90 mod_drag.F90 mod_dust.F90           \
								mod_holtbl.F90 mod_init.F90 mod_lake.F90              \
                mod_leaftemp.f90 mod_mainchem.f90 mod_main.f90        \
                mod_message.F90 mod_mppio.F90 mod_ncio.F90            \
                mod_nudge.F90 mod_o3blk.F90 mod_output.F90            \
                mod_outrad.F90 mod_param.F90 mod_pbldim.f90           \
                mod_pmoist.f90 mod_precip.f90 mod_rad.f90             \
                mod_radiation.F90 mod_runparams.F90 mod_savefile.F90  \
                mod_scenarios.F90 mod_slice.F90 mod_split.F90         \
                mod_sun.F90 mod_tendency.F90 mod_tracer.f90           \
                mod_trachem.F90 mod_tstep.f90 mod_vecbats.F90         \
                mod_vmodes.F90 mod_zengocn.F90 mod_cu_tiedtke.F90     \
								mod_memutil.f90 mod_cu_tables.f90 regcm.F90

if DO_COMPILE_LIBCLM
regcm_SOURCES += mod_clm.F90 mod_mtrxclm.F90
endif

EXTRAOBJS = mod_co2.o

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
mod_dynparam.o: mod_dynparam.F90 mod_constants.o
mod_bats_param.o: mod_bats_param.f90
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
if DO_COMPILE_LIBCLM
mod_che_semdde.o: mod_che_semdde.F90 mod_aerosol.o mod_dust.o mod_mainchem.o \
                  mod_message.o mod_dynparam.o mod_trachem.o mod_constants.o \
                  mod_ncio.o mod_mppio.o mod_bats_param.o $(CLMLIB)
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
                mod_runparams.o mod_cvaria.o
mod_cu_kuo.o: mod_cu_kuo.F90 mod_bats.o mod_constants.o mod_cvaria.o \
              mod_date.o mod_main.o mod_runparams.o mod_pmoist.o mod_rad.o \
              mod_dynparam.o mod_trachem.o mod_advection.o
mod_cu_tables.o: mod_cu_tables.f90 mod_constants.o
mod_cu_tiedke.o: mod_cu_tiedtke.F90 mod_dynparam.o mod_runparams.o \
                 mod_constants.o mod_cu_tables.o
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
mod_holtbl.o: mod_holtbl.F90 mod_bats.o mod_constants.o mod_cvaria.o \
          mod_diagnosis.o mod_main.o mod_mainchem.o mod_runparams.o \
          mod_pbldim.o mod_pmoist.o mod_dynparam.o mod_slice.o \
          mod_trachem.o mod_mppio.o
if DO_COMPILE_LIBCLM
mod_init.o: mod_init.F90 mod_bats.o mod_bdycod.o mod_constants.o mod_date.o \
        mod_diagnosis.o mod_main.o mod_mainchem.o mod_message.o \
        mod_runparams.o mod_pmoist.o mod_rad.o mod_vecbats.o mod_o3blk.o \
        mod_radiation.o mod_dynparam.o mod_mppio.o mod_trachem.o mod_lake.o \
        mod_ncio.o mod_sun.o mod_savefile.o mod_cu_bm.o $(CLMLIB)
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
mod_memutil.o: mod_memutil.f90 mod_dynparam.o mod_constants.o
mod_message.o: mod_message.F90 mod_dynparam.o
mod_mppio.o: mod_mppio.F90 mod_dynparam.o mod_message.o mod_runparams.o mod_main.o
mod_ncio.o: mod_ncio.F90 mod_dynparam.o mod_message.o mod_constants.o \
            mod_date.o mod_runparams.o mod_trachem.o
mod_nudge.o: mod_nudge.F90 mod_bdycod.o mod_main.o mod_runparams.o \
             mod_dynparam.o
mod_o3blk.o: mod_o3blk.F90 mod_dynparam.o mod_main.o mod_runparams.o mod_rad.o
if DO_COMPILE_LIBCLM
mod_output.o: mod_output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_runparams.o mod_pmoist.o mod_rad.o mod_diagnosis.o \
          mod_radiation.o mod_dynparam.o mod_split.o mod_clm.o \
          mod_trachem.o mod_mppio.o mod_ncio.o mod_lake.o mod_savefile.o \
          mod_cu_bm.o $(CLMLIB)
else
mod_output.o: mod_output.F90 mod_bats.o mod_bdycod.o mod_cvaria.o mod_date.o \
          mod_main.o mod_mainchem.o mod_message.o mod_outrad.o \
          mod_runparams.o mod_pmoist.o mod_rad.o mod_diagnosis.o \
          mod_radiation.o mod_dynparam.o mod_split.o mod_cu_bm.o \
          mod_trachem.o mod_mppio.o mod_lake.o mod_ncio.o mod_savefile.o
endif
if DO_COMPILE_LIBCLM
mod_outrad.o: mod_outrad.F90 mod_dynparam.o mod_runparams.o mod_bats.o \
              mod_rad.o mod_date.o mod_clm.o
else
mod_outrad.o: mod_outrad.F90 mod_dynparam.o mod_runparams.o mod_bats.o \
              mod_rad.o mod_date.o
endif
if DO_COMPILE_LIBCLM
mod_param.o: mod_param.F90 mod_constants.o mod_bats.o mod_cu_em.o mod_date.o \
             mod_main.o mod_message.o mod_runparams.o mod_o3blk.o \
             mod_pmoist.o mod_dynparam.o mod_trachem.o mod_mppio.o \
             mod_cu_bm.o mod_split.o mod_outrad.o mod_holtbl.o mod_chem.o \
             mod_che_semdde.o mod_ncio.o mod_leaftemp.o mod_clm.o \
             mod_aerosol.o mod_radiation.o mod_bdycod.o mod_scenarios.o \
						 mod_cu_tables.o mod_cu_tiedtke.o mod_cu_grell.o
else
mod_param.o: mod_param.F90 mod_constants.o mod_bats.o mod_cu_em.o mod_date.o \
             mod_main.o mod_message.o mod_runparams.o mod_o3blk.o \
             mod_pmoist.o mod_dynparam.o mod_trachem.o mod_mppio.o \
             mod_cu_bm.o mod_split.o mod_outrad.o mod_holtbl.o \
             mod_che_semdde.o mod_ncio.o mod_leaftemp.o mod_chem.o \
             mod_aerosol.o mod_radiation.o mod_bdycod.o mod_scenarios.o \
						 mod_cu_tables.o mod_cu_tiedtke.o mod_cu_grell.o
endif
mod_pbldim.o: mod_pbldim.f90 mod_dynparam.o
mod_pmoist.o: mod_pmoist.f90 mod_dynparam.o
mod_precip.o: mod_precip.f90 mod_bats.o mod_constants.o mod_cvaria.o \
              mod_date.o mod_main.o mod_runparams.o mod_pmoist.o \
              mod_dynparam.o mod_slice.o mod_trachem.o
mod_rad.o: mod_rad.f90 mod_dynparam.o
if DO_COMPILE_LIBCLM
mod_radiation.o: mod_radiation.F90 mod_constants.o mod_date.o mod_message.o \
                 mod_dynparam.o mod_tracer.o mod_scenarios.o mod_aerosol.o \
                 mod_bats.o mod_runparams.o mod_clm.o
else
mod_radiation.o: mod_radiation.F90 mod_constants.o mod_date.o mod_message.o \
                 mod_dynparam.o mod_tracer.o mod_scenarios.o mod_aerosol.o \
                 mod_bats.o mod_runparams.o
endif
mod_runparams.o: mod_runparams.F90 mod_message.o mod_dynparam.o mod_service.o
if DO_COMPILE_LIBCLM
mod_savefile.o: mod_savefile.F90 mod_dynparam.o mod_message.o mod_runparams.o \
                mod_bats.o mod_pmoist.o mod_main.o mod_mainchem.o mod_bdycod.o \
                mod_rad.o mod_trachem.o mod_date.o mod_radiation.o \
                mod_diagnosis.o mod_mppio.o mod_lake.o mod_cu_bm.o \
                mod_clm.o $(CLMLIB)
else
mod_savefile.o: mod_savefile.F90 mod_dynparam.o mod_message.o mod_runparams.o \
                mod_bats.o mod_pmoist.o mod_main.o mod_mainchem.o mod_bdycod.o \
                mod_rad.o mod_trachem.o mod_date.o mod_radiation.o \
                mod_diagnosis.o mod_mppio.o mod_lake.o mod_cu_bm.o
endif
mod_scenarios.o: mod_scenarios.F90 mod_message.o
mod_service.o: mod_service.F90 mod_dynparam.o
mod_slice.o: mod_slice.F90 mod_dynparam.o mod_constants.o mod_main.o \
             mod_runparams.o mod_mainchem.o mod_pbldim.o mod_pmoist.o
mod_split.o: mod_split.F90 mod_dynparam.o mod_main.o mod_runparams.o \
            mod_constants.o mod_vmodes.o mod_date.o mod_savefile.o mod_service.o
if DO_COMPILE_LIBCLM
mod_sun.o: mod_sun.F90 $(CLMLIB) mod_clm.o mod_message.o mod_date.o \
             mod_constants.o mod_dynparam.o mod_main.o
else
mod_sun.o: mod_sun.F90 mod_constants.o mod_date.o mod_message.o \
           mod_dynparam.o mod_main.o
endif
if DO_COMPILE_LIBCLM
mod_tendency.o: mod_tendency.F90 mod_aerosol.o mod_bats.o mod_bdycod.o \
                mod_constants.o mod_cvaria.o mod_date.o mod_main.o \
                mod_mainchem.o mod_message.o mod_runparams.o mod_pmoist.o \
                mod_rad.o mod_vecbats.o mod_holtbl.o mod_dynparam.o \
                mod_slice.o mod_trachem.o mod_cu_grell.o mod_cu_kuo.o \
                mod_cu_bm.o mod_cu_em.o mod_zengocn.o mod_sun.o \
                mod_precip.o mod_colmod3.o mod_diagnosis.o mod_condtq.o \
                mod_diffusion.o mod_advection.o mod_nudge.o mod_che_tend.o \
                mod_clm.o mod_mtrxclm.o $(CLMLIB)
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
              mod_constants.o mod_message.o
if DO_COMPILE_LIBCLM
mod_zengocn.o: mod_zengocn.F90 mod_bats.o mod_constants.o mod_date.o \
               mod_main.o mod_runparams.o mod_pbldim.o mod_dynparam.o \
               mod_slice.o mod_bdycod.o $(CLMLIB)
else
mod_zengocn.o: mod_zengocn.F90 mod_bats.o mod_constants.o mod_date.o \
               mod_main.o mod_runparams.o mod_pbldim.o mod_dynparam.o \
               mod_slice.o mod_bdycod.o
endif
#
# RegCM Main Program
#
if DO_COMPILE_LIBCLM
regcm.o: regcm.F90 mod_date.o mod_message.o mod_runparams.o mod_init.o \
         mod_dynparam.o mod_ncio.o mod_output.o mod_split.o mod_header.o \
         mod_param.o mod_tendency.o mod_tstep.o mod_chem.o mod_mppio.o \
         $(CLMLIB)
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
#
# Now unused, my be activated (?)
#
mod_co2.o: mod_co2.f90 mod_dynparam.o mod_runparams.o mod_bats.o mod_leaftemp.o
#
