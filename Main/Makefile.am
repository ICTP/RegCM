#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#
#    ICTP RegCM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

vpath %.o $(REGCM_ROOT)/Share:mpplib:netlib:batslib:chemlib:clmlib:\
          cumlib:radlib:pbllib:.

SUBDIRS = mpplib netlib batslib chemlib clmlib cumlib radlib pbllib

include $(REGCM_ROOT)/makeinc

MPPLIBPTH = ./mpplib
CLMLIBPTH = ./clmlib
BATLIBPTH = ./batslib
NETLIBPTH = ./netlib
CHELIBPTH = ./chemlib
CUMLIBPTH = ./cumlib
RADLIBPTH = ./radlib
PBLLIBPTH = ./pbllib
CBMZLIBPTH = chemlib/GAS_CBMZ
RRTMSWLIBPTH = radlib/RRTMG_SW

MPPLIB = $(MPPLIBPTH)/libmpp.a
CLMLIB = $(CLMLIBPTH)/libclm.a
BATLIB = $(BATLIBPTH)/libbats.a
NETLIB = $(NETLIBPTH)/libnetlib.a
CHELIB = $(CHELIBPTH)/libchem.a
CUMLIB = $(CUMLIBPTH)/libcum.a
RADLIB = $(RADLIBPTH)/librad.a
PBLLIB = $(PBLLIBPTH)/libpbl.a
CBMZLIB = $(CBMZLIBPTH)/libcbmz.a
RRTMSWLIB = $(RRTMSWLIBPTH)/librrtmgsw.a

if REGCM_PARALLEL
rcsuffix = MPI
else
rcsuffix = Serial
endif

if DO_COMPILE_LIBCLM
rcaddsuffix = _clm
endif

program_transform_name = s&$$&$(rcsuffix)$(rcaddsuffix)&

LMPPLIB = -L$(MPPLIBPTH) -lmpp
LBATLIB = -L$(BATLIBPTH) -lbats
LCLMLIB = -L$(CLMLIBPTH) -lclm
LNETLIB = -L$(NETLIBPTH) -lnetlib
LCHELIB = -L$(CHELIBPTH) -lchem
LCUMLIB = -L$(CUMLIBPTH) -lcum
LRADLIB = -L$(RADLIBPTH) -lrad
LPBLLIB = -L$(PBLLIBPTH) -lpbl
LCBMZLIB = -L$(CBMZLIBPTH) -lcbmz
LRRTMSWLIB = -L$(RRTMSWLIBPTH) -lrrtmgsw

if COMPILER_SUN
  CPPFLAGS += -M$(MPPLIBPTH) -M$(CLMLIBPTH) -M$(NETLIBPTH) -M$(CBMZLIBPTH) \
              -M$(CHELIBPTH) -M$(BATLIBPTH) -M$(CUMLIBPTH) -M$(RADLIBPTH)  \
              -M$(PBLLIBPTH) -M$(RRTMSWLIBPTH)
else
  CPPFLAGS += -I$(MPPLIBPTH) -I$(CLMLIBPTH) -I$(NETLIBPTH) -I$(CBMZLIBPTH) \
              -I$(CHELIBPTH) -I$(BATLIBPTH) -I$(CUMLIBPTH) -I$(RADLIBPTH)  \
              -I$(PBLLIBPTH) -I$(RRTMSWLIBPTH)
endif

LIBS = $(LMPPLIB) $(LCLMLIB) $(LNETLIB) $(LCHELIB) $(LCBMZLIB) $(LBATLIB) \
       $(LCUMLIB) $(LRADLIB) $(LPBLLIB) $(LRRTMSWLIB) $(SYSLIBS)

bin_PROGRAMS = regcm

regcm_SOURCES = mod_regcm_interface.F90 mod_atm_interface.F90 mod_che_interface.F90 \
  mod_lm_interface.F90 mod_cu_interface.f90 mod_rad_interface.f90 \
  mod_pbl_interface.F90 mod_mppio.F90 mod_header.F90 \
  mod_advection.F90 mod_ncio.F90 mod_bdycod.F90 mod_diffusion.F90 \
  mod_precip.F90 mod_sun.F90 mod_savefile.F90 mod_slice.F90 mod_init.F90 \
  mod_vmodes.F90 mod_split.F90 mod_output.F90 mod_cloud_s1.F90 \
  mod_tendency.F90 mod_params.F90 mod_tstep.f90 mod_slabocean.F90 \
  mod_sldepparam.f90 mod_sladvection.F90 abort.f90 regcm.f90

FCLD = $(MPIFC)

if DO_COMPILE_CPL
regcm_SOURCES += mod_update.F90
endif

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
#
# Modules
#
mod_advection.o: mod_advection.F90 mod_atm_interface.o mod_dynparam.o mod_runparams.o mod_memutil.o mod_mpmessage.o mod_service.o
mod_atm_interface.o: mod_atm_interface.F90 mod_dynparam.o mod_stdio.o mod_constants.o mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_service.o mod_memutil.o
mod_bdycod.o: mod_bdycod.F90 mod_dynparam.o mod_mppparam.o mod_memutil.o mod_atm_interface.o mod_pbl_interface.o mod_che_interface.o mod_lm_interface.o mod_mpmessage.o mod_ncio.o mod_mppio.o mod_service.o mod_slabocean.o
mod_che_interface.o: mod_che_interface.F90 mod_realkinds.o mod_atm_interface.o mod_che_common.o mod_che_cumtran.o mod_che_dust.o mod_che_indices.o mod_che_mppio.o mod_che_ncio.o mod_che_param.o mod_che_drydep.o mod_che_bdyco.o mod_che_emission.o mod_che_carbonaer.o mod_che_species.o mod_che_tend.o mod_che_start.o
mod_cloud_s1.o: mod_cloud_s1.F90 mod_realkinds.o mod_dynparam.o mod_mpmessage.o mod_memutil.o mod_atm_interface.o mod_runparams.o mod_runparams.o mod_runparams.o mod_runparams.o mod_runparams.o mod_runparams.o mod_runparams.o mod_pbl_common.o
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $< -finline-functions
mod_cu_interface.o: mod_cu_interface.f90 mod_constants.o mod_realkinds.o mod_cu_common.o mod_cu_tiedtke.o mod_cu_tables.o mod_cu_bm.o mod_cu_em.o mod_cu_kuo.o mod_cu_grell.o mod_atm_interface.o
mod_diffusion.o: mod_diffusion.F90 mod_atm_interface.o mod_runparams.o mod_mppparam.o mod_service.o
mod_header.o: mod_header.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_mppparam.o mod_date.o mod_stdio.o
mod_slabocean.o: mod_slabocean.F90 mod_runparams.o mod_memutil.o mod_atm_interface.o mod_mppparam.o mod_constants.o mod_outvars.o
if DO_COMPILE_LIBCLM
mod_init.o: mod_init.F90 mod_runparams.o mod_mppparam.o mod_lm_interface.o mod_atm_interface.o mod_che_interface.o mod_cu_interface.o mod_rad_interface.o mod_pbl_interface.o mod_slabocean.o rrtmg_sw_init.o rrtmg_lw_init.o mod_pbl_interface.o mod_precip.o mod_bdycod.o mod_mpmessage.o mod_sun.o mod_ncio.o mod_savefile.o mod_mppio.o mod_slice.o mod_constants.o mod_outvars.o mod_clm.o mod_lm_interface.o clm_varsur.o
else
mod_init.o: mod_init.F90 mod_runparams.o mod_mppparam.o mod_lm_interface.o mod_atm_interface.o mod_che_interface.o mod_cu_interface.o mod_rad_interface.o mod_pbl_interface.o mod_slabocean.o rrtmg_sw_init.o rrtmg_lw_init.o mod_pbl_interface.o mod_precip.o mod_bdycod.o mod_mpmessage.o mod_sun.o mod_ncio.o mod_savefile.o mod_mppio.o mod_slice.o mod_constants.o mod_outvars.o mod_lm_interface.o
endif
if DO_COMPILE_LIBCLM
mod_lm_interface.o: mod_lm_interface.F90 mod_bats_common.o mod_runparams.o mod_memutil.o mod_atm_interface.o mod_mtrxclm.o mod_clm.o mod_bats_mtrxbats.o mod_bats_param.o mod_bats_mppio.o mod_bats_bndry.o mod_bats_co2.o mod_bats_drag.o mod_bats_lake.o mod_bats_leaftemp.o mod_bats_mtrxbats.o mod_bats_zengocn.o
else
mod_lm_interface.o: mod_lm_interface.F90 mod_bats_common.o mod_runparams.o mod_memutil.o mod_atm_interface.o mod_bats_mtrxbats.o mod_bats_param.o mod_bats_mppio.o mod_bats_bndry.o mod_bats_co2.o mod_bats_drag.o mod_bats_lake.o mod_bats_leaftemp.o mod_bats_mtrxbats.o mod_bats_zengocn.o
endif
mod_mppio.o: mod_mppio.F90 mod_runparams.o mod_atm_interface.o mod_lm_interface.o mod_cu_interface.o mod_rad_interface.o mod_pbl_interface.o mod_mppparam.o mod_memutil.o mod_mpmessage.o
mod_ncio.o: mod_ncio.F90 mod_runparams.o mod_dynparam.o mod_ensemble.o mod_mpmessage.o mod_mppparam.o mod_memutil.o mod_nchelper.o mod_domain.o
if DO_COMPILE_LIBCLM
mod_output.o: mod_output.F90 mod_runparams.o mod_header.o mod_mpmessage.o mod_mppparam.o mod_service.o mod_atm_interface.o mod_che_interface.o mod_che_output.o mod_lm_interface.o mod_rad_interface.o mod_cu_interface.o mod_pbl_interface.o mod_ncout.o mod_bdycod.o mod_precip.o mod_split.o mod_savefile.o mod_mppio.o mod_slabocean.o mod_clm.o
else
mod_output.o: mod_output.F90 mod_runparams.o mod_header.o mod_mpmessage.o mod_mppparam.o mod_service.o mod_atm_interface.o mod_che_interface.o mod_che_output.o mod_lm_interface.o mod_rad_interface.o mod_cu_interface.o mod_pbl_interface.o mod_ncout.o mod_bdycod.o mod_precip.o mod_split.o mod_savefile.o mod_mppio.o mod_slabocean.o
endif
if DO_COMPILE_LIBCLM
mod_params.o: mod_params.F90 mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_domain.o mod_cu_interface.o mod_lm_interface.o mod_atm_interface.o mod_che_interface.o mod_rad_interface.o mod_pbl_interface.o mod_precip.o mod_cloud_s1.o mod_split.o mod_slice.o mod_bdycod.o mod_ncio.o mod_tendency.o mod_ncio.o mod_ncout.o mod_advection.o mod_mppio.o mod_sldepparam.o mod_slabocean.o mod_clm.o clm_varsur.o
else
mod_params.o: mod_params.F90 mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_domain.o mod_cu_interface.o mod_lm_interface.o mod_atm_interface.o mod_che_interface.o mod_rad_interface.o mod_pbl_interface.o mod_precip.o mod_cloud_s1.o mod_split.o mod_slice.o mod_bdycod.o mod_ncio.o mod_tendency.o mod_ncio.o mod_ncout.o mod_advection.o mod_mppio.o mod_sldepparam.o mod_slabocean.o
endif
mod_pbl_interface.o: mod_pbl_interface.F90 mod_realkinds.o mod_service.o mod_constants.o mod_dynparam.o mod_memutil.o mod_mppparam.o mod_atm_interface.o mod_pbl_common.o mod_pbl_holtbl.o mod_pbl_uwtcm.o mod_runparams.o mod_runparams.o
mod_precip.o: mod_precip.F90 mod_runparams.o mod_memutil.o mod_atm_interface.o mod_mpmessage.o
mod_rad_interface.o: mod_rad_interface.f90 mod_realkinds.o mod_atm_interface.o mod_rad_common.o mod_rad_aerosol.o mod_rad_colmod3.o mod_rrtmg_driver.o mod_rad_o3blk.o mod_rad_outrad.o mod_rad_radiation.o mod_rad_scenarios.o mod_rad_tracer.o
if DO_COMPILE_LIBCLM
mod_regcm_interface.o: mod_regcm_interface.F90 mod_memutil.o mod_service.o mod_che_interface.o mod_atm_interface.o mod_pbl_interface.o mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_ncio.o mod_ncout.o mod_output.o mod_split.o mod_bdycod.o mod_init.o mod_header.o mod_params.o mod_tendency.o mod_tstep.o mod_service.o mod_mppio.o mod_cloud_s1.o mod_sun.o perf_mod.o mod_mtrxclm.o spmdMod.o clm_varsur.o
else
if DO_COMPILE_CPL
mod_regcm_interface.o: mod_regcm_interface.F90 mod_memutil.o mod_service.o mod_che_interface.o mod_atm_interface.o mod_pbl_interface.o mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_ncio.o mod_ncout.o mod_output.o mod_split.o mod_bdycod.o mod_init.o mod_header.o mod_params.o mod_tendency.o mod_tstep.o mod_service.o mod_mppio.o mod_cloud_s1.o mod_sun.o mod_update.o
else
mod_regcm_interface.o: mod_regcm_interface.F90 mod_memutil.o mod_service.o mod_che_interface.o mod_atm_interface.o mod_pbl_interface.o mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_ncio.o mod_ncout.o mod_output.o mod_split.o mod_bdycod.o mod_init.o mod_header.o mod_params.o mod_tendency.o mod_tstep.o mod_service.o mod_mppio.o mod_cloud_s1.o mod_sun.o
endif
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
if DO_COMPILE_LIBCLM
mod_savefile.o: mod_savefile.F90 mod_runparams.o mod_mpmessage.o mod_mppparam.o mod_atm_interface.o mod_lm_interface.o mod_che_interface.o mod_rad_interface.o mod_cu_interface.o mod_pbl_interface.o mod_bdycod.o mod_mppio.o mod_clm.o restFileMod.o restFileMod.o clm_varctl.o clm_time_manager.o
else
mod_savefile.o: mod_savefile.F90 mod_runparams.o mod_mpmessage.o mod_mppparam.o mod_atm_interface.o mod_lm_interface.o mod_che_interface.o mod_rad_interface.o mod_cu_interface.o mod_pbl_interface.o mod_bdycod.o mod_mppio.o
endif
mod_sldepparam.o: mod_sldepparam.f90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o
mod_sladvection.o: mod_sladvection.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_runparams.o mod_atm_interface.o mod_sldepparam.o mod_service.o mod_mpmessage.o mod_mppparam.o
mod_slice.o: mod_slice.F90 mod_runparams.o mod_atm_interface.o mod_che_interface.o mod_pbl_interface.o mod_rad_interface.o mod_constants.o
mod_split.o: mod_split.F90 mod_dynparam.o mod_mppparam.o mod_runparams.o mod_atm_interface.o mod_vmodes.o mod_bdycod.o mod_atm_interface.o mod_memutil.o mod_service.o mod_mppio.o
mod_sun.o: mod_sun.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o mod_atm_interface.o mod_mpmessage.o mod_mppparam.o mod_service.o
if DO_COMPILE_LIBCLM
mod_tendency.o: mod_tendency.F90 mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_stdio.o mod_service.o mod_memutil.o mod_atm_interface.o mod_che_interface.o mod_cu_interface.o mod_lm_interface.o mod_rad_interface.o mod_pbl_interface.o mod_bdycod.o mod_precip.o mod_slice.o mod_sun.o mod_advection.o mod_diffusion.o mod_mppio.o mod_domain.o mod_cloud_s1.o mod_sladvection.o mod_slabocean.o mod_clm.o mod_mtrxclm.o clm_varsur.o
else
mod_tendency.o: mod_tendency.F90 mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_stdio.o mod_service.o mod_memutil.o mod_atm_interface.o mod_che_interface.o mod_cu_interface.o mod_lm_interface.o mod_rad_interface.o mod_pbl_interface.o mod_bdycod.o mod_precip.o mod_slice.o mod_sun.o mod_advection.o mod_diffusion.o mod_mppio.o mod_domain.o mod_cloud_s1.o mod_sladvection.o mod_slabocean.o
endif
mod_tstep.o: mod_tstep.f90 mod_runparams.o
mod_vmodes.o: mod_vmodes.F90 mod_runparams.o mod_memutil.o mod_mpmessage.o mod_service.o linpack.o eispack.o mod_mppparam.o
abort.o: abort.f90
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
#
# CPL model 
#
if DO_COMPILE_CPL
mod_update.o: mod_update.F90 mod_runparams.o mod_dynparam.o mod_mppparam.o \
              mod_constants.o mod_bats_common.o mod_atm_interface.o
endif
#
# RegCM Main Program
#
regcm.o: regcm.f90 mod_regcm_interface.o mod_dynparam.o \
         $(MPPLIB) $(CLMLIB) $(NETLIB) $(CHELIB) $(CBMZLIB) \
         $(BATLIB) $(CUMLIB) $(RADLIB) $(PBLLIB) $(RRTMSWLIB)
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
