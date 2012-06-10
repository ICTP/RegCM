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

vpath %.mod $(REGCM_ROOT)/Share:mpplib:netlib:batslib:chemlib:clmlib:\
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

regcm_SOURCES = mod_atm_interface.f90 mod_che_interface.F90 \
  mod_lm_interface.F90 mod_cu_interface.f90 mod_rad_interface.f90 \
  mod_pbl_interface.F90 mod_mppio.F90 mod_diagnosis.F90 mod_header.F90 \
  mod_advection.F90 mod_ncio.F90 mod_bdycod.F90 mod_diffusion.F90 \
  mod_precip.f90 mod_sun.F90 mod_savefile.F90 mod_slice.F90 mod_init.F90 \
  mod_vmodes.F90 mod_split.F90 mod_output.F90 mod_cloud_s1.F90 \
  mod_tendency.F90 mod_params.F90 mod_tstep.f90 regcm.f90 \
  mod_regcm_interface.F90

FCLD = $(MPIFC)

if DO_COMPILE_LIBESMF
LIBS += $(ESMFLIBPTH)
CPPFLAGS += $(ESMFINCPTH) -I$(ROMS_PREFIX)
regcm_SOURCES += mod_couplerr.F90 mod_esmf_atm.F90 mod_esmf_ocn.F90 \
                 mod_esmf_cpl.F90 regcm_esmf.F90
skip_files  = $(ROMS_PREFIX)/master.o
regcm_LDADD = $(filter-out $(skip_files), $(wildcard $(ROMS_PREFIX)/*.o))
else
regcm_SOURCES += regcm_nonesmf.F90
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
mod_advection.o: mod_advection.F90 mod_atm_interface.o mod_dynparam.mod mod_runparams.mod mod_memutil.mod mod_service.mod
mod_atm_interface.o: mod_atm_interface.F90 mod_dynparam.mod mod_stdio.mod mod_constants.mod mod_runparams.mod mod_mppparam.mod mod_service.mod mod_memutil.mod
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_bdycod.o: mod_bdycod.F90 mod_dynparam.mod mod_mppparam.mod mod_memutil.mod mod_atm_interface.o mod_pbl_interface.o mod_che_interface.o mod_lm_interface.o mod_mpmessage.mod mod_ncio.o mod_mppio.o mod_service.mod
mod_che_interface.o: mod_che_interface.F90 mod_realkinds.mod mod_atm_interface.o mod_che_common.mod mod_che_cumtran.mod mod_che_dust.mod mod_che_indices.mod mod_che_mppio.mod mod_che_ncio.mod mod_che_param.mod mod_che_drydep.mod mod_che_bdyco.mod mod_che_emission.mod mod_che_carbonaer.mod mod_che_species.mod mod_che_tend.mod mod_che_start.mod
mod_cloud_s1.o: mod_cloud_s1.F90 mod_runparams.mod mod_memutil.mod mod_atm_interface.o
mod_cu_interface.o: mod_cu_interface.f90 mod_constants.mod mod_realkinds.mod mod_cu_common.mod mod_cu_tiedtke.mod mod_cu_tables.mod mod_cu_bm.mod mod_cu_em.mod mod_cu_kuo.mod mod_cu_grell.mod mod_atm_interface.o
mod_diagnosis.o: mod_diagnosis.F90 mod_runparams.mod mod_atm_interface.o mod_che_interface.o mod_domain.mod mod_memutil.mod mod_mppparam.mod mod_mppio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_diffusion.o: mod_diffusion.F90 mod_atm_interface.o mod_runparams.mod mod_mppparam.mod mod_service.mod
mod_header.o: mod_header.F90 mod_constants.mod mod_realkinds.mod mod_date.mod
if DO_COMPILE_LIBCLM
mod_init.o: mod_init.F90 mod_runparams.mod mod_mppparam.mod mod_lm_interface.o mod_atm_interface.o mod_che_interface.o mod_cu_interface.o mod_rad_interface.o mod_pbl_interface.o rrtmg_sw_init.mod rrtmg_lw_init.mod mod_pbl_interface.mod mod_precip.o mod_bdycod.o mod_mpmessage.mod mod_sun.o mod_ncio.o mod_savefile.o mod_diagnosis.o mod_mppio.o mod_slice.o mod_constants.mod mod_clm.mod mod_lm_interface.mod clm_varsur.mod
else
mod_init.o: mod_init.F90 mod_runparams.mod mod_mppparam.mod mod_lm_interface.o mod_atm_interface.o mod_che_interface.o mod_cu_interface.o mod_rad_interface.o mod_pbl_interface.o rrtmg_sw_init.mod rrtmg_lw_init.mod mod_pbl_interface.mod mod_precip.o mod_bdycod.o mod_mpmessage.mod mod_sun.o mod_ncio.o mod_savefile.o mod_diagnosis.o mod_mppio.o mod_slice.o mod_constants.mod mod_lm_interface.mod
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
if DO_COMPILE_LIBCLM
mod_lm_interface.o: mod_lm_interface.F90 mod_bats_common.mod mod_runparams.mod mod_memutil.mod mod_atm_interface.o mod_mtrxclm.mod mod_clm.mod mod_bats_mtrxbats.mod mod_bats_param.mod mod_bats_mppio.mod mod_bats_bndry.mod mod_bats_co2.mod mod_bats_drag.mod mod_bats_lake.mod mod_bats_leaftemp.mod mod_bats_mtrxbats.mod mod_bats_zengocn.mod mod_bats_romsocn.mod
else
mod_lm_interface.o: mod_lm_interface.F90 mod_bats_common.mod mod_runparams.mod mod_memutil.mod mod_atm_interface.o mod_bats_mtrxbats.mod mod_bats_param.mod mod_bats_mppio.mod mod_bats_bndry.mod mod_bats_co2.mod mod_bats_drag.mod mod_bats_lake.mod mod_bats_leaftemp.mod mod_bats_mtrxbats.mod mod_bats_zengocn.mod mod_bats_romsocn.mod
endif
mod_mppio.o: mod_mppio.F90 mod_runparams.mod mod_atm_interface.o mod_lm_interface.o mod_cu_interface.o mod_rad_interface.o mod_pbl_interface.o mod_memutil.mod mod_mpmessage.mod
mod_ncio.o: mod_ncio.F90 mod_runparams.mod mod_mpmessage.mod mod_memutil.mod mod_nchelper.mod mod_domain.mod mod_che_interface.o
if DO_COMPILE_LIBCLM
mod_output.o: mod_output.F90 mod_runparams.mod mod_mpmessage.mod mod_mppparam.mod mod_service.mod mod_atm_interface.o mod_che_interface.o mod_che_output.mod mod_lm_interface.o mod_rad_interface.o mod_cu_interface.o mod_pbl_interface.o mod_ncio.o mod_bdycod.o mod_precip.o mod_split.o mod_savefile.o mod_mppio.o mod_clm.mod
else
mod_output.o: mod_output.F90 mod_runparams.mod mod_mpmessage.mod mod_mppparam.mod mod_service.mod mod_atm_interface.o mod_che_interface.o mod_che_output.mod mod_lm_interface.o mod_rad_interface.o mod_cu_interface.o mod_pbl_interface.o mod_ncio.o mod_bdycod.o mod_precip.o mod_split.o mod_savefile.o mod_mppio.o
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
if DO_COMPILE_LIBCLM
mod_params.o: mod_params.F90 mod_runparams.mod mod_mppparam.mod mod_mpmessage.mod mod_domain.mod mod_cu_interface.o mod_lm_interface.o mod_atm_interface.o mod_che_interface.o mod_rad_interface.o mod_pbl_interface.o mod_precip.o mod_cloud_s1.o mod_split.o mod_slice.o mod_bdycod.o mod_ncio.o mod_diagnosis.o mod_tendency.o mod_ncio.mod mod_advection.o mod_mppio.o mod_clm.mod clm_varsur.mod
else
mod_params.o: mod_params.F90 mod_runparams.mod mod_mppparam.mod mod_mpmessage.mod mod_domain.mod mod_cu_interface.o mod_lm_interface.o mod_atm_interface.o mod_che_interface.o mod_rad_interface.o mod_pbl_interface.o mod_precip.o mod_cloud_s1.o mod_split.o mod_slice.o mod_bdycod.o mod_ncio.o mod_diagnosis.o mod_tendency.o mod_ncio.mod mod_advection.o mod_mppio.o
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_pbl_interface.o: mod_pbl_interface.F90 mod_realkinds.mod mod_service.mod mod_constants.mod mod_dynparam.mod mod_memutil.mod mod_mppparam.mod mod_atm_interface.o mod_pbl_common.mod mod_pbl_holtbl.mod mod_pbl_uwtcm.mod mod_runparams.mod
mod_precip.o: mod_precip.f90 mod_runparams.mod mod_memutil.mod mod_atm_interface.o
mod_rad_interface.o: mod_rad_interface.f90 mod_realkinds.mod mod_atm_interface.o mod_rad_common.mod mod_rad_aerosol.mod mod_rad_colmod3.mod mod_rrtmg_driver.mod mod_rad_o3blk.mod mod_rad_outrad.mod mod_rad_radiation.mod mod_rad_scenarios.mod mod_rad_tracer.mod
if DO_COMPILE_LIBCLM
mod_regcm_interface.o: mod_regcm_interface.F90 mod_memutil.mod mod_service.mod mod_che_interface.o mod_atm_interface.o mod_runparams.mod mod_mppparam.mod mod_mpmessage.mod mod_ncio.o mod_output.o mod_split.o mod_bdycod.o mod_init.o mod_header.o mod_params.o mod_tendency.o mod_tstep.o mod_service.mod mod_mppio.o mod_cloud_s1.o mod_sun.o perf_mod.mod mod_mtrxclm.mod spmdmod.mod clm_varsur.mod
else
mod_regcm_interface.o: mod_regcm_interface.F90 mod_memutil.mod mod_service.mod mod_che_interface.o mod_atm_interface.o mod_runparams.mod mod_mppparam.mod mod_mpmessage.mod mod_ncio.o mod_output.o mod_split.o mod_bdycod.o mod_init.o mod_header.o mod_params.o mod_tendency.o mod_tstep.o mod_service.mod mod_mppio.o mod_cloud_s1.o mod_sun.o
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
if DO_COMPILE_LIBCLM
mod_savefile.o: mod_savefile.F90 mod_runparams.mod mod_mpmessage.mod mod_atm_interface.o mod_lm_interface.o mod_che_interface.o mod_rad_interface.o mod_cu_interface.o mod_pbl_interface.o mod_bdycod.o mod_diagnosis.o mod_mppio.o mod_clm.mod restfilemod.mod clm_varctl.mod clm_time_manager.mod
else
mod_savefile.o: mod_savefile.F90 mod_runparams.mod mod_mpmessage.mod mod_atm_interface.o mod_lm_interface.o mod_che_interface.o mod_rad_interface.o mod_cu_interface.o mod_pbl_interface.o mod_bdycod.o mod_diagnosis.o mod_mppio.o
endif
mod_slice.o: mod_slice.F90 mod_runparams.mod mod_atm_interface.o mod_che_interface.o mod_pbl_interface.o
mod_split.o: mod_split.F90 mod_dynparam.mod mod_mppparam.mod mod_runparams.mod mod_atm_interface.o mod_vmodes.o mod_bdycod.o mod_atm_interface.mod mod_memutil.mod mod_service.mod mod_mppio.o
mod_sun.o: mod_sun.F90 mod_constants.mod mod_dynparam.mod mod_runparams.mod mod_atm_interface.o mod_mpmessage.mod mod_service.mod
if DO_COMPILE_LIBCLM
mod_tendency.o: mod_tendency.F90 mod_runparams.mod mod_mppparam.mod mod_mpmessage.mod mod_stdio.mod mod_service.mod mod_memutil.mod mod_atm_interface.o mod_che_interface.o mod_cu_interface.o mod_lm_interface.o mod_rad_interface.o mod_pbl_interface.o mod_bdycod.o mod_che_bdyco.mod mod_precip.o mod_slice.o mod_sun.o mod_diagnosis.o mod_advection.o mod_diffusion.o mod_mppio.o mod_domain.mod mod_cloud_s1.o mod_clm.mod mod_mtrxclm.mod clm_varsur.mod
else
mod_tendency.o: mod_tendency.F90 mod_runparams.mod mod_mppparam.mod mod_mpmessage.mod mod_stdio.mod mod_service.mod mod_memutil.mod mod_atm_interface.o mod_che_interface.o mod_cu_interface.o mod_lm_interface.o mod_rad_interface.o mod_pbl_interface.o mod_bdycod.o mod_che_bdyco.mod mod_precip.o mod_slice.o mod_sun.o mod_diagnosis.o mod_advection.o mod_diffusion.o mod_mppio.o mod_domain.mod mod_cloud_s1.o
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_tstep.o: mod_tstep.f90 mod_runparams.mod
mod_vmodes.o: mod_vmodes.F90 mod_runparams.mod mod_memutil.mod mod_mpmessage.mod mod_service.mod linpack.mod eispack.mod
#
# ESMF modules
#
if DO_COMPILE_LIBESMF
regcm_esmf.o: regcm_esmf.F90 mod_couplerr.mod mod_esmf_atm.mod mod_esmf_ocn.mod mod_esmf_cpl.mod
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_couplerr.o: mod_couplerr.F90 mod_realkinds.mod mod_dynparam.mod mod_atm_interface.o
mod_esmf_atm.o: mod_esmf_atm.F90 mod_regcm_interface.o mod_atm_interface.o mod_couplerr.mod mod_runparams.mod mod_dynparam.mod mod_constants.mod mod_date.mod mod_bats_common.mod mod_bats_common.mod mod_dynparam.mod mod_bats_romsocn.mod mod_bats_romsocn.mod
mod_esmf_cpl.o: mod_esmf_cpl.F90 mod_couplerr.mod
mod_esmf_ocn.o: mod_esmf_ocn.F90 mod_regcm_interface.o mod_couplerr.mod
else
regcm_nonesmf.o: regcm_nonesmf.F90 mod_regcm_interface.o mod_runparams.mod
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
endif
#
# RegCM Main Program
#
if DO_COMPILE_LIBESMF
regcm.o: regcm.f90 mod_dynparam.mod $(MPPLIB) $(CLMLIB) $(NETLIB) $(CHELIB) $(CBMZLIB) $(BATLIB) $(CUMLIB) $(RADLIB) $(PBLLIB) $(RRTMSWLIB)
else
regcm.o: regcm.f90 mod_dynparam.mod $(MPPLIB) $(CLMLIB) $(NETLIB) $(CHELIB) $(CBMZLIB) $(BATLIB) $(CUMLIB) $(RADLIB) $(PBLLIB) $(RRTMSWLIB)
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
