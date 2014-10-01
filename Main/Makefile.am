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

vpath %.o $(REGCM_ROOT)/Share:mpplib:netlib:batslib:ocnlib:chemlib:clmlib:\
	cumlib:radlib:pbllib:clmlib/clm4.5/:.

SUBDIRS = mpplib netlib batslib ocnlib chemlib clmlib cumlib radlib pbllib

include $(REGCM_ROOT)/makeinc

if REGCM_PARALLEL
rcsuffix = MPI
else
rcsuffix = Serial
endif

MPPLIBPTH = ./mpplib
CLMLIBPTH = ./clmlib
BATLIBPTH = ./batslib
OCNLIBPTH = ./ocnlib
NETLIBPTH = ./netlib
CHELIBPTH = ./chemlib
CUMLIBPTH = ./cumlib
RADLIBPTH = ./radlib
PBLLIBPTH = ./pbllib
CBMZLIBPTH = chemlib/GAS_CBMZ_NEW
RRTMSWLIBPTH = radlib/RRTMG_SW

MPPLIB = $(MPPLIBPTH)/libmpp.a
CLMLIB = $(CLMLIBPTH)/libclm.a
BATLIB = $(BATLIBPTH)/libbats.a
OCNLIB = $(OCNLIBPTH)/libocn.a
NETLIB = $(NETLIBPTH)/libnetlib.a
CHELIB = $(CHELIBPTH)/libchem.a
CUMLIB = $(CUMLIBPTH)/libcum.a
RADLIB = $(RADLIBPTH)/librad.a
PBLLIB = $(PBLLIBPTH)/libpbl.a
CBMZLIB = $(CBMZLIBPTH)/libcbmz.a
RRTMSWLIB = $(RRTMSWLIBPTH)/librrtmgsw.a

LMPPLIB = -L$(MPPLIBPTH) -lmpp
LBATLIB = -L$(BATLIBPTH) -lbats
LOCNLIB = -L$(OCNLIBPTH) -locn
LCLMLIB = -L$(CLMLIBPTH) -lclm
if DO_COMPILE_LIBCLM45
LCLMLIB45 = -L$(CLMLIBPTH)/clm4.5 -lclm45
endif
LNETLIB = -L$(NETLIBPTH) -lnetlib
LCHELIB = -L$(CHELIBPTH) -lchem
LCUMLIB = -L$(CUMLIBPTH) -lcum
LRADLIB = -L$(RADLIBPTH) -lrad
LPBLLIB = -L$(PBLLIBPTH) -lpbl
LCBMZLIB = -L$(CBMZLIBPTH) -lcbmz
LRRTMSWLIB = -L$(RRTMSWLIBPTH) -lrrtmgsw

if COMPILER_SUN
  CPPFLAGS += -M$(MPPLIBPTH) -M$(CLMLIBPTH) -M$(NETLIBPTH) -M$(CBMZLIBPTH) \
              -M$(CHELIBPTH) -M$(BATLIBPTH) -M$(OCNLIBPTH) -M$(CUMLIBPTH)  \
              -M$(RADLIBPTH) -M$(PBLLIBPTH) -M$(RRTMSWLIBPTH)
else
  CPPFLAGS += -I$(MPPLIBPTH) -I$(CLMLIBPTH) -I$(NETLIBPTH) -I$(CBMZLIBPTH) \
              -I$(CHELIBPTH) -I$(BATLIBPTH) -I$(OCNLIBPTH) -I$(CUMLIBPTH)  \
              -I$(RADLIBPTH) -I$(PBLLIBPTH) -I$(RRTMSWLIBPTH) \
              -I$(CLMLIBPTH)/clm4.5
endif

LIBS = $(LMPPLIB) $(LCLMLIB) $(LCHELIB) $(LCBMZLIB) $(LBATLIB) \
       $(LOCNLIB) $(LCUMLIB) $(LRADLIB) $(LPBLLIB) $(LRRTMSWLIB) \
       $(LCLMLIB45) $(SYSLIBS) $(LNETLIB)

bin_PROGRAMS = regcm

regcm_SOURCES = mod_regcm_interface.F90 mod_atm_interface.F90 mod_che_interface.F90 \
  mod_lm_interface.F90 mod_cu_interface.F90 mod_rad_interface.F90 \
  mod_pbl_interface.F90 mod_header.F90 mod_cloud_variables.F90 \
  mod_advection.F90 mod_ncio.F90 mod_bdycod.F90 mod_diffusion.F90 \
  mod_precip.F90 mod_sun.F90 mod_savefile.F90 mod_slice.F90 mod_init.F90 \
  mod_vmodes.F90 mod_split.F90 mod_output.F90 mod_cloud_s1.F90 \
  mod_tendency.F90 mod_params.F90 mod_tstep.F90 mod_slabocean.F90 \
  mod_sldepparam.F90 mod_sladvection.F90 mod_oasis.F90 abort.F90 \
  mod_sound.F90 regcm.F90 

FCLD = $(MPIFC)

if DO_COMPILE_CPL
regcm_SOURCES += mod_update.F90
endif

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
#
# Modules
#
abort.o : abort.F90 mod_dynparam.o mod_intkinds.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_advection.o : mod_advection.F90 mod_dynparam.o mod_memutil.o mod_mpmessage.o mod_regcm_types.o mod_runparams.o mod_service.o
mod_atm_interface.o : mod_atm_interface.F90 mod_constants.o mod_dynparam.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_regcm_types.o mod_runparams.o mod_service.o mod_stdio.o
mod_bdycod.o : mod_bdycod.F90 mod_atm_interface.o mod_che_interface.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_pbl_interface.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_slabocean.o
mod_che_interface.o : mod_che_interface.F90 mod_che_bdyco.o mod_che_bionit.o mod_che_carbonaer.o mod_che_common.o mod_che_cumtran.o mod_che_drydep.o mod_che_dust.o mod_che_emission.o mod_che_indices.o mod_che_mppio.o mod_che_ncio.o mod_che_param.o mod_che_species.o mod_che_start.o mod_che_tend.o mod_realkinds.o mod_regcm_types.o
mod_cloud_variables.o : mod_cloud_variables.F90 mod_constants.o mod_realkinds.o
mod_cloud_s1.o : mod_cloud_s1.F90 mod_atm_interface.o mod_cloud_variables.o mod_constants.o mod_dynparam.o mod_memutil.o mod_mpmessage.o mod_precip.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o
mod_cu_interface.o : mod_cu_interface.F90 mod_atm_interface.o mod_che_interface.o mod_cu_bm.o mod_cu_common.o mod_cu_em.o mod_cu_grell.o mod_cu_kf.o mod_cu_kuo.o mod_cu_tables.o mod_cu_tiedtke.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_regcm_types.o mod_runparams.o
mod_cu_interface.o : mod_cu_interface.F90 mod_regcm_types.o mod_constants.o mod_cu_bm.o mod_cu_common.o mod_cu_em.o mod_cu_grell.o mod_cu_kuo.o mod_cu_tables.o mod_cu_tiedtke.o mod_cu_kf.o mod_realkinds.o
mod_diffusion.o : mod_diffusion.F90 mod_atm_interface.o mod_dynparam.o mod_intkinds.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_service.o
mod_header.o : mod_header.F90 mod_constants.o mod_date.o mod_intkinds.o mod_mppparam.o mod_realkinds.o mod_stdio.o
mod_init.o : mod_init.F90 mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_constants.o mod_cu_interface.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_outvars.o mod_pbl_interface.o mod_precip.o mod_rad_interface.o mod_realkinds.o mod_runparams.o mod_savefile.o mod_service.o mod_slabocean.o mod_sun.o rrtmg_lw_init.o rrtmg_sw_init.o
if DO_COMPILE_LIBCLM
mod_lm_interface.o : mod_lm_interface.F90 clm_time_manager.o clm_varctl.o clm_varsur.o mod_atm_interface.o mod_bats_common.o mod_che_common.o mod_che_interface.o mod_clm.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_mtrxclm.o mod_ocn_common.o mod_outvars.o mod_regcm_types.o mod_runparams.o mod_service.o mod_slice.o perf_mod.o restFileMod.o spmdMod.o
else
if  DO_COMPILE_LIBCLM45
mod_lm_interface.o : mod_lm_interface.F90 mod_atm_interface.o mod_bats_common.o mod_che_common.o mod_che_interface.o mod_clm_regcm.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ocn_common.o mod_outvars.o mod_regcm_types.o mod_runparams.o mod_service.o mod_slice.o
else
mod_lm_interface.o : mod_lm_interface.F90 mod_atm_interface.o mod_bats_common.o mod_che_common.o mod_che_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ocn_common.o mod_outvars.o mod_regcm_types.o mod_runparams.o mod_service.o mod_slice.o
endif
endif
mod_ncio.o : mod_ncio.F90 mod_domain.o mod_dynparam.o mod_ensemble.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_nchelper.o mod_runparams.o
mod_oasis.o : mod_oasis.F90 mod_atm_interface.o mod_realkinds.o
mod_output.o : mod_output.F90 mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_che_output.o mod_cloud_s1.o mod_cu_interface.o mod_dynparam.o mod_header.o mod_intkinds.o mod_lm_interface.o mod_mpmessage.o mod_mppparam.o mod_ncout.o mod_pbl_interface.o mod_precip.o mod_rad_interface.o mod_realkinds.o mod_runparams.o mod_savefile.o mod_service.o mod_slabocean.o mod_split.o
if DO_COMPILE_LIBCLM45
mod_params.o : mod_params.F90 mod_advection.o mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_clm_regcm.o mod_cloud_s1.o mod_cu_interface.o mod_domain.o mod_lm_interface.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_ncout.o mod_pbl_interface.o mod_precip.o mod_rad_interface.o mod_runparams.o mod_savefile.o mod_service.o mod_slabocean.o mod_sldepparam.o mod_slice.o mod_split.o mod_tendency.o
else
mod_params.o : mod_params.F90 mod_advection.o mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_cloud_s1.o mod_cu_interface.o mod_domain.o mod_lm_interface.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_ncout.o mod_pbl_interface.o mod_precip.o mod_rad_interface.o mod_runparams.o mod_savefile.o mod_service.o mod_slabocean.o mod_sldepparam.o mod_slice.o mod_split.o mod_tendency.o
endif
mod_pbl_interface.o : mod_pbl_interface.F90 mod_atm_interface.o mod_che_common.o mod_che_interface.o mod_constants.o mod_dynparam.o mod_memutil.o mod_mppparam.o mod_pbl_common.o mod_pbl_holtbl.o mod_pbl_uwtcm.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o
mod_precip.o : mod_precip.F90 mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_realkinds.o mod_regcm_types.o mod_runparams.o
mod_rad_interface.o : mod_rad_interface.F90 mod_atm_interface.o mod_che_interface.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_rad_aerosol.o mod_rad_colmod3.o mod_rad_common.o mod_rad_o3blk.o mod_rad_outrad.o mod_rad_radiation.o mod_rad_scenarios.o mod_realkinds.o mod_regcm_types.o mod_rrtmg_driver.o mod_runparams.o
if DO_COMPILE_CPL
mod_regcm_interface.o : mod_regcm_interface.F90 mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_cloud_s1.o mod_header.o mod_init.o mod_lm_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_ncout.o mod_output.o mod_params.o mod_pbl_interface.o mod_runparams.o mod_service.o mod_split.o mod_sun.o mod_tendency.o mod_tstep.o mod_update.o
else
mod_regcm_interface.o : mod_regcm_interface.F90 mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_cloud_s1.o mod_header.o mod_init.o mod_lm_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_ncout.o mod_output.o mod_params.o mod_pbl_interface.o mod_runparams.o mod_service.o mod_split.o mod_sun.o mod_tendency.o mod_tstep.o
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_savefile.o : mod_savefile.F90 mod_che_interface.o mod_che_mppio.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o
mod_slabocean.o : mod_slabocean.F90 mod_atm_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mppparam.o mod_outvars.o mod_realkinds.o mod_regcm_types.o mod_runparams.o
mod_sladvection.o : mod_sladvection.F90 mod_atm_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_service.o mod_sldepparam.o
mod_sldepparam.o : mod_sldepparam.F90 mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o
mod_slice.o : mod_slice.F90 mod_atm_interface.o mod_che_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_pbl_interface.o mod_rad_interface.o mod_realkinds.o mod_runparams.o
mod_sound.o : mod_sound.F90 mod_atm_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o
mod_split.o : mod_split.F90 mod_atm_interface.o mod_bdycod.o mod_dynparam.o mod_memutil.o mod_mppparam.o mod_runparams.o mod_service.o mod_vmodes.o
mod_sun.o : mod_sun.F90 mod_atm_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_service.o mod_sunorbit.o
mod_tendency.o : mod_tendency.F90 mod_advection.o mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_cloud_s1.o mod_constants.o mod_cu_interface.o mod_diffusion.o mod_domain.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_pbl_interface.o mod_precip.o mod_rad_interface.o mod_realkinds.o mod_runparams.o mod_service.o mod_slabocean.o mod_sladvection.o mod_slice.o mod_stdio.o mod_sun.o
mod_tstep.o : mod_tstep.F90 mod_constants.o mod_intkinds.o mod_realkinds.o mod_runparams.o
mod_vmodes.o : mod_vmodes.F90 eispack.o linpack.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_service.o
#
# CPL model 
#
if DO_COMPILE_CPL
mod_update.o : mod_update.F90 mod_atm_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_memutil.o mod_realkinds.o mod_regcm_types.o
endif
#
# RegCM Main Program
#
regcm.o : regcm.F90 mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_regcm_interface.o mod_runparams.o \
         $(MPPLIB) $(CLMLIB) $(NETLIB) $(CHELIB) $(CBMZLIB) \
         $(BATLIB) $(OCNLIB) $(CUMLIB) $(RADLIB) $(PBLLIB) $(RRTMSWLIB)
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
