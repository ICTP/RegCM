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
	microlib:cloudlib:cumlib:radlib:pbllib:clmlib/clm4.5/:.
vpath %.inc $(REGCM_ROOT)/Share:.

SUBDIRS = mpplib netlib batslib ocnlib chemlib clmlib microlib cloudlib \
          cumlib radlib pbllib

LDADD = $(OMPFLAGS)

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
MICROLIBPTH = ./microlib
CLOUDLIBPTH = ./cloudlib
CUMLIBPTH = ./cumlib
RADLIBPTH = ./radlib
PBLLIBPTH = ./pbllib
CBMZLIBPTH = chemlib/GAS_CBMZ_NEW
RRTMLIBPTH = radlib/RRTMG
PLUMELIBPTH = $(REGCM_ROOT)/external
BMILIBPTH = $(REGCM_ROOT)/external

MPPLIB = $(MPPLIBPTH)/libmpp.a
CLMLIB = $(CLMLIBPTH)/libclm.a
BATLIB = $(BATLIBPTH)/libbats.a
OCNLIB = $(OCNLIBPTH)/libocn.a
NETLIB = $(NETLIBPTH)/libnetlib.a
CHELIB = $(CHELIBPTH)/libchem.a
MICROLIB = $(MICROLIBPTH)/libmicro.a
CLOUDLIB = $(CLOUDLIBPTH)/libcloud.a
CUMLIB = $(CUMLIBPTH)/libcum.a
RADLIB = $(RADLIBPTH)/librad.a
PBLLIB = $(PBLLIBPTH)/libpbl.a
CBMZLIB = $(CBMZLIBPTH)/libcbmz.a
RRTMLIB = $(RRTMLIBPTH)/librrtmg.a
PLUMELIB = $(PLUMELIBPTH)/libaerosol.a
BMILIB = $(BMILIBPTH)/libbmi.a

LMPPLIB = -L$(MPPLIBPTH) -lmpp
LBATLIB = -L$(BATLIBPTH) -lbats
LOCNLIB = -L$(OCNLIBPTH) -locn
LCLMLIB = -L$(CLMLIBPTH) -lclm
if DO_COMPILE_LIBCLM45
LCLMLIB45 = -L$(CLMLIBPTH)/clm4.5 -lclm45
endif
LNETLIB = -L$(NETLIBPTH) -lnetlib
LCHELIB = -L$(CHELIBPTH) -lchem
LMICROLIB = -L$(MICROLIBPTH) -lmicro
LCLOUDLIB = -L$(CLOUDLIBPTH) -lcloud
LCUMLIB = -L$(CUMLIBPTH) -lcum
LRADLIB = -L$(RADLIBPTH) -lrad
LPBLLIB = -L$(PBLLIBPTH) -lpbl
LCBMZLIB = -L$(CBMZLIBPTH) -lcbmz
LRRTMLIB = -L$(RRTMLIBPTH) -lrrtmg
LPLUMELIB = -L$(PLUMELIBPTH) -laerosol
LBMILIB = -L$(BMILIBPTH) -lbmi

if COMPILER_SUN
  CPPFLAGS += -M$(MPPLIBPTH) -M$(CLMLIBPTH) -M$(NETLIBPTH) -M$(CBMZLIBPTH) \
              -M$(CHELIBPTH) -M$(BATLIBPTH) -M$(OCNLIBPTH) -M$(CUMLIBPTH)  \
              -M$(MICROLIBPTH) -M$(CLOUDLIBPTH) -M$(RADLIBPTH) -M$(PBLLIBPTH) \
              -M$(RRTMLIBPTH) -M$(CLMLIBPTH)/clm4.5
else
  CPPFLAGS += -I$(MPPLIBPTH) -I$(CLMLIBPTH) -I$(NETLIBPTH) -I$(CBMZLIBPTH) \
              -I$(CHELIBPTH) -I$(BATLIBPTH) -I$(OCNLIBPTH) -I$(CUMLIBPTH)  \
              -I$(MICROLIBPTH) -I$(CLOUDLIBPTH) -I$(RADLIBPTH) -I$(PBLLIBPTH) \
              -I$(RRTMLIBPTH) -I$(CLMLIBPTH)/clm4.5
endif

LIBS = $(LCLMLIB) $(LCHELIB) $(LCBMZLIB) $(LOCNLIB) $(LMICROLIB) \
       $(CLOUDLIB) $(LCUMLIB) $(LRADLIB) $(LPBLLIB) $(LRRTMLIB) \
       $(LCLMLIB45) $(LMPPLIB) $(SYSLIBS) $(LNETLIB) $(LBATLIB) \
       $(LPLUMELIB) $(LBMILIB)

if DO_COMPILE_LIBCLM45
bin_PROGRAMS = regcm clmsa
else
bin_PROGRAMS = regcm
endif

clmsa_SOURCES = clmsa.F90 abort.F90 mod_header.F90 mod_clm_params.F90 \
  mod_atm_stub.F90 mod_ncio.F90 mod_sun.F90

regcm_SOURCES = mod_regcm_interface.F90 mod_atm_interface.F90 \
  mod_che_interface.F90 mod_lm_interface.F90 mod_cu_interface.F90 \
  mod_rad_interface.F90 mod_pbl_interface.F90 mod_header.F90 \
  mod_advection.F90 mod_ncio.F90 mod_bdycod.F90 mod_diffusion.F90 \
  mod_micro_interface.F90 mod_sun.F90 mod_savefile.F90 mod_slice.F90 \
  mod_init.F90 mod_vmodes.F90 mod_split.F90 mod_output.F90 \
  mod_tendency.F90 mod_params.F90 mod_slabocean.F90 mod_sldepparam.F90 \
  mod_sladvection.F90 mod_oasis.F90 abort.F90 mod_sound.F90 \
  mod_timefilter.F90 mod_massck.F90 mod_moloch.F90 mod_bmiregcm.F90 regcm.F90 

FCLD = $(MPIFC)

if DO_COMPILE_CPL
regcm_SOURCES += mod_update.F90
endif

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

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
mod_advection.o : mod_advection.F90 mod_dynparam.o mod_memutil.o mod_mpmessage.o mod_regcm_types.o mod_atm_interface.o mod_runparams.o mod_service.o mod_constants.o
mod_atm_interface.o : mod_atm_interface.F90 mod_constants.o mod_dynparam.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_regcm_types.o mod_runparams.o mod_service.o mod_stdio.o
mod_bdycod.o : mod_bdycod.F90 mod_atm_interface.o mod_che_interface.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_pbl_interface.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_slabocean.o mod_zita.o mod_stdatm.o mod_date.o mod_constants.o mod_stdio.o
mod_che_interface.o : mod_che_interface.F90 mod_atm_interface.o mod_rad_interface.o mod_che_bdyco.o mod_che_bionit.o mod_che_carbonaer.o mod_che_common.o mod_che_cumtran.o mod_che_drydep.o mod_che_dust.o mod_che_emission.o mod_che_linox.o mod_che_indices.o mod_che_mppio.o mod_che_ncio.o mod_che_param.o mod_che_species.o mod_che_start.o mod_che_tend.o mod_realkinds.o mod_runparams.o mod_regcm_types.o
mod_cu_interface.o : mod_cu_interface.F90 mod_atm_interface.o mod_che_interface.o mod_cu_bm.o mod_cu_common.o mod_cu_em.o mod_cu_grell.o mod_cu_kf.o mod_cu_kuo.o mod_cu_tables.o mod_cu_tiedtke.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_regcm_types.o mod_constants.o mod_cu_bm.o mod_cu_common.o mod_cu_em.o mod_cu_grell.o mod_cu_kuo.o mod_cu_tables.o mod_cu_tiedtke.o mod_cu_kf.o mod_cu_shallow.o mod_realkinds.o mod_mppparam.o mod_constants.o mod_stdio.o clwfromt.inc
mod_diffusion.o : mod_diffusion.F90 mod_atm_interface.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_service.o mod_constants.o mod_stdio.o
mod_header.o : mod_header.F90 mod_constants.o mod_date.o mod_intkinds.o mod_dynparam.o mod_mppparam.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/version
mod_init.o : mod_init.F90 mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_constants.o mod_cu_interface.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_mpmessage.o mod_mppparam.o mod_stdio.o mod_ncio.o mod_outvars.o mod_pbl_interface.o mod_micro_interface.o mod_rad_interface.o mod_realkinds.o mod_runparams.o mod_savefile.o mod_service.o mod_slabocean.o mod_sun.o rrtmg_lw_init.o rrtmg_sw_init.o mod_slice.o mod_diffusion.o mod_sound.o mod_massck.o mod_zita.o mod_moloch.o pfwsat.inc pfesat.inc
if DO_COMPILE_LIBCLM
mod_lm_interface.o : mod_lm_interface.F90 clm_time_manager.o clm_varctl.o clm_varsur.o mod_atm_interface.o mod_bats_common.o mod_dynparam.o mod_che_interface.o mod_clm.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_mtrxclm.o mod_ocn_common.o mod_outvars.o mod_regcm_types.o mod_runparams.o mod_service.o mod_date.o perf_mod.o restFileMod.o spmdMod.o mod_stdio.o pfesat.inc pfwsat.inc wlh.inc evpt.inc
else
if  DO_COMPILE_LIBCLM45
mod_lm_interface.o : mod_lm_interface.F90 mod_atm_interface.o mod_bats_common.o mod_che_interface.o mod_dynparam.o mod_clm_regcm.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ocn_common.o mod_outvars.o mod_regcm_types.o mod_runparams.o mod_service.o mod_slice.o mod_stdio.o mod_date.o pfesat.inc pfwsat.inc wlh.inc evpt.inc pfdesatdt.inc pqderiv.inc
else
mod_lm_interface.o : mod_lm_interface.F90 mod_atm_interface.o mod_bats_common.o mod_che_interface.o mod_dynparam.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ocn_common.o mod_outvars.o mod_regcm_types.o mod_runparams.o mod_service.o mod_slice.o mod_stdio.o mod_date.o pfesat.inc pfwsat.inc wlh.inc evpt.inc pfdesatdt.inc pqderiv.inc
endif
endif
mod_moloch.o : mod_moloch.F90 mod_advection.o mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_constants.o mod_cu_interface.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_massck.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_pbl_interface.o mod_micro_interface.o mod_rad_interface.o mod_realkinds.o mod_runparams.o mod_service.o mod_slabocean.o mod_slice.o mod_stdatm.o mod_stdio.o mod_sun.o cpmf.inc wlh.inc
mod_ncio.o : mod_ncio.F90 mod_intkinds.o mod_realkinds.o mod_date.o mod_constants.o mod_stdio.o mod_domain.o mod_dynparam.o mod_ensemble.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_nchelper.o mod_runparams.o pfesat.inc
mod_oasis.o : mod_oasis.F90 mod_atm_interface.o mod_realkinds.o
mod_output.o : mod_output.F90 mod_atm_interface.o mod_bdycod.o mod_capecin.o mod_che_interface.o mod_che_output.o mod_cu_interface.o mod_dynparam.o mod_header.o mod_intkinds.o mod_lm_interface.o mod_mpmessage.o mod_mppparam.o mod_ncout.o mod_pbl_interface.o mod_micro_interface.o mod_rad_interface.o mod_realkinds.o mod_runparams.o mod_savefile.o mod_service.o mod_slabocean.o mod_moloch.o mod_split.o mod_projections.o pfwsat.inc pfesat.inc
if DO_COMPILE_LIBCLM45
mod_params.o : mod_params.F90 mod_advection.o mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_clm_regcm.o mod_cu_interface.o mod_domain.o mod_lm_interface.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_ncout.o mod_pbl_interface.o mod_micro_interface.o mod_rad_interface.o mod_runparams.o mod_savefile.o mod_service.o mod_slabocean.o mod_sldepparam.o mod_slice.o mod_split.o mod_tendency.o mod_nhinterp.o mod_diffusion.o mod_sladvection.o mod_timer.o mod_moloch.o mod_zita.o
else
mod_params.o : mod_params.F90 mod_advection.o mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_cu_interface.o mod_domain.o mod_lm_interface.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_ncout.o mod_pbl_interface.o mod_micro_interface.o mod_rad_interface.o mod_runparams.o mod_savefile.o mod_service.o mod_slabocean.o mod_sldepparam.o mod_slice.o mod_split.o mod_tendency.o mod_nhinterp.o mod_sound.o mod_diffusion.o mod_sladvection.o mod_timer.o mod_moloch.o mod_zita.o
endif
mod_pbl_interface.o : mod_pbl_interface.F90 mod_atm_interface.o mod_che_common.o mod_che_interface.o mod_constants.o mod_dynparam.o mod_memutil.o mod_mppparam.o mod_pbl_common.o mod_pbl_holtbl.o mod_pbl_uwtcm.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_pbl_gfs.o mod_pbl_myj.o
mod_micro_interface.o : mod_micro_interface.F90 mod_realkinds.o mod_service.o mod_constants.o mod_dynparam.o mod_memutil.o mod_mppparam.o mod_regcm_types.o mod_runparams.o mod_micro_nogtom.o mod_micro_subex.o mod_micro_wsm5.o mod_cloud_subex.o mod_cloud_xuran.o mod_cloud_thomp.o mod_cloud_guli2007.o mod_cloud_tompkins.o mod_cloud_echam5.o pfesat.inc pfwsat.inc clwfromt.inc mod_atm_interface.o mod_che_interface.o mod_stdio.o
mod_massck.o : mod_massck.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_runparams.o mod_mppparam.o mod_atm_interface.o mod_stdio.o
mod_rad_interface.o : mod_rad_interface.F90 mod_atm_interface.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_rad_aerosol.o mod_rad_colmod3.o mod_rad_common.o mod_rad_o3blk.o mod_rad_outrad.o mod_rad_radiation.o mod_ipcc_scenario.o mod_realkinds.o mod_regcm_types.o mod_rrtmg_driver.o mod_runparams.o mod_date.o
mod_timefilter.o : mod_timefilter.F90 mod_realkinds.o mod_intkinds.o mod_constants.o mod_dynparam.o mod_runparams.o
if DO_COMPILE_CPL
mod_regcm_interface.o : mod_regcm_interface.F90 mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_header.o mod_init.o mod_lm_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_ncout.o mod_output.o mod_params.o mod_pbl_interface.o mod_runparams.o mod_service.o mod_split.o mod_sun.o mod_tendency.o mod_update.o mod_moloch.o
else
mod_regcm_interface.o : mod_regcm_interface.F90 mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_header.o mod_init.o mod_lm_interface.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_ncio.o mod_ncout.o mod_output.o mod_params.o mod_pbl_interface.o mod_runparams.o mod_service.o mod_split.o mod_tendency.o mod_moloch.o
endif
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_savefile.o : mod_savefile.F90 mod_atm_interface.o mod_che_interface.o mod_che_mppio.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_massck.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_nchelper.o mod_date.o mod_stdio.o
mod_slabocean.o : mod_slabocean.F90 mod_atm_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mppparam.o mod_outvars.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_stdio.o
mod_sladvection.o : mod_sladvection.F90 mod_atm_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_service.o mod_sldepparam.o mod_stdio.o
mod_sldepparam.o : mod_sldepparam.F90 mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o
mod_slice.o : mod_slice.F90 mod_atm_interface.o mod_che_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_pbl_interface.o mod_rad_interface.o mod_realkinds.o mod_runparams.o pfesat.inc pfwsat.inc
mod_sound.o : mod_sound.F90 mod_atm_interface.o mod_cu_interface.o mod_timefilter.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o cpmf.inc
mod_split.o : mod_split.F90 mod_atm_interface.o mod_bdycod.o mod_dynparam.o mod_memutil.o mod_mppparam.o mod_runparams.o mod_service.o mod_vmodes.o mod_stdio.o mod_constants.o mod_intkinds.o mod_realkinds.o
mod_sun.o : mod_sun.F90 mod_atm_interface.o mod_constants.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_service.o mod_sunorbit.o mod_memutil.o mod_stdio.o
mod_tendency.o : mod_tendency.F90 mod_advection.o mod_atm_interface.o mod_bdycod.o mod_che_interface.o mod_constants.o mod_cu_interface.o mod_diffusion.o mod_domain.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_massck.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_pbl_interface.o mod_micro_interface.o mod_rad_interface.o mod_realkinds.o mod_runparams.o mod_service.o mod_slabocean.o mod_sladvection.o mod_slice.o mod_sound.o mod_stdio.o mod_timefilter.o mod_sun.o cpmf.inc
mod_vmodes.o : mod_vmodes.F90 eispack.o linpack.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_service.o mod_constants.o mod_stdio.o
#
# CPL model 
#
if DO_COMPILE_CPL
mod_update.o : mod_update.F90 mod_atm_interface.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_lm_interface.o mod_memutil.o mod_realkinds.o mod_regcm_types.o mod_runparams.o
endif
#
# RegCM Main Program
#
regcm.o : regcm.F90 mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_regcm_interface.o mod_runparams.o mod_date.o mod_constants.o mod_stdio.o \
         $(MPPLIB) $(CLMLIB) $(NETLIB) $(CHELIB) $(CBMZLIB) $(BATLIB) \
	 $(OCNLIB) $(MICROLIB) $(CLOUDLIB) $(CUMLIB) $(RADLIB) $(PBLLIB) \
         $(RRTMLIB) $(PLUMELIB) $(BMILIB)
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<

mod_atm_stub.o : mod_atm_stub.F90 mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_stdio.o mod_service.o mod_memutil.o mod_regcm_types.o mod_ncio.o mod_zita.o mod_stdatm.o mod_date.o mod_stdio.o
mod_clm_params.o : mod_clm_params.F90 mod_runparams.o mod_mppparam.o mod_mpmessage.o mod_domain.o mod_service.o mod_clm_regcm.o mod_ncio.o mod_timer.o mod_realkinds.o mod_intkinds.o mod_date.o mod_atm_stub.o mod_dynparam.o mod_ipcc_scenario.o mod_sun.o mod_zita.o mod_constants.o mod_stdio.o
clmsa.o : clmsa.F90 mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_regcm_interface.o mod_runparams.o mod_memutil.o mod_service.o mod_header.o mod_params.o mod_date.o mod_dynparam.o mod_atm_stub.o mod_clm_params.o mod_clm_regcm.o mod_ncio.o mod_date.o mod_stdio.o mod_constants.o $(MPPLIB) $(CLMLIB) $(NETLIB) $(BATLIB)
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
