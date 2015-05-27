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

include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib: \
          $(REGCM_ROOT)/Main/batslib:$(REGCM_ROOT)/Main/netlib:.

noinst_LIBRARIES = libclm45.a
libclm45_a_SOURCES = mod_clm_nchelper.F90 mod_clm_domain.F90 mod_clm_type.F90 \
  mod_clm_varpar.F90 mod_clm_varcon.F90 mod_clm_varsur.F90 mod_clm_control.F90 \
  mod_clm_time_manager.F90 mod_clm_regcm.F90 mod_clm_varctl.F90 \
  mod_clm_subgridave.F90 mod_clm_decomp.F90 mod_clm_drydep.F90 \
  mod_clm_megan.F90 mod_clm_ch4varcon.F90 mod_clm_atmlnd.F90 \
  mod_clm_activelayer.F90 mod_clm_balancecheck.F90 mod_clm_qsat.F90 \
  mod_clm_frictionvelocity.F90 mod_clm_banddiagonal.F90 \
  mod_clm_baregroundfluxes.F90 mod_clm_biogeophysics1.F90 \
  mod_clm_soiltemperature.F90 mod_clm_biogeophysics2.F90 \
  mod_clm_pftvarcon.F90 mod_clm_pft2col.F90 \
  mod_clm_urbaninput.F90 mod_clm_surfrd.F90 mod_clm_cnallocation.F90 \
  mod_clm_tridiagonal.F90 mod_clm_cnbalancecheck.F90 mod_clm_cnc14decay.F90 \
  mod_clm_cnfire.F90 mod_clm_cncisoflux.F90 mod_clm_cncstateupdate1.F90 \
  mod_clm_cnrest.F90 mod_clm_cncstateupdate2.F90 mod_clm_cncstateupdate3.F90 \
  mod_clm_cnannualupdate.F90 mod_clm_cndvecosystemdynini.F90 \
  mod_clm_cndvestablishment.F90 mod_clm_cndvlight.F90 \
  mod_clm_cnvegstructupdate.F90 mod_clm_cndecompcascadebgc.F90 \
  mod_clm_cndecompcascadecentury.F90 mod_clm_cnverticalprofile.F90 \
  mod_clm_cnnitrifdenitrif.F90 mod_clm_cndecomp.F90 mod_clm_croprest.F90 \
  mod_clm_cndv.F90 mod_clm_cnphenology.F90 mod_clm_cngresp.F90 \
  mod_clm_cngapmortality.F90 mod_clm_cnmresp.F90 mod_clm_cnndynamics.F90 \
  mod_clm_cnnstateupdate1.F90 mod_clm_cnnstateupdate2.F90 \
  mod_clm_cnnstateupdate3.F90 mod_clm_cnprecisioncontrol.F90 \
  mod_clm_cnsetvalue.F90 mod_clm_cnsoillittverttransp.F90 \
  mod_clm_cnsummary.F90 mod_clm_cnwoodproducts.F90 \
  mod_clm_cninispecial.F90 mod_clm_cninitimevar.F90 mod_clm_restfile.F90 \
  mod_clm_canopyfluxes.F90 mod_clm_dust.F90 mod_clm_drydepvelocity.F90 \
  mod_clm_dynland.F90 mod_clm_filter.F90 mod_clm_fracwet.F90 \
  mod_clm_h2osfc.F90 mod_clm_snicar.F90 mod_clm_hydrology1.F90 \
  mod_clm_slakecon.F90 mod_clm_snowhydrology.F90 mod_clm_vicmap.F90 \
  mod_clm_soilhydrology.F90 mod_clm_hydrology2.F90 mod_clm_meganfactors.F90 \
  mod_clm_slakefluxes.F90 mod_clm_slakehydrology.F90 \
  mod_clm_slaketemperature.F90 mod_clm_slakerest.F90 mod_clm_surfacealbedo.F90 \
  mod_clm_surfaceradiation.F90 mod_clm_urban.F90 mod_clm_urbaninit.F90 \
  mod_clm_vocemission.F90 mod_clm_typeinit.F90 mod_clm_driverinit.F90 \
  mod_clm_reweight.F90 mod_clm_subgrid.F90 mod_clm_decompinit.F90 \
  mod_clm_initgridcells.F90 mod_clm_initch4.F90 mod_clm_initslake.F90 \
  mod_clm_mkarbinit.F90 mod_clm_staticecosysdyn.F90 mod_clm_subgridrest.F90 \
  mod_clm_biogeophysrest.F90 mod_clm_initsurfalb.F90 mod_clm_accumul.F90 \
  mod_clm_slakerest.F90 mod_clm_initsurfalb.F90 mod_clm_accflds.F90 \
  mod_clm_histfile.F90 mod_clm_histflds.F90 mod_clm_initialize.F90 \
  mod_clm_pftdyn.F90 mod_clm_soillittverttransp.F90 mod_clm_ch4rest.F90 \
  mod_clm_cnecosystemdyn.F90 mod_clm_ch4.F90 mod_clm_driver.F90 \
  mod_clm_organicfile.F90 mod_clm_initsoilparvic.F90 mod_clm_initimeconst.F90 \
  $(REGCM_ROOT)/makeinc

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
NET = $(REGCMSRC)/netlib
NET_INCLUDE = $(NET)
BATS = $(REGCMSRC)/batslib
BATS_INCLUDE = $(BATS)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) -M$(NET_INCLUDE) -M$(BATS_INCLUDE)
else
FCFLAGS += -I$(MPP_INCLUDE) -I$(NET_INCLUDE) -I$(BATS_INCLUDE)
endif

%.o: %.f90
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_clm_accflds.o : mod_clm_accflds.F90 mod_clm_accumul.o mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_pftvarcon.o mod_clm_surfrd.o mod_clm_time_manager.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_date.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_accumul.o : mod_clm_accumul.F90 mod_clm_decomp.o mod_clm_nchelper.o mod_clm_varcon.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_activelayer.o : mod_clm_activelayer.F90 mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_constants.o mod_date.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_atmlnd.o : mod_clm_atmlnd.F90 mod_clm_ch4varcon.o mod_clm_decomp.o mod_clm_drydep.o mod_clm_megan.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_balancecheck.o : mod_clm_balancecheck.F90 mod_clm_atmlnd.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_banddiagonal.o : mod_clm_banddiagonal.F90 lapack_dgbsv.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_baregroundfluxes.o : mod_clm_baregroundfluxes.F90 mod_clm_atmlnd.o mod_clm_frictionvelocity.o mod_clm_qsat.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_biogeophysics1.o : mod_clm_biogeophysics1.F90 mod_clm_atmlnd.o mod_clm_qsat.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_constants.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_biogeophysics2.o : mod_clm_biogeophysics2.F90 mod_clm_atmlnd.o mod_clm_soiltemperature.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_biogeophysrest.o : mod_clm_biogeophysrest.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_initsurfalb.o mod_clm_mkarbinit.o mod_clm_nchelper.o mod_clm_snicar.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_canopyfluxes.o : mod_clm_canopyfluxes.F90 mod_clm_atmlnd.o mod_clm_cnallocation.o mod_clm_frictionvelocity.o mod_clm_pftvarcon.o mod_clm_qsat.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_ch4.o : mod_clm_ch4.F90 mod_clm_atmlnd.o mod_clm_ch4varcon.o mod_clm_cndecompcascadebgc.o mod_clm_cndecompcascadecentury.o mod_clm_pftvarcon.o mod_clm_subgridave.o mod_clm_tridiagonal.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_ch4rest.o : mod_clm_ch4rest.F90 mod_clm_decomp.o mod_clm_nchelper.o mod_clm_type.o mod_clm_varctl.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_ch4varcon.o : mod_clm_ch4varcon.F90 mod_dynparam.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cnallocation.o : mod_clm_cnallocation.F90 mod_clm_pftvarcon.o mod_clm_subgridave.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cnannualupdate.o : mod_clm_cnannualupdate.F90 mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cnbalancecheck.o : mod_clm_cnbalancecheck.F90 mod_clm_surfrd.o mod_clm_type.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cnc14decay.o : mod_clm_cnc14decay.F90 mod_clm_nchelper.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cncisoflux.o : mod_clm_cncisoflux.F90 mod_clm_type.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cncstateupdate1.o : mod_clm_cncstateupdate1.F90 mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cncstateupdate2.o : mod_clm_cncstateupdate2.F90 mod_clm_type.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cncstateupdate3.o : mod_clm_cncstateupdate3.F90 mod_clm_type.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cndecompcascadebgc.o : mod_clm_cndecompcascadebgc.F90 mod_clm_ch4varcon.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cndecompcascadecentury.o : mod_clm_cndecompcascadecentury.F90 mod_clm_ch4varcon.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cndecomp.o : mod_clm_cndecomp.F90 mod_clm_cnallocation.o mod_clm_cndecompcascadebgc.o mod_clm_cndecompcascadecentury.o mod_clm_cnnitrifdenitrif.o mod_clm_cnverticalprofile.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cndvecosystemdynini.o : mod_clm_cndvecosystemdynini.F90 mod_clm_decomp.o mod_clm_type.o mod_constants.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cndvestablishment.o : mod_clm_cndvestablishment.F90 mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cndv.o : mod_clm_cndv.F90 mod_clm_cndvestablishment.o mod_clm_cndvlight.o mod_clm_cnvegstructupdate.o mod_clm_decomp.o mod_clm_domain.o mod_clm_nchelper.o mod_clm_time_manager.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cndvlight.o : mod_clm_cndvlight.F90 mod_clm_type.o mod_constants.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cnecosystemdyn.o : mod_clm_cnecosystemdyn.F90 mod_clm_cnallocation.o mod_clm_cnannualupdate.o mod_clm_cnc14decay.o mod_clm_cncisoflux.o mod_clm_cncstateupdate1.o mod_clm_cncstateupdate2.o mod_clm_cncstateupdate3.o mod_clm_cndecomp.o mod_clm_cnfire.o mod_clm_cngapmortality.o mod_clm_cngresp.o mod_clm_cnmresp.o mod_clm_cnndynamics.o mod_clm_cnnstateupdate1.o mod_clm_cnnstateupdate2.o mod_clm_cnnstateupdate3.o mod_clm_cnphenology.o mod_clm_cnprecisioncontrol.o mod_clm_cnsetvalue.o mod_clm_cnsoillittverttransp.o mod_clm_cnsummary.o mod_clm_cnvegstructupdate.o mod_clm_cnwoodproducts.o mod_clm_pftdyn.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varctl.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cnfire.o : mod_clm_cnfire.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_domain.o mod_clm_histfile.o mod_clm_nchelper.o mod_clm_pftvarcon.o mod_clm_subgridave.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cngapmortality.o : mod_clm_cngapmortality.F90 mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cngresp.o : mod_clm_cngresp.F90 mod_clm_pftvarcon.o mod_clm_type.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cninispecial.o : mod_clm_cninispecial.F90 mod_clm_cnsetvalue.o mod_clm_decomp.o mod_clm_pftvarcon.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cninitimevar.o : mod_clm_cninitimevar.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_pftvarcon.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cnmresp.o : mod_clm_cnmresp.F90 mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cnndynamics.o : mod_clm_cnndynamics.F90 mod_clm_atmlnd.o mod_clm_pftvarcon.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cnnitrifdenitrif.o : mod_clm_cnnitrifdenitrif.F90 mod_clm_ch4varcon.o mod_clm_cndecompcascadebgc.o mod_clm_cndecompcascadecentury.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_date.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cnnstateupdate1.o : mod_clm_cnnstateupdate1.F90 mod_clm_pftvarcon.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cnnstateupdate2.o : mod_clm_cnnstateupdate2.F90 mod_clm_type.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cnnstateupdate3.o : mod_clm_cnnstateupdate3.F90 mod_clm_type.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cnphenology.o : mod_clm_cnphenology.F90 mod_clm_croprest.o mod_clm_pftvarcon.o mod_clm_surfrd.o mod_clm_time_manager.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_date.o mod_intkinds.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cnprecisioncontrol.o : mod_clm_cnprecisioncontrol.F90 mod_clm_pftvarcon.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cnrest.o : mod_clm_cnrest.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_nchelper.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_constants.o mod_dynparam.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cnsetvalue.o : mod_clm_cnsetvalue.F90 mod_clm_surfrd.o mod_clm_type.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_cnsoillittverttransp.o : mod_clm_cnsoillittverttransp.F90 mod_clm_tridiagonal.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cnsummary.o : mod_clm_cnsummary.F90 mod_clm_cnndynamics.o mod_clm_pftvarcon.o mod_clm_subgridave.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cnvegstructupdate.o : mod_clm_cnvegstructupdate.F90 mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varcon.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_cnverticalprofile.o : mod_clm_cnverticalprofile.F90 mod_clm_pftvarcon.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_cnwoodproducts.o : mod_clm_cnwoodproducts.F90 mod_clm_decomp.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_control.o : mod_clm_control.F90 mod_clm_canopyfluxes.o mod_clm_cnallocation.o mod_clm_cnc14decay.o mod_clm_cndecompcascadebgc.o mod_clm_cndecompcascadecentury.o mod_clm_cnndynamics.o mod_clm_cnnitrifdenitrif.o mod_clm_cnsoillittverttransp.o mod_clm_cnverticalprofile.o mod_clm_histfile.o mod_clm_hydrology1.o mod_clm_megan.o mod_clm_slakecon.o mod_clm_soilhydrology.o mod_clm_surfacealbedo.o mod_clm_urban.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_croprest.o : mod_clm_croprest.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_nchelper.o mod_clm_time_manager.o mod_clm_type.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_decomp.o : mod_clm_decomp.F90 mod_clm_domain.o mod_clm_type.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_regcm_types.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_decompinit.o : mod_clm_decompinit.F90 mod_clm_decomp.o mod_clm_subgrid.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_regcm_types.o mod_service.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_domain.o : mod_clm_domain.F90 mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_driver.o : mod_clm_driver.F90 mod_clm_accflds.o mod_clm_activelayer.o mod_clm_atmlnd.o mod_clm_balancecheck.o mod_clm_baregroundfluxes.o mod_clm_biogeophysics1.o mod_clm_biogeophysics2.o mod_clm_canopyfluxes.o mod_clm_ch4.o mod_clm_cnannualupdate.o mod_clm_cnbalancecheck.o mod_clm_cndv.o mod_clm_cnecosystemdyn.o mod_clm_cnfire.o mod_clm_cnsetvalue.o mod_clm_cnverticalprofile.o mod_clm_decomp.o mod_clm_driverinit.o mod_clm_drydep.o mod_clm_drydepvelocity.o mod_clm_dust.o mod_clm_dynland.o mod_clm_filter.o mod_clm_histfile.o mod_clm_hydrology1.o mod_clm_hydrology2.o mod_clm_pft2col.o mod_clm_pftdyn.o mod_clm_restfile.o mod_clm_reweight.o mod_clm_slakefluxes.o mod_clm_slakehydrology.o mod_clm_slaketemperature.o mod_clm_snicar.o mod_clm_staticecosysdyn.o mod_clm_surfacealbedo.o mod_clm_surfaceradiation.o mod_clm_type.o mod_clm_urban.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_vocemission.o mod_date.o mod_intkinds.o mod_mppparam.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_driverinit.o : mod_clm_driverinit.F90 mod_clm_atmlnd.o mod_clm_domain.o mod_clm_qsat.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_drydep.o : mod_clm_drydep.F90 mod_constants.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_drydepvelocity.o : mod_clm_drydepvelocity.F90 mod_clm_atmlnd.o mod_clm_drydep.o mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varcon.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_dust.o : mod_clm_dust.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_constants.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_dynland.o : mod_clm_dynland.F90 mod_clm_decomp.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_filter.o : mod_clm_filter.F90 mod_clm_decomp.o mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varcon.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_fracwet.o : mod_clm_fracwet.F90 mod_clm_type.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_frictionvelocity.o : mod_clm_frictionvelocity.F90 mod_clm_type.o mod_clm_varcon.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_h2osfc.o : mod_clm_h2osfc.F90 mod_clm_type.o mod_clm_varcon.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_histfile.o : mod_clm_histfile.F90 mod_clm_decomp.o mod_clm_domain.o mod_clm_nchelper.o mod_clm_subgridave.o mod_clm_time_manager.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_histflds.o : mod_clm_histflds.F90 mod_clm_atmlnd.o mod_clm_ch4varcon.o mod_clm_histfile.o mod_clm_megan.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_hydrology1.o : mod_clm_hydrology1.F90 mod_clm_atmlnd.o mod_clm_fracwet.o mod_clm_h2osfc.o mod_clm_snicar.o mod_clm_subgridave.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_hydrology2.o : mod_clm_hydrology2.F90 mod_clm_atmlnd.o mod_clm_snowhydrology.o mod_clm_soilhydrology.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_clm_vicmap.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_initch4.o : mod_clm_initch4.F90 mod_clm_ch4varcon.o mod_clm_decomp.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mppparam.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_initgridcells.o : mod_clm_initgridcells.F90 mod_clm_decomp.o mod_clm_domain.o mod_clm_subgrid.o mod_clm_surfrd.o mod_clm_type.o mod_clm_urbaninput.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_clm_varsur.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_initialize.o : mod_clm_initialize.F90 mod_clm_accflds.o mod_clm_atmlnd.o mod_clm_ch4varcon.o mod_clm_cndvecosystemdynini.o mod_clm_cnecosystemdyn.o mod_clm_cninitimevar.o mod_clm_control.o mod_clm_decompinit.o mod_clm_decomp.o mod_clm_domain.o mod_clm_drydep.o mod_clm_dust.o mod_clm_filter.o mod_clm_histfile.o mod_clm_histflds.o mod_clm_initch4.o mod_clm_initgridcells.o mod_clm_initimeconst.o mod_clm_initslake.o mod_clm_initsurfalb.o mod_clm_mkarbinit.o mod_clm_nchelper.o mod_clm_pftdyn.o mod_clm_pftvarcon.o mod_clm_restfile.o mod_clm_reweight.o mod_clm_staticecosysdyn.o mod_clm_surfrd.o mod_clm_time_manager.o mod_clm_typeinit.o mod_clm_urbaninit.o mod_clm_urbaninput.o mod_clm_urban.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_clm_varsur.o mod_clm_vocemission.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_stdio.o mod_sunorbit.o $(REGCM_ROOT)/makeinc
mod_clm_initimeconst.o : mod_clm_initimeconst.F90 mod_clm_atmlnd.o mod_clm_ch4varcon.o mod_clm_cndecompcascadebgc.o mod_clm_cndecompcascadecentury.o mod_clm_cninispecial.o mod_clm_decomp.o mod_clm_initsoilparvic.o mod_clm_nchelper.o mod_clm_organicfile.o mod_clm_pftvarcon.o mod_clm_snicar.o mod_clm_soilhydrology.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_clm_varsur.o mod_clm_vicmap.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_initslake.o : mod_clm_initslake.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_slakecon.o mod_clm_snicar.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_initsoilparvic.o : mod_clm_initsoilparvic.F90 mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_initsurfalb.o : mod_clm_initsurfalb.F90 mod_clm_cnecosystemdyn.o mod_clm_cnvegstructupdate.o mod_clm_decomp.o mod_clm_filter.o mod_clm_fracwet.o mod_clm_staticecosysdyn.o mod_clm_surfacealbedo.o mod_clm_type.o mod_clm_urban.o mod_clm_varcon.o mod_clm_varpar.o mod_constants.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_megan.o : mod_clm_megan.F90 mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_meganfactors.o : mod_clm_meganfactors.F90 mod_clm_nchelper.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_mkarbinit.o : mod_clm_mkarbinit.F90 mod_bats_param.o mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_snicar.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mppparam.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_nchelper.o : mod_clm_nchelper.F90 mod_clm_decomp.o mod_clm_varcon.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_message.o mod_mppparam.o mod_realkinds.o mod_regcm_types.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_organicfile.o : mod_clm_organicfile.F90 mod_clm_decomp.o mod_clm_domain.o mod_clm_nchelper.o mod_clm_type.o mod_clm_varctl.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_pft2col.o : mod_clm_pft2col.F90 mod_clm_subgridave.o mod_clm_type.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_pftdyn.o : mod_clm_pftdyn.F90 mod_clm_decomp.o mod_clm_nchelper.o mod_clm_pftvarcon.o mod_clm_time_manager.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_clm_varsur.o mod_constants.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_pftvarcon.o : mod_clm_pftvarcon.F90 mod_clm_nchelper.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_qsat.o : mod_clm_qsat.F90 mod_constants.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_regcm.o : mod_clm_regcm.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_driver.o mod_clm_initialize.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_constants.o mod_date.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_sunorbit.o $(REGCM_ROOT)/makeinc
mod_clm_restfile.o : mod_clm_restfile.F90 mod_clm_accumul.o mod_clm_biogeophysrest.o mod_clm_ch4rest.o mod_clm_cnrest.o mod_clm_croprest.o mod_clm_decomp.o mod_clm_histfile.o mod_clm_nchelper.o mod_clm_slakerest.o mod_clm_subgridrest.o mod_clm_surfrd.o mod_clm_time_manager.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_reweight.o : mod_clm_reweight.F90 mod_clm_decomp.o mod_clm_domain.o mod_clm_filter.o mod_clm_type.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_slakecon.o : mod_clm_slakecon.F90 mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_slakefluxes.o : mod_clm_slakefluxes.F90 mod_clm_atmlnd.o mod_clm_frictionvelocity.o mod_clm_qsat.o mod_clm_slakecon.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_slakehydrology.o : mod_clm_slakehydrology.F90 mod_clm_atmlnd.o mod_clm_slakecon.o mod_clm_snicar.o mod_clm_snowhydrology.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_slakerest.o : mod_clm_slakerest.F90 mod_clm_decomp.o mod_clm_nchelper.o mod_clm_type.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_slaketemperature.o : mod_clm_slaketemperature.F90 mod_clm_atmlnd.o mod_clm_qsat.o mod_clm_slakecon.o mod_clm_tridiagonal.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_snicar.o : mod_clm_snicar.F90 mod_clm_nchelper.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_snowhydrology.o : mod_clm_snowhydrology.F90 mod_clm_atmlnd.o mod_clm_slakecon.o mod_clm_snicar.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_soilhydrology.o : mod_clm_soilhydrology.F90 mod_clm_h2osfc.o mod_clm_tridiagonal.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_clm_vicmap.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_soillittverttransp.o : mod_clm_soillittverttransp.F90 mod_clm_tridiagonal.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_soiltemperature.o : mod_clm_soiltemperature.F90 mod_clm_atmlnd.o mod_clm_banddiagonal.o mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_staticecosysdyn.o : mod_clm_staticecosysdyn.F90 mod_clm_decomp.o mod_clm_domain.o mod_clm_nchelper.o mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varctl.o mod_clm_varpar.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_subgridave.o : mod_clm_subgridave.F90 mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_subgrid.o : mod_clm_subgrid.F90 mod_clm_varctl.o mod_clm_varpar.o mod_clm_varsur.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_subgridrest.o : mod_clm_subgridrest.F90 mod_clm_decomp.o mod_clm_domain.o mod_clm_nchelper.o mod_clm_type.o mod_date.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_surfacealbedo.o : mod_clm_surfacealbedo.F90 mod_clm_slakecon.o mod_clm_snicar.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o mod_sunorbit.o $(REGCM_ROOT)/makeinc
mod_clm_surfaceradiation.o : mod_clm_surfaceradiation.F90 mod_clm_atmlnd.o mod_clm_snicar.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_date.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_surfrd.o : mod_clm_surfrd.F90 mod_clm_atmlnd.o mod_clm_decomp.o mod_clm_domain.o mod_clm_nchelper.o mod_clm_pftvarcon.o mod_clm_type.o mod_clm_urbaninput.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_clm_varsur.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_time_manager.o : mod_clm_time_manager.F90 mod_clm_varctl.o mod_constants.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o $(REGCM_ROOT)/makeinc
mod_clm_tridiagonal.o : mod_clm_tridiagonal.F90 mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_type.o : mod_clm_type.F90 mod_clm_domain.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_typeinit.o : mod_clm_typeinit.F90 mod_clm_decomp.o mod_clm_drydep.o mod_clm_megan.o mod_clm_surfrd.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_urban.o : mod_clm_urban.F90 mod_clm_atmlnd.o mod_clm_filter.o mod_clm_frictionvelocity.o mod_clm_qsat.o mod_clm_type.o mod_clm_urbaninput.o mod_clm_varcon.o mod_clm_varpar.o mod_date.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_urbaninit.o : mod_clm_urbaninit.F90 mod_clm_decomp.o mod_clm_type.o mod_clm_urbaninput.o mod_clm_urban.o mod_clm_varcon.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_urbaninput.o : mod_clm_urbaninput.F90 mod_clm_decomp.o mod_clm_domain.o mod_clm_nchelper.o mod_clm_type.o mod_clm_varctl.o mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_varcon.o : mod_clm_varcon.F90 mod_clm_varpar.o mod_constants.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_varctl.o : mod_clm_varctl.F90 mod_clm_varpar.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_runparams.o mod_stdio.o $(REGCM_ROOT)/makeinc
mod_clm_varpar.o : mod_clm_varpar.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o $(REGCM_ROOT)/makeinc
mod_clm_varsur.o : mod_clm_varsur.F90 mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_vicmap.o : mod_clm_vicmap.F90 mod_clm_type.o mod_clm_varcon.o mod_clm_varpar.o mod_intkinds.o mod_realkinds.o $(REGCM_ROOT)/makeinc
mod_clm_vocemission.o : mod_clm_vocemission.F90 mod_clm_atmlnd.o mod_clm_domain.o mod_clm_meganfactors.o mod_clm_megan.o mod_clm_pftvarcon.o mod_clm_type.o mod_clm_varcon.o mod_clm_varctl.o mod_clm_varpar.o mod_intkinds.o mod_mpmessage.o mod_realkinds.o mod_stdio.o $(REGCM_ROOT)/makeinc
