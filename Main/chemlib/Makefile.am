#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#    
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:GAS_CBMZ_NEW/:.
vpath %.inc $(REGCM_ROOT)/Share:.

SUBDIRS = GAS_CBMZ_NEW

CBMZLIB = GAS_CBMZ_NEW/libcbmz.a

noinst_LIBRARIES = libchem.a
libchem_a_SOURCES = mod_che_param.F90 mod_che_common.F90 mod_che_mppio.F90 \
                    mod_che_tend.F90 mod_che_ncio.F90 mod_che_dust.F90 \
                    mod_che_drydep.F90 mod_che_species.F90 mod_che_indices.F90 \
                    mod_che_carbonaer.F90 mod_che_seasalt.F90 \
                    mod_che_start.F90 mod_che_cumtran.F90 mod_che_sox.F90 \
                    mod_che_emission.F90 mod_che_wetdep.F90 mod_che_molwg.F90 \
                    mod_che_chemistry.F90 mod_che_isorropia.F90 \
                    mod_che_common_isorropia.F90 mod_che_bdyco.F90 \
                    mod_che_output.F90 mod_che_pollen.F90 mod_che_bionit.F90 \
                    mod_che_ccn.F90 mod_che_linox.F90

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) 
else
FCFLAGS += -I$(MPP_INCLUDE)
endif

if COMPILER_SUN
  CPPFLAGS += -MGAS_CBMZ_NEW
else
  CPPFLAGS += -IGAS_CBMZ_NEW
endif

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_che_bdyco.o : mod_che_bdyco.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_nhinterp.o mod_zita.o mod_mppparam.o mod_runparams.o mod_service.o mod_mpmessage.o mod_che_common.o mod_che_mppio.o mod_che_ncio.o mod_che_species.o mod_che_indices.o mod_che_emission.o mod_mppparam.o mod_constants.o mod_date.o
mod_che_bionit.o : mod_che_bionit.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_che_common.o mod_che_dust.o
mod_che_carbonaer.o : mod_che_carbonaer.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_che_common.o mod_che_species.o mod_che_indices.o mod_runparams.o
mod_che_ccn.o : mod_che_ccn.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_che_common.o mod_che_species.o mod_che_indices.o mod_che_carbonaer.o mod_constants.o
mod_che_chemistry.o : mod_che_chemistry.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_runparams.o mod_che_common.o mod_che_indices.o mod_che_species.o  mod_cbmz_global.o mod_cbmz_parameters.o mod_cbmz_main.o mod_che_molwg.o
mod_che_common.o : mod_che_common.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_mppparam.o mod_runparams.o mod_memutil.o mod_mpmessage.o mod_che_param.o mod_che_species.o mod_che_indices.o mod_cbmz_global.o mod_stdio.o
mod_che_cumtran.o : mod_che_cumtran.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_che_common.o mod_constants.o
mod_che_drydep.o : mod_che_drydep.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_che_common.o mod_che_dust.o mod_mpmessage.o mod_service.o mod_che_ncio.o mod_che_mppio.o mod_che_indices.o mod_che_indices.o
mod_che_dust.o : mod_che_dust.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_mpmessage.o mod_memutil.o mod_mppparam.o mod_che_common.o mod_che_ncio.o mod_che_mppio.o
mod_che_emission.o : mod_che_emission.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_mpmessage.o mod_service.o mod_dynparam.o mod_che_common.o mod_che_param.o mod_che_mppio.o mod_che_dust.o mod_che_ncio.o mod_che_indices.o mod_runparams.o
mod_che_indices.o : mod_che_indices.F90 mod_intkinds.o mod_cbmz_parameters.o
mod_che_mppio.o : mod_che_mppio.F90 mod_intkinds.o mod_realkinds.o mod_mppparam.o mod_runparams.o mod_dynparam.o mod_memutil.o mod_mpmessage.o mod_che_param.o mod_che_common.o mod_che_species.o
mod_che_ncio.o : mod_che_ncio.F90 mod_intkinds.o mod_realkinds.o mod_nchelper.o mod_dynparam.o mod_mpmessage.o mod_che_indices.o mod_che_common.o mod_runparams.o mod_domain.o mod_constants.o mod_date.o
mod_che_output.o : mod_che_output.F90 mod_intkinds.o mod_realkinds.o mod_runparams.o mod_dynparam.o mod_che_param.o mod_che_common.o mod_outvars.o mod_che_indices.o mod_che_carbonaer.o mod_constants.o
mod_che_param.o : mod_che_param.F90 mod_intkinds.o mod_realkinds.o
mod_che_pollen.o : mod_che_pollen.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_che_common.o mod_che_indices.o
mod_che_seasalt.o : mod_che_seasalt.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_che_common.o
mod_che_sox.o : mod_che_sox.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_che_common.o mod_che_species.o mod_che_indices.o mod_che_bdyco.o mod_che_carbonaer.o
mod_che_species.o : mod_che_species.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o
mod_che_start.o : mod_che_start.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_mpmessage.o mod_che_common.o mod_che_indices.o mod_che_bdyco.o mod_che_wetdep.o mod_che_carbonaer.o mod_che_ncio.o mod_che_mppio.o mod_che_bdyco.o  mod_che_dust.o mod_che_sox.o mod_che_seasalt.o mod_mppparam.o mod_cbmz_hvread.o mod_che_bionit.o mod_che_chemistry.o
mod_che_tend.o : mod_che_tend.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o mod_mppparam.o mod_che_common.o mod_che_indices.o mod_che_param.o mod_che_sox.o mod_che_drydep.o mod_che_wetdep.o mod_che_emission.o mod_che_dust.o mod_che_seasalt.o mod_che_carbonaer.o mod_che_mppio.o mod_che_chemistry.o mod_che_isorropia.o mod_che_pollen.o mod_che_bionit.o mod_che_linox.o mod_che_ccn.o
mod_che_wetdep.o : mod_che_wetdep.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_che_common.o mod_che_indices.o mod_che_drydep.o
mod_che_common_isorropia.o : mod_che_common_isorropia.F90 mod_realkinds.o mod_intkinds.o
mod_che_isorropia.o : mod_che_isorropia.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_constants.o mod_che_common.o mod_che_indices.o mod_che_species.o mod_che_common_isorropia.o mod_che_molwg.o
mod_che_molwg.o: mod_che_molwg.F90 mod_che_indices.o mod_intkinds.o mod_realkinds.o
mod_che_linox.o: mod_che_linox.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_che_common.o mod_dynparam.o

