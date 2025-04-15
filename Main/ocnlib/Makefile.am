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

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:.
vpath %.inc $(REGCM_ROOT)/Share:.

noinst_LIBRARIES = libocn.a
libocn_a_SOURCES = mod_ocn_albedo.F90 mod_ocn_coare.F90 mod_ocn_internal.F90 \
                   mod_ocn_zeng.F90 mod_ocn_bats.F90 mod_ocn_common.F90  \
                   mod_ocn_lake.F90

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) 
else
FCFLAGS += -I$(MPP_INCLUDE)
endif

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_ocn_albedo.o : mod_ocn_albedo.F90 mod_constants.o mod_intkinds.o mod_ocn_internal.o mod_realkinds.o mod_runparams.o mod_service.o
mod_ocn_bats.o : mod_ocn_bats.F90 mod_constants.o mod_dynparam.o mod_intkinds.o mod_ocn_internal.o mod_realkinds.o mod_runparams.o mod_service.o pfesat.inc pfqsat.inc
mod_ocn_coare.o : mod_ocn_coare.F90 mod_constants.o mod_dynparam.o mod_intkinds.o mod_ocn_internal.o mod_realkinds.o mod_runparams.o mod_service.o cpmf.inc pfesat.inc pfqsat.inc pfdesatdt.inc pqderiv.inc wlh.inc
mod_ocn_common.o : mod_ocn_common.F90 mod_intkinds.o mod_mppparam.o mod_ocn_albedo.o mod_ocn_bats.o mod_ocn_coare.o mod_ocn_internal.o mod_ocn_lake.o mod_ocn_zeng.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_constants.o mod_mpmessage.o
mod_ocn_internal.o : mod_ocn_internal.F90 mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_regcm_types.o mod_runparams.o
mod_ocn_lake.o : mod_ocn_lake.F90 mod_dynparam.o mod_intkinds.o mod_ocn_internal.o mod_realkinds.o mod_runparams.o mod_service.o mod_constants.o pfesat.inc pfqsat.inc
mod_ocn_zeng.o : mod_ocn_zeng.F90 mod_dynparam.o mod_intkinds.o mod_ocn_internal.o mod_realkinds.o mod_runparams.o mod_service.o mod_constants.o pfesat.inc pfwsat.inc wlh.inc
