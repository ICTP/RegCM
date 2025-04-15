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

noinst_LIBRARIES = libpbl.a
libpbl_a_SOURCES = mod_pbl_common.F90 mod_pbl_thetal.F90 mod_pbl_uwtcm.F90 \
                   mod_pbl_holtbl.F90 mod_pbl_gfs.F90 mod_pbl_myj.F90 \
                   mod_pbl_shinhong.F90

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
mod_pbl_common.o : mod_pbl_common.F90 mod_intkinds.o mod_realkinds.o mod_regcm_types.o mod_constants.o mod_memutil.o
mod_pbl_holtbl.o : mod_pbl_holtbl.F90 mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_pbl_common.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_constants.o
mod_pbl_thetal.o : mod_pbl_thetal.F90 mod_constants.o mod_intkinds.o mod_realkinds.o mod_runparams.o mod_mpmessage.o mod_stdio.o pfesat.inc wlh.inc
mod_pbl_uwtcm.o : mod_pbl_uwtcm.F90 mod_constants.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_mppparam.o mod_pbl_common.o mod_pbl_thetal.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o pfesat.inc wlh.inc
mod_pbl_gfs.o : mod_pbl_gfs.F90 mod_constants.o mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_runparams.o mod_regcm_types.o mod_memutil.o
mod_pbl_myj.o : mod_pbl_myj.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o mod_regcm_types.o mod_memutil.o
mod_pbl_shinhong.o : mod_pbl_shinhong.F90 mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o mod_regcm_types.o
