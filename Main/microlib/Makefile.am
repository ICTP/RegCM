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

noinst_LIBRARIES = libmicro.a
libmicro_a_SOURCES = mod_micro_nogtom.F90 mod_micro_subex.F90 \
                     mod_micro_wsm5.F90 mod_micro_wsm7.F90 \
                     mod_micro_wdm7.F90

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
mod_micro_nogtom.o : mod_micro_nogtom.F90 mod_realkinds.o mod_dynparam.o mod_mpmessage.o mod_memutil.o mod_regcm_types.o mod_constants.o mod_service.o mod_runparams.o mod_stdio.o esatliq.inc esatice.inc
mod_micro_subex.o : mod_micro_subex.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o mod_memutil.o mod_mpmessage.o mod_regcm_types.o pfesat.inc pfwsat.inc clwfromt.inc wlh.inc
mod_micro_wsm5.o : mod_micro_wsm5.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o mod_memutil.o mod_regcm_types.o
mod_micro_wsm7.o : mod_micro_wsm7.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o mod_memutil.o mod_regcm_types.o
mod_micro_wdm7.o : mod_micro_wdm7.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o mod_memutil.o mod_regcm_types.o
