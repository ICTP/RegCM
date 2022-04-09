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

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:.
vpath %.inc $(REGCM_ROOT)/Share:.

noinst_LIBRARIES = libbats.a
libbats_a_SOURCES = mod_bats_common.F90 mod_bats_param.F90 mod_bats_bndry.F90 \
                    mod_bats_co2.F90 mod_bats_drag.F90 mod_bats_leaftemp.F90 \
                    mod_bats_albedo.F90 mod_bats_internal.F90

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
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_bats_albedo.o : mod_bats_albedo.F90 mod_bats_drag.o mod_bats_internal.o mod_bats_param.o mod_constants.o mod_intkinds.o mod_realkinds.o mod_runparams.o mod_service.o mod_bats_leaftemp.o
mod_bats_bndry.o : mod_bats_bndry.F90 mod_bats_drag.o mod_bats_internal.o mod_bats_leaftemp.o mod_bats_param.o mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_service.o pfesat.inc pfqsat.inc pfdesatdt.inc pqderiv.inc wlh.inc
mod_bats_co2.o : mod_bats_co2.F90 mod_bats_internal.o mod_bats_leaftemp.o mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_bats_param.o
mod_bats_common.o : mod_bats_common.F90 mod_bats_albedo.o mod_bats_bndry.o mod_bats_drag.o mod_bats_internal.o mod_bats_leaftemp.o mod_bats_param.o mod_constants.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_stdio.o pfesat.inc pfwsat.inc
mod_bats_drag.o : mod_bats_drag.F90 mod_bats_internal.o mod_bats_param.o mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_service.o
mod_bats_internal.o : mod_bats_internal.F90 mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_regcm_types.o mod_constants.o
mod_bats_leaftemp.o : mod_bats_leaftemp.F90 mod_bats_internal.o mod_bats_param.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_service.o pfesat.inc pfqsat.inc pfdesatdt.inc pqderiv.inc
mod_bats_param.o : mod_bats_param.F90 mod_intkinds.o mod_realkinds.o
