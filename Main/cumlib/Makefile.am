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

noinst_LIBRARIES = libcum.a
libcum_a_SOURCES = mod_cu_common.F90 mod_cu_bm.F90 mod_cu_em.F90 \
                   mod_cu_grell.F90 mod_cu_kuo.F90 mod_cu_tables.F90 \
                   mod_cu_tiedtke.F90 mod_cu_kf.F90 mod_cu_shallow.F90

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
mod_cu_bm.o : mod_cu_bm.F90 mod_cu_common.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_constants.o
mod_cu_common.o : mod_cu_common.F90 mod_constants.o mod_dynparam.o mod_intkinds.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_mppparam.o mod_memutil.o
mod_cu_em.o : mod_cu_em.F90 mod_cu_common.o mod_intkinds.o mod_memutil.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_dynparam.o mod_constants.o mod_mppparam.o pfesat.inc pfwsat.inc wlh.inc
mod_cu_grell.o : mod_cu_grell.F90 mod_cu_common.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_realkinds.o mod_dynparam.o mod_regcm_types.o mod_runparams.o mod_service.o mod_constants.o pfwsat.inc
mod_cu_kf.o : mod_cu_kf.F90 mod_constants.o mod_cu_common.o mod_intkinds.o mod_memutil.o mod_dynparam.o mod_mpmessage.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_stdio.o
mod_cu_kuo.o : mod_cu_kuo.F90 mod_cu_common.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mppparam.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_constants.o
mod_cu_shallow.o : mod_cu_shallow.F90 mod_realkinds.o mod_intkinds.o mod_constants.o mod_dynparam.o mod_cu_common.o mod_dynparam.o mod_runparams.o mod_regcm_types.o
mod_cu_tables.o : mod_cu_tables.F90 mod_constants.o mod_intkinds.o mod_realkinds.o
mod_cu_tiedtke.o : mod_cu_tiedtke.F90 mod_constants.o mod_cu_common.o mod_cu_tables.o mod_dynparam.o mod_intkinds.o mod_memutil.o mod_mpmessage.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_mppparam.o mod_stdio.o
