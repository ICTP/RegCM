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

vpath %.mod $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:.

noinst_LIBRARIES = libcum.a
libcum_a_SOURCES = mod_cu_common.f90 mod_cu_bm.F90 mod_cu_em.f90 \
                   mod_cu_grell.F90 mod_cu_kuo.F90 mod_cu_tables.f90 \
                   mod_cu_tiedtke.F90

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) 
else
FCFLAGS += -I$(MPP_INCLUDE)
endif

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

%.mod: %.f90 %.o
	@true

%.mod: %.F90 %.o
	@true

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_cu_bm.o mod_cu_bm.mod: mod_cu_bm.F90 mod_dynparam.mod mod_realkinds.mod mod_memutil.mod mod_cu_common.mod mod_runparams.mod
mod_cu_common.o mod_cu_common.mod: mod_cu_common.f90 mod_dynparam.mod mod_memutil.mod mod_realkinds.mod mod_runparams.mod
mod_cu_em.o mod_cu_em.mod: mod_cu_em.f90 mod_dynparam.mod mod_memutil.mod mod_cu_common.mod mod_runparams.mod
mod_cu_grell.o mod_cu_grell.mod: mod_cu_grell.F90 mod_dynparam.mod mod_memutil.mod mod_service.mod mod_cu_common.mod mod_runparams.mod
mod_cu_kuo.o mod_cu_kuo.mod: mod_cu_kuo.F90 mod_dynparam.mod mod_memutil.mod mod_cu_common.mod mod_runparams.mod mod_mppparam.mod
mod_cu_tables.o mod_cu_tables.mod: mod_cu_tables.f90 mod_realkinds.mod mod_constants.mod
mod_cu_tiedtke.o mod_cu_tiedtke.mod: mod_cu_tiedtke.F90 mod_realkinds.mod mod_dynparam.mod mod_mpmessage.mod mod_memutil.mod mod_constants.mod mod_cu_common.mod mod_cu_tables.mod mod_service.mod mod_runparams.mod
