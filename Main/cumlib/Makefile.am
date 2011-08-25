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

noinst_LIBRARIES = libcum.a
libcum_a_SOURCES = mod_cu_common.f90 mod_cu_bm.F90 mod_cu_em.f90 \
                   mod_cu_grell.F90 mod_cu_kuo.F90 mod_cu_tables.f90 \
                   mod_cu_tiedtke.F90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_cu_common.o mod_cu_common.mod: mod_cu_common.f90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod
mod_cu_bm.o mod_cu_bm.mod: mod_cu_bm.F90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
  $(REGCM_ROOT)/Share/mod_memutil.mod mod_cu_common.o
mod_cu_em.o mod_cu_em.mod: mod_cu_em.f90 $(REGCM_ROOT)/Share/mod_dynparam.mod \
  $(REGCM_ROOT)/Share/mod_memutil.mod mod_cu_common.o
mod_cu_grell.o mod_cu_grell.mod: mod_cu_grell.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod \
  $(REGCM_ROOT)/Share/mod_service.mod mod_cu_common.o
mod_cu_kuo.o mod_cu_kuo.mod: mod_cu_kuo.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod \
  mod_cu_common.o
mod_cu_tables.o mod_cu_tables.mod: mod_cu_tables.f90 \
  $(REGCM_ROOT)/external/m_realkinds.o $(REGCM_ROOT)/Share/mod_constants.o
mod_cu_tiedtke.o mod_cu_tiedtke.mod: mod_cu_tiedtke.F90 \
  $(REGCM_ROOT)/external/m_realkinds.o $(REGCM_ROOT)/Share/mod_constants.o \
  $(REGCM_ROOT)/Share/mod_dynparam.o $(REGCM_ROOT)/Share/mod_memutil.o \
  $(REGCM_ROOT)/Share/mod_message.o mod_cu_common.o mod_cu_tables.o
