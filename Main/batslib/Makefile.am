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

noinst_LIBRARIES = libbats.a
libbats_a_SOURCES = mod_bats_common.F90 mod_bats_param.f90 mod_bats_bndry.F90 \
                    mod_bats_co2.f90 mod_bats_drag.F90 mod_bats_lake.F90 \
                    mod_bats_leaftemp.f90 mod_bats_vecbats.F90 \
                    mod_bats_mppio.F90 mod_bats_zengocn.F90 \
                    mod_bats_romsocn.F90 mod_bats_internal.f90

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

mod_bats_mppio.o mod_bats_mppio.mod : mod_bats_mppio.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod
mod_bats_param.o mod_bats_param.mod : mod_bats_param.f90
mod_bats_common.o mod_bats_common.mod : mod_bats_common.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod \
  mod_bats_param.o mod_bats_internal.o
mod_bats_bndry.o mod_bats_bndry.mod : mod_bats_bndry.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_service.mod \
  mod_bats_common.o mod_bats_leaftemp.o mod_bats_drag.o mod_bats_internal.o \
  mod_bats_internal.o
mod_bats_co2.o mod_bats_co2.mod : mod_bats_co2.f90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod mod_bats_common.o mod_bats_leaftemp.o
mod_bats_drag.o mod_bats_drag.mod : mod_bats_drag.F90 mod_bats_internal.o \
  $(REGCM_ROOT)/Share/mod_dynparam.mod mod_bats_common.o
mod_bats_lake.o mod_bats_lake.mod : mod_bats_lake.F90 mod_bats_internal.o \
  $(REGCM_ROOT)/Share/mod_dynparam.mod mod_bats_common.o mod_bats_mppio.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_bats_leaftemp.o mod_bats_leaftemp.mod : mod_bats_leaftemp.f90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod \
  mod_bats_common.o mod_bats_internal.o
mod_bats_vecbats.o mod_bats_vecbats.mod : mod_bats_vecbats.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_message.mod \
  $(REGCM_ROOT)/Share/mod_service.mod $(REGCM_ROOT)/Share/mod_constants.o \
  mod_bats_common.o mod_bats_lake.o mod_bats_bndry.o mod_bats_drag.o \
  mod_bats_mppio.o mod_bats_zengocn.o mod_bats_romsocn.o
mod_bats_zengocn.o mod_bats_zengocn.mod : mod_bats_zengocn.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_service.mod \
  mod_bats_common.o
mod_bats_romsocn.o mod_bats_romsocn.mod : mod_bats_romsocn.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_service.mod \
  mod_bats_common.o
mod_bats_internal.o mod_bats_internal.mod: mod_bats_internal.f90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_memutil.mod
