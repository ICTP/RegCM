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

noinst_LIBRARIES = libpbl.a
libpbl_a_SOURCES = mod_pbl_common.f90 mod_pbl_thetal.f90 mod_pbl_uwtcm.F90 \
                   mod_pbl_holtbl.F90

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

mod_pbl_common.o mod_pbl_common.mod: mod_pbl_common.f90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
	$(REGCM_ROOT)/Share/mod_memutil.mod
mod_pbl_holtbl.o mod_pbl_holtbl.mod: mod_pbl_holtbl.F90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
	$(REGCM_ROOT)/Share/mod_memutil.mod mod_pbl_common.o
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_pbl_thetal.o mod_pbl_thetal.mod: mod_pbl_thetal.f90 \
  $(REGCM_ROOT)/Share/mod_realkinds.mod $(REGCM_ROOT)/Share/mod_constants.mod
mod_pbl_uwtcm.o mod_pbl_uwtcm.mod: mod_pbl_uwtcm.F90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
	$(REGCM_ROOT)/Share/mod_memutil.mod $(REGCM_ROOT)/Share/mod_message.mod \
  mod_pbl_thetal.o mod_pbl_common.o
