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

noinst_LIBRARIES = librad.a
librad_a_SOURCES = mod_rad_colmod3.F90 mod_rad_o3blk.F90 mod_rad_outrad.F90 \
                   mod_rad_common.f90 mod_rad_radiation.F90 mod_rad_tracer.f90 \
                   mod_rad_aerosol.F90 mod_rad_scenarios.F90

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

mod_rad_aerosol.o mod_rad_aerosol.mod: mod_rad_aerosol.F90 \
  $(REGCM_ROOT)/external/m_realkinds.mod $(REGCM_ROOT)/Share/mod_message.mod \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
	$(REGCM_ROOT)/Share/mod_memutil.mod mod_rad_common.o
mod_rad_common.o mod_rad_common.mod: mod_rad_common.f90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
	$(REGCM_ROOT)/Share/mod_memutil.mod
mod_rad_o3blk.o mod_rad_o3blk.mod: mod_rad_o3blk.F90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
	$(REGCM_ROOT)/Share/mod_memutil.mod mod_rad_common.o
mod_rad_colmod3.o mod_rad_colmod3.mod: mod_rad_colmod3.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod mod_rad_common.o mod_rad_radiation.o \
  mod_rad_outrad.o
mod_rad_outrad.o mod_rad_outrad.mod: mod_rad_outrad.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod mod_rad_common.o
mod_rad_scenarios.o mod_rad_scenarios.mod: mod_rad_scenarios.F90 \
  $(REGCM_ROOT)/Share/mod_dynparam.mod $(REGCM_ROOT)/Share/mod_message.mod
mod_rad_tracer.o mod_rad_tracer.mod: mod_rad_tracer.f90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod
mod_rad_radiation.o mod_rad_radiation.mod: mod_rad_radiation.F90 \
  $(REGCM_ROOT)/Share/mod_constants.mod $(REGCM_ROOT)/Share/mod_dynparam.mod \
  $(REGCM_ROOT)/Share/mod_message.mod $(REGCM_ROOT)/Share/mod_memutil.mod \
  mod_rad_tracer.o mod_rad_scenarios.o mod_rad_aerosol.o
