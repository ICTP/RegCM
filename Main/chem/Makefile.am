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

SUBDIRS = GAS_KPP_GEOS

noinst_LIBRARIES = libchem.a
libchem_a_SOURCES = mod_chem.f90 mod_mainchem.f90 mod_trachem.f90 \
		    mod_che_mppio.F90 mod_che_ncio.F90 mod_dust.F90 \
		    mod_aerosol.F90 mod_che_semdde.F90

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

mod_chem.o: mod_chem.f90 $(REGCM_ROOT)/external/m_realkinds.mod
mod_mainchem.o: mod_mainchem.f90 $(REGCM_ROOT)/external/m_realkinds.mod \
                mod_chem.o $(REGCM_ROOT)/Share/mod_dynparam.mod \
                $(REGCM_ROOT)/Share/mod_memutil.mod
mod_trachem.o: mod_trachem.f90 $(REGCM_ROOT)/external/m_realkinds.mod \
               mod_chem.o $(REGCM_ROOT)/Share/mod_dynparam.mod \
               $(REGCM_ROOT)/Share/mod_message.mod \
               $(REGCM_ROOT)/Share/mod_memutil.mod
mod_che_mppio.o: mod_che_mppio.f90 $(REGCM_ROOT)/external/m_realkinds.mod \
                 mod_chem.o $(REGCM_ROOT)/Share/mod_dynparam.mod \
                 $(REGCM_ROOT)/Share/mod_message.mod \
                 $(REGCM_ROOT)/Share/mod_memutil.mod
mod_che_ncio.o: mod_che_ncio.F90 $(REGCM_ROOT)/external/m_realkinds.mod \
                $(REGCM_ROOT)/Share/mod_dynparam.mod \
                $(REGCM_ROOT)/Share/mod_message.mod
mod_dust.o: mod_dust.F90 $(REGCM_ROOT)/external/m_realkinds.mod \
            $(REGCM_ROOT)/Share/mod_dynparam.mod \
            $(REGCM_ROOT)/Share/mod_constants.mod \
            $(REGCM_ROOT)/Share/mod_message.mod \
            $(REGCM_ROOT)/Share/mod_memutil.mod \
            mod_trachem.o mod_che_ncio.o mod_che_mppio.o
mod_aerosol.o: mod_aerosol.F90 $(REGCM_ROOT)/external/m_realkinds.mod \
               $(REGCM_ROOT)/Share/mod_dynparam.mod \
               $(REGCM_ROOT)/Share/mod_constants.mod \
               $(REGCM_ROOT)/Share/mod_memutil.mod \
               $(REGCM_ROOT)/Share/mod_message.mod \
               mod_chem.o mod_mainchem.o mod_trachem.o
mod_che_semdde.o: mod_che_semdde.F90 $(REGCM_ROOT)/external/m_realkinds.mod \
                  $(REGCM_ROOT)/Share/mod_dynparam.mod \
                  $(REGCM_ROOT)/Share/mod_message.mod \
                  mod_mainchem.o mod_trachem.o mod_dust.o mod_aerosol.o \
                  mod_che_ncio.o mod_che_mppio.o
