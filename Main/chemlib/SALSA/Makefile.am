#n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

noinst_LIBRARIES = libsalsa.a
libsalsa_a_SOURCES = aerosol_thermodynamics.f90 class_componentIndex.f90  \
        mo_constants.f90 mo_kind.f90 mo_salsa.f90 mo_salsa_cloud.f90 \
        mo_salsa_dynamics.f90 mo_salsa_init.f90 mo_salsa_driver.f90 \
        mo_salsa_nucleation.f90 mo_salsa_properties.f90 mo_salsa_sizedist.f90 \
        mo_salsa_update.f90 mo_submctl.f90


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
aerosol_thermodynamics.o: mo_kind.o
class_componentIndex.o: 
mo_constants.o: mo_kind.o
mo_kind.o: 
mo_salsa.o: class_componentIndex.o mo_kind.o mo_salsa_cloud.o \
        mo_salsa_dynamics.o mo_salsa_nucleation.o mo_salsa_properties.o \
        mo_salsa_update.o mo_submctl.o
mo_salsa_cloud.o: mo_constants.o mo_kind.o mo_submctl.o
mo_salsa_driver.o: class_componentIndex.o mo_kind.o mo_salsa.o \
        mo_salsa_properties.o mo_submctl.o
mo_salsa_dynamics.o: aerosol_thermodynamics.o class_componentIndex.o \
        mo_constants.o mo_kind.o mo_salsa_nucleation.o mo_salsa_properties.o \
        mo_submctl.o
mo_salsa_init.o: mo_kind.o mo_salsa_driver.o mo_submctl.o
mo_salsa_nucleation.o: mo_kind.o mo_submctl.o
mo_salsa_properties.o: mo_kind.o mo_submctl.o
mo_salsa_sizedist.o: mo_kind.o mo_salsa_driver.o mo_submctl.o
mo_salsa_update.o: mo_kind.o mo_submctl.o
mo_submctl.o: mo_kind.o

