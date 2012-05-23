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

vpath %.mod $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:GAS_CBMZ/:.

SUBDIRS = GAS_CBMZ

CBMZLIB = GAS_CBMZ/libcbmz.a

noinst_LIBRARIES = libchem.a
libchem_a_SOURCES = mod_che_param.f90 mod_che_common.F90 mod_che_mppio.F90 \
                    mod_che_tend.F90 mod_che_ncio.F90 mod_che_dust.F90 \
                    mod_che_drydep.F90 mod_che_species.f90 mod_che_indices.f90 \
                    mod_che_carbonaer.f90 mod_che_seasalt.f90 \
                    mod_che_start.F90 mod_che_cumtran.F90 mod_che_sox.f90 \
                    mod_che_emission.F90 mod_che_wetdep.F90 \
                    mod_che_chemistry.f90 mod_che_bdyco.F90 mod_che_output.F90

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) 
else
FCFLAGS += -I$(MPP_INCLUDE)
endif

if COMPILER_SUN
  CPPFLAGS += -MGAS_CBMZ
else
  CPPFLAGS += -IGAS_CBMZ
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
mod_che_bdyco.o mod_che_bdyco.mod: mod_che_bdyco.F90 mod_realkinds.mod mod_dynparam.mod mod_memutil.mod mod_mppparam.mod mod_service.mod mod_mpmessage.mod mod_che_common.mod mod_che_mppio.mod mod_che_ncio.mod mod_che_species.mod mod_che_indices.mod mod_che_emission.mod mod_mppparam.mod
mod_che_carbonaer.o mod_che_carbonaer.mod: mod_che_carbonaer.f90 mod_realkinds.mod mod_dynparam.mod mod_constants.mod mod_che_common.mod mod_che_species.mod mod_che_indices.mod
mod_che_chemistry.o mod_che_chemistry.mod: mod_che_chemistry.f90 mod_realkinds.mod mod_dynparam.mod mod_constants.mod mod_che_common.mod mod_che_indices.mod mod_che_species.mod mod_cbmz_boxvars.mod mod_cbmz_chemmech.mod mod_cbmz_chemvars.mod mod_cbmz_molwg.mod mod_cbmz_main1.mod
mod_che_common.o mod_che_common.mod: mod_che_common.F90 mod_realkinds.mod mod_dynparam.mod mod_mppparam.mod mod_memutil.mod mod_mpmessage.mod mod_che_param.mod mod_che_species.mod
mod_che_cumtran.o mod_che_cumtran.mod: mod_che_cumtran.F90 mod_dynparam.mod mod_realkinds.mod mod_che_common.mod
mod_che_drydep.o mod_che_drydep.mod: mod_che_drydep.F90 mod_realkinds.mod mod_constants.mod mod_dynparam.mod mod_che_common.mod mod_che_dust.mod mod_mpmessage.mod mod_service.mod mod_che_ncio.mod mod_che_mppio.mod mod_che_indices.mod mod_che_indices.mod
mod_che_dust.o mod_che_dust.mod: mod_che_dust.F90 mod_realkinds.mod mod_constants.mod mod_dynparam.mod mod_mpmessage.mod mod_memutil.mod mod_mppparam.mod mod_che_common.mod mod_che_ncio.mod mod_che_mppio.mod
	$(MPIFC) $(CPPFLAGS) $(FCFLAGS) -c $<
mod_che_emission.o mod_che_emission.mod: mod_che_emission.F90 mod_constants.mod mod_mpmessage.mod mod_service.mod mod_dynparam.mod mod_che_common.mod mod_che_param.mod mod_che_mppio.mod mod_che_dust.mod mod_che_ncio.mod mod_che_indices.mod
mod_che_indices.o mod_che_indices.mod: mod_che_indices.f90
mod_che_mppio.o mod_che_mppio.mod: mod_che_mppio.f90 mod_realkinds.mod mod_dynparam.mod mod_memutil.mod mod_mpmessage.mod mod_che_param.mod mod_che_common.mod mod_che_species.mod
mod_che_ncio.o mod_che_ncio.mod: mod_che_ncio.F90 mod_realkinds.mod mod_dynparam.mod mod_mpmessage.mod mod_che_indices.mod mod_che_common.mod
mod_che_output.o mod_che_output.mod: mod_che_output.F90 mod_constants.mod mod_mpmessage.mod mod_service.mod mod_dynparam.mod mod_che_common.mod mod_che_param.mod mod_che_mppio.mod mod_che_ncio.mod mod_che_indices.mod
mod_che_param.o mod_che_param.mod: mod_che_param.f90 mod_realkinds.mod
mod_che_seasalt.o mod_che_seasalt.mod: mod_che_seasalt.f90 mod_realkinds.mod mod_dynparam.mod mod_constants.mod mod_che_common.mod
mod_che_sox.o mod_che_sox.mod: mod_che_sox.f90 mod_constants.mod mod_dynparam.mod mod_che_common.mod mod_che_species.mod mod_che_indices.mod
mod_che_species.o mod_che_species.mod: mod_che_species.f90 mod_realkinds.mod mod_dynparam.mod mod_memutil.mod
mod_che_start.o mod_che_start.mod: mod_che_start.F90 mod_dynparam.mod mod_mpmessage.mod mod_che_common.mod mod_che_emission.mod mod_che_indices.mod mod_che_bdyco.mod mod_che_wetdep.mod mod_che_carbonaer.mod mod_che_ncio.mod mod_che_mppio.mod mod_che_bdyco.mod mod_cbmz_init1.mod mod_che_dust.mod mod_che_sox.mod mod_che_seasalt.mod mod_mppparam.mod
mod_che_tend.o mod_che_tend.mod: mod_che_tend.F90 mod_constants.mod mod_realkinds.mod mod_dynparam.mod mod_che_common.mod mod_che_indices.mod mod_che_param.mod mod_che_sox.mod mod_che_drydep.mod mod_che_wetdep.mod mod_che_emission.mod mod_che_dust.mod mod_che_seasalt.mod mod_che_carbonaer.mod mod_che_mppio.mod mod_che_chemistry.mod
mod_che_wetdep.o mod_che_wetdep.mod: mod_che_wetdep.F90 mod_realkinds.mod mod_constants.mod mod_dynparam.mod mod_che_common.mod mod_che_indices.mod mod_che_drydep.mod
