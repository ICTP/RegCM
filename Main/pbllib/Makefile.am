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

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:.

noinst_LIBRARIES = libpbl.a
libpbl_a_SOURCES = mod_pbl_common.f90 mod_pbl_thetal.f90 mod_pbl_uwtcm.F90 \
                   mod_pbl_holtbl.F90

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

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_pbl_common.o : mod_pbl_common.f90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_memutil.o mod_runparams.o
mod_pbl_holtbl.o : mod_pbl_holtbl.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_mppparam.o mod_memutil.o mod_service.o mod_mpmessage.o mod_pbl_common.o mod_runparams.o
mod_pbl_thetal.o : mod_pbl_thetal.f90 mod_intkinds.o mod_realkinds.o mod_constants.o
mod_pbl_uwtcm.o : mod_pbl_uwtcm.F90 mod_intkinds.o mod_realkinds.o mod_memutil.o mod_dynparam.o mod_constants.o mod_mpmessage.o mod_pbl_common.o mod_pbl_thetal.o mod_runparams.o mod_mppparam.o mod_runparams.o
