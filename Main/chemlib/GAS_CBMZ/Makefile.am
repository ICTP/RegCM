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

noinst_LIBRARIES = libcbmz.a
libcbmz_a_SOURCES = mod_cbmz_boxvars.F90 mod_cbmz_chemmech.F90 \
  mod_cbmz_chemvars.F90 mod_cbmz_chemlocal.F90 mod_cbmz_molwg.F90 \
  mod_cbmz_linslv.F90 mod_cbmz_jval1.F90 mod_cbmz_solve1.F90 \
  mod_cbmz_rates1.F90 mod_cbmz_main1.F90 mod_cbmz_init1.F90

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) 
else
FCFLAGS += -I$(MPP_INCLUDE)
endif

%.o: %.F90
	$(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_cbmz_boxvars.o : mod_cbmz_boxvars.F90 mod_intkinds.o mod_realkinds.o
mod_cbmz_chemlocal.o : mod_cbmz_chemlocal.F90 mod_intkinds.o mod_realkinds.o mod_cbmz_chemmech.o
mod_cbmz_chemmech.o : mod_cbmz_chemmech.F90 mod_intkinds.o mod_realkinds.o
mod_cbmz_chemvars.o : mod_cbmz_chemvars.F90 mod_intkinds.o mod_realkinds.o mod_cbmz_chemmech.o
mod_cbmz_init1.o : mod_cbmz_init1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_mpmessage.o mod_cbmz_chemmech.o mod_cbmz_chemlocal.o mod_cbmz_chemvars.o mod_cbmz_jval1.o
mod_cbmz_jval1.o : mod_cbmz_jval1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_mppparam.o mod_stdio.o
mod_cbmz_linslv.o : mod_cbmz_linslv.F90 mod_intkinds.o mod_realkinds.o
mod_cbmz_main1.o : mod_cbmz_main1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_cbmz_boxvars.o mod_cbmz_chemmech.o mod_cbmz_chemvars.o mod_cbmz_rates1.o mod_cbmz_solve1.o
mod_cbmz_molwg.o : mod_cbmz_molwg.F90 mod_intkinds.o mod_realkinds.o
mod_cbmz_rates1.o : mod_cbmz_rates1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_cbmz_chemvars.o mod_cbmz_chemmech.o mod_cbmz_chemlocal.o mod_cbmz_jval1.o
mod_cbmz_solve1.o : mod_cbmz_solve1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_cbmz_chemmech.o mod_cbmz_chemvars.o mod_cbmz_chemlocal.o mod_cbmz_linslv.o mod_cbmz_init1.o
