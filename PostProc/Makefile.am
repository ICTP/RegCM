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

vpath %.mod $(REGCM_ROOT)/Share:.

include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/PreProc/Terrain:.
vpath %.f90 $(REGCM_ROOT)/PreProc/Terrain:.

bin_PROGRAMS = sigma2p GrADSNcPlot GrADSNcPrepare average regrid
sigma2p_SOURCES = sigma2p.F90 mod_derived.f90
GrADSNcPlot_SOURCES = GrADSNcPlot.F90
GrADSNcPrepare_SOURCES = GrADSNcPrepare.F90

average: average.sh
	cp $< $@

regrid: regrid.sh
	cp $< $@

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

%.mod: %.f90 %.o
	@true

%.mod: %.F90 %.o
	@true

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
GrADSNcPlot.o GrADSNcPlot.mod: GrADSNcPlot.F90 mod_projections.mod mod_date.mod mod_message.mod mod_nchelper.mod
GrADSNcPrepare.o GrADSNcPrepare.mod: GrADSNcPrepare.F90 mod_projections.mod mod_date.mod mod_message.mod mod_nchelper.mod
mod_derived.o mod_derived.mod: mod_derived.f90 mod_constants.mod mod_realkinds.mod
sigma2p.o sigma2p.mod: sigma2p.F90 mod_realkinds.mod mod_constants.mod mod_message.mod mod_vertint.mod mod_derived.mod mod_nchelper.mod
