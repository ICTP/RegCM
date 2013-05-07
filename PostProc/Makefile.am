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

vpath %.o $(REGCM_ROOT)/Share:.

include $(REGCM_ROOT)/makeinc

bin_PROGRAMS = sigma2p sigma2z GrADSNcPlot GrADSNcPrepare average regrid
sigma2p_SOURCES = sigma2p.F90
sigma2z_SOURCES = sigma2z.F90
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

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
GrADSNcPlot.o : GrADSNcPlot.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_projections.o mod_stdio.o mod_message.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o
GrADSNcPrepare.o : GrADSNcPrepare.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_projections.o mod_stdio.o mod_message.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o
sigma2p.o : sigma2p.F90 mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o mod_vertint.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o mod_hgt.o mod_humid.o
sigma2z.o : sigma2z.F90 mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o mod_vertint.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o mod_humid.o
