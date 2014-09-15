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

if DO_COMPILE_LIBCLM45
bin_PROGRAMS = sigma2p sigma2z GrADSNcPlot GrADSNcPrepare clm45_1dto2d
else
bin_PROGRAMS = sigma2p sigma2z GrADSNcPlot GrADSNcPrepare
endif

sigma2p_SOURCES = sigma2p.F90
sigma2z_SOURCES = sigma2z.F90
GrADSNcPlot_SOURCES = GrADSNcPlot.F90
GrADSNcPrepare_SOURCES = GrADSNcPrepare.F90
clm45_1dto2d_SOURCES = clm45_1dto2d.F90

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

install-exec-hook:
	cp regrid.sh $(DESTDIR)$(bindir)/regrid
	cp average.sh $(DESTDIR)$(bindir)/average

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
GrADSNcPlot.o : GrADSNcPlot.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_projections.o mod_stdio.o mod_message.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o
GrADSNcPrepare.o : GrADSNcPrepare.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_projections.o mod_stdio.o mod_message.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o
sigma2p.o : sigma2p.F90 mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o mod_vertint.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o mod_hgt.o mod_humid.o mod_dynparam.o
sigma2z.o : sigma2z.F90 mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o mod_vertint.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o mod_humid.o mod_dynparam.o
