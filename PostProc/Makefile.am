#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP RegCM.
#    
#    Use of this source code is governed by an MIT-style license that can
#    be found in the LICENSE file or at
#
#         https://opensource.org/licenses/MIT.
#
#    ICTP RegCM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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

LDADD = $(OMPFLAGS)

%.o: %.F90
	$(FC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

install-exec-hook:
	cp regrid.sh $(DESTDIR)$(bindir)/regrid
	cp average.sh $(DESTDIR)$(bindir)/average

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
GrADSNcPlot.o : GrADSNcPlot.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_projections.o mod_stdio.o mod_message.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o mod_projections.o
GrADSNcPrepare.o : GrADSNcPrepare.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_projections.o mod_stdio.o mod_message.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o mod_projections.o
sigma2p.o : sigma2p.F90 mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o mod_vertint.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o mod_hgt.o mod_humid.o mod_dynparam.o mod_memutil.o mod_sigma.o mod_zita.o mod_stdio.o
sigma2z.o : sigma2z.F90 mod_realkinds.o mod_constants.o mod_stdio.o mod_message.o mod_vertint.o mod_date.o mod_space.o mod_memutil.o mod_dynparam.o mod_nchelper.o mod_humid.o mod_dynparam.o mod_memutil.o mod_sigma.o mod_stdio.o
clm45_1dto2d.o : clm45_1dto2d.F90 mod_dynparam.o mod_intkinds.o mod_nchelper.o mod_realkinds.o mod_stdio.o
