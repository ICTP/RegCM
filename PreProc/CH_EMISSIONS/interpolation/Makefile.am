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

SUBDIRS = tools

include $(REGCM_ROOT)/makeinc

install-exec-hook:
	cp interp_emissions.sh $(DESTDIR)$(bindir)/interp_emissions
	cp interp_pollen.sh $(DESTDIR)$(bindir)/interp_pollen
	cp interp_bionox.sh $(DESTDIR)$(bindir)/interp_bionox
	cp interp_dustparam.sh $(DESTDIR)$(bindir)/interp_dustparam
	cp interp_emissions_bioburn.sh $(DESTDIR)$(bindir)/interp_emissions_bioburn.sh

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
