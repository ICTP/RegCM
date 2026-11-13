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

# Single column tests of the physics schemes, built and run only by
# "make check". See README for the scm_tiedtke.py wrapper.

include $(REGCM_ROOT)/makeinc

MPPLIBPTH = $(REGCM_ROOT)/Main/mpplib
CUMLIBPTH = $(REGCM_ROOT)/Main/cumlib

CPPFLAGS += -I$(MPPLIBPTH) -I$(CUMLIBPTH)

check_PROGRAMS = scm_tiedtke
scm_tiedtke_SOURCES = scm_tiedtke.F90
scm_tiedtke_LDADD = $(CUMLIBPTH)/libcum.a $(MPPLIBPTH)/libmpp.a $(OMPFLAGS)

FCLD = $(MPIFC)

TESTS = run_scm_tests.sh

EXTRA_DIST = README make_profiles.py scm_tiedtke.py run_scm_tests.sh \
             scm_profiles.txt reference/scm_tiedtke.out

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt
	rm -rf scm_output

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt
	rm -rf scm_output
