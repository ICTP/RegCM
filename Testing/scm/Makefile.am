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
PBLLIBPTH = $(REGCM_ROOT)/Main/pbllib
OCNLIBPTH = $(REGCM_ROOT)/Main/ocnlib
MICLIBPTH = $(REGCM_ROOT)/Main/microlib
CLDLIBPTH = $(REGCM_ROOT)/Main/cloudlib
RADLIBPTH = $(REGCM_ROOT)/Main/radlib

CPPFLAGS += -I$(MPPLIBPTH) -I$(CUMLIBPTH) -I$(PBLLIBPTH) -I$(OCNLIBPTH) \
            -I$(MICLIBPTH) -I$(CLDLIBPTH) -I$(RADLIBPTH)

check_PROGRAMS = scm_tiedtke scm_holtslag scm_coupled scm_ocean
scm_tiedtke_SOURCES = scm_tiedtke.F90
scm_tiedtke_LDADD = $(CUMLIBPTH)/libcum.a $(MPPLIBPTH)/libmpp.a $(OMPFLAGS)
scm_holtslag_SOURCES = scm_holtslag.F90
scm_holtslag_LDADD = $(PBLLIBPTH)/libpbl.a $(MPPLIBPTH)/libmpp.a $(OMPFLAGS)
scm_coupled_SOURCES = scm_coupled.F90
scm_coupled_LDADD = $(CUMLIBPTH)/libcum.a $(PBLLIBPTH)/libpbl.a \
                    $(OCNLIBPTH)/libocn.a $(MICLIBPTH)/libmicro.a \
                    $(CLDLIBPTH)/libcloud.a $(RADLIBPTH)/librad.a \
                    $(RADLIBPTH)/RRTMG/librrtmg.a \
                    $(REGCM_ROOT)/external/libaerosol.a \
                    $(REGCM_ROOT)/external/libcosp.a \
                    $(MPPLIBPTH)/libmpp.a $(OMPFLAGS)
scm_ocean_SOURCES = scm_ocean.F90
scm_ocean_LDADD = $(OCNLIBPTH)/libocn.a $(MPPLIBPTH)/libmpp.a $(OMPFLAGS)

FCLD = $(MPIFC)

TESTS = run_scm_tests.sh

EXTRA_DIST = README make_profiles.py scm_tiedtke.py scm_holtslag.py \
             run_scm_tests.sh scm_profiles.txt reference/scm_tiedtke.out \
             reference/scm_holtslag.out testgate make_gate_profiles.py \
             scm_gate_check.txt scm_gate_check_obs.txt \
             reference/scm_tiedtke_gate.out scm_coupled.py \
             make_amazon_case.py scm_amazon.txt scm_amazon_target.txt \
             scm_amazon_era5.txt reference/scm_coupled_amazon.out \
             scm_itcz.txt scm_itcz_target.txt scm_itcz_era5.txt \
             reference/scm_coupled_itcz.out \
             reference/scm_coupled_kf_amazon.out \
             reference/scm_coupled_kf_itcz.out scm_ocean.py \
             make_dephy_case.py scm_les.py scm_bomex.txt \
             scm_bomex_target.txt scm_bomex_les.txt scm_rico.txt \
             scm_rico_target.txt scm_rico_les.txt \
             reference/scm_coupled_bomex.out reference/scm_coupled_rico.out

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt
	rm -rf scm_output

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt
	rm -rf scm_output
