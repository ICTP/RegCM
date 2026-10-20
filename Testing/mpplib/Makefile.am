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

# Regression test of Main/mpplib/mod_mppparam.F90, built and run only by
# "make check". See run_mpp_regression.sh for the options.

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:.

include $(REGCM_ROOT)/makeinc

MPPLIBPTH = $(REGCM_ROOT)/Main/mpplib

if COMPILER_SUN
  CPPFLAGS += -M$(MPPLIBPTH)
else
  CPPFLAGS += -I$(MPPLIBPTH)
endif

check_PROGRAMS = mpp_regression
mpp_regression_SOURCES = mpp_regression.F90
mpp_regression_LDADD = $(MPPLIBPTH)/libmpp.a $(OMPFLAGS)

FCLD = $(MPIFC)

TESTS = run_mpp_regression.sh

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt
	rm -rf mpp_output

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt
	rm -rf mpp_output

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mpp_regression.o : mpp_regression.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_memutil.o mod_date.o mod_mppparam.o mod_regcm_types.o
