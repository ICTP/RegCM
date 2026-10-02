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

noinst_LIBRARIES = libcmip7.a
libcmip7_a_SOURCES = mod_cmip7_helper.F90 mod_cmip7_ecea.F90

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mod_cmip7_helper.o : mod_cmip7_helper.F90 mod_realkinds.o mod_intkinds.o mod_constants.o mod_date.o mod_message.o mod_dynparam.o mod_kdinterp.o mod_stdio.o
mod_cmip7_ecea.o : mod_cmip7_ecea.F90 mod_cmip7_helper.o mod_realkinds.o mod_intkinds.o mod_message.o mod_date.o mod_stdio.o mod_dynparam.o mod_memutil.o mod_kdinterp.o
