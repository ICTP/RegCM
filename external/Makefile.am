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

include $(REGCM_ROOT)/makeinc

vpath %.f90 $(REGCM_ROOT)/external/aerosols

noinst_LIBRARIES = libaerosol.a libbmi.a

libaerosol_a_SOURCES = mo_simple_plumes.f90
libbmi_a_SOURCES = bmi.f90

%.o: %.c
	$(CC) $(CFLAGS) $(MPEUCFLAGS) -c $<

%.o: %.f90
	  $(FC) $(CPPFLAGS) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

%.o: %.F90
	  $(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt mpif.h

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt
