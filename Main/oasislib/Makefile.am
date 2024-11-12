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

include $(REGCM_ROOT)/makeinc

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:.

MPPLIBPTH = $(REGCM_ROOT)/Main/mpplib
MPPLIB = -L$(MPPLIBPTH) -lmpp
if COMPILER_SUN
  MPPINC = -M$(MPPLIBPTH)
else
  MPPINC = -I$(MPPLIBPTH)
endif

noinst_LIBRARIES = liboasis.a

if DO_COMPILE_OASIS
liboasis_a_SOURCES = mod_oasis_params.F90 mod_oasis_signature.F90 \
                     mod_oasis_generic.F90
else
liboasis_a_SOURCES = mod_oasis_void.F90
endif

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_oasis_signature.o : mod_oasis_signature.F90 mod_mppparam.o mod_stdio.o mod_dynparam.o
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) $(MPPINC) -c $<
mod_oasis_generic.o : mod_oasis_generic.F90 mod_oasis_params.o mod_oasis_signature.o mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_message.o mod_service.o mod_memutil.o mod_constants.o mod_stdio.o
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) $(MPPINC) -c $<
mod_oasis_params.o : mod_oasis_params.F90 mod_intkinds.o
