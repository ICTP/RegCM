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

vpath %.o $(REGCM_ROOT)/Share:.

bin_PROGRAMS = terrain
terrain_SOURCES = mod_header.F90 mod_fudge.F90 mod_maps.F90 \
                mod_write.F90 mod_smooth.F90 mod_maputils.F90 \
                mod_moist.F90 terrain.F90

FCLD = $(MPIFC)

LDADD = $(OMPFLAGS)

%.o: %.F90
	$(FC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_fudge.o : mod_fudge.F90 mod_intkinds.o mod_realkinds.o mod_stdio.o mod_message.o
mod_header.o : mod_header.F90 mod_intkinds.o mod_stdio.o $(REGCM_ROOT)/version
mod_maps.o : mod_maps.F90 mod_intkinds.o mod_realkinds.o mod_memutil.o
mod_maputils.o : mod_maputils.F90 mod_constants.o mod_intkinds.o mod_realkinds.o mod_projections.o
mod_smooth.o : mod_smooth.F90 mod_intkinds.o mod_realkinds.o mod_constants.o
mod_write.o : mod_write.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o mod_ncstream_types.o mod_ncstream.o mod_memutil.o
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
mod_moist.o: mod_moist.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_stdio.o
terrain.o : terrain.F90 mod_realkinds.o mod_stdio.o mod_constants.o mod_message.o mod_date.o mod_dynparam.o mod_space.o mod_memutil.o mod_maps.o mod_smooth.o mod_projections.o mod_maputils.o mod_intldtr.o mod_fudge.o mod_nchelper.o mod_rdldtr.o mod_ncstream_types.o mod_ncstream.o mod_write.o mod_header.o mod_moist.o mod_sigma.o mod_nhinterp.o mod_earth.o mod_stdatm.o mod_zita.o mod_date.o
	$(MPIFC) $(CPPFLAGS) $(OMPFLAGS) $(GITVER) $(FCFLAGS) -c $<
