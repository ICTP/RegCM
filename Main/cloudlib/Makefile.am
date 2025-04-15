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

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:.
vpath %.inc $(REGCM_ROOT)/Share:.

noinst_LIBRARIES = libcloud.a
libcloud_a_SOURCES = mod_cloud_subex.F90 mod_cloud_xuran.F90 \
                     mod_cloud_thomp.F90 mod_cloud_guli2007.F90 \
                     mod_cloud_tompkins.F90 mod_cloud_texeira.F90 \
                     mod_cloud_echam5.F90

REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE) 
else
FCFLAGS += -I$(MPP_INCLUDE)
endif

%.o: %.F90
	$(MPIFC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

clean-local:
	rm -f *.mod *__genmod.* *opt.yaml *optrpt

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_cloud_subex.o : mod_cloud_subex.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o
mod_cloud_texeira.o : mod_cloud_texeira.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o
mod_cloud_tompkins.o : mod_cloud_tompkins.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o
mod_cloud_echam5.o : mod_cloud_echam5.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o
mod_cloud_xuran.o : mod_cloud_xuran.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o
mod_cloud_thomp.o : mod_cloud_thomp.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_runparams.o
mod_cloud_guli2007.o : mod_cloud_guli2007.F90 mod_intkinds.o mod_realkinds.o mod_dynparam.o
