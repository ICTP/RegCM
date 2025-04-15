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

noinst_LIBRARIES = libcbmz.a
libcbmz_a_SOURCES = mod_cbmz_precision.F90 mod_cbmz_parameters.F90 \
   mod_cbmz_global.F90 mod_cbmz_function.F90 mod_cbmz_jacobiansp.F90 \
   mod_cbmz_jacobian.F90 mod_cbmz_jval1.F90 mod_cbmz_hvread.F90 \
   mod_cbmz_rates.F90 mod_cbmz_monitor.F90 mod_cbmz_linearalgebra.F90 \
   mod_cbmz_main.F90 mod_cbmz_integrator.F90 mod_cbmz_model.F90

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
mod_cbmz_precision.o: mod_cbmz_precision.F90 mod_realkinds.o
mod_cbmz_parameters.o: mod_cbmz_parameters.F90 mod_cbmz_precision.o
mod_cbmz_monitor.o: mod_cbmz_monitor.F90 mod_cbmz_precision.o
mod_cbmz_global.o: mod_cbmz_global.F90 mod_cbmz_parameters.o mod_cbmz_precision.o
mod_cbmz_Initialize.o: mod_cbmz_Initialize.F90
mod_cbmz_function.o: mod_cbmz_function.F90 mod_cbmz_parameters.o 
mod_cbmz_Stochastic.o: mod_cbmz_Stochastic.F90
mod_cbmz_jacobiansp.o: mod_cbmz_jacobiansp.F90 
mod_cbmz_jacobian.o: mod_cbmz_jacobian.F90 mod_cbmz_jacobiansp.o mod_cbmz_parameters.o
mod_cbmz_linearalgebra.o: mod_cbmz_linearalgebra.F90 mod_cbmz_jacobiansp.o mod_cbmz_parameters.o
mod_cbmz_rates.o: mod_cbmz_rates.F90 mod_cbmz_global.o
mod_cbmz_HessianSP.o: mod_cbmz_HessianSP.F90
mod_cbmz_Hessian.o: mod_cbmz_Hessian.F90 mod_cbmz_HessianSP.o
mod_cbmz_StoichiomSP.o: mod_cbmz_StoichiomSP.F90 
mod_cbmz_Stoichiom.o: mod_cbmz_Stoichiom.F90 mod_cbmz_StoichiomSP.o
mod_cbmz_main.o: mod_cbmz_main.F90 mod_cbmz_global.o mod_cbmz_parameters.o mod_cbmz_model.o mod_cbmz_integrator.o
mod_cbmz_model.o: mod_cbmz_model.F90 mod_cbmz_integrator.o
mod_cbmz_integrator.o: mod_cbmz_integrator.F90 mod_cbmz_global.o mod_cbmz_function.o mod_cbmz_rates.o mod_cbmz_linearalgebra.o mod_cbmz_parameters.o
mod_cbmz_jval1.o : mod_cbmz_jval1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_mppparam.o mod_stdio.o mod_cbmz_precision.o
mod_cbmz_hvread.o: mod_cbmz_hvread.F90 mod_cbmz_jval1.o mod_cbmz_precision.o
