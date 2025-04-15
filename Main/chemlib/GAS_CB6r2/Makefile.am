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

noinst_LIBRARIES = libcb6.a
libcb6_a_SOURCES = mod_cb6_Precision.F90 mod_cb6_Parameters.F90 \
   mod_cb6_Global.F90 mod_cb6_Function.F90 mod_cb6_JacobianSP.F90 \
   mod_cb6_Jacobian.F90  mod_cb6_jval2.F90 mod_che_hvread.F90 \
   mod_cb6_Rates.F90 mod_cb6_Util.F90 mod_cb6_Monitor.F90  \
   mod_cb6_LinearAlgebra.F90 mod_cb6_Main.F90  mod_cb6_Integrator.F90 \
   mod_cb6_Model.F90

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
mod_cb6_Precision.o: mod_cb6_Precision.F90
mod_cb6_Parameters.o: mod_cb6_Parameters.F90 \
                    mod_cb6_Precision.o
mod_cb6_Monitor.o: mod_cb6_Monitor.F90 \
                     mod_cb6_Precision.o
mod_cb6_Global.o: mod_cb6_Global.F90 \
                    mod_cb6_Parameters.o mod_cb6_Precision.o
mod_cb6_Initialize.o: mod_cb6_Initialize.F90  
mod_cb6_Function.o: mod_cb6_Function.F90  
mod_cb6_Stochastic.o: mod_cb6_Stochastic.F90  
mod_cb6_JacobianSP.o: mod_cb6_JacobianSP.F90 
mod_cb6_Jacobian.o: mod_cb6_Jacobian.F90   mod_cb6_JacobianSP.o
mod_cb6_LinearAlgebra.o: mod_cb6_LinearAlgebra.F90  mod_cb6_JacobianSP.o
mod_cb6_Rates.o: mod_cb6_Rates.F90  
mod_cb6_HessianSP.o: mod_cb6_HessianSP.F90  
mod_cb6_Hessian.o:  mod_cb6_Hessian.F90  mod_cb6_HessianSP.o
mod_cb6_StoichiomSP.o: mod_cb6_StoichiomSP.F90 
mod_cb6_Stoichiom.o: mod_cb6_Stoichiom.F90   mod_cb6_StoichiomSP.o
mod_cb6_Util.o: mod_cb6_Util.F90   mod_cb6_Monitor.o
mod_cb6_Main.o: mod_cb6_Main.F90   mod_cb6_Global.o mod_cb6_Parameters.o  mod_cb6_Model.o mod_cb6_Integrator.o
mod_cb6_Model.o: mod_cb6_Model.F90   mod_cb6_Integrator.o
mod_cb6_Integrator.o: mod_cb6_Integrator.F90  

mod_cb6_jval2.o : mod_cb6_jval2.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_mppparam.o mod_stdio.o
mod_che_hvread.o: mod_che_hvread.F90 mod_cb6_jval2.o 
