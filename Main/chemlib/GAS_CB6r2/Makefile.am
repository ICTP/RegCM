#n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

noinst_LIBRARIES = libcbmz.a
libcbmz_a_SOURCES = mod_cb6_Precision.F90 mod_cb6_Parameters.F90 \
   mod_cb6_Global.F90 mod_cb6_Function.F90 mod_cb6_JacobianSP.F90 \
   mod_cb6_Jacobian.F90  mod_cb6_jval1.F90 mod_che_hvread.F90 \
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
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

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

mod_cb6_jval1.o : mod_cb6_jval1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_mppparam.o mod_stdio.o
mod_che_hvread.o: mod_che_hvread.F90 mod_cb6_jval1.o 
