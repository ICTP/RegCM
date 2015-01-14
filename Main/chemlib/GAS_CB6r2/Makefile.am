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
libcbmz_a_SOURCES = cb6_Precision.F90 cb6_Parameters.F90 \
   cb6_Global.F90 cb6_Function.F90 cb6_JacobianSP.F90 \
   cb6_Jacobian.F90  cb6_jval1.F90 mod_che_hvread.F90 \
   cb6_Rates.F90 cb6_Util.F90 cb6_Monitor.F90  \
   cb6_LinearAlgebra.F90 cb6_Main.F90  cb6_Integrator.F90 \
   cb6_Model.F90

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
cb6_Precision.o: cb6_Precision.F90
cb6_Parameters.o: cb6_Parameters.F90 \
                    cb6_Precision.o
cb6_Monitor.o: cb6_Monitor.F90 \
                     cb6_Precision.o
cb6_Global.o: cb6_Global.F90 \
                    cb6_Parameters.o cb6_Precision.o
cb6_Initialize.o: cb6_Initialize.F90  
cb6_Function.o: cb6_Function.F90  
cb6_Stochastic.o: cb6_Stochastic.F90  
cb6_JacobianSP.o: cb6_JacobianSP.F90 
cb6_Jacobian.o: cb6_Jacobian.F90   cb6_JacobianSP.o
cb6_LinearAlgebra.o: cb6_LinearAlgebra.F90  cb6_JacobianSP.o
cb6_Rates.o: cb6_Rates.F90  
cb6_HessianSP.o: cb6_HessianSP.F90  
cb6_Hessian.o:  cb6_Hessian.F90  cb6_HessianSP.o
cb6_StoichiomSP.o: cb6_StoichiomSP.F90 
cb6_Stoichiom.o: cb6_Stoichiom.F90   cb6_StoichiomSP.o
cb6_Util.o: cb6_Util.F90   cb6_Monitor.o
cb6_Main.o: cb6_Main.F90   cb6_Global.o cb6_Parameters.o  cb6_Model.o cb6_Integrator.o
cb6_Model.o: cb6_Model.F90   cb6_Integrator.o
cb6_Integrator.o: cb6_Integrator.F90  

cb6_jval1.o : cb6_jval1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_mppparam.o mod_stdio.o
mod_che_hvread.o: mod_che_hvread.F90 cb6_jval1.o 
