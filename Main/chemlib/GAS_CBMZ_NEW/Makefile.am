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

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:.

noinst_LIBRARIES = libcbmz.a
libcbmz_a_SOURCES = mod_cbmz_Precision.f90 mod_cbmz_Parameters.f90 mod_cbmz_Global.f90   \
                   mod_cbmz_Function.f90 mod_cbmz_JacobianSP.f90 mod_cbmz_Jacobian.f90  \
                   mod_cbmz_jval1.f90 mod_che_hvread.f90 mod_cbmz_Rates.f90 mod_cbmz_Util.f90 mod_cbmz_Monitor.f90  \
                   mod_cbmz_LinearAlgebra.f90 mod_cbmz_Main.f90  \
                   mod_cbmz_Integrator.F90 mod_cbmz_Model.f90





REGCMSRC = $(REGCM_ROOT)/Main
MPP = $(REGCMSRC)/mpplib
MPP_INCLUDE = $(MPP)
if COMPILER_SUN
FCFLAGS += -M$(MPP_INCLUDE)
else
FCFLAGS += -I$(MPP_INCLUDE)
endif

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -cpp -c $<
%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -cpp -c $<
distclean-local:
	rm -f *.mod
clean-local:
	rm -f *.mod
# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_cbmz_Precision.o: mod_cbmz_Precision.f90
mod_cbmz_Parameters.o: mod_cbmz_Parameters.f90 \
                    mod_cbmz_Precision.o
mod_cbmz_Monitor.o: mod_cbmz_Monitor.f90 \
                     mod_cbmz_Precision.o
mod_cbmz_Global.o: mod_cbmz_Global.f90 \
                    mod_cbmz_Parameters.o mod_cbmz_Precision.o
mod_cbmz_Initialize.o: mod_cbmz_Initialize.f90  
mod_cbmz_Function.o: mod_cbmz_Function.f90  
mod_cbmz_Stochastic.o: mod_cbmz_Stochastic.f90  
mod_cbmz_JacobianSP.o: mod_cbmz_JacobianSP.f90 
mod_cbmz_Jacobian.o: mod_cbmz_Jacobian.f90   mod_cbmz_JacobianSP.o
mod_cbmz_LinearAlgebra.o: mod_cbmz_LinearAlgebra.f90  mod_cbmz_JacobianSP.o
mod_cbmz_Rates.o: mod_cbmz_Rates.f90  
mod_cbmz_HessianSP.o: mod_cbmz_HessianSP.f90  
mod_cbmz_Hessian.o:  mod_cbmz_Hessian.f90  mod_cbmz_HessianSP.o
mod_cbmz_StoichiomSP.o: mod_cbmz_StoichiomSP.f90 
mod_cbmz_Stoichiom.o: mod_cbmz_Stoichiom.f90   mod_cbmz_StoichiomSP.o
mod_cbmz_Util.o: mod_cbmz_Util.f90   mod_cbmz_Monitor.o
mod_cbmz_Main.o: mod_cbmz_Main.f90   mod_cbmz_Global.o mod_cbmz_Parameters.o  mod_cbmz_Model.o mod_cbmz_Integrator.o
mod_cbmz_Model.o: mod_cbmz_Model.f90   mod_cbmz_Integrator.o
mod_cbmz_Integrator.o: mod_cbmz_Integrator.F90  

mod_cbmz_jval1.o : mod_cbmz_jval1.f90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_mppparam.o mod_stdio.o
mod_che_hvread.o: mod_che_hvread.f90 mod_cbmz_jval1.o 
