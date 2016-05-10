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
	$(MPIFC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
mod_cbmz_precision.o: mod_cbmz_precision.F90 mod_realkinds.o
mod_cbmz_parameters.o: mod_cbmz_parameters.F90 mod_cbmz_precision.o
mod_cbmz_monitor.o: mod_cbmz_monitor.F90 mod_cbmz_precision.o
mod_cbmz_global.o: mod_cbmz_global.F90 mod_cbmz_parameters.o mod_cbmz_precision.o
mod_cbmz_Initialize.o: mod_cbmz_Initialize.F90  
mod_cbmz_function.o: mod_cbmz_function.F90  
mod_cbmz_Stochastic.o: mod_cbmz_Stochastic.F90  
mod_cbmz_jacobiansp.o: mod_cbmz_jacobiansp.F90 
mod_cbmz_jacobian.o: mod_cbmz_jacobian.F90   mod_cbmz_jacobiansp.o
mod_cbmz_linearalgebra.o: mod_cbmz_linearalgebra.F90  mod_cbmz_jacobiansp.o
mod_cbmz_rates.o: mod_cbmz_rates.F90  
mod_cbmz_HessianSP.o: mod_cbmz_HessianSP.F90  
mod_cbmz_Hessian.o:  mod_cbmz_Hessian.F90  mod_cbmz_HessianSP.o
mod_cbmz_StoichiomSP.o: mod_cbmz_StoichiomSP.F90 
mod_cbmz_Stoichiom.o: mod_cbmz_Stoichiom.F90   mod_cbmz_StoichiomSP.o
mod_cbmz_main.o: mod_cbmz_main.F90 mod_cbmz_global.o mod_cbmz_parameters.o  mod_cbmz_model.o mod_cbmz_integrator.o
mod_cbmz_model.o: mod_cbmz_model.F90 mod_cbmz_integrator.o
mod_cbmz_integrator.o: mod_cbmz_integrator.F90  

mod_cbmz_jval1.o : mod_cbmz_jval1.F90 mod_intkinds.o mod_realkinds.o mod_constants.o mod_dynparam.o mod_mppparam.o mod_stdio.o mod_cbmz_precision.o
mod_cbmz_hvread.o: mod_cbmz_hvread.F90 mod_cbmz_jval1.o mod_cbmz_precision.o
