
noinst_LIBRARIES = libkpp.a
libkpp_a_SOURCES = gasphs_precision.f90 gasphs_parameters.f90 \
                   gasphs_global.f90 gasphs_function.f90      \
                   gasphs_jacobiansp.f90 gasphs_jacobian.f90  \
                   gasphs_stoichiomsp.f90 gasphs_stoichiom.f90 \
                   gasphs_rates.f90 gasphs_util.f90 gasphs_monitor.f90 \
                   gasphs_linearalgebra.f90 gasphs_driver.f90 \
                   gasphs_initialize.f90 gasphs_integrator.F90 \
                   gasphs_model.f90

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

gasphs_precision.o: gasphs_precision.f90
gasphs_parameters.o: gasphs_parameters.f90 gasphs_precision.o
gasphs_monitor.o: gasphs_monitor.f90 gasphs_precision.o
gasphs_global.o: gasphs_global.f90 gasphs_parameters.o gasphs_precision.o
gasphs_initialize.o: gasphs_initialize.f90 gasphs_precision.o \
                     gasphs_parameters.o gasphs_global.o
gasphs_function.o: gasphs_function.f90 gasphs_precision.o \
                   gasphs_parameters.o gasphs_global.o
gasphs_jacobiansp.o: gasphs_jacobiansp.f90 gasphs_precision.o \
                     gasphs_parameters.o gasphs_global.o
gasphs_jacobian.o: gasphs_jacobian.f90 gasphs_jacobiansp.o gasphs_precision.o \
                   gasphs_parameters.o gasphs_global.o
gasphs_linearalgebra.o: gasphs_linearalgebra.f90 gasphs_jacobiansp.o \
                        gasphs_precision.o gasphs_parameters.o gasphs_global.o
gasphs_rates.o: gasphs_rates.f90 gasphs_precision.o gasphs_parameters.o \
                gasphs_global.o
gasphs_stoichiomsp.o: gasphs_stoichiomsp.f90 gasphs_precision.o \
                      gasphs_parameters.o gasphs_global.o
gasphs_stoichiom.o: gasphs_stoichiom.f90 gasphs_stoichiomsp.o \
                    gasphs_precision.o gasphs_parameters.o gasphs_global.o
gasphs_util.o: gasphs_util.f90 gasphs_monitor.o gasphs_precision.o \
               gasphs_parameters.o gasphs_global.o
gasphs_integrator.o: gasphs_integrator.F90 gasphs_util.o gasphs_initialize.o \
                     gasphs_function.o gasphs_linearalgebra.o gasphs_rates.o \
                     gasphs_stoichiom.o
gasphs_model.o: gasphs_model.f90 gasphs_integrator.o
gasphs_driver.o: gasphs_model.o
