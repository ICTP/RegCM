
all :: sptest sptest_regcm

FC = `nf-config --fc`
INCS = -I. `nf-config --fflags`
LIBS = `nf-config --flibs`

sptest : mo_simple_plumes_v1.o sp_driver_v1.f90
	$(FC) $(INCS) -o $@ sp_driver_v1.f90 mo_simple_plumes_v1.o $(LIBS)
sptest_regcm : mo_simple_plumes_v1.o sp_driver_regcm.f90
	$(FC) $(INCS) -o $@ sp_driver_regcm.f90 mo_simple_plumes_v1.o $(LIBS)

mo_simple_plumes_v1.o : mo_simple_plumes_v1.f90
	$(FC) $(INCS) -c -o $@ $<

clean :
	rm -f sptest sptest_regcm mo_simple_plumes_v1.o *.mod
