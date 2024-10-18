#FC = gfortran
FC = ifort
#FFLAGS = `nf-config --fflags` -fopenmp
FFLAGS = `nf-config --fflags` -qopenmp
LIBS = `nf-config --flibs`

all :: heatindex

mod_heatindex.o : mod_heatindex.F90
	$(FC) $(FFLAGS) -c $<

heatindex : heatindex.F90 mod_heatindex.o
	$(FC) $(FFLAGS) -o $@ heatindex.F90 mod_heatindex.o $(LIBS)

clean :
	rm -f *.o *.mod heatindex

