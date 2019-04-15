
FC = gfortran
FCFLAGS = -O2 -g -Wall -pedantic `nf-config --fflags`
LIBS = `nf-config --flibs`

all:: plot

plot: cumclouds
	./cumclouds
	paste gnuplot_*.dat > dat.dat
	gnuplot plotclf.gp

cumclouds: cumclouds.F90
	$(FC) $(FCFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f cumclouds profile_*.txt gnuplot_*.dat dat.dat *.nc *.ps
