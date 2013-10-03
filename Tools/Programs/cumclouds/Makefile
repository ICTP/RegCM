
FC = gfortran
FCFLAGS = -O2 -g -Wall -pedantic

all:: plot

plot: cumclouds
	./cumclouds
	paste gnuplot_*.dat > dat.dat
	gnuplot plotclf.gp

cumclouds: cumclouds.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f cumclouds profile_*.txt gnuplot_*.dat dat.dat
