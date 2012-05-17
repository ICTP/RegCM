
FC = gfortran
FCFLAGS = -O2 -g -I`nf-config --includedir`
LIBS = `nf-config --flibs`

all:: changeland

changeland: changeland.f90
	$(FC) $(FCFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f *.o *.mod changeland
