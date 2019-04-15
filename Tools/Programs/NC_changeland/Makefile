
FC = gfortran
FCFLAGS = -O2 -g -I`nf-config --includedir`
LIBS = `nf-config --flibs`

all:: changeland

changeland: changeland.F90
	$(FC) $(FCFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f *.o *.mod changeland
