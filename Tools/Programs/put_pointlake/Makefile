
FC = gfortran
FCFLAGS = -O2 -g -I`nf-config --includedir`
LIBS = `nf-config --flibs`

all:: put_pointlake

put_pointlake: put_pointlake.F90
	$(FC) $(FCFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f *.o *.mod put_pointlake
