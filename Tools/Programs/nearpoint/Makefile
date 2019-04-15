
FC = `nf-config --fc`
FCFLAGS = -O2 `nf-config --fflags` -g -I../../../Share

all:: nearpoint

nearpoint: nearpoint.F90
	$(FC) $(FCFLAGS) -o $@ $< ../../../Share/librcmlib.a `nf-config --flibs`

clean:
	rm -f *.o *.mod nearpoint
