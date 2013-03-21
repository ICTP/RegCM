
FC = gfortran
FCFLAGS = -O2 -g -I../../../Share

all:: nearpoint

nearpoint: nearpoint.f90
	$(FC) $(FCFLAGS) -o $@ $< ../../../Share/librcmlib.a

clean:
	rm -f *.o *.mod nearpoint
