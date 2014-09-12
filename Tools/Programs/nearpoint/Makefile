
FC = gfortran
FCFLAGS = -O2 -g -I../../../Share

all:: nearpoint

nearpoint: nearpoint.F90
	$(FC) $(FCFLAGS) -o $@ $< ../../../Share/librcmlib.a

clean:
	rm -f *.o *.mod nearpoint
