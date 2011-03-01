
FC = gfortran
FCFLAGS = -O2 -g

all:: lakeprog

OBJS = mod_constants.o mod_lake.o

lakeprog: lakeprog.f90 $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $< $(OBJS)

mod_constants.o: mod_constants.f90
	$(FC) $(FCFLAGS) -c $<

mod_lake.o: mod_lake.f90 mod_constants.o
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod lakeprog
