
FC = gfortran
FCFLAGS = -O2 -g

all:: lakeprog

OBJS = mod_constants.o mod_lake.o

lakeprog: lakeprog.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $< $(OBJS)

mod_constants.o: mod_constants.F90
	$(FC) $(FCFLAGS) -c $<

mod_lake.o: mod_lake.F90 mod_constants.o
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod lakeprog
