
FC = gfortran
FCFLAGS = -I`nf-config --includedir` -I. -O3
LDFLAGS = 
LIBS = `nf-config --flibs`
OBJS = mod_date.o

all: int_dur

mod_date.o: mod_date.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

int_dur: int_dur.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $< $(LDFLAGS) $(LIBS) $(OBJS)

clean:
	rm -f int_dur *.o *.mod
