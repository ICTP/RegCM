
FC = ifort
FCFLAGS = -I${NETCDF}/include -I. -O3
LDFLAGS = -L${NETCDF}/lib
LIBS = -lnetcdff -lnetcdf
OBJS = mod_date.o

all: int_dur

mod_date.o: mod_date.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

int_dur: int_dur.f90 $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $< $(LDFLAGS) $(LIBS) $(OBJS)

clean:
	rm -f int_dur *.o *.mod
