
FC = ifort
FCFLAGS = -O2 -g -I${NETCDF}/include
LIBS = -L${NETCDF}/lib -lnetcdff -lnetcdf

all:: put_pointlake

put_pointlake: put_pointlake.f90
	$(FC) $(FCFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f *.o *.mod put_pointlake
