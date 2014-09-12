
FC = gfortran
FCFLAGS = -I`nf-config --includedir` -I. -O3
LDFLAGS = 
LIBS = `nf-config --flibs`

all: trend

trend: trend.F90 $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $< $(LDFLAGS) $(LIBS)

clean:
	rm -f trend
