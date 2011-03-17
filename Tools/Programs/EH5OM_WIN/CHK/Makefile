
FC = gfortran
FCFLAGS = -O2 -g -fconvert=big-endian

all: toGRADS1 toGRADS2

toGRADS1: toGRADS1.f90
	$(FC) $(FCFLAGS) -o $@ $<

toGRADS2: toGRADS2.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f toGRADS1 toGRADS2
