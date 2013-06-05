
FC = gfortran
FCFLAGS = -O2 -g -Wall -pedantic

all:: cumclouds

cumclouds: cumclouds.f90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f cumclouds cumcloud_*.txt
