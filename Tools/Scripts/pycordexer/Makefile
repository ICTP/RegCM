# gfortran
export F2PYC=f2py --opt="-fopenmp" -lgomp
# intel
#INTELPATH = $(shell dirname `which ifort` | sed -e 's!/bin/!/lib/!')
#export F2PYC=f2py3 --fcompiler=intelem --opt="-qopenmp" -L$(INTELPATH) -liomp5 -lsvml
#export F2PYC=f2py3

all:
	$(MAKE) $(MFLAGS) -C utilities

clean:
	$(MAKE) $(MFLAGS) -C utilities clean
