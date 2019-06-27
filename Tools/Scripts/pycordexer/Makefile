# gfortran
export F2PYC=f2py3 --opt="-fopenmp" -lgomp
# intel
#export F2PYC=f2py3 --opt="-qopenmp" -liomp5
#export F2PYC=f2py3

all:
	$(MAKE) $(MFLAGS) -C utilities

clean:
	$(MAKE) $(MFLAGS) -C utilities clean
