# Makefile for netlib

include ../../Makefile.inc

# check if we need the MPI wrapper 
ifneq ($(origin PARALLEL), undefined)
  F90 = $(MPIF90)
endif  

all:: libnetlib.a

NETLIBOBJS = linpack.o eispack.o

libnetlib.a: $(NETLIBOBJS)
	$(AR) $(ARFLAGS) $@ $(NETLIBOBJS)

%.o: %.F90
	$(F90) -c $(NETCDFINC) $(CPPFLAGS) $(F90FLAGS)  -c $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS)  -c $<

clean:
	rm -f *.o libnetlib.a
#
# Netlib libraries
#
linpack.o: linpack.f90
eispack.o: eispack.f90
