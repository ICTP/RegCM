# Makefile for RegCM4
  
include ../../../Makefile.inc

all: ReadRegcmio

OBJS = mod_date.o

mod_date.o: ../../../PreProc/ICBC/mod_date.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c -o $@ $<

ReadRegcmio: ReadRegcmio.F90 $(OBJS)
	$(F90) $(LDFLAGS) $(CPPFLAGS) $(NETCDFINC) $(F90FLAGS) -o $@ $< $(OBJS) \
               $(NETCDFLIB)

%.o: %.F90
	$(F90) $(CPPFLAGS) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o ReadRegcmio

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
