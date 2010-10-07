include ../../../Makefile.inc

all: CheckSun

OBJS = mod_date.o mod_constants.o

mod_constants.o: ../../../Config/mod_constants.F90
	$(F90) -c $(CPPFLAGS) $(NETCDFINC) $(F90FLAGS) $(USCORING) -c -o $@ $<

mod_date.o: ../../../PreProc/ICBC/mod_date.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c -o $@ $<

CheckSun: CheckSun.F90 $(OBJS)
	$(F90) $(LDFLAGS) $(CPPFLAGS) $(NETCDFINC) $(F90FLAGS) -o $@ $< $(OBJS) \
               $(NETCDFLIB)

%.o: %.F90
	$(F90) $(CPPFLAGS) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o CheckSun

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
