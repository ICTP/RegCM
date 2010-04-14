
include ../../Makefile.inc

all: small3D small smallHT smallPS

small3D: small3D.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

small: small.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

smallHT: smallHT.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

smallPS: smallPS.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

clean:
	rm -f small3D small smallHT smallPS
