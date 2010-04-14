
include ../../../Makefile.inc

all: day_ATM day_RAD day_SRF day_SRF2

day_ATM: day_ATM.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

day_RAD: day_RAD.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

day_SRF: day_SRF.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

day_SRF2: day_SRF2.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

clean:
	rm -f day_ATM day_RAD day_SRF day_SRF2
