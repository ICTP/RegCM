
include ../../../../Makefile.inc

all: toGRADS1 toGRADS2

toGRADS1: toGRADS1.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

toGRADS2: toGRADS2.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $<

clean:
	rm -f toGRADS1 toGRADS2
