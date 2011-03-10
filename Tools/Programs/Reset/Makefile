
include ../../../Makefile.inc

vpath %.F90 $(REGCM_BASE_DIR)/Config:.

all: reSAV

%.o: %.F90
	  $(MPIF90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) -c $<

reSAV: reSAV.F90 mod_constants.o mod_dynparam.o
	$(MPIF90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $< mod_constants.o mod_dynparam.o

clean:
	rm -f *.o *.mod reSAV
