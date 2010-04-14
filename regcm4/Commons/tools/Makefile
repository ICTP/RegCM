
include ../../Makefile.inc

all: readCRU reSAV resetSST

mod_regcm_param.o: ../../Config/mod_regcm_param.F90
	$(F90) $(CPPFLAGS) $(F90FLAGS) $(USCORING) -c -o $@ $<

readCRU: readCRU.f90 mod_regcm_param.o
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $< mod_regcm_param.o

reSAV: reSAV.f90 mod_regcm_param.o
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $< mod_regcm_param.o

resetSST: resetSST.f90 mod_regcm_param.o
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
           $(USCORING) -o $@ $< mod_regcm_param.o

clean:
	rm -f *.o *.mod readCRU reSAV resetSST
