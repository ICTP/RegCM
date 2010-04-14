# Makefile for RegCM4
# for Input tools code.
  
include ../../Makefile.inc

all: day_ICBC mon_ICBC sigmatop

SIGTOPOBJ = height.o htsig.o humid1.o humid2.o intlin.o intlog.o slpres.o \
            sigmatop.o

SHRMOD = mod_regcm_param.o mod_preproc_param.o mod_constants.o

mod_regcm_param.o: ../../Config/mod_regcm_param.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

mod_preproc_param.o: ../../Config/mod_preproc_param.f90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

mod_constants.o: ../../Config/mod_constants.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) \
            $(USCORING) -c -o $@ $<

header_ICBC.o: header_ICBC.F90
	$(F90) $(CPPFLAGS) $(SVNDEF) $(F90FLAGS) $(INTELINC) $(USCORING) -c $<

day_ICBC: $(SHRMOD) day_ICBC.o
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ day_ICBC.o $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

mon_ICBC: $(SHRMOD) mon_ICBC.o
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ mon_ICBC.o $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

sigmatop: $(SHRMOD) $(SIGTOPOBJ)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(SIGTOPOBJ) $(SHRMOD) \
               $(LIB) $(NETCDFLIB)

%.o: %.f90
	$(F90) -c $(NETCDFINC) $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o day_ICBC mon_ICBC sigmatop

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

day_ICBC.o: day_ICBC.f90
height.o: height.f90 mod_constants.o mod_preproc_param.o
htsig.o: htsig.f90 mod_constants.o mod_preproc_param.o
humid1.o: humid1.f90 mod_constants.o mod_preproc_param.o
humid2.o: humid2.f90 mod_constants.o mod_preproc_param.o
intlin.o: intlin.f90 mod_constants.o mod_preproc_param.o
intlog.o: intlog.f90 mod_constants.o mod_preproc_param.o
mon_ICBC.o: mon_ICBC.f90
sigmatop.o: sigmatop.f90 mod_regcm_param.o
slpres.o: slpres.f90 mod_constants.o mod_preproc_param.o
