
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/external/mpi-serial
vpath %.f90 $(REGCM_ROOT)/external/aerosols
vpath %.c $(REGCM_ROOT)/external/mpi-serial
vpath %.h $(REGCM_ROOT)/external/mpi-serial

if REGCM_PARALLEL
noinst_LIBRARIES = libaerosol.a libbmi.a
else
noinst_LIBRARIES = libmpiser.a libaerosol.a libbmi.a
endif
libmpiser_a_SOURCES = mpi.c send.c recv.c collective.c req.c list.c \
                      handles.c comm.c group.c time.c pack.c fort.F90 \
                      mpistub.F90
libaerosol_a_SOURCES = mo_simple_plumes.f90
libbmi_a_SOURCES = bmi.f90

MPIFH= mpif.$(FORT_SIZE).h

mpif.h:
	cp -f $(REGCM_ROOT)/external/mpi-serial/$(MPIFH) mpif.h

%.o: %.c
	$(CC) $(CFLAGS) $(MPEUCFLAGS) -c $<

%.o: %.f90
	  $(FC) $(CPPFLAGS) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

%.o: %.F90
	  $(FC) $(CPPFLAGS) $(GITVER) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

distclean-local:
	rm -f *.mod mpif.h

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mpi.o: mpif.h
fort.o: fort.F90
mpistub.o: mpistub.F90
