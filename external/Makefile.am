
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/external/mpi-serial
vpath %.f90 $(REGCM_ROOT)/external/mpi-serial
vpath %.c $(REGCM_ROOT)/external/mpi-serial
vpath %.h $(REGCM_ROOT)/external/MCT/mpi-serial

if REGCM_PARALLEL
noinst_LIBRARIES =
else
noinst_LIBRARIES = libmpiser.a
endif
libmpiser_a_SOURCES = mpi.c send.c recv.c collective.c req.c list.c \
                      handles.c comm.c group.c time.c pack.c fort.F90 \
                      mpistub.f90

MPIFH= mpif.$(FORT_SIZE).h

mpif.h:
	cp -f $(REGCM_ROOT)/external/MCT/mpi-serial/$(MPIFH) mpif.h

%.o: %.c
	$(CC) $(CFLAGS) $(MPEUCFLAGS) -c $<

%.o: %.f90
	  $(FC) $(CPPFLAGS) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

%.o: %.F90
	  $(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mpi.o: mpif.h
fort.o: fort.F90
mpistub.o: mpistub.f90
