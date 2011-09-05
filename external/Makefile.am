
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/external/mpeu:$(REGCM_ROOT)/external/MCT/mpi-serial
vpath %.f90 $(REGCM_ROOT)/external/mpeu:$(REGCM_ROOT)/external/MCT/mpi-serial
vpath %.c $(REGCM_ROOT)/external/mpeu:$(REGCM_ROOT)/external/MCT/mpi-serial
vpath %.h $(REGCM_ROOT)/external/MCT/mpi-serial

if REGCM_PARALLEL
noinst_LIBRARIES = libmpeu.a
else
noinst_LIBRARIES = libmpiser.a libmpeu.a
endif
libmpiser_a_SOURCES = mpi.c send.c recv.c collective.c req.c list.c \
                      handles.c comm.c group.c time.c pack.c fort.F90 \
                      mpistub.f90
libmpeu_a_SOURCES = get_zeits.c m_chars.F90 m_die.F90 m_dropdead.F90      \
                    m_FcComms.F90 m_Filename.F90 m_FileResolv.F90         \
                    m_flow.F90 m_IndexBin_char.F90 m_IndexBin_integer.F90 \
                    m_IndexBin_logical.F90 m_inpak90.F90 m_ioutil.F90     \
                    m_List.F90 m_mall.F90 m_MergeSorts.F90 m_mpif90.F90   \
                    m_mpif.F90 m_mpout.F90 m_Permuter.F90 m_rankMerge.F90 \
                    m_realkinds.F90 m_SortingTools.F90 m_stdio.F90        \
                    m_String.F90 m_StringLinkedList.F90 m_StrTemplate.F90 \
                    m_TraceBack.F90 m_zeit.F90

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
m_IndexBin_char.o: m_die.o m_stdio.o
m_IndexBin_integer.o: m_die.o m_stdio.o
m_IndexBin_logical.o: m_die.o m_stdio.o
m_List.o: m_String.o m_die.o m_mall.o
m_MergeSorts.o: m_die.o m_realkinds.o m_stdio.o
m_Filename.o:
m_FcComms.o: m_die.o
m_StringLinkedList.o: m_String.o
m_Permuter.o: m_die.o m_realkinds.o
m_SortingTools.o: m_IndexBin_char.o m_IndexBin_integer.o m_IndexBin_logical.o m_MergeSorts.o m_Permuter.o m_rankMerge.o
m_String.o: m_die.o m_mall.o m_mpif90.o
m_StrTemplate.o: m_chars.o m_die.o m_stdio.o
m_chars.o:
m_die.o: m_dropdead.o m_flow.o m_mpif90.o m_mpout.o m_stdio.o
m_dropdead.o: m_mpif90.o m_stdio.o
m_flow.o: m_chars.o
m_inpak90.o: m_die.o m_ioutil.o m_mall.o m_mpif90.o m_realkinds.o m_stdio.o
m_ioutil.o: m_stdio.o
m_mall.o: m_chars.o m_die.o m_ioutil.o m_realkinds.o m_stdio.o
m_mpif.o:
m_mpif90.o: m_mpif.o m_realkinds.o m_stdio.o
m_mpout.o: m_dropdead.o m_ioutil.o m_mpif90.o m_stdio.o
m_rankMerge.o:
m_realkinds.o:
m_stdio.o:
m_zeit.o: m_SortingTools.o m_die.o m_ioutil.o m_mpif90.o m_stdio.o get_zeits.o
get_zeits.o:
m_FileResolv.o: m_die.o m_StrTemplate.o
m_TraceBack.o:	m_die.o m_stdio.o m_String.o
