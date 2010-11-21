
include $(REGCM_ROOT)/makeinc

vpath %.F90 $(REGCM_ROOT)/external/MCT/mpeu:$(REGCM_ROOT)/external/MCT/mct:\
	          $(REGCM_ROOT)/external/MCT/mpi-serial
vpath %.f90 $(REGCM_ROOT)/external/MCT/mpeu:$(REGCM_ROOT)/external/MCT/mct:\
	          $(REGCM_ROOT)/external/MCT/mpi-serial
vpath %.c $(REGCM_ROOT)/external/MCT/mpeu:$(REGCM_ROOT)/external/MCT/mct:\
	        $(REGCM_ROOT)/external/MCT/mpi-serial
vpath %.h $(REGCM_ROOT)/external/MCT/mpi-serial

if MPISERIAL
noinst_LIBRARIES = libmpeu.a libmct.a
else
noinst_LIBRARIES = libmpiser.a libmpeu.a libmct.a
endif
libmpiser_a_SOURCES = mpi.c send.c recv.c collective.c req.c list.c \
	                    handles.c comm.c group.c time.c pack.c fort.F90
libmpeu_a_SOURCES = get_zeits.c m_chars.F90 m_die.F90 m_dropdead.F90      \
										m_FcComms.F90 m_Filename.F90 m_FileResolv.F90         \
										m_flow.F90 m_IndexBin_char.F90 m_IndexBin_integer.F90 \
                    m_IndexBin_logical.F90 m_inpak90.F90 m_ioutil.F90     \
										m_List.F90 m_mall.F90 m_MergeSorts.F90 m_mpif90.F90   \
										m_mpif.F90 m_mpout.F90 m_Permuter.F90 m_rankMerge.F90 \
										m_realkinds.F90 m_SortingTools.F90 m_stdio.F90        \
										m_String.F90 m_StringLinkedList.F90 m_StrTemplate.F90 \
										m_TraceBack.F90 m_zeit.F90
libmct_a_SOURCES = m_AccumulatorComms.F90 m_Accumulator.F90              \
									 m_AttrVectComms.F90 m_AttrVect.F90                    \
									 m_AttrVectReduce.F90 m_ConvertMaps.F90                \
									 m_ExchangeMaps.F90 m_GeneralGridComms.F90             \
                   m_GeneralGrid.F90 m_GlobalMap.F90                     \
									 m_GlobalSegMapComms.F90 m_GlobalSegMap.F90            \
                   m_GlobalToLocal.F90 m_MatAttrVectMul.F90              \
                   m_MCTWorld.F90 m_Merge.F90 m_Navigator.F90            \
									 m_Rearranger.F90 m_Router.F90 m_SparseMatrixComms.F90 \
                   m_SparseMatrixDecomp.F90 m_SparseMatrix.F90           \
                   m_SparseMatrixPlus.F90 m_SparseMatrixToMaps.F90       \
									 m_SpatialIntegral.F90 m_SpatialIntegralV.F90          \
									 m_Transfer.F90

MPIFH= mpif.$(FORT_SIZE).h

mpif.h:
	cp -f $(REGCM_ROOT)/external/MCT/mpi-serial/$(MPIFH) mpif.h

%.o: %.c
	$(CC) $(CFLAGS) $(MPEUCFLAGS) -c $<

%.o: %.f90
	  $(FC) $(CPPFLAGS) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

%.o: %.F90
	  $(FC) $(CPPFLAGS) $(SVNDEF) $(FCFLAGS) $(MPEUFCFLAGS) -c $<

clean:
	    @echo -n Removing objects, mods and $(noinst_LIBRARIES) ... 
			@rm -f *.o *.mod $(noinst_LIBRARIES) mpif.h
			@echo Done

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.

mpi.o: mpif.h
m_IndexBin_char.o: m_die.o m_stdio.o
m_IndexBin_integer.o: m_die.o m_stdio.o
m_IndexBin_logical.o: m_die.o m_stdio.o
m_List.o: m_String.o m_die.o m_mall.o
m_MergeSorts.o: m_die.o m_realkinds.o m_stdio.o
m_Filename.o:
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
m_AttrVect.o: 
m_Accumulator.o: m_AttrVect.o
m_GlobalMap.o:
m_GlobalSegMap.o: m_MCTWorld.o
m_GlobalSegMapComms.o: m_GlobalSegMap.o
m_Navigator.o:
m_AttrVectComms.o: m_AttrVect.o m_GlobalMap.o m_GlobalSegMap.o
m_AttrVectReduce.o: m_AttrVect.o
m_AccumulatorComms.o: m_AttrVect.o m_GlobalMap.o m_AttrVectComms.o m_Accumulator.o
m_SparseMatrix.o: m_AttrVect.o m_GlobalMap.o m_AttrVectComms.o
m_GeneralGrid.o: m_AttrVect.o
m_GeneralGridComms.o: m_AttrVect.o m_GeneralGrid.o m_AttrVectComms.o m_GlobalMap.o m_GlobalSegMap.o
m_MatAttrVectMul.o: m_AttrVect.o m_SparseMatrix.o m_GlobalMap.o m_GlobalSegMap.o m_SparseMatrixPlus.o m_Rearranger.o
m_Merge.o: m_AttrVect.o m_GeneralGrid.o
m_Router.o: m_GlobalToLocal.o m_MCTWorld.o m_GlobalSegMap.o m_ExchangeMaps.o
m_Rearranger.o: m_Router.o m_MCTWorld.o m_GlobalSegMap.o m_AttrVect.o
m_GlobalToLocal.o: m_GlobalSegMap.o m_Navigator.o m_SparseMatrix.o
m_ExchangeMaps.o: m_GlobalMap.o m_GlobalSegMap.o m_MCTWorld.o m_ConvertMaps.o
m_ConvertMaps.o: m_GlobalMap.o m_GlobalSegMap.o m_MCTWorld.o
m_SparseMatrixDecomp.o: m_SparseMatrix.o m_GlobalSegMap.o
m_SparseMatrixToMaps.o: m_SparseMatrix.o m_GlobalSegMap.o
m_SparseMatrixComms.o:	m_SparseMatrix.o m_SparseMatrixDecomp.o m_GlobalSegMap.o m_AttrVectComms.o
accumulate.o: m_AttrVect.o m_Accumulator.o
m_SpatialIntegral.o: m_SpatialIntegralV.o m_GeneralGrid.o m_AttrVect.o m_AttrVectReduce.o
m_SpatialIntegralV.o: m_AttrVect.o m_AttrVectReduce.o
m_Transfer.o: m_AttrVect.o m_Router.o m_MCTWorld.o
m_SparseMatrixPlus.o: m_GlobalSegMap.o m_Rearranger.o m_SparseMatrix.o m_SparseMatrixComms.o m_SparseMatrixToMaps.o m_GlobalToLocal.o
