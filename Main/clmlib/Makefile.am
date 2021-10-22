include $(REGCM_ROOT)/makeinc

if DO_COMPILE_LIBCLM45
SUBDIRS = clm4.5
endif

REGCMSRC = $(REGCM_ROOT)/Main
BATS     = $(REGCMSRC)/batslib
BSP      = $(REGCMSRC)/clmlib
CLMSRC   = $(BSP)/clm3.5
MPP      = $(REGCMSRC)/mpplib

MPP_INCLUDE  = $(MPP)
BATS_INCLUDE = $(BATS)
ESMF_INCLUDE = $(CLMSRC)/src/utils/esmf_wrf_timemgr
UTIL_INCLUDE = $(CLMSRC)/src/utils/timing/
GLOB_INCLUDE = .

vpath %.o $(REGCM_ROOT)/Share:$(REGCM_ROOT)/Main/mpplib:\
          $(REGCM_ROOT)/Main/batslib:.

vpath %.F90 $(BSP):$(CLMSRC)/src/csm_share/shr:$(CLMSRC)/src/csm_share/eshr:\
            $(CLMSRC)/src/utils/timing:$(CLMSRC)/src/utils/mct/mpeu:\
            $(CLMSRC)/src/utils/mct/mct:$(CLMSRC)/src/utils/esmf_wrf_timemgr:\
            $(CLMSRC)/src/main:$(CLMSRC)/src/biogeophys:\
            $(CLMSRC)/src/biogeochem:$(CLMSRC)/src/riverroute:\
            $(REGCM_ROOT)/Config:$(REGCM_ROOT)/Main

vpath %.f90 $(BSP):$(CLMSRC)/src/csm_share/shr:$(CLMSRC)/src/csm_share/eshr:\
            $(CLMSRC)/src/utils/timing:$(CLMSRC)/src/utils/mct/mpeu:\
            $(CLMSRC)/src/utils/mct/mct:$(CLMSRC)/src/utils/esmf_wrf_timemgr:\
            $(CLMSRC)/src/main:$(CLMSRC)/src/biogeophys:\
            $(CLMSRC)/src/biogeochem:$(CLMSRC)/src/riverroute:\
            $(REGCM_ROOT)/Config:$(REGCM_ROOT)/Main

vpath %.inc $(BSP):$(CLMSRC)/src/csm_share/shr:$(CLMSRC)/src/csm_share/eshr:\
            $(CLMSRC)/src/utils/timing:$(CLMSRC)/src/utils/mct/mpeu:\
            $(CLMSRC)/src/utils/mct/mct:$(CLMSRC)/src/utils/esmf_wrf_timemgr:\
            $(CLMSRC)/src/main:$(CLMSRC)/src/biogeophys:\
            $(CLMSRC)/src/biogeochem:$(CLMSRC)/src/riverroute

vpath %.c $(BSP):$(CLMSRC)/src/csm_share/shr:$(CLMSRC)/src/csm_share/eshr:\
          $(CLMSRC)/src/utils/timing:$(CLMSRC)/src/utils/mct/mpeu:\
          $(CLMSRC)/src/utils/mct/mct:$(CLMSRC)/src/utils/esmf_wrf_timemgr:\
          $(CLMSRC)/src/main:$(CLMSRC)/src/biogeophys:\
          $(CLMSRC)/src/biogeochem:$(CLMSRC)/src/riverroute

vpath %.h $(BSP):$(CLMSRC)/src/csm_share/shr:$(CLMSRC)/src/csm_share/eshr:\
          $(CLMSRC)/src/utils/timing:$(CLMSRC)/src/utils/mct/mpeu:\
          $(CLMSRC)/src/utils/mct/mct:$(CLMSRC)/src/utils/esmf_wrf_timemgr:\
          $(CLMSRC)/src/main:$(CLMSRC)/src/biogeophys:\
          $(CLMSRC)/src/biogeochem:$(CLMSRC)/src/riverroute

noinst_LIBRARIES = libclm.a

if DO_COMPILE_LIBCLM
libclm_a_SOURCES =  mod_clm.F90 BalanceCheckMod.F90 BiogeophysRestMod.F90 \
Biogeophysics1Mod.F90 Biogeophysics2Mod.F90 \
C13SummaryMod.F90 CASAMod.F90 CASAPhenologyMod.F90 CNAllocationMod.F90 \
CNAnnualUpdateMod.F90 CNBalanceCheckMod.F90 CNC13FluxMod.F90 \
CNC13StateUpdate1Mod.F90 CNC13StateUpdate2Mod.F90 CNC13StateUpdate3Mod.F90 \
CNCStateUpdate1Mod.F90 CNCStateUpdate2Mod.F90 CNCStateUpdate3Mod.F90 \
CNDecompMod.F90 CNEcosystemDynMod.F90 CNFireMod.F90 CNGRespMod.F90 \
CNGapMortalityMod.F90 CNMRespMod.F90 CNNDynamicsMod.F90 CNNStateUpdate1Mod.F90 \
CNNStateUpdate2Mod.F90 CNNStateUpdate3Mod.F90 CNPhenologyMod.F90 \
CNPrecisionControlMod.F90 CNSetValueMod.F90 CNSummaryMod.F90 \
CNVegStructUpdateMod.F90 CNiniSpecial.F90 CNiniTimeVar.F90 CNrestMod.F90 \
DGVMAllocationMod.F90 DGVMEcosystemDynMod.F90 \
DGVMEstablishmentMod.F90 DGVMFireMod.F90 DGVMKillMod.F90 DGVMLightMod.F90 \
DGVMMod.F90 DGVMMortalityMod.F90 DGVMReproductionMod.F90 DGVMRestMod.F90 \
DGVMTurnoverMod.F90 DUSTMod.F90 DriverInitMod.F90 ESMF_AlarmClockMod.F90 \
ESMF_AlarmMod.F90 ESMF_BaseMod.F90 ESMF_BaseTimeMod.F90 ESMF_CalendarMod.F90 \
ESMF_ClockMod.F90 ESMF_FractionMod.F90 ESMF_Mod.F90 ESMF_Stubs.F90 \
ESMF_TimeIntervalMod.F90 ESMF_TimeMod.F90 FracWetMod.F90 VOCRestMod.F90 \
FrictionVelocityMod.F90 Hydrology1Mod.F90 Hydrology2Mod.F90 \
HydrologyLakeMod.F90 Meat.F90 QSatMod.F90 RtmMod.F90 RunoffMod.F90 \
SnowHydrologyMod.F90 SoilHydrologyMod.F90 SoilTemperatureMod.F90 \
SurfaceRadiationMod.F90 TridiagonalMod.F90 abortutils.F90 \
accumulMod.F90 areaMod.F90 clm_csmMod.F90 clm_mct_mod.F90 \
clm_varcon.F90 clm_varorb.F90 do_close_dispose.F90 \
domainMod.F90 eshr_estate_mod.F90 eshr_inputinfo_mod.F90 \
eshr_rc_mod.F90 eshr_timemgr_mod.F90 fileutils.F90 \
filterMod.F90 get_zeits.c getdatetime.F90 gptl.c gptl_papi.c gptlutil.c \
clm_varvoc.o inicFileMod.F90 initGridCellsMod.F90 initSurfAlbMod.F90 \
iobinary.F90 lnd_comp_mct.F90 m_Accumulator.F90 m_AccumulatorComms.F90 \
m_AttrVect.F90 m_AttrVectComms.F90  m_AttrVectReduce.F90 m_ConvertMaps.F90 \
m_ExchangeMaps.F90 m_FileResolv.F90 m_Filename.F90 m_GeneralGrid.F90 \
m_GeneralGridComms.F90 m_GlobalMap.F90 m_GlobalSegMap.F90 \
m_GlobalSegMapComms.F90 m_GlobalToLocal.F90 m_IndexBin_char.F90 \
m_IndexBin_integer.F90 m_IndexBin_logical.F90 m_List.F90 m_MCTWorld.F90 \
m_MatAttrVectMul.F90 m_Merge.F90 m_MergeSorts.F90 m_Navigator.F90 \
m_Permuter.F90 m_Rearranger.F90 m_Router.F90 m_SortingTools.F90 \
m_SparseMatrix.F90 m_SparseMatrixComms.F90 m_SparseMatrixDecomp.F90 \
m_SparseMatrixPlus.F90 m_SparseMatrixToMaps.F90 m_SpatialIntegral.F90 \
m_SpatialIntegralV.F90 m_StrTemplate.F90 m_String.F90 m_TraceBack.F90 \
m_Transfer.F90 m_chars.F90 m_die.F90 m_dropdead.F90 m_flow.F90 \
m_inpak90.F90 m_ioutil.F90 m_mall.F90 m_mpif.F90 m_mpif90.F90 \
m_mpout.F90 m_rankMerge.F90 m_realkinds.F90 m_stdio.F90 m_zeit.F90 \
mct_mod.F90 nanMod.F90 ndepFileMod.F90 pft2colMod.F90 pftdynMod.F90 \
pftvarcon.F90 program_csm.F90 program_offMod.F90 rtcrate.F90 \
scam_setlatlonidx.F90 shr_alarm_mod.F90 shr_cal_mod.F90 shr_const_mod.F90 \
shr_date_mod.F90 shr_file_mod.F90 shr_flux_mod.F90 shr_inputinfo_mod.F90 \
shr_jlcp.c shr_kind_mod.F90 shr_map_mod.F90 shr_mpi_mod.F90 shr_msg_mod.F90 \
shr_ncio_mod.F90 shr_ncread_mod.F90 shr_stream_mod.F90 \
shr_string_mod.F90 shr_tInterp_mod.F90 shr_timer_mod.F90 \
shr_vmath_fwrap.c shr_vmath_mod.F90 snowdp2lev.F90 spmdGathScatMod.F90 \
spmdMod.F90 subgridAveMod.F90 subgridRestMod.F90 \
system_messages.F90 threadutil.c wrf_error_fatal.F90 \
wrf_message.F90 clm_varvoc.F90 program_offMod.F90 VOCRestMod.F90 \
BalanceCheckMod.F90 initGridCellsMod.F90 clm_time_manager.F90 \
abortutils.F90 Meat.F90 wrf_error_fatal.F90 wrf_message.F90 \
scam_setlatlonidx.F90 gptl.c gptl_papi.c gptlutil.c f_wrappers.c \
threadutil.c getdatetime.F90 snowdp2lev.F90 do_close_dispose.F90 \
accFldsMod.F90 atmdrvMod.F90 BareGroundFluxesMod.F90 BiogeophysicsLakeMod.F90 \
CanopyFluxesMod.F90 clm_atmlnd.F90 clm_comp.F90 clmtype.F90 clmtypeInitMod.F90 \
clm_varctl.F90 clm_varpar.F90 clm_varsur.F90 controlMod.F90 \
decompMod.F90 driver.F90 histFileMod.F90 histFldsMod.F90 initializeMod.F90 \
iniTimeConst.F90 mkarbinitMod.F90 perf_mod.F90 restFileMod.F90 shr_orb_mod.F90 \
STATICEcosysDynMod.F90 subgridMod.F90 SurfaceAlbedoMod.F90 surfrdMod.F90 \
VOCEmissionMod.F90  ncdio.F90 shr_sys_mod.F90 mod_mtrxclm.F90 DryDepVelocityMod.F90 clm_drydep.F90
else
libclm_a_SOURCES = mod_clm_void.F90
endif

# Set up special characters
null  :=

FPPFLAGS = $(CPPFLAGS) $(DEFINE)MAXPATCH_PFT=4 $(DEFINE)OFFLINE

if COMPILER_IBM
if OSNAME_AIX
  FPPFLAGS += $(DEFINE)AIX $(DEFINE)FORTRAN_SAME
  DEFFLAGS = -DAIX -DFORTRAN_SAME
else
  FPPFLAGS += $(DEFINE)LINUX $(DEFINE)FORTRAN_SAME
  DEFFLAGS = -DFORTRAN_SAME -DLINUX
endif
else
  FPPFLAGS += $(DEFINE)LINUX $(DEFINE)FORTRANUNDERSCORE
  DEFFLAGS = -DLINUX -DFORTRANUNDERSCORE
endif

if COMPILER_GNU
FCFLAGS += -fno-range-check
FPPFLAGS += $(DEFINE)G95
DEFFLAGS += -DFORTRANUNDERSCORE -DLINUX
endif

if COMPILER_SUN
  FCFLAGS += -M$(MPP_INCLUDE) -M$(ESMF_INCLUDE) -M$(UTIL_INCLUDE) \
             -M$(GLOB_INCLUDE) -M$(BATS_INCLUDE) $(NETCDFINC)
else
  FCFLAGS += -I$(MPP_INCLUDE) -I$(ESMF_INCLUDE) -I$(UTIL_INCLUDE) \
             -I$(GLOB_INCLUDE) -I$(BATS_INCLUDE) $(NETCDFINC)
endif

%.o: %.f90
	$(MPIFC) $(FPPFLAGS) $(FCFLAGS) -c $<

%.o: %.F90
	$(MPIFC) $(FPPFLAGS) $(GITVER) $(FCFLAGS) -c $<

%.o: %.c
	$(CC) -o $@ $(CFLAGS) $(DEFFLAGS) -c $<

distclean-local:
	rm -f *.mod

clean-local:
	rm -f *.mod

# --------------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend depends on it.
#

BalanceCheckMod.o:  BalanceCheckMod.F90 shr_kind_mod.o subgridAveMod.o clmtype.o clm_varpar.o clm_time_manager.o clm_atmlnd.o abortutils.o misc.h preproc.h
BareGroundFluxesMod.o:  BareGroundFluxesMod.F90 shr_kind_mod.o clm_varcon.o clmtype.o shr_const_mod.o clm_varpar.o clm_atmlnd.o FrictionVelocityMod.o misc.h preproc.h
BiogeophysRestMod.o:  BiogeophysRestMod.F90 shr_kind_mod.o clm_varcon.o clmtype.o decompMod.o clm_varpar.o clm_varctl.o abortutils.o initSurfAlbMod.o ncdio.o misc.h preproc.h
Biogeophysics1Mod.o:  Biogeophysics1Mod.F90 shr_kind_mod.o clm_varcon.o clmtype.o clm_varpar.o clm_atmlnd.o QSatMod.o misc.h preproc.h
Biogeophysics2Mod.o:  Biogeophysics2Mod.F90 shr_kind_mod.o subgridAveMod.o clm_varcon.o clmtype.o clm_varpar.o SoilTemperatureMod.o clm_time_manager.o clm_atmlnd.o misc.h preproc.h
BiogeophysicsLakeMod.o:  BiogeophysicsLakeMod.F90 shr_kind_mod.o clm_varcon.o clmtype.o clm_varpar.o TridiagonalMod.o clm_time_manager.o clm_atmlnd.o QSatMod.o FrictionVelocityMod.o misc.h preproc.h
C13SummaryMod.o:  C13SummaryMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o pft2colMod.o misc.h preproc.h
CASAMod.o:  CASAMod.F90 shr_kind_mod.o subgridAveMod.o spmdGathScatMod.o spmdMod.o clm_varcon.o clmtype.o shr_const_mod.o decompMod.o clm_varpar.o domainMod.o clm_time_manager.o clm_varctl.o clm_atmlnd.o abortutils.o fileutils.o ncdio.o misc.h preproc.h
CASAPhenologyMod.o:  CASAPhenologyMod.F90 shr_kind_mod.o clm_varcon.o clmtype.o decompMod.o clm_varpar.o clm_time_manager.o clm_varctl.o misc.h preproc.h
CNAllocationMod.o:  CNAllocationMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_time_manager.o clm_varctl.o pft2colMod.o misc.h preproc.h
CNAnnualUpdateMod.o:  CNAnnualUpdateMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_time_manager.o clm_varctl.o pft2colMod.o misc.h preproc.h
CNBalanceCheckMod.o:  CNBalanceCheckMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_varctl.o clm_time_manager.o abortutils.o misc.h preproc.h
CNC13FluxMod.o:  CNC13FluxMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_varctl.o clm_time_manager.o pft2colMod.o misc.h preproc.h
CNC13StateUpdate1Mod.o:  CNC13StateUpdate1Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_varctl.o clm_time_manager.o misc.h preproc.h
CNC13StateUpdate2Mod.o:  CNC13StateUpdate2Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_varctl.o clm_time_manager.o misc.h preproc.h
CNC13StateUpdate3Mod.o:  CNC13StateUpdate3Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_varctl.o clm_time_manager.o misc.h preproc.h
CNCStateUpdate1Mod.o:  CNCStateUpdate1Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_varctl.o clm_time_manager.o misc.h preproc.h
CNCStateUpdate2Mod.o:  CNCStateUpdate2Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_varctl.o clm_time_manager.o misc.h preproc.h
CNCStateUpdate3Mod.o:  CNCStateUpdate3Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_varctl.o clm_time_manager.o misc.h preproc.h
CNDecompMod.o:  CNDecompMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o shr_const_mod.o clm_varpar.o clm_time_manager.o clm_varctl.o CNAllocationMod.o pft2colMod.o misc.h preproc.h
CNEcosystemDynMod.o:  CNEcosystemDynMod.F90 CNSummaryMod.o CNAnnualUpdateMod.o CNNStateUpdate3Mod.o shr_kind_mod.o CNC13StateUpdate3Mod.o CNNStateUpdate2Mod.o spmdMod.o CNCStateUpdate1Mod.o CNPrecisionControlMod.o CNDecompMod.o CNCStateUpdate3Mod.o CNNDynamicsMod.o clmtype.o CNPhenologyMod.o CNC13FluxMod.o CNSetValueMod.o CNC13StateUpdate1Mod.o CNBalanceCheckMod.o CNNStateUpdate1Mod.o C13SummaryMod.o CNFireMod.o CNVegStructUpdateMod.o CNCStateUpdate2Mod.o CNC13StateUpdate2Mod.o CNGapMortalityMod.o CNGRespMod.o CNMRespMod.o misc.h preproc.h
CNFireMod.o:  CNFireMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o shr_const_mod.o clm_varpar.o clm_varctl.o clm_time_manager.o pft2colMod.o misc.h preproc.h
CNGRespMod.o:  CNGRespMod.F90 shr_kind_mod.o clmtype.o clm_varcon.o misc.h preproc.h
CNGapMortalityMod.o:  CNGapMortalityMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o misc.h preproc.h
CNMRespMod.o:  CNMRespMod.F90 shr_kind_mod.o clmtype.o clm_varcon.o shr_const_mod.o clm_varpar.o clm_time_manager.o clm_varctl.o misc.h preproc.h
CNNDynamicsMod.o:  CNNDynamicsMod.F90 shr_kind_mod.o clmtype.o clm_varcon.o clm_varpar.o clm_time_manager.o clm_varctl.o clm_atmlnd.o misc.h preproc.h
CNNStateUpdate1Mod.o:  CNNStateUpdate1Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_time_manager.o clm_varctl.o misc.h preproc.h
CNNStateUpdate2Mod.o:  CNNStateUpdate2Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_time_manager.o clm_varctl.o misc.h preproc.h
CNNStateUpdate3Mod.o:  CNNStateUpdate3Mod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o clm_time_manager.o clm_varctl.o misc.h preproc.h
CNPhenologyMod.o:  CNPhenologyMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o shr_const_mod.o clm_varpar.o clm_time_manager.o clm_varctl.o misc.h preproc.h
CNPrecisionControlMod.o:  CNPrecisionControlMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o misc.h preproc.h
CNSetValueMod.o:  CNSetValueMod.F90 shr_kind_mod.o clmtype.o clm_varpar.o misc.h preproc.h
CNSummaryMod.o:  CNSummaryMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o pft2colMod.o misc.h preproc.h
CNVegStructUpdateMod.o:  CNVegStructUpdateMod.F90 shr_kind_mod.o spmdMod.o clmtype.o clm_varcon.o shr_const_mod.o clm_varpar.o clm_atmlnd.o pftvarcon.o misc.h preproc.h
CNiniSpecial.o:  CNiniSpecial.F90 shr_kind_mod.o clmtype.o clm_varcon.o CNSetValueMod.o decompMod.o pftvarcon.o misc.h preproc.h
CNiniTimeVar.o:  CNiniTimeVar.F90 shr_kind_mod.o clm_varcon.o clmtype.o decompMod.o clm_atmlnd.o pftvarcon.o misc.h preproc.h
CNrestMod.o:  CNrestMod.F90 shr_kind_mod.o spmdMod.o clmtype.o decompMod.o clm_varpar.o clm_varctl.o clm_atmlnd.o ncdio.o misc.h preproc.h
CanopyFluxesMod.o:  CanopyFluxesMod.F90 shr_kind_mod.o spmdMod.o clm_varcon.o clmtype.o shr_const_mod.o clm_varpar.o clm_time_manager.o clm_atmlnd.o abortutils.o QSatMod.o FrictionVelocityMod.o misc.h preproc.h
DGVMAllocationMod.o:  DGVMAllocationMod.F90 shr_kind_mod.o clmtype.o shr_const_mod.o misc.h preproc.h
DGVMEcosystemDynMod.o:  DGVMEcosystemDynMod.F90 shr_kind_mod.o clmtype.o shr_const_mod.o decompMod.o clm_varpar.o clm_time_manager.o pftvarcon.o nanMod.o misc.h preproc.h
DGVMEstablishmentMod.o:  DGVMEstablishmentMod.F90 shr_kind_mod.o clm_varcon.o clmtype.o shr_const_mod.o clm_varpar.o abortutils.o pftvarcon.o misc.h preproc.h
DGVMFireMod.o:  DGVMFireMod.F90 shr_kind_mod.o clmtype.o misc.h preproc.h
DGVMKillMod.o:  DGVMKillMod.F90 shr_kind_mod.o clmtype.o misc.h preproc.h
DGVMLightMod.o:  DGVMLightMod.F90 shr_kind_mod.o clmtype.o misc.h preproc.h
DGVMMod.o:  DGVMMod.F90 shr_kind_mod.o spmdGathScatMod.o spmdMod.o SurfaceAlbedoMod.o DGVMMortalityMod.o DGVMAllocationMod.o DGVMEcosystemDynMod.o DGVMKillMod.o DGVMEstablishmentMod.o DGVMReproductionMod.o clm_varcon.o clmtype.o shr_const_mod.o DGVMLightMod.o decompMod.o clm_varpar.o shr_sys_mod.o domainMod.o DGVMTurnoverMod.o clm_time_manager.o clm_varctl.o abortutils.o pftvarcon.o fileutils.o DGVMFireMod.o ncdio.o misc.h preproc.h
DGVMMortalityMod.o:  DGVMMortalityMod.F90 shr_kind_mod.o clmtype.o misc.h preproc.h
DGVMReproductionMod.o:  DGVMReproductionMod.F90 shr_kind_mod.o clmtype.o misc.h preproc.h
DGVMRestMod.o:  DGVMRestMod.F90 shr_kind_mod.o spmdMod.o clmtype.o decompMod.o clm_varctl.o abortutils.o ncdio.o misc.h preproc.h
DGVMTurnoverMod.o:  DGVMTurnoverMod.F90 shr_kind_mod.o clmtype.o misc.h preproc.h
DUSTMod.o:  DUSTMod.F90 subgridAveMod.o clm_varcon.o clmtype.o clm_varpar.o abortutils.o shr_kind_mod.o shr_const_mod.o decompMod.o clm_atmlnd.o misc.h preproc.h
DriverInitMod.o:  DriverInitMod.F90 subgridAveMod.o clmtype.o clm_varpar.o shr_kind_mod.o misc.h preproc.h
ESMF_AlarmClockMod.o:  ESMF_AlarmClockMod.F90 ESMF_TimeIntervalMod.o ESMF_TimeMod.o ESMF_AlarmMod.o ESMF_ClockMod.o ESMF_Macros.inc ESMF_TimeMgr.inc
ESMF_AlarmMod.o:  ESMF_AlarmMod.F90 ESMF_TimeIntervalMod.o ESMF_TimeMod.o ESMF_BaseMod.o ESMF_Macros.inc ESMF_TimeMgr.inc
ESMF_BaseMod.o:  ESMF_BaseMod.F90  
ESMF_BaseTimeMod.o:  ESMF_BaseTimeMod.F90 ESMF_BaseMod.o ESMF_Macros.inc ESMF_TimeMgr.inc
ESMF_CalendarMod.o:  ESMF_CalendarMod.F90 ESMF_BaseMod.o ESMF_BaseTimeMod.o ESMF_Macros.inc ESMF_TimeMgr.inc ESMF.inc
ESMF_ClockMod.o:  ESMF_ClockMod.F90 ESMF_TimeIntervalMod.o ESMF_TimeMod.o ESMF_BaseMod.o ESMF_AlarmMod.o ESMF_Macros.inc ESMF_TimeMgr.inc
ESMF_FractionMod.o:  ESMF_FractionMod.F90  
ESMF_Mod.o:  ESMF_Mod.F90 ESMF_TimeIntervalMod.o ESMF_FractionMod.o ESMF_TimeMod.o ESMF_BaseMod.o ESMF_AlarmMod.o ESMF_CalendarMod.o ESMF_AlarmClockMod.o ESMF_ClockMod.o ESMF_Stubs.o ESMF_BaseTimeMod.o ESMF_Macros.inc ESMF_TimeMgr.inc
ESMF_Stubs.o:  ESMF_Stubs.F90 ESMF_TimeMod.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF.inc
ESMF_TimeIntervalMod.o:  ESMF_TimeIntervalMod.F90 ESMF_FractionMod.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_BaseTimeMod.o ESMF_Macros.inc ESMF_TimeMgr.inc
ESMF_TimeMod.o:  ESMF_TimeMod.F90 ESMF_TimeIntervalMod.o ESMF_BaseMod.o ESMF_CalendarMod.o ESMF_BaseTimeMod.o ESMF_Macros.inc ESMF_TimeMgr.inc
FracWetMod.o:  FracWetMod.F90 clmtype.o shr_kind_mod.o misc.h preproc.h
FrictionVelocityMod.o:  FrictionVelocityMod.F90 clm_varcon.o clmtype.o shr_kind_mod.o shr_const_mod.o clm_atmlnd.o misc.h preproc.h
Hydrology1Mod.o:  Hydrology1Mod.F90 subgridAveMod.o clm_varcon.o clmtype.o shr_kind_mod.o FracWetMod.o clm_time_manager.o clm_atmlnd.o misc.h preproc.h
Hydrology2Mod.o:  Hydrology2Mod.F90 clm_varcon.o clmtype.o clm_varpar.o shr_kind_mod.o SoilHydrologyMod.o clm_time_manager.o clm_atmlnd.o SnowHydrologyMod.o misc.h preproc.h
HydrologyLakeMod.o:  HydrologyLakeMod.F90 clm_varcon.o clmtype.o shr_kind_mod.o clm_time_manager.o clm_atmlnd.o misc.h preproc.h
Meat.o:  Meat.F90 ESMF_TimeIntervalMod.o ESMF_FractionMod.o ESMF_TimeMod.o ESMF_BaseMod.o ESMF_AlarmMod.o ESMF_CalendarMod.o ESMF_ClockMod.o ESMF_BaseTimeMod.o ESMF_Macros.inc ESMF_TimeMgr.inc ESMF.inc
QSatMod.o:  QSatMod.F90 shr_kind_mod.o shr_const_mod.o misc.h preproc.h
RtmMod.o:  RtmMod.F90 clm_mct_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o shr_sys_mod.o clm_varctl.o RunoffMod.o abortutils.o shr_kind_mod.o spmdGathScatMod.o perf_mod.o shr_const_mod.o decompMod.o domainMod.o clm_time_manager.o ncdio.o areaMod.o misc.h preproc.h
RunoffMod.o:  RunoffMod.F90 spmdMod.o clm_mct_mod.o clmtype.o abortutils.o shr_kind_mod.o misc.h preproc.h
STATICEcosysDynMod.o:  STATICEcosysDynMod.F90 spmdMod.o clmtype.o clm_varpar.o clm_varctl.o abortutils.o fileutils.o shr_kind_mod.o spmdGathScatMod.o decompMod.o clm_time_manager.o pftvarcon.o ncdio.o nanMod.o misc.h preproc.h
SnowHydrologyMod.o:  SnowHydrologyMod.F90 clm_varcon.o clmtype.o clm_varpar.o shr_kind_mod.o clm_time_manager.o misc.h preproc.h
SoilHydrologyMod.o:  SoilHydrologyMod.F90 clm_varcon.o clmtype.o clm_varpar.o shr_kind_mod.o shr_const_mod.o clm_time_manager.o TridiagonalMod.o misc.h preproc.h
SoilTemperatureMod.o:  SoilTemperatureMod.F90 clm_varcon.o clmtype.o clm_varpar.o shr_kind_mod.o TridiagonalMod.o clm_time_manager.o clm_atmlnd.o misc.h preproc.h
SurfaceAlbedoMod.o:  SurfaceAlbedoMod.F90 shr_orb_mod.o spmdMod.o clmtype.o clm_varcon.o clm_varpar.o mod_clm.o shr_kind_mod.o decompMod.o clm_time_manager.o misc.h preproc.h
SurfaceRadiationMod.o:  SurfaceRadiationMod.F90 clm_varcon.o clmtype.o clm_varpar.o shr_kind_mod.o clm_time_manager.o clm_atmlnd.o misc.h preproc.h
TridiagonalMod.o:  TridiagonalMod.F90 shr_kind_mod.o misc.h preproc.h
VOCEmissionMod.o:  VOCEmissionMod.F90 clm_varvoc.o clmtype.o clm_varpar.o clm_varctl.o abortutils.o shr_kind_mod.o shr_const_mod.o clm_time_manager.o clm_atmlnd.o misc.h preproc.h
VOCRestMod.o:  VOCRestMod.F90 clm_varvoc.o clm_varcon.o clmtype.o clm_varctl.o abortutils.o shr_kind_mod.o decompMod.o ncdio.o misc.h preproc.h
abortutils.o:  abortutils.F90 shr_sys_mod.o spmdMod.o misc.h preproc.h
accFldsMod.o:  accFldsMod.F90 shr_kind_mod.o clm_varctl.o clm_time_manager.o clm_atmlnd.o abortutils.o pftvarcon.o clm_varcon.o mod_clm.o accumulMod.o clmtype.o shr_const_mod.o nanMod.o decompMod.o misc.h preproc.h
accumulMod.o:  accumulMod.F90 shr_kind_mod.o spmdMod.o clm_varctl.o clm_time_manager.o abortutils.o clm_varcon.o ncdio.o nanMod.o shr_const_mod.o decompMod.o misc.h preproc.h
areaMod.o:  areaMod.F90 clm_varpar.o clm_varsur.o shr_kind_mod.o clm_mct_mod.o spmdMod.o domainMod.o abortutils.o clm_varcon.o nanMod.o shr_const_mod.o decompMod.o misc.h preproc.h
atmdrvMod.o:  atmdrvMod.F90 shr_kind_mod.o perf_mod.o clm_mct_mod.o spmdMod.o clm_varcon.o shr_const_mod.o decompMod.o clm_varpar.o domainMod.o clm_varctl.o clm_time_manager.o clm_atmlnd.o abortutils.o mod_clm.o fileutils.o ncdio.o areaMod.o surfrdMod.o nanMod.o misc.h preproc.h
clm_atmlnd.o:  clm_atmlnd.F90 clm_varsur.o shr_kind_mod.o subgridAveMod.o spmdMod.o clmtype.o clm_varcon.o decompMod.o clm_varpar.o domainMod.o abortutils.o mod_clm.o QSatMod.o areaMod.o nanMod.o misc.h preproc.h clm_drydep.o
clm_comp.o:  clm_comp.F90 shr_kind_mod.o shr_orb_mod.o eshr_timemgr_mod.o perf_mod.o driver.o initializeMod.o shr_sys_mod.o clm_varorb.o clm_time_manager.o clm_varctl.o clm_atmlnd.o shr_inputinfo_mod.o mod_clm.o initSurfAlbMod.o misc.h preproc.h
clm_csmMod.o:  clm_csmMod.F90 shr_kind_mod.o subgridAveMod.o perf_mod.o spmdGathScatMod.o spmdMod.o clmtype.o clm_varcon.o shr_const_mod.o decompMod.o clm_varpar.o shr_sys_mod.o domainMod.o clm_time_manager.o clm_varctl.o RunoffMod.o clm_atmlnd.o abortutils.o ncdio.o nanMod.o misc.h preproc.h
clm_mct_mod.o:  clm_mct_mod.F90 shr_kind_mod.o m_stdio.o m_AttrVect.o m_SparseMatrixPlus.o m_String.o m_SparseMatrixComms.o m_SparseMatrix.o m_List.o m_Rearranger.o m_GlobalSegMap.o shr_const_mod.o m_MatAttrVectMul.o m_GeneralGridComms.o shr_sys_mod.o m_AttrVectComms.o m_GlobalToLocal.o shr_mpi_mod.o m_GeneralGrid.o m_SparseMatrixToMaps.o m_MCTWorld.o m_die.o m_Permuter.o m_MergeSorts.o 
clm_time_manager.o:  clm_time_manager.F90 shr_kind_mod.o spmdMod.o ESMF_Mod.o clm_varctl.o abortutils.o ncdio.o misc.h preproc.h mod_date.o mod_mppparam.o mod_mpmessage.o clm_varsur.o
clm_varcon.o:  clm_varcon.F90 shr_kind_mod.o shr_const_mod.o clm_varpar.o misc.h preproc.h
clm_varctl.o:  clm_varctl.F90 shr_kind_mod.o misc.h preproc.h
clm_varorb.o:  clm_varorb.F90 shr_kind_mod.o misc.h preproc.h
clm_varpar.o:  clm_varpar.F90 shr_kind_mod.o preproc.h misc.h
clm_varsur.o:  clm_varsur.F90 shr_kind_mod.o misc.h preproc.h
clm_varvoc.o:  clm_varvoc.F90 shr_kind_mod.o clm_varpar.o misc.h preproc.h
clmtype.o:  clmtype.F90 shr_kind_mod.o clm_varpar.o domainMod.o misc.h preproc.h
clmtypeInitMod.o:  clmtypeInitMod.F90 shr_kind_mod.o CASAMod.o clm_varcon.o clm_varvoc.o clmtype.o decompMod.o clm_varpar.o nanMod.o misc.h preproc.h
controlMod.o:  controlMod.F90 clm_varsur.o shr_kind_mod.o CASAMod.o spmdMod.o shr_string_mod.o shr_const_mod.o decompMod.o restFileMod.o clm_varpar.o shr_sys_mod.o histFileMod.o clm_time_manager.o clm_varctl.o shr_inputinfo_mod.o abortutils.o mod_clm.o fileutils.o misc.h preproc.h
decompMod.o:  decompMod.F90 clm_varsur.o shr_kind_mod.o subgridMod.o clm_mct_mod.o spmdMod.o clmtype.o shr_sys_mod.o domainMod.o clm_varctl.o abortutils.o mod_clm.o misc.h preproc.h
do_close_dispose.o:  do_close_dispose.F90 shr_kind_mod.o clm_csmMod.o clm_time_manager.o misc.h preproc.h
domainMod.o:  domainMod.F90 shr_kind_mod.o spmdMod.o abortutils.o nanMod.o misc.h preproc.h
driver.o:  driver.F90 CASAPhenologyMod.o SurfaceAlbedoMod.o pftdynMod.o spmdMod.o clm_varcon.o clmtype.o Biogeophysics1Mod.o CanopyFluxesMod.o VOCEmissionMod.o DGVMMod.o Hydrology2Mod.o restFileMod.o shr_sys_mod.o clm_csmMod.o histFileMod.o RtmMod.o CNBalanceCheckMod.o clm_varctl.o abortutils.o DUSTMod.o mod_clm.o Hydrology1Mod.o BalanceCheckMod.o SurfaceRadiationMod.o shr_kind_mod.o STATICEcosysDynMod.o CASAMod.o perf_mod.o BiogeophysicsLakeMod.o DGVMEcosystemDynMod.o BareGroundFluxesMod.o decompMod.o filterMod.o accFldsMod.o inicFileMod.o DriverInitMod.o clm_time_manager.o HydrologyLakeMod.o clm_atmlnd.o pft2colMod.o ndepFileMod.o Biogeophysics2Mod.o CNEcosystemDynMod.o misc.h preproc.h DryDepVelocityMod.o
clm_drydep.o : clm_drydep.F90 shr_kind_mod.o spmdMod.o shr_const_mod.o shr_sys_mod.o
DryDepVelocityMod.o : DryDepVelocityMod.F90 shr_kind_mod.o spmdMod.o clm_varcon.o clmtype.o shr_const_mod.o clm_drydep.o clm_time_manager.o clm_atmlnd.o abortutils.o pftvarcon.o misc.h preproc.h
eshr_estate_mod.o:  eshr_estate_mod.F90 ESMF_Mod.o shr_string_mod.o shr_sys_mod.o shr_kind_mod.o shr_const_mod.o shr_mpi_mod.o eshr_rc_mod.o 
eshr_inputinfo_mod.o:  eshr_inputinfo_mod.F90 ESMF_Mod.o eshr_estate_mod.o shr_inputinfo_mod.o shr_kind_mod.o eshr_rc_mod.o 
eshr_rc_mod.o:  eshr_rc_mod.F90 ESMF_Mod.o shr_sys_mod.o 
eshr_timemgr_mod.o:  eshr_timemgr_mod.F90 shr_orb_mod.o ESMF_Mod.o shr_string_mod.o shr_sys_mod.o eshr_estate_mod.o shr_file_mod.o shr_kind_mod.o shr_mpi_mod.o eshr_rc_mod.o shr_ncio_mod.o 
f_wrappers.o:  f_wrappers.c  gptl.h private.h
fileutils.o:  fileutils.F90 shr_file_mod.o abortutils.o misc.h preproc.h
filterMod.o:  filterMod.F90 clm_varcon.o clmtype.o abortutils.o shr_kind_mod.o decompMod.o misc.h preproc.h
get_zeits.o:  get_zeits.c  
getdatetime.o:  getdatetime.F90  misc.h preproc.h
gptl.o:  gptl.c  gptl.h private.h
gptl_papi.o:  gptl_papi.c  gptl.h private.h
gptlutil.o:  gptlutil.c  gptl.h private.h
histFileMod.o:  histFileMod.F90 clm_varsur.o subgridAveMod.o spmdMod.o clmtype.o clm_varcon.o iobinary.o clm_varpar.o shr_sys_mod.o clm_varctl.o RunoffMod.o abortutils.o fileutils.o shr_kind_mod.o CASAMod.o spmdGathScatMod.o shr_const_mod.o decompMod.o domainMod.o clm_time_manager.o ncdio.o misc.h preproc.h clm_varsur.o mod_runparams.o mod_dynparam.o
histFldsMod.o:  histFldsMod.F90 clm_varcon.o clmtype.o histFileMod.o clm_varctl.o RunoffMod.o shr_kind_mod.o CASAMod.o clm_atmlnd.o misc.h preproc.h
iniTimeConst.o:  iniTimeConst.F90 clm_varsur.o spmdMod.o clm_varvoc.o clm_varcon.o clmtype.o clm_varpar.o clm_varctl.o abortutils.o mod_clm.o fileutils.o shr_kind_mod.o decompMod.o domainMod.o clm_time_manager.o clm_atmlnd.o pftvarcon.o ncdio.o ndepFileMod.o nanMod.o misc.h preproc.h
inicFileMod.o:  inicFileMod.F90 spmdMod.o clm_varcon.o clmtype.o clm_varpar.o shr_sys_mod.o clm_varctl.o abortutils.o fileutils.o shr_kind_mod.o decompMod.o ncdio.o misc.h preproc.h
initGridCellsMod.o:  initGridCellsMod.F90 clm_varsur.o spmdMod.o clm_varcon.o clmtype.o clm_varpar.o clm_varctl.o abortutils.o shr_kind_mod.o subgridMod.o spmdGathScatMod.o shr_const_mod.o decompMod.o domainMod.o misc.h preproc.h
initSurfAlbMod.o:  initSurfAlbMod.F90 SurfaceAlbedoMod.o spmdMod.o shr_orb_mod.o clm_varcon.o clmtype.o CNSetValueMod.o clm_varpar.o abortutils.o STATICEcosysDynMod.o shr_kind_mod.o DGVMEcosystemDynMod.o shr_const_mod.o FracWetMod.o decompMod.o filterMod.o clm_time_manager.o CNVegStructUpdateMod.o CNEcosystemDynMod.o misc.h preproc.h
initializeMod.o:  initializeMod.F90 clm_varsur.o CASAPhenologyMod.o pftdynMod.o spmdMod.o clm_varvoc.o initGridCellsMod.o DGVMMod.o restFileMod.o clm_varpar.o shr_sys_mod.o histFileMod.o RtmMod.o controlMod.o clm_varctl.o shr_inputinfo_mod.o DUSTMod.o abortutils.o mod_clm.o fileutils.o histFldsMod.o STATICEcosysDynMod.o shr_kind_mod.o CASAMod.o mkarbinitMod.o eshr_timemgr_mod.o DGVMEcosystemDynMod.o decompMod.o accFldsMod.o filterMod.o domainMod.o clm_time_manager.o clm_atmlnd.o clmtypeInitMod.o areaMod.o ndepFileMod.o surfrdMod.o misc.h preproc.h
iobinary.o:  iobinary.F90 spmdMod.o clmtype.o RunoffMod.o abortutils.o shr_kind_mod.o spmdGathScatMod.o decompMod.o misc.h preproc.h
lnd_comp_mct.o:  lnd_comp_mct.F90 shr_orb_mod.o spmdMod.o clm_varcon.o shr_sys_mod.o controlMod.o clm_varorb.o clm_varctl.o shr_inputinfo_mod.o shr_kind_mod.o eshr_timemgr_mod.o perf_mod.o clm_comp.o decompMod.o domainMod.o clm_time_manager.o clm_atmlnd.o mct_mod.o misc.h preproc.h
m_Accumulator.o:  m_Accumulator.F90 m_AttrVect.o m_String.o m_List.o m_mall.o m_die.o m_stdio.o m_realkinds.o 
m_AccumulatorComms.o:  m_AccumulatorComms.F90 m_AttrVect.o m_GlobalMap.o m_AttrVectComms.o m_die.o m_stdio.o m_GlobalSegMap.o m_Accumulator.o m_mpif90.o 
m_AttrVect.o:  m_AttrVect.F90 m_String.o m_SortingTools.o m_List.o m_mall.o m_die.o m_TraceBack.o m_stdio.o m_realkinds.o 
m_AttrVectComms.o:  m_AttrVectComms.F90 m_AttrVect.o m_String.o m_List.o m_GlobalMap.o m_die.o m_stdio.o m_GlobalSegMap.o m_mpif90.o 
m_AttrVectReduce.o:  m_AttrVectReduce.F90 m_AttrVect.o m_List.o m_die.o m_stdio.o m_mpif90.o m_realkinds.o 
m_ConvertMaps.o:  m_ConvertMaps.F90 m_SortingTools.o m_GlobalMap.o m_die.o m_stdio.o m_GlobalSegMap.o m_MCTWorld.o 
m_ExchangeMaps.o:  m_ExchangeMaps.F90 m_GlobalMap.o m_ConvertMaps.o m_die.o m_stdio.o m_GlobalSegMap.o m_MCTWorld.o m_mpif90.o 
m_FileResolv.o:  m_FileResolv.F90 m_StrTemplate.o m_die.o 
m_Filename.o:  m_Filename.F90  
m_GeneralGrid.o:  m_GeneralGrid.F90 m_AttrVect.o m_String.o m_List.o m_die.o m_TraceBack.o m_stdio.o m_realkinds.o 
m_GeneralGridComms.o:  m_GeneralGridComms.F90 m_List.o m_GlobalMap.o m_AttrVectComms.o m_GeneralGrid.o m_die.o m_stdio.o m_GlobalSegMap.o m_MCTWorld.o m_mpif90.o 
m_GlobalMap.o:  m_GlobalMap.F90 m_die.o m_stdio.o m_mpif90.o 
m_GlobalSegMap.o:  m_GlobalSegMap.F90 m_SortingTools.o m_die.o m_stdio.o m_mpif90.o 
m_GlobalSegMapComms.o:  m_GlobalSegMapComms.F90 m_die.o m_stdio.o m_GlobalSegMap.o m_MCTWorld.o m_mpif90.o 
m_GlobalToLocal.o:  m_GlobalToLocal.F90 m_GlobalMap.o m_die.o m_Navigator.o m_stdio.o m_SparseMatrix.o m_GlobalSegMap.o m_mpif90.o 
m_IndexBin_char.o:  m_IndexBin_char.F90 m_die.o m_stdio.o 
m_IndexBin_integer.o:  m_IndexBin_integer.F90 m_die.o m_stdio.o 
m_IndexBin_logical.o:  m_IndexBin_logical.F90 m_die.o m_stdio.o 
m_List.o:  m_List.F90 m_String.o m_mall.o m_die.o m_stdio.o m_mpif90.o 
m_MCTWorld.o:  m_MCTWorld.F90 m_List.o m_die.o m_stdio.o m_mpif90.o 
m_MatAttrVectMul.o:  m_MatAttrVectMul.F90 m_AttrVect.o m_SparseMatrixPlus.o m_String.o m_List.o m_die.o m_stdio.o m_SparseMatrix.o m_Rearranger.o m_mpif90.o m_realkinds.o 
m_Merge.o:  m_Merge.F90 m_AttrVect.o m_String.o m_List.o m_GeneralGrid.o m_die.o m_stdio.o m_realkinds.o 
m_MergeSorts.o:  m_MergeSorts.F90 m_die.o m_stdio.o m_realkinds.o 
m_Navigator.o:  m_Navigator.F90 m_mall.o m_die.o m_stdio.o 
m_Permuter.o:  m_Permuter.F90 m_die.o m_realkinds.o 
m_Rearranger.o:  m_Rearranger.F90 m_AttrVect.o m_Router.o m_die.o m_stdio.o m_GlobalSegMap.o m_MCTWorld.o m_mpif90.o m_realkinds.o 
m_Router.o:  m_Router.F90 m_GlobalToLocal.o m_ExchangeMaps.o m_die.o m_stdio.o m_GlobalSegMap.o m_MCTWorld.o m_mpif90.o m_realkinds.o 
m_SortingTools.o:  m_SortingTools.F90 m_rankMerge.o m_IndexBin_logical.o m_Permuter.o m_IndexBin_char.o m_MergeSorts.o m_IndexBin_integer.o 
m_SparseMatrix.o:  m_SparseMatrix.F90 m_AttrVect.o m_String.o m_List.o m_die.o m_TraceBack.o m_stdio.o m_mpif90.o m_realkinds.o 
m_SparseMatrixComms.o:  m_SparseMatrixComms.F90 m_List.o m_GlobalMap.o m_AttrVectComms.o m_die.o m_SparseMatrixDecomp.o m_stdio.o m_SparseMatrix.o m_GlobalSegMap.o m_mpif90.o 
m_SparseMatrixDecomp.o:  m_SparseMatrixDecomp.F90 m_AttrVect.o m_List.o m_GlobalMap.o m_AttrVectComms.o m_die.o m_SparseMatrix.o m_GlobalSegMap.o 
m_SparseMatrixPlus.o:  m_SparseMatrixPlus.F90 m_SparseMatrixComms.o m_String.o m_List.o m_GlobalToLocal.o m_SparseMatrixToMaps.o m_die.o m_stdio.o m_SparseMatrix.o m_GlobalSegMap.o m_Rearranger.o m_mpif90.o 
m_SparseMatrixToMaps.o:  m_SparseMatrixToMaps.F90 m_List.o m_die.o m_stdio.o m_SparseMatrix.o m_GlobalSegMap.o m_mpif90.o 
m_SpatialIntegral.o:  m_SpatialIntegral.F90 m_AttrVect.o m_SpatialIntegralV.o m_String.o m_List.o m_GeneralGrid.o m_die.o m_stdio.o m_AttrVectReduce.o m_mpif90.o m_realkinds.o 
m_SpatialIntegralV.o:  m_SpatialIntegralV.F90 m_AttrVect.o m_List.o m_die.o m_stdio.o m_AttrVectReduce.o m_mpif90.o m_realkinds.o 
m_StrTemplate.o:  m_StrTemplate.F90 m_die.o m_stdio.o m_chars.o 
m_String.o:  m_String.F90 m_mall.o m_die.o m_mpif90.o 
m_TraceBack.o:  m_TraceBack.F90 m_String.o m_die.o m_stdio.o 
m_Transfer.o:  m_Transfer.F90 m_AttrVect.o m_Router.o m_die.o m_stdio.o m_MCTWorld.o m_mpif90.o 
m_chars.o:  m_chars.F90  
m_die.o:  m_die.F90 m_flow.o m_dropdead.o m_stdio.o m_mpout.o m_mpif90.o 
m_dropdead.o:  m_dropdead.F90 m_stdio.o m_mpif90.o 
m_flow.o:  m_flow.F90 m_chars.o 
m_inpak90.o:  m_inpak90.F90 m_ioutil.o m_mall.o m_die.o m_stdio.o m_mpif90.o m_realkinds.o 
m_ioutil.o:  m_ioutil.F90 m_stdio.o 
m_mall.o:  m_mall.F90 m_ioutil.o m_die.o m_stdio.o m_chars.o m_realkinds.o 
m_mpif.o:  m_mpif.F90  
m_mpif90.o:  m_mpif90.F90 m_stdio.o m_mpif.o m_realkinds.o 
m_mpout.o:  m_mpout.F90 m_ioutil.o m_dropdead.o m_stdio.o m_mpif90.o 
m_rankMerge.o:  m_rankMerge.F90 m_realkinds.o 
m_realkinds.o:  m_realkinds.F90  
m_stdio.o:  m_stdio.F90  
m_zeit.o:  m_zeit.F90 m_SortingTools.o m_ioutil.o m_die.o m_stdio.o m_mpif90.o 
mct_mod.o:  mct_mod.F90 m_AttrVect.o m_SparseMatrixPlus.o m_String.o m_SparseMatrixComms.o m_List.o shr_sys_mod.o m_AttrVectComms.o m_GlobalToLocal.o m_Router.o m_GeneralGrid.o m_SparseMatrixToMaps.o m_die.o m_Transfer.o shr_kind_mod.o m_stdio.o m_SparseMatrix.o m_Rearranger.o m_GlobalSegMap.o m_inpak90.o shr_const_mod.o m_MatAttrVectMul.o m_GeneralGridComms.o shr_mpi_mod.o m_MCTWorld.o m_Permuter.o m_MergeSorts.o 
mkarbinitMod.o:  mkarbinitMod.F90 clm_varsur.o spmdMod.o clm_varcon.o clmtype.o clm_varpar.o shr_kind_mod.o shr_const_mod.o decompMod.o pftvarcon.o misc.h preproc.h mod_bats_param.o
mod_clm.o:  mod_clm.F90 mod_runparams.o mod_dynparam.o
mod_mtrxclm.o : mod_mtrxclm.F90 atmdrvMod.o clm_comp.o clm_drydep.o clm_time_manager.o clmtype.o clm_varpar.o clm_varsur.o initializeMod.o mod_clm.o mod_date.o mod_dynparam.o mod_intkinds.o mod_mpmessage.o mod_mppparam.o mod_realkinds.o mod_regcm_types.o mod_runparams.o mod_service.o mod_sunorbit.o perf_mod.o program_offMod.o shr_orb_mod.o mod_stdio.o mod_constants.o
nanMod.o:  nanMod.F90 shr_kind_mod.o 
ncdio.o:  ncdio.F90 clm_mct_mod.o spmdMod.o clm_varcon.o clm_varpar.o shr_sys_mod.o clm_varctl.o RunoffMod.o abortutils.o shr_kind_mod.o spmdGathScatMod.o decompMod.o nanMod.o misc.h preproc.h
ndepFileMod.o:  ndepFileMod.F90 spmdMod.o clm_varcon.o clmtype.o clm_varpar.o clm_varctl.o abortutils.o fileutils.o shr_kind_mod.o decompMod.o clm_time_manager.o clm_atmlnd.o ncdio.o misc.h preproc.h
perf_mod.o:  perf_mod.F90 shr_sys_mod.o shr_file_mod.o shr_kind_mod.o shr_mpi_mod.o gptl.inc
pft2colMod.o:  pft2colMod.F90 subgridAveMod.o clmtype.o shr_kind_mod.o misc.h preproc.h
pftdynMod.o:  pftdynMod.F90 clm_varsur.o spmdMod.o clm_varcon.o clmtype.o clm_varpar.o clm_varctl.o abortutils.o fileutils.o shr_kind_mod.o spmdGathScatMod.o decompMod.o clm_time_manager.o ncdio.o misc.h preproc.h
pftvarcon.o:  pftvarcon.F90 spmdMod.o clm_varpar.o clm_varctl.o abortutils.o fileutils.o shr_kind_mod.o misc.h preproc.h
program_csm.o:  program_csm.F90 spmdMod.o shr_orb_mod.o ESMF_Mod.o clm_csmMod.o controlMod.o clm_varorb.o clm_varctl.o shr_file_mod.o abortutils.o shr_kind_mod.o perf_mod.o clm_comp.o clm_time_manager.o misc.h preproc.h
program_offMod.o:  program_offMod.F90 spmdMod.o clm_mct_mod.o shr_orb_mod.o ESMF_Mod.o controlMod.o atmdrvMod.o abortutils.o mod_clm.o shr_kind_mod.o perf_mod.o clm_comp.o clm_time_manager.o misc.h preproc.h
restFileMod.o:  restFileMod.F90 spmdMod.o clm_varpar.o shr_sys_mod.o histFileMod.o clm_csmMod.o RtmMod.o clm_varctl.o RunoffMod.o abortutils.o fileutils.o subgridRestMod.o shr_kind_mod.o CASAMod.o decompMod.o CNrestMod.o DGVMRestMod.o clm_time_manager.o accumulMod.o BiogeophysRestMod.o ncdio.o misc.h preproc.h mod_dynparam.o
rtcrate.o:  rtcrate.F90  
scam_setlatlonidx.o:  scam_setlatlonidx.F90 abortutils.o shr_kind_mod.o 
shr_alarm_mod.o:  shr_alarm_mod.F90 shr_sys_mod.o shr_date_mod.o shr_kind_mod.o shr_cal_mod.o 
shr_cal_mod.o:  shr_cal_mod.F90 shr_sys_mod.o shr_kind_mod.o shr_const_mod.o 
shr_const_mod.o:  shr_const_mod.F90 shr_kind_mod.o 
shr_date_mod.o:  shr_date_mod.F90 shr_sys_mod.o shr_kind_mod.o shr_cal_mod.o 
shr_file_mod.o:  shr_file_mod.F90 shr_sys_mod.o shr_kind_mod.o 
shr_flux_mod.o:  shr_flux_mod.F90 shr_sys_mod.o shr_kind_mod.o shr_const_mod.o 
shr_inputinfo_mod.o:  shr_inputinfo_mod.F90 shr_string_mod.o shr_sys_mod.o shr_file_mod.o shr_kind_mod.o shr_mpi_mod.o shr_ncio_mod.o 
shr_jlcp.o:  shr_jlcp.c  
shr_kind_mod.o:  shr_kind_mod.F90  
shr_map_mod.o:  shr_map_mod.F90 shr_sys_mod.o shr_kind_mod.o shr_const_mod.o 
shr_mpi_mod.o:  shr_mpi_mod.F90 shr_kind_mod.o 
shr_msg_mod.o:  shr_msg_mod.F90 shr_sys_mod.o shr_file_mod.o 
shr_ncio_mod.o:  shr_ncio_mod.F90 shr_string_mod.o shr_sys_mod.o shr_kind_mod.o shr_const_mod.o shr_mpi_mod.o 
shr_ncread_mod.o:  shr_ncread_mod.F90 shr_string_mod.o shr_sys_mod.o shr_file_mod.o shr_kind_mod.o shr_cal_mod.o 
shr_orb_mod.o:  shr_orb_mod.F90 clm_varsur.o shr_sys_mod.o shr_kind_mod.o shr_const_mod.o 
shr_stream_mod.o:  shr_stream_mod.F90 shr_string_mod.o shr_sys_mod.o shr_kind_mod.o shr_const_mod.o shr_ncread_mod.o 
shr_string_mod.o:  shr_string_mod.F90 shr_sys_mod.o shr_kind_mod.o shr_timer_mod.o shr_cal_mod.o 
shr_sys_mod.o:  shr_sys_mod.F90 shr_kind_mod.o shr_mpi_mod.o 
shr_tInterp_mod.o:  shr_tInterp_mod.F90 shr_sys_mod.o shr_kind_mod.o shr_const_mod.o shr_cal_mod.o 
shr_timer_mod.o:  shr_timer_mod.F90 shr_sys_mod.o shr_kind_mod.o 
shr_vmath_fwrap.o:  shr_vmath_fwrap.c  
shr_vmath_mod.o:  shr_vmath_mod.F90 shr_kind_mod.o 
snowdp2lev.o:  snowdp2lev.F90 clmtype.o clm_varpar.o shr_kind_mod.o misc.h preproc.h
spmdGathScatMod.o:  spmdGathScatMod.F90 clm_mct_mod.o spmdMod.o clmtype.o RunoffMod.o abortutils.o shr_kind_mod.o decompMod.o misc.h preproc.h
spmdMod.o:  spmdMod.F90 shr_kind_mod.o misc.h preproc.h
subgridAveMod.o:  subgridAveMod.F90 clm_varcon.o clmtype.o clm_varpar.o abortutils.o shr_kind_mod.o misc.h preproc.h
subgridMod.o:  subgridMod.F90 clm_varsur.o spmdMod.o clmtype.o clm_varpar.o clm_varctl.o abortutils.o shr_kind_mod.o nanMod.o 
subgridRestMod.o:  subgridRestMod.F90 spmdMod.o clmtype.o shr_kind_mod.o decompMod.o clm_time_manager.o ncdio.o 
surfrdMod.o:  surfrdMod.F90 clm_varsur.o spmdMod.o clm_varvoc.o clm_varcon.o clmtype.o clm_varpar.o clm_varctl.o abortutils.o mod_clm.o fileutils.o shr_kind_mod.o decompMod.o domainMod.o clm_time_manager.o pftvarcon.o areaMod.o ncdio.o misc.h preproc.h
system_messages.o:  system_messages.F90 abortutils.o 
threadutil.o:  threadutil.c  gptl.h private.h
wrf_error_fatal.o:  wrf_error_fatal.F90 shr_sys_mod.o 
wrf_message.o:  wrf_message.F90  
