include ../../Makefile.inc

CLMSRC = ../clm3.5

VPATH = $(CLMSRC)/src/csm_share/shr:$(CLMSRC)/src/csm_share/eshr:\
        $(CLMSRC)/src/utils/timing:$(CLMSRC)/src/utils/mct/mpeu:\
        $(CLMSRC)/src/utils/mct/mct:$(CLMSRC)/src/utils/esmf_wrf_timemgr:\
        $(CLMSRC)/src/main:$(CLMSRC)/src/biogeophys:$(CLMSRC)/src/biogeochem:\
        $(CLMSRC)/src/riverroute

SRCS =  \
BalanceCheckMod.F90 BareGroundFluxesMod.F90 BiogeophysRestMod.F90 \
Biogeophysics1Mod.F90 Biogeophysics2Mod.F90 BiogeophysicsLakeMod.F90 \
C13SummaryMod.F90 CASAMod.F90 CASAPhenologyMod.F90 CNAllocationMod.F90 \
CNAnnualUpdateMod.F90 CNBalanceCheckMod.F90 CNC13FluxMod.F90 \
CNC13StateUpdate1Mod.F90 CNC13StateUpdate2Mod.F90 CNC13StateUpdate3Mod.F90 \
CNCStateUpdate1Mod.F90 CNCStateUpdate2Mod.F90 CNCStateUpdate3Mod.F90 \
CNDecompMod.F90 CNEcosystemDynMod.F90 CNFireMod.F90 CNGRespMod.F90 \
CNGapMortalityMod.F90 CNMRespMod.F90 CNNDynamicsMod.F90 CNNStateUpdate1Mod.F90 \
CNNStateUpdate2Mod.F90 CNNStateUpdate3Mod.F90 CNPhenologyMod.F90 \
CNPrecisionControlMod.F90 CNSetValueMod.F90 CNSummaryMod.F90 \
CNVegStructUpdateMod.F90 CNiniSpecial.F90 CNiniTimeVar.F90 CNrestMod.F90 \
CanopyFluxesMod.F90 DGVMAllocationMod.F90 DGVMEcosystemDynMod.F90 \
DGVMEstablishmentMod.F90 DGVMFireMod.F90 DGVMKillMod.F90 DGVMLightMod.F90 \
DGVMMod.F90 DGVMMortalityMod.F90 DGVMReproductionMod.F90 DGVMRestMod.F90 \
DGVMTurnoverMod.F90 DUSTMod.F90 DriverInitMod.F90 ESMF_AlarmClockMod.F90 \
ESMF_AlarmMod.F90 ESMF_BaseMod.F90 ESMF_BaseTimeMod.F90 ESMF_CalendarMod.F90 \
ESMF_ClockMod.F90 ESMF_FractionMod.F90 ESMF_Mod.F90 ESMF_Stubs.F90 \
ESMF_TimeIntervalMod.F90 ESMF_TimeMod.F90 FracWetMod.F90 \
FrictionVelocityMod.F90 Hydrology1Mod.F90 Hydrology2Mod.F90 \
HydrologyLakeMod.F90 Meat.F90 QSatMod.F90 RtmMod.F90 RunoffMod.F90 \
STATICEcosysDynMod.F90 SnowHydrologyMod.F90 SoilHydrologyMod.F90 \
SoilTemperatureMod.F90 SurfaceAlbedoMod.F90 SurfaceRadiationMod.F90 \
TridiagonalMod.F90 VOCEmissionMod.F90 abortutils.F90 accFldsMod.F90 \
accumulMod.F90 areaMod.F90 atmdrvMod.F90 clm_atmlnd.F90 clm_comp.F90 \
clm_csmMod.F90 clm_mct_mod.F90 clm_time_manager.F90 clm_varcon.F90 \
clm_varctl.F90 clm_varorb.F90 clm_varpar.F90 clm_varsur.F90 clmtype.F90 \
clmtypeInitMod.F90 controlMod.F90 decompMod.F90 do_close_dispose.F90 \
domainMod.F90 driver.F90 eshr_estate_mod.F90 eshr_inputinfo_mod.F90 \
eshr_rc_mod.F90 eshr_timemgr_mod.F90 f_wrappers.c fileutils.F90 \
filterMod.F90 get_zeits.c getdatetime.F90 gptl.c gptl_papi.c gptlutil.c \
histFileMod.F90 histFldsMod.F90 iniTimeConst.F90 \
inicFileMod.F90 initGridCellsMod.F90 initSurfAlbMod.F90 initializeMod.F90 \
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
mct_mod.F90 mkarbinitMod.F90 nanMod.F90 ncdio.F90 ndepFileMod.F90 \
perf_mod.F90 pft2colMod.F90 pftdynMod.F90 pftvarcon.F90 program_csm.F90 \
program_off.F90 restFileMod.F90 rtcrate.F90 scam_setlatlonidx.F90 \
shr_alarm_mod.F90 shr_cal_mod.F90 shr_const_mod.F90 shr_date_mod.F90 \
shr_file_mod.F90 shr_flux_mod.F90 shr_inputinfo_mod.F90 shr_jlcp.c \
shr_kind_mod.F90 shr_map_mod.F90 shr_mpi_mod.F90 shr_msg_mod.F90 \
shr_ncio_mod.F90 shr_ncread_mod.F90 shr_orb_mod.F90 shr_stream_mod.F90 \
shr_string_mod.F90 shr_sys_mod.F90 shr_tInterp_mod.F90 shr_timer_mod.F90 \
shr_vmath_fwrap.c shr_vmath_mod.F90 snowdp2lev.F90 spmdGathScatMod.F90 \
spmdMod.F90 subgridAveMod.F90 subgridMod.F90 subgridRestMod.F90 \
surfrdMod.F90 system_messages.F90 threadutil.c wrf_error_fatal.F90 \
wrf_message.F90

OBJS = $(addsuffix .o, $(basename $(SRCS)))

all: libclm.a

libclm.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

ESMF_INCLUDE = $(CLMSRC)/src/utils/esmf_wrf_timemgr
UTIL_INCLUDE = $(CLMSRC)/src/utils/timing/
GLOB_INCLUDE = .

CPPFLAGS += -DMAXPATCH_PFT=4 -DOFFLINE
F90FLAGS += -I$(ESMF_INCLUDE) -I$(UTIL_INCLUDE) -I$(GLOB_INCLUDE) \
            $(NETCDFINC)

%.o: %.F90
	$(F90) -c -o $@ $(CPPFLAGS) $(F90FLAGS) $(USCORING) -c $<

%.o: %.f90
	$(F90) -c -o $@ $(F90FLAGS) $(USCORING) -c $<

clean:
	rm -f *.mod *.o libclm.a

include clmdeps
