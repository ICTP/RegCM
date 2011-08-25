program test_eshr_timemgr
use ESMF_Mod
use shr_kind_mod,       only: SHR_KIND_CL, r8=>SHR_KIND_R8
use shr_orb_mod,        only: SHR_ORB_UNDEF_REAL
use eshr_timemgr_mod,   only: &
                              eshr_timemgr_NMLinfoType,       &
                              eshr_timemgr_NMLinfoRead,       &
                              eshr_timemgr_NMLinfoSetDefault, &
                              eshr_timemgr_NMLinfoPutData,    &
                              eshr_timemgr_NMLinfoGetData,    &
                              eshr_timemgr_NMLinfoRestRead,   &
                              eshr_timemgr_clockType,         &
                              eshr_timemgr_clockInitNMLinfo,  &
                              eshr_timemgr_clockIsOnLastStep, &
                              eshr_timeMgr_curTimeLEstopTime, &
                              eshr_timemgr_clockGet,          &
                              eshr_timemgr_clockInfoType,     &
                              eshr_timemgr_clockInfoGet,      &
                              eshr_timemgr_clockInfoPutData,  &
                              eshr_timemgr_clockPutData,      &
                              eshr_timemgr_clockGetRunDays,   &
                              eshr_timemgr_clockRestWrite,    &
                              eshr_timemgr_clockAlarmOffRest, &
                              eshr_timemgr_clockPrint,        &
                              eshr_timemgr_clockAdvance,      &
                              eshr_timemgr_clockAlarmOnRest,  &
                              eshr_timemgr_clockAlarmIsOnRes, &
                              eshr_timemgr_clockGetCalDay,    &
                              eshr_timemgr_clockClocksInSync, &
                              eshr_timemgr_clockDateInSync,   &
                              eshr_timemgr_clockInitClock,    &
                              eshr_timemgr_clockInfoType,     &
                              eshr_timemgr_clockStopAtDayEnd, &
                              eshr_timemgr_clockInfoGet,      &
                              eshr_timemgr_clockGetPerpYMD
use shr_inputInfo_mod, only:  shr_inputInfo_initType,         &
                              shr_inputInfo_initSetDefaults,  &
                              shr_inputInfo_initRead,         &
                              shr_inputInfo_initPutData,      &
                              shr_inputInfo_initGetData,      &
                              shr_inputInfo_initIsBranch,     &
                              shr_inputInfo_initRPointerWrite,&
                              shr_inputInfo_initRestWrite,    &
                              shr_inputInfo_initRPointerRead, &
                              shr_inputInfo_initRestRead

use shr_const_mod,      only: SHR_CONST_PI
use shr_sys_mod,        only: shr_sys_abort, shr_sys_system
use shr_orb_mod,        only: SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
use shr_mpi_mod,        only: shr_mpi_init, shr_mpi_commrank, shr_mpi_finalize

implicit none

type(shr_inputInfo_InitType) :: CCSMInit
type(eshr_timemgr_clockType) :: clock
type(eshr_timemgr_clockType) :: clock2
type(eshr_timemgr_clockType) :: ClockRest
type(eshr_timemgr_clockInfoType) :: ClockInfo
type(eshr_timemgr_NMLinfoType) :: ClockNMLinit, ClockNMLinit2
character(len=SHR_KIND_CL) :: NLFilename
character(len=SHR_KIND_CL) :: rest_file, rest_file2
type(ESMF_Alarm) :: restart
type(ESMF_Clock) :: esmfClock
type(ESMF_Time) :: RingTime
integer :: MPICom = 1
integer :: year = 0, month = 12, day = 1, tod = 0
integer :: ymd, sec
integer :: YearRest = 0, MonthRest = 12, DayRest = 1, TODRest = 0
logical :: MasterTask
integer :: MPIRank
logical :: perpetual_run
integer :: perpetual_ymd
integer :: test
integer :: iType
integer :: CurrentYMD, CurrentTOD, CurrentYear, CurrentMonth, CurrentDay
integer :: PrevYMD, PrevTOD, PrevDays, PrevSecs
integer :: CurrentDays, CurrentSecs, StepNo
integer, parameter :: DTime = 1200    ! Time-step
integer :: TimeStep
! Stop times to use for all the different ways of setting stop-time
integer, parameter :: NstepsType     = 1
integer, parameter :: NdaysType      = NstepsType     + 1
integer, parameter :: MonthlyType    = NdaysType      + 1
integer, parameter :: NMonthsType    = MonthlyType    + 1
integer, parameter :: YearlyType     = NMonthsType    + 1
integer, parameter :: DateType       = YearlyType     + 1
integer, parameter :: StopEndDayType = DateType       + 1
integer, parameter :: DateTypeII     = StopEndDayType + 1
integer, parameter :: nTestTypes     = DateTypeII     ! Number of test types
integer, parameter :: nTests = 5                  ! Number of tests for each type
integer, parameter :: stop_nsteps(nTests)    = (/   1,    5,   10,    20,    60 /)
integer, parameter :: stop_ndays(nTests)     = (/   1,    5,   10,    20,    50 /)
integer, parameter :: stop_nmonths(nTests)   = (/   1,    5,    6,    12,    24 /)
integer, parameter :: stop_nyears(nTests)    = (/   1,    2,    3,     4,     5 /)
! Note: dayend MUST NOT be the end of a month!
integer, parameter :: stop_dayend_ymd(nTests)= (/   530,    718,   726,    1030,   1204 /)
integer, parameter :: stop_ymd(nTests)       = (/   527,    603,   616,    1023,   1104 /)
integer, parameter :: stop_ymdII(nTests)     = (/990927, 991103, 991216,1000123,1000104 /)
integer, parameter :: start_tod(nTests)      = (/   0, DTime, DTime*3, DTime*35, DTime*70/)
! Start time for each of the above tests
integer, parameter :: start_ymd(nTests,nTestTypes)    = reshape( (/                  &
                            (/    101,     228,    331,    430,     531 /), & ! nsteps
                            (/    102,     223,    322, 990415,20000531 /), & ! ndays
                            (/    102,     223,    322, 990415,20000531 /), & ! nmonths
                            (/    102,     223,    322, 990415,20000531 /), & ! nmonths/end
                            (/    102,     223,    322, 990415,20000531 /), & ! years
                            (/    115,     213,    312,    415,     507 /), & ! ymd
                            (/    115,     213,    312,    415,     507 /), & ! dayend
                            (/ 990615,  990713, 990812, 990915,  991007 /)  & ! ymd II
                            /), (/nTests,nTestTypes/) )
! End date for each of the above tests
integer, parameter :: EndYMDTable(nTests,4)    = reshape( (/                &
                            (/    101,     228,    331,    430,     601 /), & ! nsteps
                            (/    103,     228,    401, 990505,20000720 /), & ! ndays
                            (/    202,     723,    922,1000415,20020531 /), & ! nmonths
                            (/    202,     723,    922,1000415,20020531 /)  & ! nmonths type
                            /), (/nTests,4/) )
integer :: restart_nsteps(nTests) = (/   1,    2,    5,    10,     6 /)
integer :: restart_ndays(nTests)  = (/   1,    2,    3,     4,     5 /)
integer :: EndYMD
integer :: EndTOD
integer :: FirstYMD
integer :: FirstTOD
integer :: iStep
integer, parameter :: perp_ymd = 20060321
integer :: idx
integer :: orb_iyear_AD
real(r8) :: orb_eccen, orb_mvelp, orb_obliq
integer :: orb_iyear_AD_in(2) = (/ SHR_ORB_UNDEF_INT, 1992 /)
real(r8) :: orb_eccen_in(2) = (/ 0.05_r8, SHR_ORB_UNDEF_REAL /)
real(r8) :: orb_mvelp_in(2) = (/ 272._r8, SHR_ORB_UNDEF_REAL /)
real(r8) :: orb_obliq_in(2) = (/ 89.9_r8, SHR_ORB_UNDEF_REAL /)
real(r8) :: orb_lambm0_in(1) = (/ 0.05_r8*SHR_CONST_PI/180.0_r8 /)
real(r8) :: orb_mvelpp_in(1) = (/ 272._r8*SHR_CONST_PI/180.0_r8 /)
real(r8) :: orb_obliqr_in(1) = (/ 89.9_r8*SHR_CONST_PI/180.0_r8 /)
real(r8) :: orb_obliqr
real(r8) :: orb_lambm0
real(r8) :: orb_mvelpp
real(r8) :: calday
integer  :: rc, rcode
logical :: LastRestart = .false.
logical :: fromclock
integer :: restart_n, restart_n2, restart_n3
character(len=256) :: restart_option, restart_option2, restart_option3
character(len=*), parameter :: F00 = &
"('Advance clock: current date: ', i2.2, '/', i2.2, '/', i4.4, ' : ', i6.6, ' (Sec) CalDay:', f17.12 )"
#include <mpif.h>

  call shr_mpi_init( )
  MPICom = MPI_COMM_WORLD
  call shr_mpi_commrank( MPICom, MPIRank )
  MasterTask = (MPIRank == 0)

call ESMF_Initialize( rc=rc )

NLFilename = "timemgr.namelist"
call shr_inputInfo_initSetDefaults( CCSMInit )
call shr_inputInfo_initPutData( CCSMInit, start_type="startup", case_name="csmrun" )
call shr_inputInfo_initRead( NLFilename, LogPrint=MasterTask, MPICom=MPICom, MasterTask=MasterTask, &
                           CCSMInitOut=CCSMInit )
call shr_inputInfo_initGetData( CCSMInit, perpetual_run=perpetual_run, perpetual_ymd=perpetual_ymd )
call eshr_timemgr_NMLinfoSetDefault(   perpetual_run, perpetual_ymd, ClockNMLinit )
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, restart_n=1 )
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, desc="master test synchronization clock" )
call eshr_timemgr_NMLinfoRead( NLFilename, logPrint=MasterTask, MPICom=MPICom, MasterTask=MasterTask, &
                         clockNMLinfoOut=ClockNMLinit )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, &
                                   clockOut=clock )

if ( MasterTask ) call eshr_timemgr_clockPrint( clock )
if ( MasterTask ) print *, 'Loop until stop_time'
iStep = 0
do while ( eshr_timeMgr_curTimeLEstopTime(   clock ) )
    ! restart alarm rings if restart time OR if next time-step is stop_time
    if ( eshr_timemgr_clockAlarmIsOnRes( clock ) )then
      if ( MasterTask ) print *, 'Is restart time: '
      call eshr_timemgr_clockGet(     clock, CurrentYMD=ymd, CurrentTOD=tod )
      call shr_inputInfo_initRPointerWrite(  ymd, tod, MPICom, MasterTask, &
                                        CCSMInit=CCSMInit, restart_file=rest_file )
      call shr_inputInfo_initRestWrite(   rest_file, MPICom, MasterTask, &
                                         CCSMInit=CCSMInit  )
      call eshr_timemgr_clockRestWrite(     rest_file, MPICom, MasterTask, clock )
      call eshr_timemgr_clockAlarmOffRest( clock )
    end if
    call eshr_timemgr_clockAdvance(  clock )
    call eshr_timemgr_clockGet(     clock, year=year, month=month, day=day, CurrentTOD=tod )
    calday = eshr_timemgr_clockGetCalDay( clock )
    iStep = iStep + 1
    if ( MasterTask ) print F00, month, day, year, tod, calday
    if ( iStep > 70*24*3600/DTime )then
       call shr_sys_abort( 'Simulation proceeding for too long--1' )
    end if
end do
! Now check that can get clockInfo orbital info and set it to something different
call eshr_timemgr_clockGet(     clock, info=clockInfo )
call eshr_timemgr_clockInfoGet( clockInfo, orb_eccen=orb_eccen, orb_obliq=orb_obliq, &
                                          orb_mvelp=orb_mvelp, orb_obliqr=orb_obliqr, &
                                          orb_mvelpp=orb_mvelpp, orb_lambm0=orb_lambm0 )
if ( orb_eccen == SHR_ORB_UNDEF_REAL )then
  call shr_sys_abort( 'orb_eccen NOT set' )
end if
if ( orb_obliq == SHR_ORB_UNDEF_REAL )then
  call shr_sys_abort( 'orb_obliq NOT set' )
end if
if ( orb_mvelp == SHR_ORB_UNDEF_REAL )then
  call shr_sys_abort( 'orb_mvelp NOT set' )
end if
if ( orb_obliqr == SHR_ORB_UNDEF_REAL )then
  call shr_sys_abort( 'orb_obliqr NOT set' )
end if
if ( orb_mvelpp == SHR_ORB_UNDEF_REAL )then
  call shr_sys_abort( 'orb_mvelpp NOT set' )
end if
if ( orb_lambm0 == SHR_ORB_UNDEF_REAL )then
  call shr_sys_abort( 'orb_lambm0 NOT set' )
end if
idx = 1
call eshr_timemgr_clockInfoPutData( clockInfo, orb_eccen=orb_eccen_in(idx),  &
                                              orb_obliq=orb_obliq_in(idx),   &
                                              orb_mvelp=orb_mvelp_in(idx),   &
                                              orb_obliqr=orb_obliqr_in(idx), &
                                              orb_mvelpp=orb_mvelpp_in(idx), &
                                              orb_lambm0=orb_lambm0_in(idx) )
call eshr_timemgr_clockPutData(     clock, clockInfo=clockInfo )
call verify_orbit( clock, orb_eccen_in=orb_eccen_in(idx),   &
                   orb_obliq_in=orb_obliq_in(idx),   &
                   orb_mvelp_in=orb_mvelp_in(idx),   &
                   orb_obliqr_in=orb_obliqr_in(idx), &
                   orb_mvelpp_in=orb_mvelpp_in(idx), &
                   orb_lambm0_in=orb_lambm0_in(idx) )

! Test that if do not set orbit, can get it set properly
call eshr_timemgr_NMLinfoSetDefault(   ClockNMLinfo=ClockNMLinit2 )
call eshr_timemgr_NMLinfoPutData( ClockNMLinit2, orb_notSet=.true., &
                     stop_option="ndays", stop_n=999999, start_ymd=20060527, &
                     atm_cpl_dt=3600, lnd_cpl_dt=3600, ocn_cpl_dt=3600, &
                     ice_cpl_dt=3600 )
call eshr_timemgr_NMLinfoPutData( ClockNMLinit2, desc="test sync clock orbit set later" )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit2, logPrint=MasterTask, &
                                   clockOut=clock2 )
call verify_orbit( clock2, orb_eccen_in=SHR_ORB_UNDEF_REAL, &
                   orb_obliq_in=SHR_ORB_UNDEF_REAL,  &
                   orb_mvelp_in=SHR_ORB_UNDEF_REAL,  &
                   orb_obliqr_in=SHR_ORB_UNDEF_REAL, &
                   orb_mvelpp_in=SHR_ORB_UNDEF_REAL, &
                   orb_lambm0_in=SHR_ORB_UNDEF_REAL )
call eshr_timemgr_clockGet(     clock2, Info=clockInfo )
idx = 1
call eshr_timemgr_clockInfoPutData( clockInfo, orb_eccen=orb_eccen_in(idx),   &
                                              orb_obliq=orb_obliq_in(idx),   &
                                              orb_mvelp=orb_mvelp_in(idx),   &
                                              orb_obliqr=orb_obliqr_in(idx), &
                                              orb_mvelpp=orb_mvelpp_in(idx), &
                                              orb_lambm0=orb_lambm0_in(idx) )
call eshr_timemgr_clockPutData(     clock2, clockInfo=clockInfo, orb_set=.true. )
call verify_orbit( clock2, orb_eccen_in=orb_eccen_in(idx),   &
                   orb_obliq_in=orb_obliq_in(idx),   &
                   orb_mvelp_in=orb_mvelp_in(idx),   &
                   orb_obliqr_in=orb_obliqr_in(idx), &
                   orb_mvelpp_in=orb_mvelpp_in(idx), &
                   orb_lambm0_in=orb_lambm0_in(idx) )
rest_file2 = "tmp."//trim(rest_file)
call eshr_timemgr_clockRestWrite(     rest_file2, MPICom, MasterTask, clock2 )

! Get current date and print out clock
call eshr_timemgr_clockGet(     clock, CurrentYMD=ymd, CurrentTOD=tod )
if ( MasterTask ) call eshr_timemgr_clockPrint( clock )

if ( MasterTask ) print *, 'Now read in the restart file for a continue run: '
call shr_inputInfo_initSetDefaults( CCSMInit )
call shr_inputInfo_initPutData( CCSMInit, start_type="continue", case_name="csmrun" )
call shr_inputInfo_initRead( NLFilename, LogPrint=MasterTask, MPICom=MPICom, MasterTask=MasterTask, &
                           CCSMInitOut=CCSMInit )
call shr_inputInfo_initGetData( CCSMInit, perpetual_run=perpetual_run, perpetual_ymd=perpetual_ymd )
call eshr_timemgr_NMLinfoSetDefault(   perpetual_run, perpetual_ymd, ClockNMLinit )
call shr_inputInfo_initRPointerRead( MPICom, MasterTask, CCSMInit, rest_file )
call shr_inputInfo_initRestRead(  rest_file, MPICom=MPICom, &
                                 MasterTask=MasterTask, CCSMInitOut=CCSMInit )
call eshr_timemgr_NMLinfoRead( NLFilename, LogPrint=MasterTask, MPICom=MPICom, &
                         MasterTask=MasterTask, ClockNMLinfoOut=ClockNMLinit )
call eshr_timemgr_NMLinfoRestRead(    rest_file, MPICom=MPICom, &
                                     MasterTask=MasterTask, clockNMLinfo=ClockNMLinit )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, clockOut=clockRest )
if ( MasterTask ) call eshr_timemgr_clockPrint( ClockRest )
call eshr_timemgr_clockGet(     ClockRest, year=YearRest, month=MonthRest, day=DayRest, &
                               CurrentTOD=TODRest, CurrentYMD=ymd )
if ( (YearRest /= year) .or. (MonthRest /= month) .or. (DayRest /= day) .or. &
     (TODRest /= tod) )then
  write(6,*) "yearRest, MonthRest, DayRest, TODRest = ",  &
              yearRest, MonthRest, DayRest, TODRest
  write(6,*) "year, Month, Day, TOD = ",  &
              year, Month, Day, TOD
  call shr_sys_abort( 'Restart time not same as time before' )
end if
! Test that if create new clock from this one that they will be in sync
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, CurrentYMD=ymd, CurrentTOD=TODRest )
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, desc="test clock 2",    &
                                 MasterSyncClock=.false., restart=.false. )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, clockOut=clock2 )
if ( .not. eshr_timemgr_clockClocksInSync( clock2, clockRest ) )then
  call shr_sys_abort( 'restart clock and clock created with current time NOT in sync' )
end if
call eshr_timemgr_clockGet(     ClockRest, CurrentTOD=TODRest, CurrentYMD=ymd )
if ( .not. eshr_timemgr_clockDateInSync( clock2, ymd, TODRest ) )then
  call shr_sys_abort( 'restart clock and clock created with current time NOT in sync by Date' )
end if
! Now if advance second clock -- they won't be
call eshr_timemgr_clockAdvance( clock2 )
if ( eshr_timemgr_clockClocksInSync( clock2, clockRest ) )then
  call shr_sys_abort( 'restart clock and clock created with current time '// &
                      'are in sync and should NOT be' )
end if
call eshr_timemgr_clockGet(     ClockRest, CurrentTOD=TODRest, CurrentYMD=ymd )
if ( .not. eshr_timemgr_clockDateInSync( clock2, ymd, TODRest, prev=.true. ) )then
  call shr_sys_abort( 'restart clock and clock created with current time after advanced NOT in sync by Date with previous' )
end if
! Test that if create a clock in perpetual mode that will NOT be in sync with this one
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, perpetual_run=.true., perpetual_ymd=perp_ymd)
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, CurrentYMD=ymd, CurrentTOD=TODRest, &
                                 restart=.false. )
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, desc="test clock",    &
                                 MasterSyncClock=.false. )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, clockOut=clock2 )
if ( eshr_timemgr_clockClocksInSync( clock2, clockRest ) )then
  call shr_sys_abort( 'restart clock and perpetual clock created with current '// &
                      'time in sync and should NOT be' )
end if
! Also test eshr_timemgr_clockGetPerpYMD
call eshr_timemgr_clockGetPerpYMD( clock2, ymd, sec, offset=0 )
if ( ymd /= perp_ymd .and. sec /= TODRest )then
  call shr_sys_abort( 'Bad results from clockGetPerpYMD method' )
end if
call eshr_timemgr_clockGetPerpYMD( clock2, ymd, sec, offset=3600 )
if ( ymd /= perp_ymd .and. sec /= (TODRest+3600) )then
  call shr_sys_abort( 'Bad results from clockGetPerpYMD method' )
end if
! Test with offset
! Test that if create another perpetual mode clock that they will be in sync
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, perpetual_run=.true., perpetual_ymd=perp_ymd)
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, CurrentYMD=ymd, CurrentTOD=TODRest )
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, desc="test perpetual clock",    &
                                 MasterSyncClock=.false. )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, clockOut=clock )
if ( .not. eshr_timemgr_clockClocksInSync( clock2, clock ) )then
  call shr_sys_abort( 'perpetual clock and perpetual clock created with current '// &
                      'time are NOT in sync' )
end if
! Now if advance second clock -- they won't be in sync
call eshr_timemgr_clockAdvance( clock )
if ( eshr_timemgr_clockClocksInSync( clock2, clock ) )then
  call shr_sys_abort( 'perpetual clocks '// &
                      'are in sync and should NOT be' )
end if
! Test that can setup a new clock from RestClock and that it will be in sync
call eshr_timemgr_clockGet( clockRest, dtime=TimeStep )
call eshr_timemgr_clockInitClock( ClockRest, TimeStep, LogPrint=MasterTask, &
                                 desc="test init from clock", ClockOut=clock2 )
if ( .not. eshr_timemgr_clockClocksInSync( clock2, clockRest ) )then
  call shr_sys_abort( 'clock from and restart clock and restart clock '// &
                      'are NOT in sync' )
end if
! Now if advance second clock -- they won't be
call eshr_timemgr_clockAdvance( clock2 )
if ( eshr_timemgr_clockClocksInSync( clock2, clockRest ) )then
  call shr_sys_abort( 'clock and restart clock '// &
                      'are in sync and should NOT be' )
end if
! Test if can override restart frequency
call eshr_timemgr_NMLinfoGetData( ClockNMLinit, restart_option=restart_option, restart_n=restart_n )
restart_n3      = restart_n
restart_option3 = restart_option
restart_n       = 18
restart_option  = "nmonths"
restart_n2      = restart_n
restart_option2 = restart_option
call eshr_timemgr_NMLinfoRestRead(    rest_file, MPICom=MPICom, &
                                     MasterTask=MasterTask, clockNMLinfo=ClockNMLinit )
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, &
                           restart_option=restart_option, restart_n=restart_n )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, clockOut=clockRest )
! Make sure both clock and ClockNMLinit get changed properly
call eshr_timemgr_clockGet(     ClockRest, RestartIntervalMonths=restart_n )
call eshr_timemgr_NMLinfoGetData(     ClockNMLinit, restart_option=restart_option )
if ( (restart_n2 /= restart_n) .and. (restart_option /= restart_option2) )then
  call shr_sys_abort( 'Restart option and n not correctly overridden' )
end if
call eshr_timemgr_NMLinfoGetData(     ClockNMLinit, restart_n=restart_n )
if ( restart_n2 /= restart_n )then
  call shr_sys_abort( 'ClockNMLinit not properly updated after readrestart for override case' )
end if
! Test if won't override if override option not set
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, restart_file_TMGoverride="" )
call eshr_timemgr_NMLinfoRestRead(    rest_file, MPICom=MPICom, &
                                     MasterTask=MasterTask, clockNMLinfo=ClockNMLinit )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, clockOut=clockRest )
call eshr_timemgr_NMLinfoGetData( ClockNMLinit, restart_option=restart_option )
call eshr_timemgr_clockGet(     ClockRest, RestartIntervalYears=restart_n )
if ( (restart_n3 /= restart_n) .and. (restart_option /= restart_option3) )then
  call shr_sys_abort( 'Restart n and restart_option were overridden and should not have been' )
end if
call eshr_timemgr_NMLinfoGetData(     ClockNMLinit, restart_n=restart_n )
if ( restart_n3 /= restart_n )then
  call shr_sys_abort( 'ClockNMLinit not properly updated after readrestart for non-override case' )
end if
restart_n2      = restart_n
restart_option2 = restart_option
restart_n       = 18
restart_option  = "nmonths"
call eshr_timemgr_NMLinfoPutData( ClockNMLinit, restart_option=restart_option, restart_n=restart_n )
call eshr_timemgr_NMLinfoRestRead(    rest_file, MPICom=MPICom, &
                                     MasterTask=MasterTask, clockNMLinfo=ClockNMLinit )
call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, clockOut=clockRest )
call eshr_timemgr_NMLinfoGetData( ClockNMLinit, restart_option=restart_option )
call eshr_timemgr_clockGet(     ClockRest, RestartIntervalYears=restart_n )
if ( (restart_n2 /= restart_n) .or. (restart_option /= restart_option2) )then
  call shr_sys_abort( 'Restart option and n not correctly saved from restart file' )
end if

if ( MasterTask ) print *, 'Loop until stop_time'
iStep = 0
do while ( eshr_timeMgr_curTimeLEstopTime(   clockRest ) )
    ! restart alarm rings if restart time OR if next time-step is stop_time
    if ( eshr_timemgr_clockAlarmIsOnRes( ClockRest ) )then
      if ( MasterTask ) print *, 'Is restart time: '
         call shr_inputInfo_initRPointerWrite(  ymd, tod, MPICom, MasterTask, &
                                            CCSMInit=CCSMInit, restart_file=rest_file )
         call shr_inputInfo_initRestWrite(   rest_file, MPICom, MasterTask, &
                                            CCSMInit=CCSMInit  )
         call eshr_timemgr_clockRestWrite(     rest_file, MPICom, MasterTask, ClockRest )
         call eshr_timemgr_clockAlarmOffRest(  ClockRest )
    end if
    call eshr_timemgr_clockAdvance(  ClockRest )
    call eshr_timemgr_clockGet(     ClockRest, CurrentYMD=ymd, year=year, month=month, &
                                   day=day, CurrentTOD=tod )
    calday = eshr_timemgr_clockGetCalDay( ClockRest )
    iStep = iStep + 1
    if ( MasterTask ) print F00, month, day, year, tod, calday
    if ( iStep > 66*24*3600/DTime )then
       call shr_sys_abort( 'Simulation proceeding for too long-2' )
    end if
end do

if ( MasterTask ) print *, 'Now do more extensive testing with start and stop times'
NLFilename = "timemgr.minimumnamelist"
! Loop over test types
do iType = 1, nTestTypes
   ! Loop over seperate tests for each test type
   do test = 1, nTests
      call shr_inputInfo_initSetDefaults( CCSMInit )
      call shr_inputInfo_initPutData( CCSMInit, start_type="startup", case_name="csmruntest" )
      call shr_inputInfo_initRead( NLFilename, LogPrint=MasterTask, MPICom=MPICom, &
                                 MasterTask=MasterTask, CCSMInitOut=CCSMInit )
      if ( test == 4 )then
         call shr_inputInfo_initPutData( CCSMInit, aqua_planet = .true. )
      end if
      call shr_inputInfo_initGetData( CCSMInit, perpetual_run=perpetual_run, &
                              perpetual_ymd=perpetual_ymd )
      call eshr_timemgr_NMLinfoSetDefault( perpetual_run, perpetual_ymd, ClockNMLinit )
      if ( mod(test,2) == 0 )then
        idx = 1
      else
        idx = 2
      end if
      call eshr_timemgr_NMLinfoPutData( ClockNMLinit, &
                               orb_iyear_AD=orb_iyear_AD_in(idx), orb_mvelp=orb_mvelp_in(idx),   &
                               orb_eccen=orb_eccen_in(idx), orb_obliq=orb_obliq_in(idx) )
      if (      iType == NstepsType )then
         call eshr_timemgr_NMLinfoPutData( ClockNMLinit, stop_option="nsteps", stop_n= stop_nsteps(test), &
                                  restart_option="nsteps", restart_n= restart_nsteps(test) )
         EndYMD = EndYMDTable(test,iType)
         EndTOD = start_tod(test) + (stop_nsteps(test)+1)*DTime
         if ( EndTOD > 24*3600 ) EndTOD = EndTOD - 24*3600
      else if ( iType == NdaysType )then
         call eshr_timemgr_NMLinfoPutData( ClockNMLinit, stop_option="ndays", stop_n = stop_ndays(test), &
                                  restart_option="ndays", restart_n = 1 )
         EndYMD = EndYMDTable(test,iType)
         EndTOD = start_tod(test) + DTime
      else if ( iType == MonthlyType )then
         call eshr_timemgr_NMLinfoPutData( ClockNMLinit, stop_option = "nmonths", stop_n = stop_nmonths(test), &
                                  restart_option = 'monthly' )
         EndYMD = EndYMDTable(test,iType)
         EndTOD = start_tod(test) + DTime
      else if ( iType == NMonthsType )then
         call eshr_timemgr_NMLinfoPutData( ClockNMLinit, stop_option = "nmonths", stop_n = stop_nmonths(test), &
                                          restart_option = "end" )
         EndYMD = EndYMDTable(test,iType)
         EndTOD = start_tod(test) + DTime
      else if ( iType == YearlyType )then
         call eshr_timemgr_NMLinfoPutData( ClockNMLinit, stop_option = "nyears", stop_n = stop_nyears(test), &
                                  restart_option = "yearly" )
         EndYMD = start_ymd(test,iType) + stop_nyears(test)*10000
         EndTOD = start_tod(test) + DTime
      else if ( iType == DateType )then
         call eshr_timemgr_NMLinfoPutData( ClockNMLinit, restart_option = "end", stop_option = "date", &
                                  stop_ymd    = stop_ymd(test) )
         EndYMD = stop_ymd(test)
         EndTOD = DTime
      else if ( iType == StopEndDayType )then
         call eshr_timemgr_NMLinfoPutData( ClockNMLinit, restart_option = "none",   &
                                          stop_option = "ndays", stop_n      = 9999 )
         EndYMD = stop_dayend_ymd(test) + 1
         EndTOD = DTime
      else if ( iType == DateTypeII )then
         call eshr_timemgr_NMLinfoPutData( ClockNMLinit, restart_option = "ndays", restart_n = 3, &
                                           stop_option = "date", stop_ymd = stop_ymdII(test) )
         EndYMD = stop_ymdII(test)
         EndTOD = DTime
      else
         call shr_sys_abort( 'ERROR:: type out of bounds' )
      end if
      call eshr_timemgr_NMLinfoPutData( ClockNMLinit, atm_cpl_dt = DTime,     &
                                  start_ymd=start_ymd(test,iType), &
                                  start_tod=start_tod(test) )
      call eshr_timemgr_NMLinfoRead( NLFilename, logPrint=MasterTask, MPICom=MPICom, &
                               MasterTask=MasterTask, clockNMLinfoOut=ClockNMLinit )
      call eshr_timemgr_clockInitNMLinfo( ClockNMLinit, logPrint=MasterTask, &
                                         ClockOut=clock )
      fromclock = .false.
      call eshr_timemgr_clockGet(   clock, CurrentYMD=CurrentYMD, CurrentTOD=tod, eClock=esmfClock, &
                                   info=clockInfo )
      if ( CurrentYMD /= start_ymd(test,iType) .or. (tod /= start_tod(test)) )then
         call shr_sys_abort( 'Clock startup current time NOT equal to start time given in setup' )
      end if
      call eshr_timemgr_clockInfoGet( clockInfo, orb_iyear_AD=orb_iyear_AD,     &
                                     orb_mvelp=orb_mvelp, orb_eccen=orb_eccen, &
                                     orb_obliq=orb_obliq, perpetual_ymd=perpetual_ymd )
      if ( orb_iyear_AD_in(idx) == SHR_ORB_UNDEF_INT )then
      ! Orbital info set NOT based on year
         if ( (orb_iyear_AD /= orb_iyear_AD_in(idx)) .or. (orb_mvelp /= orb_mvelp_in(idx)) .or. &
              (orb_eccen /= orb_eccen_in(idx))       .or. (orb_obliq /= orb_obliq_in(idx)) )then
            call shr_sys_abort( 'orbital info got was different that what was set' )
         end if
      else
      ! Orbital info set based on year, so specific values should NOT correspond to original setting
         if ( (orb_iyear_AD /= orb_iyear_AD_in(idx)) .or. (orb_mvelp == orb_mvelp_in(idx)) .or. &
              (orb_eccen == orb_eccen_in(idx))       .or. (orb_obliq == orb_obliq_in(idx)) )then
            call shr_sys_abort( 'orbital info got was different that what was set' )
         end if
      end if
      if ( MasterTask ) call eshr_timemgr_clockPrint( clock )
      !call ESMF_clockGetAlarm( esmfClock, "Restart alarm", restart )
      !call ESMF_AlarmGet( restart, RingTime=RingTime )
      !if ( MasterTask ) print *, 'Restart alarm ring time: '
      !call ESMF_TimePrint( RingTime )
      if ( MasterTask ) print *, 'Loop until stop_time'
      iStep = 0
      ! Loop until stop time....
      do while ( eshr_timeMgr_curTimeLEstopTime(   clock ) )
          ! restart alarm rings if restart time OR if next time-step is stop_time
          if ( eshr_timemgr_clockAlarmIsOnRes( clock ) )then
            if ( MasterTask ) print *, 'Is restart time: '
               call eshr_timemgr_clockGet(     clock, year=year, month=month, CurrentYMD=ymd, &
                                              day=day, CurrentTOD=tod )
               ! If not last step, check if restarts are at the right frequency
               if ( .not. eshr_timemgr_clockIsOnLastStep( clock, nextStep=.true. ) )then
                  if ( iType == NstepsType .and. (mod(iStep,restart_nsteps(test)) /= 0) )then
                     call shr_sys_abort( 'Restarts should occur on restart_nsteps or just before stopping' )
                  end if
                  if ( iType == NdaysType .and. (tod /= start_tod(test) ) )then
                     call shr_sys_abort( 'Restarts should be daily' )
                  end if
                  if ( iType == MonthlyType .and.  ((tod /= 0) .or. (day /= 1)) )then
                     call shr_sys_abort( 'Restarts should be monthly' )
                  end if
                  if ( iType == YearlyType .and. ((tod /= 0) .or. (day /= 1) .or. (month /= 1)) ) then
                     call shr_sys_abort( 'Restarts should be yearly' )
                  end if
                  if ( iType == DateType )then
                     call shr_sys_abort( 'no restarts are done on stop_ymd case, except end' )
                  end if
                  if ( iType == DateTypeII .and. (tod /= start_tod(test) ) )then
                     call shr_sys_abort( 'restarts should be every 3 days' )
                  end if
               end if
               if ( iType == StopEndDayType )then
                  call shr_sys_abort( 'no restarts are done at all on stop at day end case' )
               end if
               call shr_inputInfo_initRPointerWrite(  ymd, tod, MPICom, &
                                                  MasterTask, CCSMInit=CCSMInit, &
                                                  restart_file=rest_file )
               call shr_inputInfo_initRestWrite(   rest_file, MPICom, MasterTask, &
                                                  CCSMInit=CCSMInit  )
               call eshr_timemgr_clockRestWrite(     rest_file, MPICom, MasterTask, clock )
               call eshr_timemgr_clockAlarmOffRest( clock )
               call eshr_timemgr_clockGet(     clock, year=year, month=month, day=day, CurrentTOD=tod )
               LastRestart = .true.
               ! Test that can setup a new clock from RestClock and that it will be in sync
               call eshr_timemgr_clockGet( clock, dtime=TimeStep )
               call eshr_timemgr_clockInitClock( Clock, TimeStep, LogPrint=MasterTask, &
                                                desc="from clock", ClockOut=clock2 )
               fromclock = .true.
               ! Remove file after writing so don't create tons of useless files
               call shr_sys_system( "/bin/rm -f "//trim(rest_file), rcode )
          else
               call eshr_timemgr_clockGet(     clock, year=year, month=month, day=day, CurrentTOD=tod )
               if ( iType == NdaysType .and. (iStep /= 0) .and. (tod == start_tod(test) ) )then
                  call shr_sys_abort( 'Restarts should be daily, did NOT restart' )
               end if
               if ( iType == MonthlyType .and. ((iStep /= 0) .and. (tod == 0) .and. (day == 1)) )then
                  call shr_sys_abort( 'Restarts should be monthly, did NOT restart' )
               end if
               if ( iType == YearlyType .and. ((iStep /= 0) .and. (tod == 0) .and. (day == 1) .and. month == 1) )then
                  call shr_sys_abort( 'Restarts should be yearly, did NOT restart' )
               end if
               LastRestart = .false.
               if ( fromclock )then
                  if ( .not. eshr_timemgr_clockClocksInSync( clock2, clock ) )then
                    call shr_sys_abort( 'clock from '// &
                                        'are NOT in sync' )
                  end if
               end if
          end if
          call eshr_timemgr_clockGet(     clock, year=year, month=month, day=day, CurrentTOD=tod )
          calday = eshr_timemgr_clockGetCalDay( clock )
          if ( MasterTask ) print F00, month, day, year, tod, calday
          call eshr_timemgr_clockGet( clock, CurrentYMD=CurrentYMD, &
                                     CurrentTOD=CurrentTOD )
          if ( iType == StopEndDayType )then
             if ( CurrentYMD == stop_dayend_ymd(test) )then
                call eshr_timemgr_clockStopAtDayEnd( LogPrint=MasterTask, clock=clock )
             end if
          end if
          call eshr_timemgr_clockGetRunDays( clock, CurrentDays=CurrentDays, &
                                            CurrentSecs=CurrentSecs )
          call eshr_timemgr_clockAdvance(  clock )
          if ( fromclock )then
             call eshr_timemgr_clockAdvance(  clock2 )
          end if
          if ( istep > 1 )then
             call eshr_timemgr_clockGet( clock, PrevYMD=PrevYMD, PrevTOD=PrevTOD )
             if ( (CurrentYMD /= PrevYMD) .or. (CurrentTOD /= PrevTOD) )then
                call shr_sys_abort( 'Error comparing current time to previous '//&
                                    'time after advance' )
             end if
             call eshr_timemgr_clockGetRunDays( clock, PrevDays=PrevDays, &
                                               PrevSecs=PrevSecs )
             if ( (CurrentDays /= PrevDays) .or. (CurrentSecs /= PrevSecs) )then
                call shr_sys_abort( 'Error comparing current run days to previous '//&
                                    'run days after advance' )
             end if
          end if
          iStep = iStep + 1
          call eshr_timemgr_clockGet( clock, StepNo=StepNo )
          if ( iStep /= StepNo )then
             call shr_sys_abort( 'Time-steps out of sync' )
          end if
          if ( iStep > maxval(stop_nyears)*366*24*3600/DTime )then
             call shr_sys_abort( 'Simulation proceeding for too long' )
          end if
      end do
      ! Reached stop-time now do end testing...
      if ( .not. LastRestart .and. (itype /= StopEndDayType) )then
         call shr_sys_abort( 'Did not do a restart on the last time-step' )
      end if
      ! Check that current time equals predicted end time
      call eshr_timemgr_clockGet(     clock, CurrentYMD=ymd, CurrentTOD=tod )
      if ( iType == NstepsType .and. (iStep .ne. (stop_nsteps(test)+1)) )then
        call shr_sys_abort( 'Did not run the correct number of nsteps' )
      end if
      if ( (ymd /= EndYMD) .or. (tod /= EndTOD) )then
        call shr_sys_abort( 'Ending date for simulation not at predicted date/time' )
      end if
   end do   ! End of test number
end do      ! End of test type

if ( MasterTask ) print *, 'Done!'
if ( MasterTask ) print *, 'Testing passed!'

  call shr_mpi_finalize( )

contains

subroutine verify_orbit( clock, orb_eccen_in,      &
                   orb_obliq_in,   &
                   orb_mvelp_in,   &
                   orb_obliqr_in, &
                   orb_mvelpp_in, &
                   orb_lambm0_in )
type(eshr_timemgr_clockType), intent(IN) :: clock
real(r8), intent(IN) :: orb_eccen_in
real(r8), intent(IN) :: orb_mvelp_in
real(r8), intent(IN) :: orb_obliq_in
real(r8), intent(IN) :: orb_lambm0_in
real(r8), intent(IN) :: orb_mvelpp_in
real(r8), intent(IN) :: orb_obliqr_in

type(eshr_timemgr_clockInfoType) :: ClockInfo
integer :: orb_iyear_AD
integer :: orb_iyear_AD_in
real(r8) :: orb_eccen, orb_mvelp, orb_obliq
real(r8) :: orb_obliqr
real(r8) :: orb_lambm0
real(r8) :: orb_mvelpp

call eshr_timemgr_clockGet(     clock, Info=clockInfo )
call eshr_timemgr_clockInfoGet( clockInfo, orb_eccen=orb_eccen, orb_obliq=orb_obliq, &
                                          orb_mvelp=orb_mvelp, orb_obliqr=orb_obliqr, &
                                          orb_mvelpp=orb_mvelpp, orb_lambm0=orb_lambm0 )
if ( orb_eccen /= orb_eccen_in )then
  call shr_sys_abort( 'orb_eccen did NOT get set properly' )
end if
if ( orb_mvelp /= orb_mvelp_in )then
  call shr_sys_abort( 'orb_mvelp did NOT get set properly' )
end if
if ( orb_obliq /= orb_obliq_in )then
  call shr_sys_abort( 'orb_obliq did NOT get set properly' )
end if
if ( orb_obliqr /= orb_obliqr_in )then
  call shr_sys_abort( 'orb_obliqr did NOT get set properly' )
end if
if ( orb_lambm0 /= orb_lambm0_in )then
  call shr_sys_abort( 'orb_lambm0 did NOT get set properly' )
end if
end subroutine verify_orbit

end program test_eshr_timemgr

