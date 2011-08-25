!===============================================================================
! SVN $Id$
! SVN $URL$
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: eshr_timemgr_mod --- Time-manager module
!
! !DESCRIPTION:
!
!     A module to create derived types to manage time and clock information 
!     for use with CCSM drivers and models.
!
! Typical usage:
!
! ! Set the default values for the timemgr clockNMLinfo object
! call eshr_timemgr_NMLinfoSetDefault( clockNMLinfo=clockNMLinfo )
! ! Change some values in the clockNMLinfo object
! call eshr_timemgr_NMLinfoPutData( clockNMLinfo, start_ymd=my_start, desc="my_clock" )
! ! Read in the timemgr namelist
! call eshr_timemgr_NMLinfoRead( nlfilename, LogPrint, MPICom, MasterTask, &
!                                clockNMLinfo )
! ! Get info from the clockNMLinfo object
! call eshr_timemgr_NMLinfoGetData( clockNMLinfo, stop_option=stop_option )
! ! Setup the clock
! call eshr_timemgr_clockInitNMLinfo( clockNMLinfo, LogPrint=.true., clockOut=clock )
! Print out the clock
! call eshr_timemgr_clockPrint( clock )
! ! Get info from the clock
! call eshr_timemgr_clockGet( clock, CurrentYMD=ymd, currentTOD=tod )
! ! Advance the clock
! call eshr_timemgr_clockAdvance( clock )
! ! Write restart file out
! if ( eshr_timemgr_clockAlarmIsOnRes( clock ) ) then
!     ! ... figure out restart filename...
!     call eshr_timemgr_clockRestWrite( restart_file, MPICom, MasterTask, clock )
!     ! .. write any other date to driver restart file and archive it
!     ! turn the restart alarm off
!     call eshr_timemgr_clockAlarmOffRest( clock )
! end if
! ! Check if this is the last time-step
! if ( eshr_timemgr_clockIsOnLastStep( clock ) ) ....
!
! ! Now startup a new clock based on a restart file...
! call eshr_timemgr_NMLinfoSetDefault( clockNMLinfo=clockNMLinfo ) ! Set to defaults
! call eshr_timemgr_NMLinfoRead( nlfilename, LogPrint, MPICom, MasterTask, &
!                                clockNMLinfo )     ! Read in namelist
! ! Get restart_file from driver restart pointer file (inputinfo object)
! call eshr_timemgr_NMLinfo( restart_file, clock, LogPrint, MPICom, &
!                            MasterTask, clockNMLinfo )
! ! Setup the clock based on NMLinfo object modified by both Namelist read and 
! ! restart-file
! call eshr_timemgr_clockInitNMLinfo( clockNMLinfo, LogPrint=.true., ClockOut=clock )
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2005-Nov-11 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------

module eshr_timemgr_mod

! !USES:
   use ESMF_Mod
   use SHR_KIND_mod,     only: SHR_KIND_IN, SHR_KIND_R8, SHR_KIND_CS, &
                               SHR_KIND_CL, SHR_KIND_I8
   use shr_orb_mod,      only: SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
   use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
   use shr_ncio_mod,     only: shr_ncio_descripType, shr_ncio_open,   &
                               shr_ncio_close

   implicit none

   private    ! default private

! ! PUBLIC TYPES:

   public :: eshr_timemgr_clockType         ! Wrapped clock object
   public :: eshr_timemgr_NMLinfoType       ! Namelist Information to setup a clock
   public :: eshr_timemgr_clockInfoType     ! Miscellaneous clock information

! !PUBLIC MEMBER FUNCTIONS

   ! --- Namelist Information object methods -----------------------------------
   public :: eshr_timemgr_NMLinfoSetDefault ! Set default values for clockNMLinfo object
   public :: eshr_timemgr_NMLinfoRead       ! Read in the namelist for clock information
   public :: eshr_timemgr_NMLinfoCheck      ! Check clockNMLinfo object for valid settings
   public :: eshr_timemgr_NMLinfoPutData    ! Change clock clockNMLinfo object values
   public :: eshr_timemgr_NMLinfoGetData    ! Get clock clockNMLinfo object values
   public :: eshr_timemgr_NMLinfoRestRead   ! Read time-manager restart information
   ! --- Clock Information object methods --------------------------------------
   public :: eshr_timemgr_clockInfoPutData  ! Change data in the clock-Info object
   public :: eshr_timemgr_clockInfoGet      ! Get values from clock Info object
   public :: eshr_timemgr_clockInfoIsSame   ! Check if two clock Info objects are same
#ifdef SEQ_ESMF
   ! --- Get or put clock info to ESMF State object ----------------------------
   public :: eshr_timemgr_info2EState       ! Put clock-info object on ESMF State
   public :: eshr_timemgr_EState2Info       ! Get clock-info object from ESMF State
#endif
   ! --- Clock object methods --------------------------------------------------
   public :: eshr_timemgr_clockGet          ! Get data from the clock
   public :: eshr_timemgr_clockGetCalDay    ! Get current calendar day
   public :: eshr_timemgr_clockGetRunDays   ! Get running days from reference date
   public :: eshr_timemgr_clockGetPerpYMD   ! Get perpetual date from clock
   public :: eshr_timemgr_clockInitNMLinfo  ! Setup the clock based on InitNML object
   public :: eshr_timemgr_clockInitClock    ! Setup the clock based on an input clock
   public :: eshr_timemgr_clockAdvance      ! Advance the clock
   public :: eshr_timemgr_clockStopAtDayEnd ! Set the clock to stop at end of this day
   public :: eshr_timemgr_clockAlarmOnRest  ! Turn the restart alarm on
   public :: eshr_timemgr_clockAlarmIsOnRes ! Return true if the restart alarm is on
   public :: eshr_timemgr_clockIsOnLastStep ! Return true if this is the last time-step
   public :: eshr_timeMgr_curTimeLEstopTime ! Return true if current time is <= stop time
   public :: eshr_timemgr_clockRestWrite    ! Write restart information out to file
   public :: eshr_timemgr_clockAlarmOffRest ! Turn the restart alarm off
   public :: eshr_timemgr_clockPrint        ! Print clock information 
   public :: eshr_timemgr_clockClocksInSync ! Check that two clocks are in sync...
   public :: eshr_timemgr_clockDateInSync   ! Check that input date and clock in sync
   public :: eshr_timemgr_clockPutData      ! Put object data into the clock object
   ! --- ESMF Time object methods ----------------------------------------------
   public :: eshr_timemgr_ETimeInit         ! Set a ESMF time-object
   public :: eshr_timemgr_ETimeGetDataYMD   ! Get YYYYMMDD from an ESMF time instance
   public :: eshr_timemgr_ETimeGetCalDay    ! Get calendar day from an ESMF time instance
   public :: eshr_timemgr_ETimeGetDiff      ! Get difference between two dates

! !PUBLIC DATA MEMBERS:

   character(len=*), public, parameter :: eshr_timemgr_noLeap    = "NO_LEAP"
   character(len=*), public, parameter :: eshr_timemgr_gregorian = "GREGORIAN"
   character(len=*), public, parameter  :: &
              eshr_timemgr_rest_optNone      = "none"
   character(len=*), public, parameter  :: &
              eshr_timemgr_rest_optNSteps    = "nsteps"
   character(len=*), public, parameter  :: &
              eshr_timemgr_rest_optNDays     = "ndays"
   character(len=*), public, parameter  :: &
              eshr_timemgr_rest_optNMonths   = "nmonths"
   character(len=*), public, parameter  :: &
              eshr_timemgr_rest_optMonthly   = "monthly"
   character(len=*), public, parameter  :: &
              eshr_timemgr_rest_optNYears    = "nyears"
   character(len=*), public, parameter  :: &
              eshr_timemgr_rest_optYearly    = "yearly"
   character(len=*), public, parameter  :: &
              eshr_timemgr_rest_optEnd       = "end"
   character(len=*), public, parameter  :: &
              eshr_timemgr_stop_optionNSteps = "nsteps"
   character(len=*), public, parameter  :: &
              eshr_timemgr_stop_optionNDays  = "ndays"
   character(len=*), public, parameter  :: &
              eshr_timemgr_stop_optionNMons  = "nmonths"
   character(len=*), public, parameter  :: &
              eshr_timemgr_stop_optionNYears = "nyears"
   character(len=*), public, parameter  :: &
              eshr_timemgr_stop_optionDate   = "date"
   character(len=*), public, parameter  :: eshr_timemgr_restartAlarmName = &
              "CCSM_Restart-Alarm"

!EOP

   ! Private member functions:

   ! --- Namelist Information object methods -----------------------------------
   private :: eshr_timemgr_NMLinfoPrint      ! Print the clockNMLinfo object
   private :: eshr_timemgr_NMLinfoOrbReconcl ! Reconcile orbit mode in clockNMLinfo
   ! --- Clock information object methods --------------------------------------
   private :: eshr_timemgr_clockInfoPrint    ! Print the info object
   private :: eshr_timemgr_clockInfoIsMaster ! Return true if this is Master Sync clock
   private :: eshr_timemgr_clockInfoNoRest   ! Return true if no restarts should be done
   private :: eshr_timemgr_clockInfoIsPerpet ! Return true if in perpetual mode
   private :: eshr_timemgr_clockInfoRestRead ! Read clockInfo restart information
   private :: eshr_timemgr_clockInfoRestWrit ! Write clockInfo restart information out
   private :: eshr_timemgr_clockInfoInitRest ! Setup clock-info restart file info
   private :: eshr_timemgr_clockInfoSetOrb   ! Set the orital mode
   private :: eshr_timemgr_clockInfoInit     ! Initialize the clockInfo object
   ! --- Clock object methods --------------------------------------------------
   private :: eshr_timemgr_clockAlarmRstInit ! Set restart alarm from clockNMLinfo object
   private :: eshr_timemgr_clockSetupRest    ! Setup restart file info
   private :: eshr_timemgr_clockGetESMFClock ! Get the ESMF Clock
   private :: eshr_timemgr_clockGetInfo      ! Get the info object
   ! --- ESMF Clock object methods ---------------------------------------------
   private :: eshr_timemgr_EClockInit        ! Initialize the ESMF clock
   ! --- ESMF error return code object methods ---------------------------------
   private :: eshr_timemgr_ErCodeCheck       ! Check the ESMF return code
   ! --- timemgr methods (not dependent on input object) -----------------------
   private :: eshr_timemgr_initCalendar      ! Initialize the calendar
   !
   ! Information included in a CCSM clock
   !
   integer, parameter :: NInfoVars = 10   ! Variables to save on the restart file
   type eshr_timemgr_clockInfoType
      ! CCSM time information: perpetual mode info and orbital information
      private            ! This is an opaque type
      ! ------------------------------------------------------------------------
      ! ----------- Information saved on restart -------------------------------
      ! ------------------------------------------------------------------------
      ! Description of this clock
      character(len=SHR_KIND_CL) :: desc
      ! If this is a master synchronization clock or not
      logical                    :: MasterSyncClock
      ! Calendar type
      character(len=SHR_KIND_CS) :: calendar
      ! Perpetual date information
      logical              :: perpetual_run = .false.   ! If running in perpetual mode
      type(ESMF_Time)      :: perpetual_time            ! Time/date of perpetual mode
      ! Orbital information to set: either year or all of rest
      integer(SHR_KIND_IN) :: orb_mode     = SHR_ORB_UNDEF_INT  ! Orbit mode
      integer(SHR_KIND_IN) :: orb_iyear_AD = SHR_ORB_UNDEF_INT  ! year for orbit
      real(SHR_KIND_R8)    :: orb_obliq    = SHR_ORB_UNDEF_REAL ! Obliquity of orbit
      real(SHR_KIND_R8)    :: orb_eccen    = SHR_ORB_UNDEF_REAL ! Eccentricity of orbit
      real(SHR_KIND_R8)    :: orb_mvelp    = SHR_ORB_UNDEF_REAL ! Locatn of vernal equinox
      ! ------------------------------------------------------------------------
      ! ----------- Information not saved on restart ---------------------------
      ! ------------------------------------------------------------------------
      ! If wish to turn restarts off 
      ! (only do restarts if explicitly call clockAlarmOnRest)
      logical              :: NoRestarts   = .false.
      ! Orbital information derived from above
      real(SHR_KIND_R8)    :: orb_obliqr   = SHR_ORB_UNDEF_REAL ! Obliquity in radians
      real(SHR_KIND_R8)    :: orb_lambm0   = SHR_ORB_UNDEF_REAL ! Long of perh. radians
      real(SHR_KIND_R8)    :: orb_mvelpp   = SHR_ORB_UNDEF_REAL ! mvelp in radians
      !------------------------------------------------------------------------
      ! Variables that are written out to restart file from above
      !-------------------------------------------------------------------------
      type(shr_ncio_descripType) :: var(NInfoVars)
   end type eshr_timemgr_clockInfoType
   !
   ! ---- List of ClockInfo variables to save to restart file ------
   !
   character(len=*), parameter :: ClockInfoSave(NInfoVars) = (/ &
          "desc           ",      &
          "MasterSyncClock",      &
          "calendar       ",      &
          "perpetual_run  ",      &
          "perpetual_ymd  ",      &
          "orb_mode       ",      &
          "orb_iyear_AD   ",      &
          "orb_obliq      ",      &
          "orb_eccen      ",      &
          "orb_mvelp      " /)
   ! Description of orbit modes 
   character(len=*), private, parameter  :: &
         eshr_timemgr_orb_mode_fixed_yr = &
        'Orb_mode_calculated_based_on_fixed_year'  ! Calculate based on fixed year
   character(len=*), private, parameter  :: &
         eshr_timemgr_orb_mode_fixed = &
        'Orb_mode_fixed_by_orbit_parameters'       ! Use parameters given from user
   integer, parameter :: NOrbList = 2
   character(len=SHR_KIND_CL), private, parameter :: &
        eshr_timemgr_orbModeDescrips = &
        eshr_timemgr_orb_mode_fixed_yr//':'// &
        eshr_timemgr_orb_mode_fixed
   integer, save, private :: &
           eshr_timemgr_orbModeList(NOrbList) = (/1,2/) ! Index values of above list
   !
   ! CCSM Clock, information, time and restart alarm
   !
   integer, parameter :: NClockVars = 12 ! Number of variables to save on restart file
   type eshr_timemgr_clockType
      private            ! This is an opaque type
      ! ------------------------------------------------------------------------
      ! ----------- Information saved on restart -------------------------------
      ! ------------------------------------------------------------------------
      type(ESMF_Clock)                 :: EClock    ! Time information
      type(ESMF_Alarm)                 :: restart   ! Restart alarm
      type(eshr_timemgr_clockInfoType) :: info      ! Misc. Info: perpetual mode, orbit
      ! ------------------------------------------------------------------------
      ! Variables that are written out to restart file from above
      ! ------------------------------------------------------------------------
      type(shr_ncio_descripType) :: var(NClockVars)
   end type eshr_timemgr_clockType
   !
   ! ---- List of Clock variables to save to restart file ------
   !
   character(len=*), parameter :: ClockSave(NClockVars) = (/ &
          "start_ymd            ", &
          "start_tod            ", &
          "ref_ymd              ", &
          "ref_tod              ", &
          "CurrentYMD           ", &
          "CurrentTOD           ", &
          "DTime                ", &
          "RestartIntervalSec   ", &
          "RestartIntervalMonths", &
          "RestartIntervalYears ", &
          "RestartNextAlarmYMD  ", &
          "RestartNextAlarmTOD  " /)
 
   character(len=*), parameter :: &
              prefix = "eshr_timemgr_"  ! Prefix to put on variable on restart file

   !
   ! Namelist values for setting up a CCSM clock
   !
   type eshr_timemgr_NMLinfoType
      private            ! This is an opaque type
      ! Description of this clock
      character(len=SHR_KIND_CL) :: desc
      ! If restart or not
      logical                 :: restart
      ! If this is a master synchronization clock or not
      logical                 :: MasterSyncClock
      ! Calendar to use: NO_LEAP or GREGORIAN
      character(SHR_KIND_CS)  :: calendar
      ! Frequency of restarts
      character(SHR_KIND_CS)  :: restart_option  ! Restart type units
      integer(SHR_KIND_IN)    :: restart_n       ! number to next restart
      ! Time-step for this clock
      integer(SHR_KIND_IN)    :: dtime           ! clock time-step
      ! Coupling intervals of each model component
      integer(SHR_KIND_IN)    :: atm_cpl_dt      ! Atmosphere coupling interval
      integer(SHR_KIND_IN)    :: lnd_cpl_dt      ! Land coupling interval
      integer(SHR_KIND_IN)    :: ice_cpl_dt      ! Sea-ice coupling interval
      integer(SHR_KIND_IN)    :: ocn_cpl_dt      ! Ocean coupling interval
      ! Simulation stop time
      character(SHR_KIND_CS)  :: stop_option     ! Stop type units
      integer(SHR_KIND_IN)    :: stop_n          ! Number until stop
      integer(SHR_KIND_IN)    :: stop_ymd        ! YYYYMMDD date to stop on
      integer(SHR_KIND_IN)    :: stop_tod        ! Time of day to stop on (if date option)
      integer(SHR_KIND_IN)    :: stop_final_ymd  ! Final date (in YYYYMMDD) to run to
      ! Simulation start time
      integer(SHR_KIND_IN)    :: start_ymd       ! Start date (YYYYMMDD)
      integer(SHR_KIND_IN)    :: start_tod       ! Start time of day (seconds)
      ! Simulation reference time
      integer(SHR_KIND_IN)    :: ref_ymd         ! Reference date (YYYYMMDD)
      integer(SHR_KIND_IN)    :: ref_tod         ! Reference time of day (seconds)
      ! Simulation current time
      integer(SHR_KIND_IN)    :: CurrentYMD      ! Current date (YYYYMMDD)
      integer(SHR_KIND_IN)    :: CurrentTOD      ! Current time of day (seconds)
      ! Simulation current time
      integer(SHR_KIND_IN)    :: NextRestYMD     ! Next restart date (YYYYMMDD)
      integer(SHR_KIND_IN)    :: NextRestTOD     ! Next restart time of day (seconds)
      ! Perpetual date to hold to
      logical                 :: perpetual_run   ! Flag if this is a perpetual date or not
      integer(SHR_KIND_IN)    :: perpetual_ymd   ! Perpetual date (YYYYMMDD)
      ! Orbital information to set: either year or all of rest
      logical                 :: orb_notSet      ! Flag if orbit is not set here
      integer(SHR_KIND_IN)    :: orb_iyear_AD    ! Orbit year
      real(SHR_KIND_R8)       :: orb_obliq       ! Orbital obliquity
      real(SHR_KIND_R8)       :: orb_eccen       ! Orbital eccentricity
      real(SHR_KIND_R8)       :: orb_mvelp       ! Orbital location of vernal equinox
      ! Items to overRide restart_file on namelist
      character(SHR_KIND_CL)  :: restart_file_TMGoverRide
   end type eshr_timemgr_NMLinfoType
   ! List of values read in from restart_file that can be over-ridden on namelist
   character(len=*), private, parameter :: restart_file_overRideList = &
              "restart_option" &  ! Restart frequency type
           //":restart_n"         ! Restart frequency number
   !
   ! Private local data
   !
   type(ESMF_Calendar), save, target :: eshr_timemgr_cal        ! calendar
   logical, save :: eshr_timemgr_setCalendar = .false.          ! if calendar has been set
   integer, parameter :: eshr_timemgr_unInitInteger = -99999999 ! value if uninitialized
   integer, parameter :: eshr_timemgr_secPerDay = 86400         ! Seconds per day
   character(len=*), private, parameter :: ClockInfoEStateName = "ClockInfoEState"
   character(len=*), private, parameter :: LogPrefix = "(eshr_timemgr) "
   character(len=*), private, parameter :: FA  = "( '"//LogPrefix//"', A)"
   character(len=*), private, parameter :: FAA = "( '"// &
     LogPrefix//"', A, T45, ' = ', A )"
    character(len=*), parameter :: F0Date = &
       "('"//Logprefix//"', A, T45, ' = ', I8.8, ' (YYYYMMDD) ', I5.5, ' (sec)' )"
    character(len=*), parameter :: F0I = &
       "('"//LogPrefix//"', A, T45, ' = ', I8, ' ', A )"
    character(len=*), parameter :: F0I8 = &
       "('"//LogPrefix//"', A, T45, ' = ', I8.8 )"
    character(len=*), parameter :: F0A    = "( '"//LogPrefix//"',A,T45,' = ', A)"
    character(len=*), parameter :: F0F    = &
       "( '"//LogPrefix//"',A,T45,' = ', F15.10, ' ', A)"

!===============================================================================
CONTAINS
!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_NMLinfoSetDefault -- Set defaults for clockNMLinfo object.
!   
! !DESCRIPTION:
!   
!     Set the defaults for the input clockNMLinfo object.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_NMLinfoSetDefault( perpetual_run, perpetual_ymd, &
                                           clockNMLinfo )
! !USES:

   USE shr_orb_mod,  ONLY: SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   logical, intent(IN), optional :: perpetual_run    ! Flag for perpetual mode or not
   integer, intent(IN), optional :: perpetual_ymd    ! perpetual mode default date to use
   type(eshr_timemgr_NMLinfoType), intent(OUT) :: clockNMLinfo ! Output clock clockNMLinfo

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_NMLinfoDefaults) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    clockNMLinfo%desc             = ' '
    clockNMLinfo%restart          = .false.
    clockNMLinfo%MasterSyncClock  = .true.
    clockNMLinfo%restart_n        = 1
    clockNMLinfo%restart_option   = eshr_timemgr_rest_optYearly
    clockNMLinfo%dtime            = 0
    clockNMLinfo%atm_cpl_dt       = 0
    clockNMLinfo%lnd_cpl_dt       = 0
    clockNMLinfo%ice_cpl_dt       = 0
    clockNMLinfo%ocn_cpl_dt       = 0
    clockNMLinfo%calendar         = eshr_timemgr_noLeap
    clockNMLinfo%stop_option      = ' '
    clockNMLinfo%stop_n           = 0
    clockNMLinfo%stop_ymd         = 0
    clockNMLinfo%stop_final_ymd   = 99991231
    clockNMLinfo%stop_tod         = 0
    clockNMLinfo%start_ymd        = 0
    clockNMLinfo%start_tod        = 0
    clockNMLinfo%CurrentYMD       = 0
    clockNMLinfo%CurrentTOD       = 0
    clockNMLinfo%NextRestYMD      = 0
    clockNMLinfo%NextRestTOD      = 0
    clockNMLinfo%ref_ymd          = 0
    clockNMLinfo%ref_tod          = 0
    if ( present(perpetual_run) )then
       if ( .not. present(perpetual_ymd) )then
          call shr_sys_abort( subname//':: perpetual_run AND perpetual_ymd '//&
                              'should be given' )
       end if
       clockNMLinfo%perpetual_run      = perpetual_run
    else
       clockNMLinfo%perpetual_run      = .false.
    end if
    if ( present(perpetual_ymd) )then
       if ( .not. present(perpetual_run) )then
          call shr_sys_abort( subname//':: perpetual_run AND perpetual_ymd '//&
                              'should be given' )
       end if
       clockNMLinfo%perpetual_ymd         = perpetual_ymd
    else
       clockNMLinfo%perpetual_ymd         = 0
    end if
    clockNMLinfo%orb_notSet               = .false.
    clockNMLinfo%orb_obliq                = SHR_ORB_UNDEF_REAL
    clockNMLinfo%orb_eccen                = SHR_ORB_UNDEF_REAL
    clockNMLinfo%orb_mvelp                = SHR_ORB_UNDEF_REAL
    clockNMLinfo%orb_iyear_AD             = SHR_ORB_UNDEF_INT
    clockNMLinfo%restart_file_TMGoverRide = ' '

END SUBROUTINE eshr_timemgr_NMLinfoSetDefault

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_NMLinfoRead -- Set clock clockNMLinfo values from namelist
!   
! !DESCRIPTION:
!   
!     Read in the clock clockNMLinfo values from the namelist.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_NMLinfoRead( nlfilename, LogPrint, MPICom, MasterTask, &
                                     clockNMLinfoOut )
! !USES:

   USE shr_file_mod, ONLY : shr_file_getUnit, shr_file_freeUnit
   USE shr_orb_mod,   ONLY : SHR_ORB_UNDEF_REAL, SHR_ORB_UNDEF_INT
   USE shr_mpi_mod,   ONLY : shr_mpi_bcast

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(len=*),            intent(IN)    :: nlfilename  ! Namelist filename
   logical,                     intent(IN)    :: LogPrint    ! If should print info to log
   integer, optional,           intent(IN)    :: MPICom      ! MPI communicator
   logical, optional,           intent(IN)    :: MasterTask  ! If root MPI task or not
   type(eshr_timemgr_NMLinfoType), intent(INOUT) :: clockNMLinfoOut ! Output clock setup

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_NMLinfoRead) '
    logical                 :: MasterTask2           ! If root MPI task or not
    character(SHR_KIND_CS)  :: restart_option        ! Restart option units
    integer(SHR_KIND_IN)    :: restart_n             ! Number until restart interval
    integer(SHR_KIND_IN)    :: dtime                 ! Time-step of this clock
    integer(SHR_KIND_IN)    :: atm_cpl_dt            ! Atmosphere coupling interval
    integer(SHR_KIND_IN)    :: lnd_cpl_dt            ! Land coupling interval
    integer(SHR_KIND_IN)    :: ice_cpl_dt            ! Sea-Ice coupling interval
    integer(SHR_KIND_IN)    :: ocn_cpl_dt            ! Ocean coupling interval
    character(SHR_KIND_CS)  :: calendar              ! Calendar type
    character(SHR_KIND_CL)  :: restart_file_TMGoverRide ! Override items from restart file
    character(SHR_KIND_CS)  :: stop_option           ! Stop option units
    integer(SHR_KIND_IN)    :: stop_n                ! Number until stop
    integer(SHR_KIND_IN)    :: stop_ymd              ! Stop date (YYYYMMDD)
    integer(SHR_KIND_IN)    :: stop_tod              ! Stop time-of-day
    integer(SHR_KIND_IN)    :: stop_final_ymd        ! Final stop date (YYYYMMDD)
    integer(SHR_KIND_IN)    :: start_ymd             ! Start date (YYYYMMDD)
    integer(SHR_KIND_IN)    :: start_tod             ! Start time of day (seconds)
    integer(SHR_KIND_IN)    :: ref_ymd               ! Reference date (YYYYMMDD)
    integer(SHR_KIND_IN)    :: ref_tod               ! Reference time of day (seconds)
    logical                 :: perpetual_run         ! Flag if perpetual date
    integer(SHR_KIND_IN)    :: perpetual_ymd         ! Perpetual date (YYYYMMDD)
    real(SHR_KIND_R8)       :: orb_obliq             ! Orbital obliquity
    real(SHR_KIND_R8)       :: orb_eccen             ! Orbital eccentricity
    real(SHR_KIND_R8)       :: orb_mvelp             ! Location of vernal equinox
    integer(SHR_KIND_IN)    :: orb_iyear_AD          ! Orbital year
    integer(SHR_KIND_IN)    :: nlUnit                ! Namelist unit number
    integer(SHR_KIND_IN)    :: ierr                  ! Return code

    namelist /timemgr_inparm/  &   
         calendar, stop_option, stop_n, stop_ymd, stop_tod, stop_final_ymd, &
         restart_option, restart_n, start_ymd, start_tod, ref_ymd, ref_tod, &
         perpetual_run, perpetual_ymd, atm_cpl_dt, orb_obliq, orb_eccen,    &
         orb_mvelp, orb_iyear_AD, restart_file_TMGoverRide

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(MasterTask) )then
       MasterTask2 = MasterTask
    else
       MasterTask2 = .true.
    end if
    !---------------------------------------------------------------------------
    ! Set initial values to input object values
    !---------------------------------------------------------------------------
    restart_n             = clockNMLinfoOut%restart_n
    restart_option        = clockNMLinfoOut%restart_option
    stop_n                = clockNMLinfoOut%stop_n
    stop_option           = clockNMLinfoOut%stop_option
    stop_ymd              = clockNMLinfoOut%stop_ymd
    stop_tod              = clockNMLinfoOut%stop_tod
    stop_final_ymd        = clockNMLinfoOut%stop_final_ymd
    dtime                 = clockNMLinfoOut%dtime
    atm_cpl_dt            = clockNMLinfoOut%atm_cpl_dt
    lnd_cpl_dt            = clockNMLinfoOut%lnd_cpl_dt
    ice_cpl_dt            = clockNMLinfoOut%ice_cpl_dt
    ocn_cpl_dt            = clockNMLinfoOut%ocn_cpl_dt
    perpetual_ymd         = clockNMLinfoOut%perpetual_ymd
    start_ymd             = clockNMLinfoOut%start_ymd
    start_tod             = clockNMLinfoOut%start_tod
    ref_ymd               = clockNMLinfoOut%ref_ymd
    ref_tod               = clockNMLinfoOut%ref_tod
    perpetual_run         = clockNMLinfoOut%perpetual_run
    orb_iyear_AD          = clockNMLinfoOut%orb_iyear_AD
    orb_obliq             = clockNMLinfoOut%orb_obliq
    orb_eccen             = clockNMLinfoOut%orb_eccen
    orb_mvelp             = clockNMLinfoOut%orb_mvelp
    calendar              = clockNMLinfoOut%calendar
    restart_file_TMGoverRide = clockNMLinfoOut%restart_file_TMGoverRide
    !---------------------------------------------------------------------------
    ! Read in the namelist
    !---------------------------------------------------------------------------
    if ( MasterTask2 )then
       nlUnit = shr_file_getUnit()
       if ( LogPrint ) write(6,FA) 'Read in timemgr_inparm namelist from: '// &
                                    trim(nlfilename)
       open( nlUnit, file=trim(nlfilename), status='old' )
       ierr = 1
       do while ( ierr /= 0 )
          read(nlUnit,timemgr_inparm,iostat=ierr)
          if (ierr < 0) then
             call shr_sys_abort( subname//':: namelist read returns an'// &
                                 ' end of file or end of record condition' )
          end if
       end do
       close(nlUnit)            
       call shr_file_freeUnit( nlUnit )
    end if
    !---------------------------------------------------------------------------
    ! If MPI broadcast values to all tasks
    !---------------------------------------------------------------------------
    if ( present(MPICOM) )then
       call shr_mpi_bcast( restart_n,                MPICom )
       call shr_mpi_bcast( restart_option,           MPICom )
       call shr_mpi_bcast( stop_n,                   MPICom )
       call shr_mpi_bcast( stop_option,              MPICom )
       call shr_mpi_bcast( stop_ymd,                 MPICom )
       call shr_mpi_bcast( stop_tod,                 MPICom )
       call shr_mpi_bcast( stop_final_ymd,           MPICom )
       call shr_mpi_bcast( atm_cpl_dt,               MPICom )
       call shr_mpi_bcast( lnd_cpl_dt,               MPICom )
       call shr_mpi_bcast( ice_cpl_dt,               MPICom )
       call shr_mpi_bcast( ocn_cpl_dt,               MPICom )
       call shr_mpi_bcast( perpetual_ymd,            MPICom )
       call shr_mpi_bcast( start_ymd,                MPICom )
       call shr_mpi_bcast( start_tod,                MPICom )
       call shr_mpi_bcast( ref_ymd,                  MPICom )
       call shr_mpi_bcast( ref_tod,                  MPICom )
       call shr_mpi_bcast( perpetual_run,            MPICom )
       call shr_mpi_bcast( orb_iyear_AD,             MPICom )
       call shr_mpi_bcast( orb_obliq,                MPICom )
       call shr_mpi_bcast( orb_eccen,                MPICom )
       call shr_mpi_bcast( orb_mvelp,                MPICom )
       call shr_mpi_bcast( calendar,                 MPICom )
       call shr_mpi_bcast( restart_file_TMGoverRide, MPICom )
    end if

    ! --- Error check that perpetual_run not changed if it's already true by default ----
    if ( clockNMLinfoOut%perpetual_run .and. (.not. perpetual_run) )then
       call shr_sys_abort( subname//': perpetual_run can NOT be changed on '//&
                           'namelist if it is set to  true by default' )
    end if
    !---------------------------------------------------------------------------
    ! Time-manager clockNMLinfoOut derived type
    !---------------------------------------------------------------------------
    clockNMLinfoOut%restart_n             = restart_n
    clockNMLinfoOut%restart_option        = restart_option
    clockNMLinfoOut%stop_n                = stop_n
    clockNMLinfoOut%stop_option           = stop_option
    clockNMLinfoOut%stop_ymd              = stop_ymd
    clockNMLinfoOut%stop_tod              = stop_tod
    clockNMLinfoOut%stop_final_ymd        = stop_final_ymd
    clockNMLinfoOut%atm_cpl_dt            = atm_cpl_dt
    clockNMLinfoOut%lnd_cpl_dt            = atm_cpl_dt ! Copy atm coupling time into land
    clockNMLinfoOut%ice_cpl_dt            = atm_cpl_dt ! Copy atm coupling time into ice
    clockNMLinfoOut%ocn_cpl_dt            = atm_cpl_dt ! Copy atm coupling time into ocean
    clockNMLinfoOut%perpetual_ymd         = perpetual_ymd
    clockNMLinfoOut%start_ymd             = start_ymd
    clockNMLinfoOut%start_tod             = start_tod
    clockNMLinfoOut%ref_ymd               = ref_ymd
    clockNMLinfoOut%ref_tod               = ref_tod
    clockNMLinfoOut%perpetual_run         = perpetual_run
    clockNMLinfoOut%orb_iyear_AD          = orb_iyear_AD
    clockNMLinfoOut%orb_obliq             = orb_obliq
    clockNMLinfoOut%orb_eccen             = orb_eccen
    clockNMLinfoOut%orb_mvelp             = orb_mvelp
    clockNMLinfoOut%calendar              = calendar
    clockNMLinfoOut%restart_file_TMGoverRide = restart_file_TMGoverRide

END SUBROUTINE eshr_timemgr_NMLinfoRead

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_NMLinfoCheck -- Check the clock clockNMLinfoOut values
!   
! !DESCRIPTION:
!   
!     Check the clock clockNMLinfoOut values to make sure they are correct for both
! validity and for consistency.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_NMLinfoCheck( clockNMLinfo )

! !USES:

   use shr_string_mod,   only: shr_string_listIntersect, &
                               shr_string_listGetIndexF

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_NMLinfoType),  intent(IN) :: clockNMLinfo ! clock clockNMLinfo

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_NMLinfoCheck) '

    ! ---- Intersection of over-ride list to master overRide list ------
    character(len=SHR_KIND_CL) :: OverrideList 
    integer :: rc                              ! Return code
    integer :: list                            ! List index into overRide list

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- Override values ------------------------------------------------------
    if ( clockNMLinfo%restart .and. &
        (len_trim(clockNMLinfo%restart_file_TMGoverRide) > 0) )then
       call shr_string_listIntersect(trim(clockNMLinfo%restart_file_TMGoverRide), &
                                     restart_file_overRideList, &
                                     OverRideList, rc )
       if ( trim(OverRideList) /= trim(clockNMLinfo%restart_file_TMGoverRide) )then
          write(6,FA)'ERROR: the only values that can be overridden are: '// &
                    trim(restart_file_overRideList)
          write(6,FA)'List from namelist: '// &
                      trim(clockNMLinfo%restart_file_TMGoverRide)
          call shr_sys_abort( subname//': list of values to overRide '//&
                              'includes values that can NOT' &
                              //" be overridden" )
       end if
       ! --- If overRide restart_option, must also over-ride restart_n ---------
       list  = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                         "restart_option" )
       if ( (list /= 0) .and. &
       ((trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optNone)   .or. &
        (trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optEnd)    .or. &
        (trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optYearly) .or. &
         trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optMonthly) )then
          list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                           "restart_n" )
          if ( list == 0 ) &
             call shr_sys_abort( subname//': if restart_option overridden, '//&
                                 'so must restart_n' )
       end if
    end if
    if ( .not. clockNMLinfo%restart .and. &
         (len_trim(clockNMLinfo%restart_file_TMGoverRide) > 0) )then
       call shr_sys_abort( subname//': if not restart there is no need to '//&
                           'set the restart_file_TMGoverRide option' )
    end if
    ! --- Restart frequency ----------------------------------------------------
    if ( clockNMLinfo%restart_n < 0 )then
       call shr_sys_abort( subname//': restart_n can not be negative' )
    end if
    if (  trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optNone    &
    .and. trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optEnd     &
    .and. trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optNSteps  &
    .and. trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optNDays   &
    .and. trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optNMonths &
    .and. trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optNYears  &
    .and. trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optMonthly &
    .and. trim(clockNMLinfo%restart_option) /= eshr_timemgr_rest_optYearly  &
    )then
       call shr_sys_abort( subname//': restart_option can only be set '// &
                           'to: '//&
                            eshr_timemgr_rest_optNone   //', '//&
                            eshr_timemgr_rest_optNSteps //', '//&
                            eshr_timemgr_rest_optNDays  //', '//&
                            eshr_timemgr_rest_optNMonths//', '//&
                            eshr_timemgr_rest_optNYears //', '//&
                            eshr_timemgr_rest_optMonthly//', '//&
                            eshr_timemgr_rest_optYearly //', or '//&
                            eshr_timemgr_rest_optEnd )
    end if
    if ( (trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNSteps  &
    .or.  trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNDays   &
    .or.  trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNMonths &
    .or.  trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNYears) &
    .and. (clockNMLinfo%restart_n <= 0) )then
       call shr_sys_abort( subname//': restart_n required to be set '//&
                           'if restart_option: ' //                  &
                            eshr_timemgr_rest_optNSteps //', '//&
                            eshr_timemgr_rest_optNDays  //', '//&
                            eshr_timemgr_rest_optNMonths//', or '//&
                            eshr_timemgr_rest_optNYears )
    end if
    ! --- Stop time ------------------------------------------------------------
    if ( trim(clockNMLinfo%stop_option) /= eshr_timemgr_stop_optionNSteps    .and. &
         trim(clockNMLinfo%stop_option) /= eshr_timemgr_stop_optionNDays     .and. &
         trim(clockNMLinfo%stop_option) /= eshr_timemgr_stop_optionNMons   .and. &
         trim(clockNMLinfo%stop_option) /= eshr_timemgr_stop_optionNYears    .and. &
         trim(clockNMLinfo%stop_option) /= eshr_timemgr_stop_optionDate )then
       call shr_sys_abort( subname//': stop_option can only be set to: '// &
                            eshr_timemgr_stop_optionNSteps //', '//&
                            eshr_timemgr_stop_optionNDays  //', '//&
                            eshr_timemgr_stop_optionNMons//', '//&
                            eshr_timemgr_stop_optionNYears //', or '//&
                            eshr_timemgr_stop_optionDate   )
    end if
    if ( (trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNSteps    .or. &
          trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNDays     .or. &
          trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNMons   .or. &
          trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNYears) &
          .and. (clockNMLinfo%stop_n <= 0) )then
       call shr_sys_abort( subname//': stop_n required to be set if '// &
                           'stop_option: '//&
                            eshr_timemgr_stop_optionNSteps //', '//&
                            eshr_timemgr_stop_optionNDays  //', '//&
                            eshr_timemgr_stop_optionNMons//', or '//&
                            eshr_timemgr_stop_optionNYears )
    end if
    if ( (trim(clockNMLinfo%stop_option) /= eshr_timemgr_stop_optionDate) &
    .and. (clockNMLinfo%stop_ymd /= 0) )then
       call shr_sys_abort( subname//': can only set stop_ymd if '// &
                           'stop_option set to '//eshr_timemgr_stop_optionDate )
    end if
    if ( clockNMLinfo%stop_n < 0 )then
       call shr_sys_abort( subname//': stop_n can not be negative' )
    end if
    if ( (trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionDate) .and. &
         ((clockNMLinfo%stop_ymd < 101) .or. (clockNMLinfo%stop_ymd > 99991231))&
    )then
       call shr_sys_abort( subname//': bad date for stop_ymd' )
    end if
    ! --- Calendar  ------------------------------------------------------------
    if ( trim(clockNMLinfo%calendar) /= eshr_timemgr_noLeap .and. &
    trim(clockNMLinfo%calendar) /= eshr_timemgr_Gregorian )then
       call shr_sys_abort( subname//': calendar must be: '// &
                           eshr_timemgr_noLeap//' or '// &
                           eshr_timemgr_Gregorian )
    end if
    ! --------------------------------------------------------------------------
    ! Check these if NOT restart type
    ! --------------------------------------------------------------------------
    if ( .not. clockNMLinfo%restart )then
       ! --- Start time date ---------------------------------------------------
       if ( clockNMLinfo%start_ymd == 0 )then
          call shr_sys_abort( subname//': Must set start_ymd' )
       end if
       if ( (clockNMLinfo%start_ymd < 101) .or. (clockNMLinfo%start_ymd > 99991231)&
       )then
          call shr_sys_abort( subname//': bad date for start_ymd' )
       end if
       ! --- Coupling intervals ------------------------------------------------
       if ( clockNMLinfo%dtime == 0 )then
          if ( clockNMLinfo%atm_cpl_dt <= 0 )then
             call shr_sys_abort( subname//': Required to set atm_cpl_dt' )
          end if
          if ( clockNMLinfo%atm_cpl_dt /= clockNMLinfo%lnd_cpl_dt )then
             call shr_sys_abort( subname//': Required to set lnd_cpl_dt '//&
                                 'to atm_cpl_dt value' )
          end if
          if ( clockNMLinfo%atm_cpl_dt /= clockNMLinfo%ice_cpl_dt )then
             call shr_sys_abort( subname//': Required to set ice_cpl_dt '//&
                                 'to atm_cpl_dt value' )
          end if
          if ( clockNMLinfo%atm_cpl_dt /= clockNMLinfo%ocn_cpl_dt )then
             call shr_sys_abort( subname//': Required to set ocn_cpl_dt '//&
                                 'to atm_cpl_dt value' )
          end if
       end if
       ! --- Perpetual mode ----------------------------------------------------
       if ( clockNMLinfo%perpetual_run .and. (clockNMLinfo%perpetual_ymd == 0) &
       )then
          call shr_sys_abort( subname//': If set perpetual_run, must '//&
                              'also set perpetual_ymd' )
       end if
       if ( (.not. clockNMLinfo%perpetual_run) .and. &
       (clockNMLinfo%perpetual_ymd /= 0) )then
          call shr_sys_abort( subname//': If set perpetual_ymd, must '//&
                              'also set perpetual_run to true' )
       end if
       if ( (clockNMLinfo%perpetual_ymd /= 0) .and. &
       ((clockNMLinfo%perpetual_ymd < 101) .or. &
       (clockNMLinfo%perpetual_ymd > 99991231)) )then
          call shr_sys_abort( subname//': bad date for perpetual_ymd' )
       end if
       ! --- Orbital information -----------------------------------------------
       ! --- Check that orbit information set correctly ----
       if ( .not. clockNMLinfo%orb_notSet )then
          if ( clockNMLinfo%orb_iyear_AD == SHR_ORB_UNDEF_INT ) then
             if ( (clockNMLinfo%orb_obliq == SHR_ORB_UNDEF_REAL) &
             .or. (clockNMLinfo%orb_eccen == SHR_ORB_UNDEF_REAL) .or. &
                  (clockNMLinfo%orb_mvelp == SHR_ORB_UNDEF_REAL) ) then
                 call shr_sys_abort( subname//': Need to set iyear_AD, or '// &
                                     'all of obliq, mvelp, AND eccen' )
             end if
          end if
          if ( (clockNMLinfo%orb_iyear_AD /= SHR_ORB_UNDEF_INT) .and. &
               ((clockNMLinfo%orb_obliq /= SHR_ORB_UNDEF_REAL) .or.   &
               (clockNMLinfo%orb_eccen /= SHR_ORB_UNDEF_REAL) .or. &
                (clockNMLinfo%orb_mvelp /= SHR_ORB_UNDEF_REAL)) ) then
                 call shr_sys_abort( subname//': Need to set iyear_AD, or '// &
                                     'all of obliq, mvelp, AND eccen NOT both' )
          end if
       ! --- If orbit information NOT set now -- make sure set to UNDEF
       else
          if ( (clockNMLinfo%orb_iyear_AD /= SHR_ORB_UNDEF_INT) .or. &
               (clockNMLinfo%orb_obliq /= SHR_ORB_UNDEF_REAL) .or.   &
               (clockNMLinfo%orb_eccen /= SHR_ORB_UNDEF_REAL) .or. &
               (clockNMLinfo%orb_mvelp /= SHR_ORB_UNDEF_REAL) ) then
                 call shr_sys_abort( subname//': iyear_AD, '// &
                                     'obliq, mvelp, and eccen should NOT be set' )
          end if
       end if
    end if

END SUBROUTINE eshr_timemgr_NMLinfoCheck

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_NMLinfoPutData -- Change clock clockNMLinfo values
!   
! !DESCRIPTION:
!   
!     Change given values in the clock clockNMLinfo. Currently only has options to
!  set a small fraction of the clockNMLinfo information, could be expanded to change all
!  clockNMLinfo information. But, we are only adding methods to change that which we need
!  at this point.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_NMLinfoPutData( clockNMLinfo, orb_iyear_AD, orb_obliq,    &
                                        orb_eccen, orb_mvelp, atm_cpl_dt,         &
                                        lnd_cpl_dt, ocn_cpl_dt, ice_cpl_dt,       &
                                        start_ymd, start_tod, stop_option,        &
                                        stop_n, stop_ymd, restart_n,              &
                                        restart_option, restart_file_TMGoverRide, &
                                        calendar, perpetual_run, orb_notSet,      &
                                        perpetual_ymd, dtime, ref_ymd,            &
                                        ref_tod, stop_tod, CurrentYMD,            &
                                        CurrentTOD, desc, MasterSyncClock,        &
                                        restart )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_NMLinfoType), intent(INOUT) :: clockNMLinfo     ! clock setup
    integer,              intent(IN), optional :: orb_iyear_AD! Orbital year
    real(SHR_KIND_R8),    intent(IN), optional :: orb_obliq  ! Orbital obliquity
    real(SHR_KIND_R8),    intent(IN), optional :: orb_eccen  ! Orbital eccentricity
    real(SHR_KIND_R8),    intent(IN), optional :: orb_mvelp  ! Orbital vernal equinox pos.
    logical,              intent(IN), optional :: orb_notSet ! Orbital values not set
    integer(SHR_KIND_IN), intent(IN), optional :: dtime      ! Clock time-step
    integer(SHR_KIND_IN), intent(IN), optional :: atm_cpl_dt ! Atmos. couple interval (s)
    integer(SHR_KIND_IN), intent(IN), optional :: lnd_cpl_dt ! Land couple interval (s)
    integer(SHR_KIND_IN), intent(IN), optional :: ice_cpl_dt ! Sea-Ice couple interval (s)
    integer(SHR_KIND_IN), intent(IN), optional :: ocn_cpl_dt ! Ocean couple interval (s)
    integer(SHR_KIND_IN), intent(IN), optional :: start_ymd  ! Start time (YYYYMMDD)
    integer(SHR_KIND_IN), intent(IN), optional :: start_tod  ! Start time of day (sec)
    integer(SHR_KIND_IN), intent(IN), optional :: ref_ymd    ! Reference time (YYYYMMDD)
    integer(SHR_KIND_IN), intent(IN), optional :: ref_tod    ! Reference time of day (sec)
    integer(SHR_KIND_IN), intent(IN), optional :: CurrentYMD ! Current time (YYYYMMDD)
    integer(SHR_KIND_IN), intent(IN), optional :: CurrentTOD ! Current time of day (sec)
    character(len=*),     intent(IN), optional :: stop_option! Frequency type to stop at
    integer(SHR_KIND_IN), intent(IN), optional :: stop_n     ! No. of units to run
    integer(SHR_KIND_IN), intent(IN), optional :: stop_ymd   ! Date to stop at
    integer(SHR_KIND_IN), intent(IN), optional :: stop_tod   ! time of day to stop at
    character(len=*),     intent(IN), optional :: restart_option ! Frequency type to rest
    integer(SHR_KIND_IN), intent(IN), optional :: restart_n  ! Restart frequency number
    character(len=*),     intent(IN), optional :: calendar   ! Calendar type
    character(len=*),     intent(IN), optional :: desc       ! Clock description
    logical,              intent(IN), optional :: perpetual_run   ! Perpetual mode
    logical,              intent(IN), optional :: MasterSyncClock ! If Master clock
    logical,              intent(IN), optional :: restart         ! If restart
    integer(SHR_KIND_IN), intent(IN), optional :: perpetual_ymd   ! Perpetual ymd

    ! ----- Items to overRide from namelist ------
    character(len=*), intent(IN), optional :: restart_file_TMGoverRide 

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_NMLinfoPutData) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(orb_iyear_AD) )   clockNMLinfo%orb_iyear_AD   = orb_iyear_AD
    if ( present(orb_obliq) )      clockNMLinfo%orb_obliq      = orb_obliq
    if ( present(orb_eccen) )      clockNMLinfo%orb_eccen      = orb_eccen
    if ( present(orb_mvelp) )      clockNMLinfo%orb_mvelp      = orb_mvelp
    if ( present(orb_notSet) )     clockNMLinfo%orb_notSet     = orb_notSet
    if ( present(dtime) )          clockNMLinfo%dtime          = dtime
    if ( present(atm_cpl_dt) )     clockNMLinfo%atm_cpl_dt     = atm_cpl_dt
    if ( present(ocn_cpl_dt) )     clockNMLinfo%ocn_cpl_dt     = ocn_cpl_dt
    if ( present(ice_cpl_dt) )     clockNMLinfo%ice_cpl_dt     = ice_cpl_dt
    if ( present(lnd_cpl_dt) )     clockNMLinfo%lnd_cpl_dt     = lnd_cpl_dt
    if ( present(start_ymd) )      clockNMLinfo%start_ymd      = start_ymd
    if ( present(start_tod) )      clockNMLinfo%start_tod      = start_tod
    if ( present(CurrentYMD) )     clockNMLinfo%CurrentYMD     = CurrentYMD
    if ( present(CurrentTOD) )     clockNMLinfo%CurrentTOD     = CurrentTOD
    if ( present(ref_ymd) )        clockNMLinfo%ref_ymd        = ref_ymd
    if ( present(ref_tod) )        clockNMLinfo%ref_tod        = ref_tod
    if ( present(stop_option) )    clockNMLinfo%stop_option    = trim(stop_option)
    if ( present(stop_n) )         clockNMLinfo%stop_n         = stop_n
    if ( present(stop_ymd) )       clockNMLinfo%stop_ymd       = stop_ymd
    if ( present(stop_tod) )       clockNMLinfo%stop_tod       = stop_tod
    if ( present(desc) )           clockNMLinfo%desc           = trim(desc)
    if ( present(perpetual_run) )  clockNMLinfo%perpetual_run  = perpetual_run
    if ( present(perpetual_ymd) )  clockNMLinfo%perpetual_ymd  = perpetual_ymd
    if ( present(restart_option) ) clockNMLinfo%restart_option = trim(restart_option)
    if ( present(MasterSyncClock)) clockNMLinfo%MasterSyncClock= MasterSyncClock
    if ( present(restart) )        clockNMLinfo%restart        = restart
    if ( present(restart_n) )      clockNMLinfo%restart_n      = restart_n
    if ( present(restart_file_TMGoverRide) ) &
         clockNMLinfo%restart_file_TMGoverRide  = trim(restart_file_TMGoverRide)

END SUBROUTINE eshr_timemgr_NMLinfoPutData

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_NMLinfoGetData -- Get clock clockNMLinfo values
!   
! !DESCRIPTION:
!   
!  Gets values out of the clock clockNMLinfo object. Only a subset of the values in the
!  object can be recovered. In the future as more values are needed they can 
!  easily be added.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_NMLinfoGetData( clockNMLinfo, orb_iyear_AD, orb_obliq,    &
                                        orb_eccen, orb_mvelp, atm_cpl_dt,         &
                                        lnd_cpl_dt, ocn_cpl_dt, ice_cpl_dt,       &
                                        start_ymd, start_tod, stop_option,        &
                                        stop_n, stop_ymd, restart_n,              &
                                        restart_option, restart_file_TMGoverRide, &
                                        dtime   )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_NMLinfoType), intent(IN)  :: clockNMLinfo  ! clock clockNMLinfo
    integer,              intent(OUT), optional :: orb_iyear_AD  ! Orbital year
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_obliq     ! Orbital obliquity
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_eccen     ! Orbital eccentricity
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_mvelp     ! Orbital vernal equinox 
    integer(SHR_KIND_IN), intent(OUT), optional :: dtime         ! Time-step
    integer(SHR_KIND_IN), intent(OUT), optional :: atm_cpl_dt    ! Atmos. interval (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: lnd_cpl_dt    ! Land interval (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: ice_cpl_dt    ! Ice interval (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: ocn_cpl_dt    ! Ocean interval (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: start_ymd     ! Start time (YYYYMMDD)
    integer(SHR_KIND_IN), intent(OUT), optional :: start_tod     ! Start time of day (sec)
    character(len=*),     intent(OUT), optional :: stop_option   ! Type to stop at
    integer(SHR_KIND_IN), intent(OUT), optional :: stop_n        ! No. of units to run
    integer(SHR_KIND_IN), intent(OUT), optional :: stop_ymd      ! Date to stop at
    character(len=*),     intent(OUT), optional :: restart_option! Restart frequency type
    integer(SHR_KIND_IN), intent(OUT), optional :: restart_n     ! Restart frequency no.
    ! ------ Items to overRide from namelist -------
    character(len=*),     intent(OUT), optional :: restart_file_TMGoverRide 

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_NMLinfoGetData) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(orb_iyear_AD) )    orb_iyear_AD    = clockNMLinfo%orb_iyear_AD
    if ( present(orb_obliq) )       orb_obliq       = clockNMLinfo%orb_obliq
    if ( present(orb_eccen) )       orb_eccen       = clockNMLinfo%orb_eccen
    if ( present(orb_mvelp) )       orb_mvelp       = clockNMLinfo%orb_mvelp
    if ( present(dtime) )           dtime           = clockNMLinfo%dtime
    if ( present(atm_cpl_dt) )      atm_cpl_dt      = clockNMLinfo%atm_cpl_dt
    if ( present(ocn_cpl_dt) )      ocn_cpl_dt      = clockNMLinfo%ocn_cpl_dt
    if ( present(ice_cpl_dt) )      ice_cpl_dt      = clockNMLinfo%ice_cpl_dt
    if ( present(lnd_cpl_dt) )      lnd_cpl_dt      = clockNMLinfo%lnd_cpl_dt
    if ( present(start_ymd) )       start_ymd       = clockNMLinfo%start_ymd
    if ( present(start_tod) )       start_tod       = clockNMLinfo%start_tod
    if ( present(stop_option) )     stop_option     = clockNMLinfo%stop_option
    if ( present(stop_n) )          stop_n          = clockNMLinfo%stop_n
    if ( present(stop_ymd) )        stop_ymd        = clockNMLinfo%stop_ymd
    if ( present(restart_option) )  restart_option  = clockNMLinfo%restart_option
    if ( present(restart_n) )       restart_n       = clockNMLinfo%restart_n
    if ( present(restart_file_TMGoverRide) ) &
         restart_file_TMGoverRide  = clockNMLinfo%restart_file_TMGoverRide

END SUBROUTINE eshr_timemgr_NMLinfoGetData

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_NMLinfoRestRead -- Read in the clock restart file
!   
! !DESCRIPTION:
!   
! Read in clock restart information from given input netCDF file into the
! clockNMLinfo object. This sets up the object so that it then can be used
! to initialize a clock object with a time that will correspond to the restart 
! file.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_NMLinfoRestRead( restart_file, MPICom, &
                                         MasterTask, clockNMLinfo )
! !USES:

   USE shr_ncio_mod,   ONLY : shr_ncio_descripRead, shr_ncio_descripName,     &
                              shr_ncio_descripGetInteger,                     &
                              shr_ncio_descripGetRealR8,                      &
                              shr_ncio_descripGetString,                      &
                              shr_ncio_descripGetLogical
   USE shr_string_mod, ONLY : shr_string_listGetIndexF

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(SHR_KIND_CL),      intent(IN) :: restart_file! Restart local filename
   integer, optional,           intent(IN) :: MPICom      ! MPI communicator
   logical, optional,           intent(IN) :: MasterTask  ! If root MPI task or not
   type(eshr_timemgr_NMLinfoType), intent(INOUT) :: clockNMLinfo ! Input Clock NML info

!EOP

   !----- local -----
   character(len=*), parameter :: subname = '(eshr_timemgr_NMLinfoRestRead) '
   logical :: MasterTask2                   ! If root MPI task or not
   integer :: i                             ! Index
   integer :: RestartIntervalSec            ! Restart interval in seconds
   integer :: RestartIntervalMonths         ! Restart interval in months
   integer :: RestartIntervalYears          ! Restart interval in years
   integer :: RestartNextAlarmYMD           ! Date of next restart alarm (YYYYMMDD)
   integer :: RestartNextAlarmTOD           ! Time of day of next restart alarm (seconds)
   integer :: nSet                          ! Number of values set
   integer :: list                          ! List index of item in overRide list
   integer :: ncId                          ! NetCDF file id
   logical :: exists                        ! If file exists or not
   type(eshr_timemgr_clockType) :: clock    ! Local clock to setup restart information

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( present(MasterTask) )then
      MasterTask2 = MasterTask
   else
      MasterTask2 = .true.
   end if

   clockNMLinfo%restart = .true.

   call eshr_timemgr_NMLinfoCheck( clockNMLinfo )
   ! --- Setup restart information so can read from restart file ---------------
   call eshr_timemgr_clockSetupRest( Clock )
   ! --- Read in restart file -------
   call shr_ncio_open( restart_file, MasterTask2, FileType=prefix//"restart_file", &
                       ncId=ncId, exists=exists )
   if ( present(MPICom) )then
      call shr_ncio_descripRead( ncId, NClockVars, prefix//"clock_", MPICom, &
                                 MasterTask2, var=Clock%var )
   else
      call shr_ncio_descripRead( ncId, NClockVars, prefix//"clock_", &
                                 var=Clock%var )
   end if
   call shr_ncio_close( ncId, MasterTask2, type=prefix//"restart_file", &
                        NCFileName=restart_file )
   !----------------------------------------------------------------------------
   ! Put information read into derived type for clock variables
   ! Check if item is on over-ride list (if list_index != 0)
   ! If on over-ride list use value from clockNMLinfo else use the value from the restart
   !----------------------------------------------------------------------------
   nSet = 0
   do i = 1, NClockVars
      if (      trim(shr_ncio_descripName(Clock%var(i))) == "start_ymd"   )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                      trim(shr_ncio_descripName(Clock%var(i)) ) )
         if ( list == 0 )then
            clockNMLinfo%start_ymd = shr_ncio_descripGetInteger( Clock%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "start_tod"   )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                     trim(shr_ncio_descripName(Clock%var(i)) ) )
         if ( list == 0 )then
            clockNMLinfo%start_tod = shr_ncio_descripGetInteger( Clock%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "ref_ymd"     )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                     trim(shr_ncio_descripName(Clock%var(i)) ) )
         if ( list == 0 )then
            clockNMLinfo%ref_ymd = shr_ncio_descripGetInteger( Clock%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "ref_tod"     )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                     trim(shr_ncio_descripName(Clock%var(i)) ) )
         if ( list == 0 )then
            clockNMLinfo%ref_tod = shr_ncio_descripGetInteger( Clock%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "CurrentYMD"    )then
         clockNMLinfo%CurrentYMD  = shr_ncio_descripGetInteger( Clock%var(i) )
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "CurrentTOD"    )then
         clockNMLinfo%CurrentTOD = shr_ncio_descripGetInteger( Clock%var(i) )
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "DTime"       )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                     trim(shr_ncio_descripName(Clock%var(i)) ) )
         if ( list == 0 )then
            clockNMLinfo%Dtime = shr_ncio_descripGetInteger( Clock%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "RestartIntervalSec" &
      )then
         RestartIntervalSec = shr_ncio_descripGetInteger( Clock%var(i) )
         if ( RestartIntervalSec > 0 ) nSet = nSet + 1
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "RestartIntervalMonths" &
      )then
         RestartIntervalMonths = shr_ncio_descripGetInteger( Clock%var(i) )
         if ( RestartIntervalMonths > 0 ) nSet = nSet + 1
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "RestartIntervalYears"&
      )then
         RestartIntervalYears = shr_ncio_descripGetInteger( Clock%var(i) )
         if ( RestartIntervalYears > 0 ) nSet = nSet + 1
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "RestartNextAlarmYMD" &
      )then
         RestartNextAlarmYMD = shr_ncio_descripGetInteger( Clock%var(i) )
      else if ( trim(shr_ncio_descripName(Clock%var(i))) == "RestartNextAlarmTOD" &
      )then
         RestartNextAlarmTOD = shr_ncio_descripGetInteger( Clock%var(i) )
      else
         call shr_sys_abort( subname//'ERROR: unknown clock variable name: '// &
                             trim(shr_ncio_descripName(Clock%var(i) ) ) )
      end if
   end do
   if ( nSet > 1 )then
      call shr_sys_abort( subname//'ERROR: more than one restart interval set' )
   end if
   !----------------------------------------------------------------------------
   ! Read the clockInfo object restart information
   !----------------------------------------------------------------------------
   if ( present(MPICom) )then
      call eshr_timemgr_clockInfoRestRead( restart_file, MPICom,     &
                                           MasterTask2,              &
                                           clockInfo=Clock%info,     &
                                           clockNMLinfo=clockNMLinfo )
   else
      call eshr_timemgr_clockInfoRestRead( restart_file,              &
                                           clockInfo=Clock%info,      &
                                           clockNMLinfo=clockNMLinfo )
   end if

   !----------------------------------------------------------------------------
   ! Check overRide list to see if should use restart frequency set on input namelist 
   !----------------------------------------------------------------------------
   list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                    "restart_option" )
   !----------------------------------------------------------------------------
   ! Use the restart frequency from the restart file
   !----------------------------------------------------------------------------
   if ( list == 0 )then
      ! --- Change clockNMLinfo to correspond to restart file -------------------------
      if (      RestartIntervalSec    > 0 )then
        if ( mod(RestartIntervalSec,eshr_timemgr_secPerDay) == 0 )then
           clockNMLinfo%restart_option = eshr_timemgr_rest_optNDays
           clockNMLinfo%restart_n      = RestartIntervalSec / eshr_timemgr_secPerDay
        else if ( mod(RestartIntervalSec,clockNMLinfo%DTime) == 0 )then
           clockNMLinfo%restart_option = eshr_timemgr_rest_optNSteps
           clockNMLinfo%restart_n      = RestartIntervalSec / clockNMLinfo%DTime
        else
            call shr_sys_abort( subname//'ERROR: RestartIntervalSec on '// &
                                'restart file not evenly divisible by ' // &
                                'seconds per day or time-step' )
        end if
      else if ( RestartIntervalMonths > 0 )then
        clockNMLinfo%restart_option = eshr_timemgr_rest_optNMonths
        clockNMLinfo%restart_n      = RestartIntervalMonths
      else if ( RestartIntervalYears  > 0 )then
        clockNMLinfo%restart_option = eshr_timemgr_rest_optNYears
        clockNMLinfo%restart_n      = RestartIntervalYears
      end if
      clockNMLinfo%NextRestYMD = RestartNextAlarmYMD
      clockNMLinfo%NextRestTOD = RestartNextAlarmTOD
   end if

END SUBROUTINE eshr_timemgr_NMLinfoRestRead       

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoPutData -- Change data in the clock-Info object
!   
! !DESCRIPTION:
!   
!     Put input data into the clock-info object
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInfoPutData( Info, NoRestarts, orb_eccen,      &
                                          orb_mvelp, orb_obliq, orb_obliqr, &
                                          orb_lambm0, orb_mvelpp )

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockInfoType), intent(INOUT) :: Info   ! Input clockInfo object
   logical,             optional,   intent(IN) :: NoRestarts ! Do NOT write restarts
   real(SHR_KIND_R8),   optional,   intent(IN) :: orb_eccen  ! Orbital eccentricity
   real(SHR_KIND_R8),   optional,   intent(IN) :: orb_mvelp  ! Orbital vernal equinox pos.
   real(SHR_KIND_R8),   optional,   intent(IN) :: orb_obliq  ! Orbital obliquity
   real(SHR_KIND_R8),   optional,   intent(IN) :: orb_obliqr ! Orbital obliq. in radians
   real(SHR_KIND_R8),   optional,   intent(IN) :: orb_lambm0 ! Longitude of perihelion rad
   real(SHR_KIND_R8),   optional,   intent(IN) :: orb_mvelpp ! mvelp in radians

!EOP
    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockInfoPutData) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(NoRestarts) ) Info%NoRestarts = NoRestarts
    if ( present(orb_obliq) )  Info%orb_obliq  = orb_obliq
    if ( present(orb_eccen) )  Info%orb_eccen  = orb_eccen
    if ( present(orb_mvelp) )  Info%orb_mvelp  = orb_mvelp
    if ( present(orb_obliqr) ) Info%orb_obliqr = orb_obliqr
    if ( present(orb_lambm0) ) Info%orb_lambm0 = orb_lambm0
    if ( present(orb_mvelpp) ) Info%orb_mvelpp = orb_mvelpp

END SUBROUTINE eshr_timemgr_clockInfoPutData

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoGet -- Get values from Clock-Info object
!   
! !DESCRIPTION:
!   
!  Gets values out of the clock info object.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInfoGet( ClockInfo, orb_iyear_AD, orb_eccen,    &
                                      orb_mvelp, orb_obliq, orb_obliqr,      &
                                      orb_lambm0, orb_mvelpp, perpetual_ymd, &
                                      perpetual_run, calendar, desc,         &
                                      ETimePerpetual )
    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_clockInfoType), intent(IN) :: ClockInfo    ! input object
    ! ---- Perpetual date YYYYMMDD --------
    logical,              intent(OUT), optional :: perpetual_run ! If in perpetual mode
    integer(SHR_KIND_IN), intent(OUT), optional :: perpetual_ymd ! Perpetual date
    integer,              intent(OUT), optional :: orb_iyear_AD  ! Orbital year
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_obliq     ! Obliquity of orbit
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_eccen     ! Eccentricity of orbit
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_mvelp     ! vernal equinox location
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_obliqr    ! Orbital obliq. in rad
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_lambm0    ! Longitude of perh. rad
    real(SHR_KIND_R8),    intent(OUT), optional :: orb_mvelpp    ! mvelp in radians
    type(ESMF_Time),      intent(OUT), optional :: ETimePerpetual! Perpetual date in ESMF
    character(len=SHR_KIND_CS), intent(OUT), optional :: calendar! Calendar type
    character(len=SHR_KIND_CL), intent(OUT), optional :: desc    ! Clock description

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(orb_iyear_AD) )  orb_iyear_AD = ClockInfo%orb_iyear_AD
    if ( present(orb_obliq)    )  orb_obliq    = ClockInfo%orb_obliq
    if ( present(orb_eccen)    )  orb_eccen    = ClockInfo%orb_eccen
    if ( present(orb_mvelp)    )  orb_mvelp    = ClockInfo%orb_mvelp
    if ( present(orb_obliqr)   )  orb_obliqr   = ClockInfo%orb_obliqr
    if ( present(orb_mvelpp)   )  orb_mvelpp   = ClockInfo%orb_mvelpp
    if ( present(orb_lambm0)   )  orb_lambm0   = ClockInfo%orb_lambm0
    ! --- If want perpetual date -----------------------------------------------
    if ( present(perpetual_ymd) )then
        if ( ClockInfo%perpetual_run ) then
           perpetual_ymd = eshr_timemgr_ETimeGetDataYMD( ClockInfo%perpetual_time )
        else
           perpetual_ymd = 0
        end if
    end if
    if ( present(perpetual_run) ) perpetual_run = ClockInfo%perpetual_run
    ! --- Calendar -------------------------------------------------------------
    if ( present(calendar) )      calendar = ClockInfo%calendar
    if ( present(desc) )          desc     = ClockInfo%desc
    if ( present(ETimePerpetual) )then
       ETimePerpetual = ClockInfo%perpetual_time
    end if

END SUBROUTINE eshr_timemgr_clockInfoGet

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoIsSame -- Check if two clock-info objects are same
!   
! !DESCRIPTION:
!   
!  Returns true if the two input clock-info objects have the same values.
!      
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION eshr_timemgr_clockInfoIsSame( ClockInfo, ClockInfo2 )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(eshr_timemgr_clockInfoType), intent(IN) :: ClockInfo    ! input object
    type(eshr_timemgr_clockInfoType), intent(IN) :: ClockInfo2   ! 2nd input object
!EOP
    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
    eshr_timemgr_clockInfoIsSame = .false.
    if ( clockInfo%NoRestarts       .neqv. clockInfo2%NoRestarts      ) return
    if ( clockInfo%perpetual_run    .neqv. clockInfo2%perpetual_run   ) return
    if ( clockInfo%MasterSyncClock  .neqv. clockInfo2%MasterSyncClock ) return
    if ( clockInfo%perpetual_run )then
       if ( clockInfo%perpetual_time   /=     clockInfo2%perpetual_time  ) return
    end if
    if ( clockInfo%calendar         /=     clockInfo2%calendar        ) return
    if ( clockInfo%desc             /=     clockInfo2%desc            ) return
    if ( clockInfo%orb_iyear_AD     /=     clockInfo2%orb_iyear_AD    ) return
    if ( clockInfo%orb_obliq        /=     clockInfo2%orb_obliq       ) return
    if ( clockInfo%orb_eccen        /=     clockInfo2%orb_eccen       ) return
    if ( clockInfo%orb_mvelp        /=     clockInfo2%orb_mvelp       ) return
    if ( clockInfo%orb_obliqr       /=     clockInfo2%orb_obliqr      ) return
    if ( clockInfo%orb_mvelpp       /=     clockInfo2%orb_mvelpp      ) return
    if ( clockInfo%orb_mode         /=     clockInfo2%orb_mode        ) return
    eshr_timemgr_clockInfoIsSame = .true.
    return

END FUNCTION eshr_timemgr_clockInfoIsSame

!===============================================================================
#ifdef SEQ_ESMF

!===============================================================================
! !IROUTINE: eshr_timemgr_info2EState -- put input clock-info values on ESMF State
!   
! !DESCRIPTION:
!   
!  Put the values on the input clock-info object onto the input ESMF State.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_timemgr_info2EState( ClockInfo, eState, ClockInfoEState, print )
    use eshr_rc_mod,      only: eshr_rc_check
    use eshr_estate_mod,  only: eshr_estate_FLog2ELog,         &
                                eshr_estate_printAttributes
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(eshr_timemgr_clockInfoType), intent(IN)    :: ClockInfo       ! input object
    type(ESMF_State),                 intent(INOUT) :: eState          ! Output ESMF State
    type(ESMF_State),                 intent(OUT)   :: ClockInfoEState ! ClockInfo E-State
    logical,                optional, intent(IN)    :: print           ! If print or not
!EOP
    !----- local -----
    integer                    :: rc              ! ESMF error return code
    logical                    :: printIt         ! If should print state when done
    character(len=SHR_KIND_CL) :: desc            ! Description of this clock
    character(len=SHR_KIND_CS) :: calendar        ! Type of calendar
    logical                    :: perpetual_run   ! If running in perpetual mode
    integer(SHR_KIND_IN)       :: perpetual_ymd   ! date of perpetual mode (YYYYMMDD)
    integer(SHR_KIND_IN)       :: orb_iyear_AD    ! year for orbit
    real(SHR_KIND_R8)          :: orb_obliq       ! Obliquity of orbit
    real(SHR_KIND_R8)          :: orb_eccen       ! Eccentricity of orbit
    real(SHR_KIND_R8)          :: orb_mvelp       ! Locatn of vernal equinox
    logical                    :: NoRestarts      ! Do NOT do any restarts
    real(SHR_KIND_R8)          :: orb_obliqr      ! Obliquity in radians
    real(SHR_KIND_R8)          :: orb_lambm0      ! Long of perh. radians
    real(SHR_KIND_R8)          :: orb_mvelpp      ! Locatn of vernal equinox at perh.

!-------------------------------------------------------------------------------
! Notes: Don't pass MasterSyncClock and orb_mode as not needed for
!        subcomponents. These are only needed at the driver level.
!-------------------------------------------------------------------------------

    if ( .not. present(print) )then
       printIt = .false.
    else
       printIt = print
    end if

    !------- Create a seperate state for the ClockInfo object ------------------

    ClockInfoEState = ESMF_StateCreate( ClockInfoEStateName, rc=rc )
    call eshr_rc_check( rc, "Error in creation of ClockInfo ESMF State" )

    !------  Add a long_name attribute to describe this state ------------------

    call ESMF_StateSetAttribute( ClockInfoEState, name="long_name", &
                                 value="eshr_timemgr_clockInfo object ESMF State", &
                                 rc=rc )
    call eshr_rc_check( rc, "Error adding long_name to ClockInfo EState" )

    !------ Get the data from the ClockInfo object -----------------------------

    NoRestarts = clockInfo%NoRestarts
    call eshr_timemgr_clockInfoGet( ClockInfo, orb_iyear_AD=orb_iyear_AD,         &
                                    orb_eccen=orb_eccen, orb_mvelp=orb_mvelp,     &
                                    orb_obliq=orb_obliq, orb_obliqr=orb_obliqr,   &
                                    orb_lambm0=orb_lambm0, orb_mvelpp=orb_mvelpp, &
                                    perpetual_ymd=perpetual_ymd,                  &
                                    perpetual_run=perpetual_run,                  &
                                    calendar=calendar, desc=desc )
    !
    !------ Add ClockInfo object data as attributes to this state
    !

    call ESMF_StateSetAttribute( ClockInfoEState, name="orb_iyear_AD", &
                                 value=orb_iyear_AD, rc=rc )
    call eshr_rc_check( rc, "Error adding orb_iyear_AD to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="orb_eccen", &
                                 value=orb_eccen, rc=rc )
    call eshr_rc_check( rc, "Error adding orb_eccen to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="orb_mvelp", &
                                 value=orb_mvelp, rc=rc )
    call eshr_rc_check( rc, "Error adding orb_mvelp to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="orb_obliq", &
                                 value=orb_obliq, rc=rc )
    call eshr_rc_check( rc, "Error adding orb_obliq to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="orb_obliqr", &
                                 value=orb_obliqr, rc=rc )
    call eshr_rc_check( rc, "Error adding orb_obliqr to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="orb_lambm0", &
                                 value=orb_lambm0, rc=rc )
    call eshr_rc_check( rc, "Error adding orb_lambm0 to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="orb_mvelpp", &
                                 value=orb_mvelpp, rc=rc )
    call eshr_rc_check( rc, "Error adding orb_mvelpp to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="calendar", &
                                 value=calendar, &
                                 rc=rc )
    call eshr_rc_check( rc, "Error adding calendar to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="desc", &
                                 value=desc, &
                                 rc=rc )
    call eshr_rc_check( rc, "Error adding calendar to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="perpetual_run", &
                                 value=eshr_estate_FLog2ELog(perpetual_run), rc=rc )
    call eshr_rc_check( rc, "Error adding perpetual_run to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="NoRestarts", &
                                 value=eshr_estate_FLog2ELog(NoRestarts), rc=rc )
    call eshr_rc_check( rc, "Error adding NoRestarts to clockInfo State" )
    call ESMF_StateSetAttribute( ClockInfoEState, name="perpetual_ymd", &
                                 value=perpetual_ymd, rc=rc )
    call eshr_rc_check( rc, "Error adding perpetual_ymd to clockInfo State" )

    if ( printIt ) call eshr_estate_printAttributes( ClockInfoEState )
    !
    !------ Add the ClockInfo state to the input ESMF state ---------------------
    !

    call ESMF_StateAddState( eState, ClockInfoEState, rc=rc )
    call eshr_rc_check( rc, "Error in adding of ClockInfo  ESMF State to input eState" )

END SUBROUTINE eshr_timemgr_info2EState

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_timemgr_EState2Info -- get input clock-info values from ESMF State
!   
! !DESCRIPTION:
!   
!  Get the values from the input ESMF State and put on the input clock-info object.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_timemgr_eState2Info( eState, ClockInfo, print )
    use eshr_rc_mod,      only: eshr_rc_check
    use eshr_estate_mod,  only: eshr_estate_ELog2FLog,         &
                                eshr_estate_printAttributes
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),                 intent(IN)  :: eState      ! Input ESMF State
    type(eshr_timemgr_clockInfoType), intent(OUT) :: ClockInfo   ! Output object
    logical,                optional, intent(IN)  :: print       ! If print or not
!EOP
    !----- local -----
    type(ESMF_State)               :: ClockInfoEState ! ClockInfo E-State
    integer                        :: rc              ! ESMF error return code
    logical                        :: printIt         ! If should print state when done
    character(len=SHR_KIND_CL)     :: string          ! String to read from ESMF State
    type(ESMF_Logical)             :: ELogValue       ! ESMF logical value
    type(eshr_timemgr_NMLinfoType) :: ClockNMLInfo    ! Namelist object to create Info
    character(len=SHR_KIND_CL)     :: desc            ! Description of this clock
    character(len=SHR_KIND_CS)     :: calendar        ! Type of calendar
    logical                        :: perpetual_run   ! If running in perpetual mode
    integer(SHR_KIND_IN)           :: perpetual_ymd   ! date of perpetual mode (YYYYMMDD)
    integer(SHR_KIND_IN)           :: orb_iyear_AD    ! year for orbit
    real(SHR_KIND_R8)              :: orb_obliq       ! Obliquity of orbit
    real(SHR_KIND_R8)              :: orb_eccen       ! Eccentricity of orbit
    real(SHR_KIND_R8)              :: orb_mvelp       ! Locatn of vernal equinox
    logical                        :: NoRestarts      ! Do NOT do any restarts
    real(SHR_KIND_R8)              :: orb_obliqr      ! Obliquity in radians
    real(SHR_KIND_R8)              :: orb_lambm0      ! Long of perh. radians
    real(SHR_KIND_R8)              :: orb_mvelpp      ! Locatn of vernal equinox at perh.

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( .not. present(print) )then
       printIt = .false.
    else
       printIt = print
    end if

    !------- Get the state for the ClockInfo object from the input EState -------

    call ESMF_StateGetState( eState, ClockInfoEStateName, ClockInfoEState, rc=rc)
    call eshr_rc_check( rc, "Error in getting ClockInfo ESMF State from input state" )

    if ( printIt ) call eshr_estate_printAttributes( ClockInfoEState  )
    !
    !------ Get ClockInfo object data as attributes from this state
    !

    call ESMF_StateGetAttribute( ClockInfoEState, name="calendar", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting calendar from ClockInfo EState" )
    calendar  = string
    call ESMF_StateGetAttribute( ClockInfoEState, name="desc", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting desc from ClockInfo EState" )
    desc  = string
    call ESMF_StateGetAttribute( ClockInfoEState, name="perpetual_run", &
                                 value=ELogValue, rc=rc )
    call eshr_rc_check( rc, "Error getting perpetual_run from ClockInfo EState" )
    perpetual_run = eshr_estate_ELog2FLog(ELogValue)
    call ESMF_StateGetAttribute( ClockInfoEState, name="perpetual_ymd", &
                                 value=perpetual_ymd, rc=rc )
    call eshr_rc_check( rc, "Error getting perpetual_ymd from ClockInfo EState" )
    call ESMF_StateGetAttribute( ClockInfoEState, name="NoRestarts", &
                                 value=ELogValue, rc=rc )
    call eshr_rc_check( rc, "Error getting NoRestarts from ClockInfo EState" )
    NoRestarts = eshr_estate_ELog2FLog(ELogValue)
    call ESMF_StateGetAttribute( ClockInfoEState, name="orb_iyear_AD", &
                                 value=orb_iyear_AD, rc=rc )
    call eshr_rc_check( rc, "Error getting orb_iyear_AD from ClockInfo EState" )
    call ESMF_StateGetAttribute( ClockInfoEState, name="orb_mvelp", &
                                 value=orb_mvelp, rc=rc )
    call eshr_rc_check( rc, "Error getting orb_mvelp from ClockInfo EState" )
    call ESMF_StateGetAttribute( ClockInfoEState, name="orb_obliq", &
                                 value=orb_obliq, rc=rc )
    call eshr_rc_check( rc, "Error getting orb_obliq from ClockInfo EState" )
    call ESMF_StateGetAttribute( ClockInfoEState, name="orb_eccen", &
                                 value=orb_eccen, rc=rc )
    call eshr_rc_check( rc, "Error getting orb_eccen from ClockInfo EState" )
    call ESMF_StateGetAttribute( ClockInfoEState, name="orb_lambm0", &
                                 value=orb_lambm0, rc=rc )
    call eshr_rc_check( rc, "Error getting orb_lambm0 from ClockInfo EState" )
    call ESMF_StateGetAttribute( ClockInfoEState, name="orb_obliqr", &
                                 value=orb_obliqr, rc=rc )
    call eshr_rc_check( rc, "Error getting orb_obliqr from ClockInfo EState" )
    call ESMF_StateGetAttribute( ClockInfoEState, name="orb_mvelpp", &
                                 value=orb_mvelpp, rc=rc )
    call eshr_rc_check( rc, "Error getting orb_mvelpp from ClockInfo EState" )

    ! --- Create a new clockInfo object with the above data on it

    call eshr_timemgr_NMLinfoSetDefault( clockNMLinfo=clockNMLinfo )
    call eshr_timemgr_NMLinfoPutData( clockNMLinfo, orb_iyear_AD=orb_iyear_AD,  &
                                      orb_obliq=orb_obliq, orb_eccen=orb_eccen, &
                                      orb_mvelp=orb_mvelp, calendar=calendar,   &
                                      perpetual_run=perpetual_run,              &
                                      perpetual_ymd=perpetual_ymd,              &
                                      orb_notSet=.false.,                       &
                                      desc=desc )
    call eshr_timemgr_clockInfoInit( clockNMLinfo, ClockInfoOut=ClockInfo )
    clockInfo%NoRestarts = NoRestarts

END SUBROUTINE eshr_timemgr_eState2Info

!===============================================================================
#endif

!===============================================================================
! !IROUTINE: eshr_timemgr_clockGet -- Get information from the clock
!   
! !DESCRIPTION:
!   
!     Get various values from the clock.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockGet( clock, year, month, day, CurrentYMD,         &
                                  CurrentTOD, PrevYMD, PrevTOD, start_ymd,     &
                                  start_tod, StepNo, ref_ymd, ref_tod,         &
                                  stop_ymd, stop_tod, DTime, perpetual_ymd,    &
                                  RestartNextAlarmYMD, RestartNextAlarmTOD,    &
                                  RestartIntervalSec, RestartIntervalMonths,   &
                                  RestartIntervalYears, EClock, ECurrTime,     &
                                  info, calendar, perpetual_run )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_clockType),  intent(IN) :: clock        ! Input clock object
    integer(SHR_KIND_IN), intent(OUT), optional :: year       ! Current year
    integer(SHR_KIND_IN), intent(OUT), optional :: month      ! Current month
    integer(SHR_KIND_IN), intent(OUT), optional :: day        ! Current day in month
    integer(SHR_KIND_IN), intent(OUT), optional :: CurrentYMD ! Current date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: CurrentTOD ! Current time of day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: PrevYMD    ! Previous date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: PrevTOD    ! Previous time of day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: StepNo     ! Number of steps taken
    integer(SHR_KIND_IN), intent(OUT), optional :: start_ymd  ! Starting date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: start_tod  ! Starting time-of-day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: stop_ymd   ! Stop date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: stop_tod   ! Stop time-of-day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: ref_ymd    ! Reference date YYYYMMDD
    integer(SHR_KIND_IN), intent(OUT), optional :: ref_tod    ! Reference time-of-day (s)
    integer(SHR_KIND_IN), intent(OUT), optional :: DTime      ! Time-step (seconds)
    ! ---- Perpetual date YYYYMMDD --------
    logical,              intent(OUT), optional :: perpetual_run ! If in perpetual mode
    integer(SHR_KIND_IN), intent(OUT), optional :: perpetual_ymd        
    ! ---- Next restart alarm date (YYYYMMDD) -----
    integer(SHR_KIND_IN), intent(OUT), optional :: RestartNextAlarmYMD  
    ! ---- Next restart alarm time-of-day (sec) -----
    integer(SHR_KIND_IN), intent(OUT), optional :: RestartNextAlarmTOD  
    ! ---- Restart interval for seconds -----
    integer(SHR_KIND_IN), intent(OUT), optional :: RestartIntervalSec   
    ! ---- Restart interval for months -----
    integer(SHR_KIND_IN), intent(OUT), optional :: RestartIntervalMonths
    ! ---- Restart interval for years -----
    integer(SHR_KIND_IN), intent(OUT), optional :: RestartIntervalYears 
    character(len=SHR_KIND_CS), intent(OUT), optional :: calendar   ! Calendar type
    type(ESMF_Clock),           intent(OUT), optional :: EClock     ! Output ESMF clock
    type(ESMF_Time),            intent(OUT), optional :: ECurrTime  ! Current ESMF time
    type(eshr_timemgr_clockInfoType), intent(OUT), optional :: info ! Clock info

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockGet) '
    type(ESMF_Time) :: CurrentTime            ! Current time
    type(ESMF_Time) :: PreviousTime           ! Previous time
    type(ESMF_Time) :: StartTime              ! Start time
    type(ESMF_Time) :: StopTime               ! Stop time
    type(ESMF_Time) :: RefTime                ! Reference time
    type(ESMF_TimeInterval) :: timeStep       ! Clock, time-step
    integer :: rc                             ! Return code
    integer(SHR_KIND_I8) :: advSteps          ! Number of time-steps that have advanced
    integer :: yy, mm, dd, sec                ! Return time values
    integer :: ymd                            ! Date (YYYYMMDD)
    type(ESMF_TimeInterval) :: alarmInterval  ! Alarm interval
    type(ESMF_Time) :: ringTime               ! Next alarm ring time

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- Get current time from clock ------------------------------------------
    call ESMF_ClockGet( clock%EClock, currTime=CurrentTime, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname//"Error from ESMF_ClockGet" )
    call ESMF_TimeGet( CurrentTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname//"Error from ESMF_TimeGet" )
    if ( present(year) )     year     = yy
    if ( present(month) )    month    = mm
    if ( present(day) )      day      = dd
    if ( present(CurrentTOD) ) CurrentTOD = sec
    if ( present(CurrentYMD) ) CurrentYMD = eshr_timemgr_ETimeGetDataYMD( CurrentTime )
    if ( present(ECurrTime)  ) ECurrTime  = CurrentTime
    ! --- Number of advance steps ----------------------------------------------
    if ( present(StepNo)  )then
       call ESMF_ClockGet(clock%EClock, advanceCount=advSteps, rc=rc)
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                 "Error from ESMF_ClockGet, advanceCount" )
       StepNo = advSteps
    end if
    ! --- Previous time --------------------------------------------------------
    if ( present(PrevYMD) .or. present(PrevTOD) )then
       call ESMF_ClockGet( clock%EClock, prevTime=previousTime, rc=rc )
       ymd = eshr_timemgr_ETimeGetDataYMD( PreviousTime, tod=sec )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname//"Error from ESMF_ClockGet" )
       if ( present(PrevTOD) ) PrevTOD = sec
       if ( present(PrevYMD) ) PrevYMD = ymd
    end if
    !---------------------------------------------------------------------------
    ! If want information about restart alarm
    !---------------------------------------------------------------------------
    if ( present(RestartNextAlarmTOD) .or. present(RestartNextAlarmYMD)  .or. &
         present(RestartIntervalSec)  .or. present(RestartIntervalMonths) .or.&
         present(RestartIntervalYears) )then
       call ESMF_AlarmGet( clock%restart, RingTime=RingTime, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_AlarmGet" )
       call ESMF_TimeGet( ringTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeGet" )
       ! --- If want restart next alarm information ----------------------------
       if ( present(RestartNextAlarmTOD) ) RestartNextAlarmTOD = sec
       if ( present(RestartNextAlarmYMD) ) RestartNextAlarmYMD = &
                                                 eshr_timemgr_ETimeGetDataYMD( RingTime )
       call ESMF_AlarmGet( clock%restart, RingInterval=AlarmInterval, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_AlarmGet" )
       call ESMF_TimeIntervalGet( alarmInterval, yy=yy, mm=mm, d=dd, s=sec, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalGet" )
       sec = sec + dd*86400.0_SHR_KIND_R8
       ! --- If want restart next interval information -------------------------
       if ( present(RestartIntervalSec)    ) RestartIntervalSec    = sec
       if ( present(RestartIntervalMonths) ) RestartIntervalMonths = mm
       if ( present(RestartIntervalYears)  ) RestartIntervalYears  = yy
    end if
    ! --- If want start date ---------------------------------------------------
    if ( present(start_ymd) .or. present(start_tod) )then
       call ESMF_ClockGet( clock%EClock, startTime=StartTime, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_ClockGet" )
       call ESMF_TimeGet( StartTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeGet" )
       if ( present(start_ymd) ) start_ymd = eshr_timemgr_ETimeGetDataYMD( StartTime )
       if ( present(start_tod) ) start_tod = sec
    end if
    ! --- If want stop date ----------------------------------------------------
    if ( present(stop_ymd) .or. present(stop_tod) )then
       call ESMF_ClockGet( clock%EClock, stopTime=StopTime, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_ClockGet" )
       call ESMF_TimeGet( StopTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeGet" )
       if ( present(stop_ymd) ) stop_ymd = eshr_timemgr_ETimeGetDataYMD( StopTime )
       if ( present(stop_tod) ) stop_tod = sec
    end if
    ! --- If want reference date -----------------------------------------------
    if ( present(ref_ymd) .or. present(ref_tod) )then
       call ESMF_ClockGet( clock%EClock, refTime=RefTime, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_ClockGet" )
       call ESMF_TimeGet( RefTime, yy=yy, mm=mm, dd=dd, s=sec, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeGet" )
       if ( present(ref_ymd) ) ref_ymd = eshr_timemgr_ETimeGetDataYMD( RefTime )
       if ( present(ref_tod) ) ref_tod = sec
    end if
    ! --- If want time step ----------------------------------------------------
    if ( present(DTime) )then
       call ESMF_ClockGet( clock%EClock, TimeStep=timeStep, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_ClockGet" )
       call ESMF_TimeIntervalGet( timeStep, s=DTime, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalGet" )
    end if
    ! --- If want ESMF clock ---------------------------------------------------
    if ( present(EClock) )then
       call eshr_timemgr_clockGetESMFClock( clock, EClock )
    end if
    ! --- If want eshr_clock info object ----------------------------------------
    if ( present(info) )then
       call eshr_timemgr_clockGetInfo( clock, info )
    end if
    !---------------------------------------------------------------------------
    ! --- If want data from the eshr_clock info object --------------------------
    !---------------------------------------------------------------------------
    ! --- If want perpetual date -----------------------------------------------
    if ( present(perpetual_ymd) )then
       call eshr_timemgr_clockInfoGet( clock%info, perpetual_ymd=perpetual_ymd )
    end if
    if ( present(perpetual_run) )then
       call eshr_timemgr_clockInfoGet( clock%info, perpetual_run=perpetual_run )
    end if
    ! --- Calendar -------------------------------------------------------------
    if ( present(calendar) )then
       call eshr_timemgr_clockInfoGet( clock%info, calendar=calendar )
    end if

END SUBROUTINE eshr_timemgr_clockGet

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockGetCalDay -- Get current day in year
!   
! !DESCRIPTION:
!   
!  Get the current day into the year
!      
! !INTERFACE: ------------------------------------------------------------------

real(SHR_KIND_R8) FUNCTION eshr_timemgr_clockGetCalDay( Clock, Offset )

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_clockType), intent(IN) :: Clock    ! input object
    integer,            optional, intent(IN) :: Offset   ! Offset time (sec)

!EOP

    !----- local -----
   character(len=*), parameter :: subname = '(eshr_timemgr_clockGetCalDay) '
   integer :: rc                          ! Return code
   type(ESMF_Time) :: date                ! Current date
   type(ESMF_Time) :: PerpDate            ! Perpetual date
   type(ESMF_TimeInterval) :: off         ! Offset
   type(ESMF_TimeInterval) :: diurnal     ! Diurnal cycle
   integer :: year, month, day, tod       ! Year, month day and time-of-day

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
   call eshr_timemgr_clockGet( clock, ECurrTime=date )

   ! --- If offset exists add it to current date -------------------------------
   if (present(offset)) then
      if (offset > 0) then
         call ESMF_TimeIntervalSet( off, s=offset, rc=rc )
         call eshr_timemgr_ErCodeCheck(rc, subname// &
                                       ': error return from ESMF_TimeIntervalSet')
         date = date + off
      else if (offset < 0) then
         call ESMF_TimeIntervalSet( off, s=-offset, rc=rc )
         call eshr_timemgr_ErCodeCheck(rc, subname// &
                                       ': error return from ESMF_TimeIntervalSet')
         date = date - off
      end if
   end if
   !
   ! --- For perpetual date -- add time of day to perpetual date ---------------
   !
   if ( eshr_timemgr_clockInfoIsPerpet( Clock%info ) ) then
      call eshr_timemgr_clockInfoGet( Clock%info, ETimePerpetual=PerpDate )
      ! -- Get current time-of-day ---------------------------------------------
      call ESMF_TimeGet(date, yy=year, mm=month, dd=day, s=tod, rc=rc)
      call eshr_timemgr_ErCodeCheck(rc, subname//&
                                    ': error return from ESMF_TimeGet')
      ! --- Get date from perpetual date add time-of-day to it -----------------
      call ESMF_TimeIntervalSet( diurnal, s=tod, rc=rc )
      call eshr_timemgr_ErCodeCheck(rc, subname// &
                        ': error return from ESMF_TimeIntervalSet')
      date = PerpDate + diurnal
   end if
   !
   ! --- Now get the interval from start of year to current date ---------------
   !
   eshr_timemgR_clockGetCalDay = eshr_timemgr_ETimeGetCalDay( date )

END FUNCTION eshr_timemgr_clockGetCalDay

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockGetRunDays -- Get running days from ref date
!   
! !DESCRIPTION:
!   
!  Get the running days from the reference date.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockGetRunDays( Clock, CurrentDays, CurrentSecs, &
                                         PrevDays, PrevSecs )
    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_clockType), intent(IN) :: Clock       ! input object
    integer, optional,           intent(OUT) :: CurrentDays ! Days from current time
    integer, optional,           intent(OUT) :: CurrentSecs ! Left over days
    integer, optional,           intent(OUT) :: PrevDays    ! Days from previous time
    integer, optional,           intent(OUT) :: PrevSecs    ! Left over days

!EOP

    !----- local -----
   character(len=*), parameter :: subname = '(eshr_timemgr_clockGetRunDays) '
   integer :: rc                          ! Return code
   type(ESMF_Time) :: date                ! Current or previous date
   type(ESMF_Time) :: RefDate             ! Reference date
   type(ESMF_TimeInterval) :: diff        ! Time interval from reference date
   integer :: days                        ! Days from reference date
   integer :: secs                        ! Left over seconds

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( (present(CurrentDays) .or. present(CurrentSecs)) .and. &
        (present(PrevDays)    .or. present(PrevSecs)) )then
      call shr_sys_abort(subname//"CurrentDays/Secs and PrevDays/Secs"// &
                         "can NOT both be entered to this subroutine" )
   end if
   if ( present(CurrentDays) .or. present(CurrentSecs) )then
      call eshr_timemgr_clockGet( clock, ECurrTime=date )
   else if ( present(PrevDays) .or. present(PrevSecs) )then
      call ESMF_ClockGet( clock%EClock, prevTime=date, rc=rc )
      call eshr_timemgr_ErCodeCheck(rc, subname// &
                           ': error return from ESMF_ClockGet')
   end if
   call ESMF_ClockGet( clock%EClock, refTime=RefDate, rc=rc )
   call eshr_timemgr_ErCodeCheck(rc, subname// &
                        ': error return from ESMF_ClockGet')
   diff = date - RefDate
   call ESMF_TimeIntervalGet( diff, d=days, s=secs, rc=rc )
   call eshr_timemgr_ErCodeCheck(rc, subname// &
                        ': error return from ESMF_TimeIntervalGet')
   if ( present(CurrentDays) ) CurrentDays = days
   if ( present(CurrentSecs) ) CurrentSecs = secs
   if ( present(PrevDays)    ) PrevDays    = days
   if ( present(PrevSecs)    ) PrevSecs    = secs

END SUBROUTINE eshr_timemgr_clockGetRunDays

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_timemgr_clockGetPerpYMD  -- Get the perpetual date in terms of YYYYMMDD
!   
! !DESCRIPTION:
!   
! Get the perpetual date as an integer (optionally send an offset time and get 
! seconds in day as well).
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_timemgr_clockGetPerpYMD( clock, ymd, tod, offset )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(IN)  :: clock  ! Input clock to query
   integer,                      intent(OUT) :: ymd    ! date (YYYYMMDD)
   integer, optional,            intent(OUT) :: tod    ! time of day (sec)
   integer, optional,            intent(IN)  :: offset ! Offset in seconds

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockGetPerpYMD) '
    integer :: yy, mm, dd               ! Return time values
    type(ESMF_Time) :: CurrentTime      ! Current time
    integer :: rc                       ! Return code
    type(ESMF_TimeInterval) :: DelTime  ! Delta time to add to time

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- If perpetual mode get perpetual date -------
    if ( eshr_timemgr_clockInfoIsPerpet( clock%info ) )then
      CurrentTime = clock%info%Perpetual_Time
    ! --- Otherwise get current date ----
    else
      call ESMF_ClockGet( clock%EClock, currTime=CurrentTime, rc=rc )
      call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                     "Error from ESMF_ClockGet" )
    end if
    if ( present(offset) )then
       call ESMF_TimeIntervalSet( DelTime, s=offset, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalSet" )
       CurrentTime = CurrentTime + DelTime
    end if
    ymd = eshr_timemgr_ETimeGetDataYMD( CurrentTime )
    if ( present(tod) )then
       call ESMF_TimeGet( CurrentTime, yy=yy, mm=mm, dd=dd, s=tod, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname//"Error from ESMF_TimeGet" )
    end if

END SUBROUTINE eshr_timemgr_clockGetPerpYMD

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_timemgr_clockInitNMLinfo -- Setup the clock based on clockNMLinfo
!   
! !DESCRIPTION:
!   
!     Setup the clock based on the input clockNMLinfo information.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInitNMLinfo( clockNMLinfo, LogPrint, ClockOut )

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_NMLinfoType), intent(IN) :: clockNMLinfo ! clock clockNMLinfo
    logical,                        intent(IN) :: LogPrint     ! Whether to print log info
    type(eshr_timemgr_clockType),   intent(OUT):: ClockOut     ! Output clock

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockInitNMLinfo) '
    type(ESMF_Time) :: StartTime              ! Start time 
    type(ESMF_Time) :: RefTime                ! Reference time
    type(ESMF_Time) :: CurrentTime            ! Current time
    type(ESMF_TimeInterval) :: TimeStep       ! Clock time-step
    type(ESMF_TimeInterval) :: AlarmInterval  ! Alarm interval
    integer :: rc                             ! Return code
    integer :: min_cpl_dt                     ! Minimum coupling interval
    integer :: dtime                          ! time-step to use

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- Check the clock clockNMLinfo values ----------------------------------
    call eshr_timemgr_NMLinfoCheck( clockNMLinfo )
    if( LogPrint ) call eshr_timemgr_NMLinfoPrint( clockNMLinfo )
    call eshr_timemgr_initCalendar( clockNMLinfo%calendar )
    StartTime = eshr_timemgr_ETimeInit( clockNMLinfo%start_ymd, &
                                       clockNMLinfo%start_tod, "Start date" )
    ! --- Initialize reference date for time coordinate. -----------------------

    if ( clockNMLinfo%ref_ymd /= 0 ) then
       RefTime = eshr_timemgr_ETimeInit( clockNMLinfo%ref_ymd, &
                                        clockNMLinfo%ref_tod, &
                                        "Reference date" )
    else
       RefTime = StartTime
    end if

    if ( clockNMLinfo%currentYMD == 0 )then
       CurrentTime = StartTime
    else
       CurrentTime = eshr_timemgr_ETimeInit( clockNMLinfo%CurrentYMD, &
                                            clockNMLinfo%CurrentTOD, &
                                            "Current date" )
    end if

    !
    ! --- Figure out what CCSM time-stepping interval should be. ---------------
    !
    min_cpl_dt = min(clockNMLinfo%atm_cpl_dt, clockNMLinfo%lnd_cpl_dt)
    min_cpl_dt = min(min_cpl_dt, clockNMLinfo%ocn_cpl_dt)
    min_cpl_dt = min(min_cpl_dt, clockNMLinfo%ice_cpl_dt)
    if ( clockNMLinfo%dtime == 0 )then
       if ( mod(min_cpl_dt,clockNMLinfo%atm_cpl_dt) /= 0 )then
          call shr_sys_abort( subname//' : error atm coupling interval not '// &
                              'compatible' )
       end if
       if ( mod(min_cpl_dt,clockNMLinfo%lnd_cpl_dt) /= 0 )then
          call shr_sys_abort( subname//' : error lnd coupling interval not '// &
                              'compatible' )
       end if
       if ( mod(min_cpl_dt,clockNMLinfo%ice_cpl_dt) /= 0 )then
          call shr_sys_abort( subname//' : error ice coupling interval not '// &
                              'compatible' )
       end if
       if ( mod(min_cpl_dt,clockNMLinfo%ocn_cpl_dt) /= 0 )then
          call shr_sys_abort( subname//' : error ocn coupling interval not '// &
                              'compatible' )
       end if
       dtime = min_cpl_dt
    else
       dtime = clockNMLinfo%dtime
       if ( (min_cpl_dt > 0) .and. mod(min_cpl_dt,clockNMLinfo%dtime) /= 0 )then
          write(6,F0I) 'Dtime:', clockNMLinfo%dtime, "(sec)"
          write(6,F0I) 'Minimum coupling time-step:', min_cpl_dt, "(sec)"
          call shr_sys_abort( subname// &
                              ' : clockNMLinfo dtime not compatible with '//&
                              'coupling intervals' )
        end if
    end if
    call ESMF_TimeIntervalSet( TimeStep, s=dtime, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, subname//': error return from '// &
                                   'ESMF_TimeIntervalSet' )

    ! --- Initialize clock and stop date ---------------------------------------

    call eshr_timemgr_EClockInit( TimeStep, StartTime, RefTime, CurrentTime, &
                                  clockNMLinfo%desc, clockNMLinfo,           &
                                  ESMFClockOut=ClockOut%EClock )
    ! --- Initialize clockInfo object ------------------------------------------
    call eshr_timemgr_clockInfoInit( clockNMLinfo=clockNMLinfo, &
                                     ClockInfoOut=ClockOut%info )
    !---------------------------------------------------------------------------
    ! Restart alarm
    !---------------------------------------------------------------------------
    call eshr_timemgr_clockAlarmRstInit( clockNMLinfo, ClockOut )
    if ( .not. eshr_timemgr_clockInfoNoRest( ClockOut%info ) )then
       call ESMF_AlarmGet( ClockOut%restart, RingInterval=AlarmInterval, rc=rc )
    end if
    ! --- Setup restart information --------------------------------------------
    call eshr_timemgr_clockSetupRest( ClockOut )
    ! --- Advance clock if restart --------
    if ( clockNMLinfo%restart )then
       call eshr_timemgr_clockAdvance( ClockOut )
    end if
    ! --- Print clock out ------
    if ( LogPrint ) call eshr_timemgr_clockPrint( ClockOut )

END SUBROUTINE eshr_timemgr_clockInitNMLinfo

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInitClock -- Setup the clock based on another input
!		      		           clock and a input time-step (in seconds)
!   
! !DESCRIPTION:
!   
!     Setup the clock based on the input clock.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInitClock( ClockIn, DTime, LogPrint, desc, &
                                        ClockOut )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(IN) :: ClockIn   ! Input clock to create new clock
   integer,                      intent(IN) :: DTime     ! Time-step to use (seconds)
   logical,                      intent(IN) :: LogPrint  ! Whether to print log info out
   character(len=*),             intent(IN) :: desc      ! Description of this clock
   type(eshr_timemgr_clockType), intent(OUT) :: ClockOut ! Output clock

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockInitClock) '
    type(ESMF_TimeInterval) :: TimeStep      ! Time-step
    type(ESMF_Time)         :: RefTime       ! Reference time
    type(ESMF_Time)         :: StopTime      ! Stop time
    type(ESMF_Time)         :: CurrentTime   ! Current time
    type(ESMF_Time)         :: StartTime     ! Start time
    type(ESMF_Time)         :: FirstAlarm    ! First restart alarm time
    type(ESMF_TimeInterval) :: AlarmInterval ! Restart alarm interval
    integer                 :: rc            ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- Get information from input clock -------------------------------------
    call ESMF_ClockGet( ClockIn%EClock, starttime = StartTime, &
                        stoptime=StopTime, reftime=RefTime, &
                        currtime=CurrentTime, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname//"Error from ESMF_ClockGet" )
    call ESMF_TimeIntervalSet( TimeStep, s=DTime, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   "Error from ESMF_TimeIntervalSet" )
    
    ! --- Initialize clock and stop date ---------------------------------------

    call eshr_timemgr_EClockInit( TimeStep=TimeStep, StartTime=StartTime,   &
                                  RefTime=RefTime, CurrentTime=CurrentTime, &
                                  desc=desc, StopTime=StopTime,             &
                                  ESMFClockOut=ClockOut%EClock )
    ! --- Restart alarm --------------------------------------------------------

    call ESMF_AlarmGet( ClockIn%restart, RingTime=FirstAlarm, &
                        RingInterval=AlarmInterval, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname//"Error from ESMF_AlarmGet" )
    ClockOut%restart = ESMF_AlarmCreate( name=eshr_timemgr_restartAlarmName, &
                                         clock=ClockOut%EClock,             &
                                         ringTime=FirstAlarm,               &
					 ringInterval=AlarmInterval, rc=rc)
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   "Error from ESMF_AlarmCreate" )
    ! --- Initialize the clockInfo object -----
    call eshr_timemgr_clockInfoInit( desc=desc, MasterSyncClock=.false., &
                                     ClockInfoIn=ClockIn%info,           &
                                     ClockInfoOut=ClockOut%info )
    ! --- Setup restart info -----
    call eshr_timemgr_clockSetupRest( ClockOut )
    ! --- Print out resultant clock -----
    if ( LogPrint ) call eshr_timemgr_clockPrint( ClockOut )

END SUBROUTINE eshr_timemgr_clockInitClock

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockAdvance  -- Advance the clock
!   
! !DESCRIPTION:
!   
! Advance this clock
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockAdvance( clock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(INOUT) :: clock    ! Input clock to advance

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockAdvance) '
    integer :: rc    ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   call ESMF_ClockAdvance( clock%EClock, rc=rc )
   call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                  "Error from ESMF_ClockAdvance" )

END SUBROUTINE eshr_timemgr_clockAdvance

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockStopAtDayEnd -- Set stop time to end of day.
!   
! !DESCRIPTION:
!   
! Set the stop-time of this clock to stop at the end of the current day.
! This is designed for clocks that listen to the coupler for the signal that
! need to stop at the end of this day.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockStopAtDayEnd( LogPrint, clock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   logical,                      intent(IN)    :: LogPrint ! If should print log
   type(eshr_timemgr_clockType), intent(INOUT) :: clock    ! Input clock

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockStopAtDayEnd) '
    integer                 :: rc         ! Return code
    integer                 :: ymd        ! Date of time to stop
    type(ESMF_Time)         :: StopTime   ! Current stop time
    type(ESMF_Time)         :: time2Stop  ! Time that should stop now
    type(ESMF_TimeInterval) :: oneDay     ! One day
    integer                 :: CurrentYMD ! Current date (YYYYMMDD)
    character(len=*), parameter :: F0I = "( '"//LogPrefix//"',A, I8.8)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   !
   ! --- Figure out what the end of the current day is -----
   !
   call eshr_timemgr_clockGet( clock, CurrentYMD=CurrentYMD )
   time2Stop = eshr_timemgr_ETimeInit( CurrentYMD, tod=0, desc="Time to Stop" )
   call ESMF_TimeIntervalSet( oneDay, d=1, rc=rc )
   call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                  "Error from ESMF_TimeIntervalSet" )
   time2Stop = time2Stop + oneDay
   ! --- Get the current stop time so we can see if it differs from above
   call ESMF_ClockGet( clock%EClock, stopTime=StopTime, rc=rc )
   call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                  "Error from ESMF_ClockGet" )
   !
   ! --- Now check if current stop time differs -- if so reset it
   !
   if ( StopTime /= time2Stop )then
      if ( LogPrint )then
         ymd = eshr_timemgr_ETimeGetDataYMD( time2Stop )
         write(6,F0I) 'Stopping at end of the current day at 0Z on date:', ymd
      end if
      call ESMF_ClockSet( clock%EClock, stopTime=time2Stop, rc=rc )
      call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                     "Error from ESMF_ClockSet" )
   end if

END SUBROUTINE eshr_timemgr_clockStopAtDayEnd

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockAlarmOnRest  -- Turn the restart alarm on
!   
! !DESCRIPTION:
!   
! Turn restart alarm for given clock to on.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockAlarmOnRest( clock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(INOUT) :: clock    ! Input clock to set restart

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockAlarmOnRest) '
    integer :: rc    ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call ESMF_AlarmRingerOn( clock%restart, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   "Error from ESMF_AlarmRingerOn" )

END SUBROUTINE eshr_timemgr_clockAlarmOnRest

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockAlarmIsOnRes  -- Test if restart alarm is on or not
!   
! !DESCRIPTION:
!   
! Test if it is time to write restart file out
!      
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION eshr_timemgr_clockAlarmIsOnRes( clock, SyncClock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(IN) :: clock               ! Input clock to query
   type(eshr_timemgr_clockType), intent(IN), optional :: SyncClock ! Synchronization clock

! If SyncClock entered use it to get the next restart alarm time. Then use input clock
! to test if the current time equals the next restart alarm time -- if it is, return true.

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockAlarmIsOnRes) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- If no restarts, restart alarm is never on -----
    if ( eshr_timemgr_clockInfoNoRest( clock%info ) )then
       eshr_timemgr_clockAlarmIsOnRes = .false.
    ! --- Otherwise check more carefully -----
    else
       eshr_timemgr_clockAlarmIsOnRes = .false.
       ! --- If last step, restart alarm is on ------
       if ( eshr_timemgr_clockIsOnLastStep( clock ) ) &
                   eshr_timemgr_clockAlarmIsOnRes = .true.
       ! --- If restart alarm is on, then return true ------
       if ( ESMF_AlarmIsRinging( clock%restart ) ) &
                   eshr_timemgr_clockAlarmIsOnRes = .true.
       if ( present(SyncClock) )then
          call shr_sys_abort( subname// &
                          'ERROR: can NOT handle syncClock option right now' )
       end if
    end if

END FUNCTION eshr_timemgr_clockAlarmIsOnRes

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockIsOnLastStep -- Test if last step or not
!   
! !DESCRIPTION:
!   
! Test if it's the last step or not
!      
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION eshr_timemgr_clockIsOnLastStep( clock, NextStep )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(IN) :: clock    ! Input clock to query
   logical, optional,            intent(IN) :: NextStep ! If should test on next-time step

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockIsOnLastStep) '
    integer :: rc                           ! Return code
    type(ESMF_Time) :: CurrentTime          ! Current time
    type(ESMF_Time) :: StopTime             ! Stop time
    type(ESMF_TimeInterval) :: TimeStep     ! Time-step

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call ESMF_ClockGet( clock%EClock, TimeStep=TimeStep,  &
                        currtime=CurrentTime,             &
                        stoptime=StopTime, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   "Error from ESMF_ClockGet" )
    CurrentTime = CurrentTime
    if ( CurrentTime >= StopTime ) then
       eshr_timemgr_clockIsOnLastStep = .true.
    else
       eshr_timemgr_clockIsOnLastStep = .false.
    end if
    if ( present(NextStep) )then
       if ( NextStep )then
          if ( (CurrentTime+TimeStep) >= StopTime ) &
                       eshr_timemgr_clockIsOnLastStep = .true.
       end if
    end if

END FUNCTION eshr_timemgr_clockIsOnLastStep

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_curTimeLEstopTime -- true if current time <= stop time
!   
! !DESCRIPTION:
!   
! Return true if the current time is less than or equal to the stop time
!      
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION eshr_timemgr_curTimeLEstopTime(clock)

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   type(eshr_timemgr_clockType), intent(IN) :: clock    ! Input clock to query
!EOP
    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_curTimeLEstopTime) '
    integer :: rc                           ! Return code
    type(ESMF_Time) :: CurrentTime          ! Current time
    type(ESMF_Time) :: StopTime             ! Stop time

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call ESMF_ClockGet( clock%EClock, currtime=CurrentTime, &
                        stoptime=StopTime, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   "Error from ESMF_ClockGet" )

    eshr_timemgr_curTimeLEstopTime = .false.
    if ( CurrentTime <= StopTime ) eshr_timemgr_curTimeLEstopTime = .true.

END FUNCTION eshr_timemgr_curTimeLEstopTime

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_timemgr_clockRestWrite -- Write out the clock restart file
!   
! !DESCRIPTION:
!   
! Write time-manager clock information out to input netCDF file
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockRestWrite( restart_file, MPICom, MasterTask, clock )

! !USES:

   use shr_ncio_mod,  only: shr_ncio_descripName, shr_ncio_descripPutData, &
                            shr_ncio_descripWrite

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(SHR_KIND_CL), intent(IN) :: restart_file   ! Restart local filename
   integer, optional,      intent(IN) :: MPICom         ! MPI communicator
   logical, optional,      intent(IN) :: MasterTask     ! If root MPI task or not
   type(eshr_timemgr_clockType), intent(INOUT) :: clock ! Input clock

!EOP

  !----- local -----
  character(len=*), parameter :: subname = '(eshr_timemgr_clockRestWrite) '
  logical :: MasterTask2                ! If root MPI task or not
  integer :: start_ymd                  ! Start date (YYYYMMDD)
  integer :: start_tod                  ! Start time of day (seconds)
  integer :: CurrentYMD                 ! Current date (YYYYMMDD)
  integer :: CurrentTOD                 ! Current time of day (seconds)
  integer :: ref_ymd                    ! Reference date (YYYYMMDD)
  integer :: ref_tod                    ! Reference time of day (seconds)
  integer :: DTime                      ! Time-step
  integer :: RestartNextAlarmYMD        ! Next restart alarm date (YYYYMMDD)
  integer :: RestartNextAlarmTOD        ! Next restart alarm time-of-day (seconds)
  integer :: RestartIntervalSec         ! Restart interval for seconds
  integer :: RestartIntervalMonths      ! Restart interval for months
  integer :: RestartIntervalYears       ! Restart interval for years
  integer :: i                          ! Index
  integer :: ncId                       ! NetCDF file id
  logical :: exists                     ! If file exists or not

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if ( present(MasterTask) )then
     MasterTask2 = MasterTask
  else
     MasterTask2 = .true.
  end if

  call eshr_timemgr_clockGet( clock, RestartIntervalSec = RestartIntervalSec, &
                              RestartIntervalMonths = RestartIntervalMonths,  &
                              RestartIntervalYears = RestartIntervalYears,    &
                              CurrentYMD=CurrentYMD, CurrentTOD=CurrentTOD )
  !-----------------------------------------------------------------------------
  ! Loop over clock variables, set the data to the NetCDF descriptor
  !-----------------------------------------------------------------------------
  do i = 1, NClockVars
      if (      trim(shr_ncio_descripName( clock%var(i) ) ) == "start_ymd" )then
         call eshr_timemgr_clockGet( clock, start_ymd = start_ymd )
         call shr_ncio_descripPutData( clock%var(i), "start_ymd", &
                              IntegerData=start_ymd )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "start_tod" )then
         call eshr_timemgr_clockGet( clock, start_tod = start_tod )
         call shr_ncio_descripPutData( clock%var(i), "start_tod", &
                              IntegerData=start_tod )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "ref_ymd" )then
         call eshr_timemgr_clockGet( clock, ref_ymd = ref_ymd )
         call shr_ncio_descripPutData( clock%var(i), "ref_ymd", IntegerData=ref_ymd )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "ref_tod" )then
         call eshr_timemgr_clockGet( clock, ref_tod = ref_tod )
         call shr_ncio_descripPutData( clock%var(i), "ref_tod", IntegerData=ref_tod )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "CurrentYMD" )then
         call eshr_timemgr_clockGet( clock, CurrentYMD = CurrentYMD )
         call shr_ncio_descripPutData( clock%var(i), "CurrentYMD", &
                                   IntegerData=CurrentYMD )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "CurrentTOD" )then
         call eshr_timemgr_clockGet( clock, CurrentTOD = CurrentTOD )
         call shr_ncio_descripPutData( clock%var(i), "CurrentTOD", &
                                   IntegerData=CurrentTOD )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "DTime" )then
         call eshr_timemgr_clockGet( clock, DTime = DTime )
         call shr_ncio_descripPutData( clock%var(i), "DTime", IntegerData=DTime )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "RestartNextAlarmYMD" &
      )then
         call eshr_timemgr_clockGet( clock, RestartNextAlarmYMD = &
                                     RestartNextAlarmYMD )
         call shr_ncio_descripPutData( clock%var(i), "RestartNextAlarmYMD", &
                                       IntegerData=RestartNextAlarmYMD )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "RestartNextAlarmTOD" &
      )then
         call eshr_timemgr_clockGet( clock, RestartNextAlarmTOD = &
                                     RestartNextAlarmTOD )
         call shr_ncio_descripPutData( clock%var(i), "RestartNextAlarmTOD", &
                                       IntegerData=RestartNextAlarmTOD )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "RestartIntervalSec" &
      )then
         call shr_ncio_descripPutData( clock%var(i), "RestartIntervalSec",       &
                                       IntegerData=RestartIntervalSec )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "RestartIntervalMonths" &
      )then
         call shr_ncio_descripPutData( clock%var(i), "RestartIntervalMonths",    &
                                       IntegerData=RestartIntervalMonths )
      else if ( trim(shr_ncio_descripName( clock%var(i) ) ) == "RestartIntervalYears" &
      )then
         call shr_ncio_descripPutData( clock%var(i), "RestartIntervalYears",     &
                                       IntegerData=RestartIntervalYears )
      else
         call shr_sys_abort( subname//': unknown clock variable name: '// &
                             trim(shr_ncio_descripName( clock%var(i) ) ) )
      end if
  end do
  ! --- Write the clock variables to the restart file ------
  call shr_ncio_open( NCFileName=restart_file, MasterTask=masterTask2,  &
                      FileType=trim(prefix)//"restart_file",            &
                      ncId=ncId, exists=exists, writing=.true. )
  if ( present(MPICom) )then
     call shr_ncio_descripWrite( ncId, NClockVars, prefix//"clock_", mpicom, &
                                 MasterTask2, exists, clock%var )
  else
     call shr_ncio_descripWrite( ncId, NClockVars, prefix//"clock_", &
                                 exists=exists, var=clock%var )
  end if
  call shr_ncio_close( ncId, MasterTask2, type=prefix//"restart_file",  &
                       NCFileName=restart_file )
  ! --- Write the info variables out to the restart file ------
  if ( present(MPICom) )then
     call eshr_timemgr_clockInfoRestWrit( restart_file, MPICom, MasterTask2, &
                                          clockInfo=clock%info )
  else
     call eshr_timemgr_clockInfoRestWrit( restart_file, clockInfo=clock%info )
  end if

END SUBROUTINE eshr_timemgr_clockRestWrite      

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockAlarmOffRest -- Turn the restart alarm off
!   
! !DESCRIPTION:
!   
! Turn the restart alarm off on the given input clock
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockAlarmOffRest( clock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(INOUT) :: clock ! Input clock to turn restart off

!EOP

   !----- local -----
   character(len=*), parameter :: subname = '(eshr_timemgr_clockAlarmOffRest) '
   integer :: rc    ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call ESMF_AlarmRingerOff( clock%restart, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   "Error from ESMF_AlarmRingerOff" )

END SUBROUTINE eshr_timemgr_clockAlarmOffRest   

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockPrint -- Print clock information out
!   
! !DESCRIPTION:
!   
!      Print clock information out.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockPrint( clock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(eshr_timemgr_clockType), intent(in) :: clock   ! Input clock to print

!EOP

    integer(SHR_KIND_IN) :: CurrentYMD ! Current date YYYYMMDD
    integer(SHR_KIND_IN) :: CurrentTOD ! Current time of day (s)
    integer(SHR_KIND_IN) :: StepNo     ! Number of steps taken
    integer(SHR_KIND_IN) :: start_ymd  ! Starting date YYYYMMDD
    integer(SHR_KIND_IN) :: start_tod  ! Starting time-of-day (s)
    integer(SHR_KIND_IN) :: stop_ymd   ! Stop date YYYYMMDD
    integer(SHR_KIND_IN) :: stop_tod   ! Stop time-of-day (s)
    integer(SHR_KIND_IN) :: ref_ymd    ! Reference date YYYYMMDD
    integer(SHR_KIND_IN) :: ref_tod    ! Reference time-of-day (s)
    integer(SHR_KIND_IN) :: DTime      ! Time-step (seconds)
    integer(SHR_KIND_IN) :: RestartNextAlarmYMD   ! Next restart alarm date (YYYYMMDD)
    integer(SHR_KIND_IN) :: RestartNextAlarmTOD   ! Next restart alarm time-of-day (sec)
    integer(SHR_KIND_IN) :: RestartIntervalSec    ! Restart interval for seconds
    integer(SHR_KIND_IN) :: RestartIntervalMonths ! Restart interval for months
    integer(SHR_KIND_IN) :: RestartIntervalYears  ! Restart interval for years

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  call eshr_timemgr_clockGet( clock, CurrentYMD=CurrentYMD,                  &
                              CurrentTOD=CurrentTOD, start_ymd=start_ymd,    &
                              start_tod=start_tod, StepNo=StepNo,            &
                              ref_ymd=ref_ymd, ref_tod=ref_tod,              &
                              stop_ymd=stop_ymd, stop_tod=stop_tod,          &
                              DTime=DTime )
  write(6,F0Date) "Start",       start_ymd, start_tod
  write(6,F0Date) "Current",     CurrentYMD, CurrentTOD
  write(6,F0Date) "Reference",   ref_ymd, ref_tod
  write(6,F0Date) "Stop",        stop_ymd, stop_tod
  write(6,F0I8)   "Step number", StepNo
  write(6,F0I)    "TimeStep",    DTime, " (sec)"
                                
  if ( .not. eshr_timemgr_clockInfoNoRest( clock%info ) )then
     write(6,FA) "Restart alarm    "
     call eshr_timemgr_clockGet( clock, &
                                 RestartNextAlarmYMD=RestartNextAlarmYMD,     &
                                 RestartNextAlarmTOD=RestartNextAlarmTOD,     &
                                 RestartIntervalSec=RestartIntervalSec,       &
                                 RestartIntervalMonths=RestartIntervalMonths, &
                                 RestartIntervalYears=RestartIntervalYears )
      write(6,F0Date) "Next Restart", RestartNextAlarmYMD, RestartNextAlarmTOD
      if (      RestartIntervalSec    > 0 )then
         write(6,F0I) "Restart interval", RestartIntervalSec, " (sec)"
      else if ( RestartIntervalMonths > 0 )then
         RestartIntervalMonths = RestartIntervalMonths + RestartIntervalYears*12
         write(6,F0I) "Restart interval", RestartIntervalMonths, " (months)"
      else if ( RestartIntervalYears  > 0 )then
         write(6,F0I) "Restart interval", RestartIntervalYears, " (years)"
      end if
  end if
  call eshr_timemgr_clockInfoPrint( clock%info )

END SUBROUTINE eshr_timemgr_clockPrint

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockClocksInSync -- Check that two clocks are in sync
!   
! !DESCRIPTION:
!   
!     Check that the two input clocks are in sync
!      
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION eshr_timemgr_clockClocksInSync( clock, SyncClock )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(IN) :: clock
   type(eshr_timemgr_clockType), intent(IN) :: SyncClock

!EOP

   !----- local -----
   character(len=*), parameter :: subname = "(eshr_timemgr_clockClocksInSync) "
   integer         :: rc         ! return code from ESMF
   type(ESMF_Time) :: curr_date1 ! Current time
   type(ESMF_Time) :: curr_date2 ! Current time
   integer         :: tod1       ! Time of day
   integer         :: tod2       ! Time of day

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  call eshr_timemgr_clockGet( clock, ECurrTime=curr_date1, CurrentTOD=tod1 )
  call eshr_timemgr_clockGet( SyncClock, ECurrTime=curr_date2, CurrentTOD=tod2 )
  !
  ! --- If neither clock is in perpetual mode ----------------------------------
  !
  if ( (.not. eshr_timemgr_clockInfoIsPerpet( clock%info ) ) .and. &
       (.not. eshr_timemgr_clockInfoIsPerpet( SyncClock%info ) ) )then
     ! --- If current dates agree return true -- else false
     if ( curr_date1 == curr_date2 )then
        eshr_timemgr_clockClocksInSync = .true.
     else
        eshr_timemgr_clockClocksInSync = .false.
     end if
  !
  ! --- If both clocks are in perpetual mode -----------------------------------
  !
  else if ( eshr_timemgr_clockInfoIsPerpet( clock%info ) .and. &
            eshr_timemgr_clockInfoIsPerpet( SyncClock%info ) )then
      ! --- If perpetual dates and time-of-day agree -- return true -- else false
      call eshr_timemgr_clockInfoGet( clock%info, ETimePerpetual=curr_date1 )
      call eshr_timemgr_clockInfoGet( clock%info, ETimePerpetual=curr_date2 )
      if ( (curr_date1 == curr_date2) .and.  (tod1 == tod2) )then
        eshr_timemgr_clockClocksInSync = .true.
      else
        eshr_timemgr_clockClocksInSync = .false.
      end if
  !
  ! ---- Otherwise clocks do NOT agree on perpetual mode -----------------------
  !
  else
    eshr_timemgr_clockClocksInSync = .false.
  end if

END FUNCTION eshr_timemgr_clockClocksInSync

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockDateInSync -- Check that input date in sync with clock
!   
! !DESCRIPTION:
!   
!     Check that the given input date/time is in sync with clock time
!      
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION eshr_timemgr_clockDateInSync( clock, ymd, tod, prev )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(IN) :: clock   ! Input clock to compare
   integer,                      intent(IN) :: ymd     ! Date (YYYYMMDD)
   integer,                      intent(IN) :: tod     ! Time of day (sec)
   logical, optional,            intent(IN) :: prev    ! If should get previous time

!EOP

   !----- local -----
   character(len=*), parameter :: subname = "(eshr_timemgr_clockDateInSync) "
   integer :: ymd1     ! Date (YYYYMMDD)
   integer :: tod1     ! Time of day
   logical :: previous ! If need to get previous time for comparison

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if ( present(prev) )then
    previous = prev
  else
    previous = .false.
  end if
  if ( .not. previous )then
     call eshr_timemgr_clockGet( clock, CurrentYMD=ymd1, CurrentTOD=tod1 )
  else
     call eshr_timemgr_clockGet( clock, PrevYMD=ymd1, PrevTOD=tod1 )
  end if
  !
  ! --- If current dates agree return true -- else false
  !
  if ( (ymd == ymd1) .and. (tod == tod1) )then
     eshr_timemgr_clockDateInSync = .true.
  else
     eshr_timemgr_clockDateInSync = .false.
  end if

END FUNCTION eshr_timemgr_clockDateInSync

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockPutData  -- Put data objects into the clock object.
!   
! !DESCRIPTION:
!   
! Put data objects into the clock object.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockPutData( clock, EClock, EAlarmRest, ClockInfo, &
                                      orb_set )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType),     intent(INOUT) :: clock  ! Output share clock
   ! clock information object
   type(eshr_timemgr_clockInfoType), intent(IN), optional :: ClockInfo
   ! If should set the orbital information
   logical,                          intent(IN), optional :: orb_set
   ! ESMF clock
   type(ESMF_Clock),                 intent(IN), optional :: EClock
   ! ESMF restart alarm
   type(ESMF_Alarm),                 intent(IN), optional :: EAlarmRest

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockPutData) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( present(ClockInfo) )then
      clock%info    = ClockInfo
      if ( present(orb_set) )then
         if ( orb_set )then
            call eshr_timemgr_clockInfoSetOrb( clock%info )
         else
            call shr_sys_abort( subname//'ERROR: if orb_set present -- it'// &
                                'should only be set to true' )
         end if
      end if
   else
      if ( present(orb_set) ) &
         call shr_sys_abort( subname//'ERROR: orb_set can ONLY be given '// &
                             'if clockinfo input' )
   end if
   if ( present(EClock) )     clock%EClock  = EClock
   if ( present(EAlarmRest) ) clock%restart = EAlarmRest

END SUBROUTINE eshr_timemgr_clockPutData

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_ETimeInit -- Return ESMF_Time object set to input YMD values
!   
! !DESCRIPTION:
!   
!     Return the ESMF_Time object corresponding to the given input time, given in
!  YMD (Year Month Day) and TOD (Time-of-day) format.
! Set the time by an integer as YYYYMMDD and integer seconds in the day
!      
! !INTERFACE: ------------------------------------------------------------------

FUNCTION eshr_timemgr_ETimeInit( ymd, tod, desc )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer, intent(in) :: ymd                  ! Year, month, day YYYYMMDD
   integer, intent(in) :: tod                  ! Time of day in seconds
   character(len=*), intent(in) :: desc        ! Description of time to set
   type(ESMF_Time) :: eshr_timemgr_ETimeInit   ! Return value

!EOP

   !----- local -----
   character(len=*), parameter :: subname = '(eshr_timemgr_ETimeInit) '
   integer :: yr, mon, day          ! Year, month, day as integers
   integer :: rc                    ! return code
   character(len=*), parameter :: F00 = "('"//LogPrefix// &
                                  "',A,T60, I8.8, ' (YYYYMMDD)', I8.8,' (sec)')"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( (ymd < 0) .or. (tod < 0) .or. (tod > eshr_timemgr_secPerDay) )then
      write(6,F00) subname//': error yymmdd is a negative number or '// &
                          'time-of-day out of bounds', ymd, tod
      call shr_sys_abort( subname//'ERROR: Bad input' )
   end if
   yr  = ymd / 10000
   mon = (ymd - yr*10000) / 100
   day =  ymd - yr*10000 - mon*100
   call ESMF_TimeSet( eshr_timemgr_ETimeInit, yy=yr, mm=mon, dd=day, s=tod, &
                      calendar=eshr_timemgr_cal, rc=rc )
   call eshr_timemgr_ErCodeCheck(rc, subname//': error return from '// &
                                 'ESMF_TimeSet: setting '//trim(desc))

END FUNCTION eshr_timemgr_ETimeInit

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_ETimeGetDataYMD -- Get the date in YYYYMMDD form from ESMF Time
!   
! !DESCRIPTION:
!   
!     Get the date in YYYYMMDD format from a ESMF time object.
!      
! !INTERFACE: ------------------------------------------------------------------

INTEGER FUNCTION eshr_timemgr_ETimeGetDataYMD( ETime, offset, tod )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Time),   intent(IN)  :: ETime   ! Input ESMF time
    integer, optional, intent(IN)  :: offset  ! Offset from input time (sec)
    integer, optional, intent(OUT) :: tod     ! Time of day

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_ETimeGetDataYMD) '
    type(ESMF_Time)         :: ETimeAdd ! ESMF time + offset
    type(ESMF_TimeInterval) :: ETimeOff ! ESMF offset time-interval
    integer                 :: year     ! Year
    integer                 :: month    ! Month
    integer                 :: day      ! Day in month
    integer                 :: sec      ! Day in month
    integer                 :: rc       ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ETimeAdd = ETime
    if ( present(offset) )then
        if ( offset > 0 )then
           call ESMF_TimeIntervalSet( ETimeOff, s=offset, rc=rc )
           call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                          ": Error from ESMF_TimeIntervalSet" )
           ETimeAdd = ETime + ETimeOff
        else if ( offset < 0 )then
           call ESMF_TimeIntervalSet( ETimeOff, s=-offset, rc=rc )
           call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                          ": Error from ESMF_TimeIntervalSet" )
           ETimeAdd = ETime - ETimeOff
        end if
    end if
    call ESMF_TimeGet( ETimeAdd, yy=year, mm=month, dd=day, s=sec, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   ": Error from ESMF_TimeGet" )
    eshr_timemgr_ETimeGetDataYMD = year*10000 + month*100 + day
    if ( present(tod) ) tod = sec

END FUNCTION eshr_timemgr_ETimeGetDataYMD

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_ETimeGetCalDay -- Get the calendar day from ESMF Time
!   
! !DESCRIPTION:
!   
!     Get the calendar day from an ESMF time-instant
!      
! !INTERFACE: ------------------------------------------------------------------

real(SHR_KIND_R8) FUNCTION eshr_timemgr_ETimeGetCalDay( ETime )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Time), intent(IN)  :: ETime     ! Input ESMF time

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_ETimeGetCalDay) '
    type(ESMF_TimeInterval) :: ETimeInt ! ESMF time interval from beg of year
    integer                 :: rc       ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call ESMF_TimeGet( ETime, dayOfYear_intvl=ETimeInt, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   ": Error from ESMF_TimeGet" )
    call ESMF_TimeIntervalGet( ETimeInt, d_r8=eshr_timemgr_ETimeGetCalDay, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   ": Error from ESMF_TimeIntervalGet" )
    eshr_timemgr_ETimeGetCalDay = eshr_timemgr_ETimeGetCalDay + 1.0_SHR_KIND_R8
    if ( (eshr_timemgr_ETimeGetCalDay < 1.0) .or.  &
         (eshr_timemgr_ETimeGetCalDay > 366.0) )then
        call shr_sys_abort( subname//": Error calday out of range" )
    end if

END FUNCTION eshr_timemgr_ETimeGetCalDay

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_ETimeGetDiff -- Get difference between two times
!   
! !DESCRIPTION:
!   
!     Get the difference between two ESMF time instants.
!      
! !INTERFACE: ------------------------------------------------------------------

real(SHR_KIND_R8) FUNCTION eshr_timemgr_ETimeGetDiff( ETime, ETime1 )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Time), intent(IN)  :: ETime     ! Input ESMF time
    type(ESMF_Time), intent(IN)  :: ETime1    ! 2nd Input ESMF time

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_ETimeGetDiff) '
    type(ESMF_TimeInterval) :: ETimeInt ! ESMF time interval between the two dates
    integer                 :: rc       ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ETimeInt = ETime - ETime1
    call ESMF_TimeIntervalGet( ETimeInt, d_r8=eshr_timemgr_ETimeGetDiff, rc=rc )
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   ": Error from ESMF_TimeIntervalGet" )

END FUNCTION eshr_timemgr_ETimeGetDiff

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_NMLinfoOrbReconcl  -- Reconcile orbit infor in clockNMLinfo
!   
! !DESCRIPTION:
!   
!     Reconcile the orbit information in the clockNMLinfo object.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_NMLinfoOrbReconcl( clockNMLinfo )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_NMLinfoType), intent(INOUT) :: clockNMLinfo

!EOP

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( clockNMLinfo%orb_iyear_AD /= SHR_ORB_UNDEF_INT )then
      clockNMLinfo%orb_eccen = SHR_ORB_UNDEF_REAL
      clockNMLinfo%orb_mvelp = SHR_ORB_UNDEF_REAL
      clockNMLinfo%orb_obliq = SHR_ORB_UNDEF_REAL
   end if

END SUBROUTINE eshr_timemgr_NMLinfoOrbReconcl

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_NMLinfoPrint -- Private routine to print the clock clockNMLinfo
!   
! !DESCRIPTION:
!   
!     Print the clock clockNMLinfo values to make sure they are correct.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_NMLinfoPrint( clockNMLinfo )

! !USES:

    USE shr_orb_mod,   ONLY : SHR_ORB_UNDEF_INT
    USE shr_orb_mod,   ONLY : shr_orb_print

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_NMLinfoType),  intent(IN) :: clockNMLinfo ! clock clockNMLinfo

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    write(6,F0A)   'Calendar                      ', trim(clockNMLinfo%calendar)
    write(6,F0Date)'Start date                    ', clockNMLinfo%start_ymd, &
                  clockNMLinfo%start_tod
    if ( clockNMLinfo%ref_ymd > 0 )then
       write(6,F0Date) 'Reference date                ', clockNMLinfo%ref_ymd, &
                  clockNMLinfo%ref_tod
    end if
    if (      trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNSteps &
    )then
       write(6,F0I) 'Stop NSteps                   ', clockNMLinfo%stop_n
    else if ( trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNDays  &
    )then
       write(6,F0I) 'Stop NDays                    ', clockNMLinfo%stop_n
    else if ( trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNMons &
    )then
       write(6,F0I) 'Stop nmonths                  ', clockNMLinfo%stop_n
    else if ( trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNYears &
    )then
       write(6,F0I) 'Stop NYears                   ', clockNMLinfo%stop_n
    else if ( trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionDate &
    )then
       write(6,F0Date) 'Stop date                     ', clockNMLinfo%stop_ymd, &
                  clockNMLinfo%stop_tod
    end if
    write(6,F0Date)    'Stop final date               ', &
                       clockNMLinfo%stop_final_ymd, 0
    if (      trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNSteps   &
    )then
       write(6,F0I) 'Restart NSteps                ', clockNMLinfo%restart_n
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNDays    &
    )then
       write(6,F0I) 'Restart NDays  ', clockNMLinfo%restart_n
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optMonthly  &
    )then
       write(6,FA) 'Restart monthly at first of month'
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNMonths  &
    )then
       write(6,F0I) 'Restart nmonths  ', clockNMLinfo%restart_n, ' from current day'
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optYearly   &
    )then
       write(6,FA) 'Restart yearly at first of year'
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNYears   &
    )then
       write(6,F0I) 'Restart NYears   ', clockNMLinfo%restart_n, ' from current day'
    end if
    write(6,F0I)    'Atmosphere coupling interval  ', clockNMLinfo%atm_cpl_dt
    write(6,F0I)    'Land coupling interval        ', clockNMLinfo%lnd_cpl_dt
    write(6,F0I)    'Sea-ice coupling interval     ', clockNMLinfo%ice_cpl_dt
    write(6,F0I)    'Ocean coupling interval       ', clockNMLinfo%ocn_cpl_dt
    if ( clockNMLinfo%perpetual_run )then
       write(6,F0Date) 'Perpetual date                ', &
                       clockNMLinfo%perpetual_ymd, 0
    end if
    if ( .not. clockNMLinfo%orb_notSet )then
       call shr_orb_print( ClockNMLinfo%orb_iyear_AD, ClockNMLinfo%orb_eccen, &
                           ClockNMLinfo%orb_obliq, ClockNMLinfo%orb_mvelp )
    end if
    if ( clockNMLinfo%restart )then
       if ( len_trim(clockNMLinfo%restart_file_TMGoverRide) == 0 )then
          write(6,F0A) 'Do NOT over-ride any settings from the restart_file'
       else
          write(6,F0A) 'Override the following list of items on the '// &
                     'restart_file from the namelist: '
          write(6,F0A) 'Override list                 ', &
                     trim(clockNMLinfo%restart_file_TMGoverRide)
       end if
    end if

END SUBROUTINE eshr_timemgr_NMLinfoPrint

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoPrint -- Print clock information out
!   
! !DESCRIPTION:
!   
!      Print clock information out.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInfoPrint( ClockInfo )

  use shr_orb_mod,   only: shr_orb_print

   implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(eshr_timemgr_clockInfoType), intent(in) :: ClockInfo  ! Input clock info

!EOP

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   integer :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
   character(len=*), parameter :: F0I = "('"//LogPrefix//"',A,T45,' = ', I8.8,A)"

   if (ClockInfo%MasterSyncClock )then
      write(6,FA) 'This is a Master Synchronization clock'
   else
      write(6,FA) 'This is an internal component clock'
   end if
   write(6,FA) trim(clockInfo%desc)
   if ( ClockInfo%perpetual_run )then
      call eshr_timemgr_clockInfoGet( clockInfo, perpetual_ymd=perpetual_ymd )
      write(6,F0I) 'Perpetual date: ', perpetual_ymd, " (YYYYMMDD)"
   end if
   if ( ClockInfo%NoRestarts )then
      write(6,FA) 'No restarts will be written (unless explicitly triggered)'
   end if
   write(6,FAA) 'Calendar type', trim(ClockInfo%calendar)
   if ( ClockInfo%orb_eccen /= SHR_ORB_UNDEF_REAL .and. &
        ClockInfo%orb_obliq /= SHR_ORB_UNDEF_REAL .and. &
        ClockInfo%orb_mvelp /= SHR_ORB_UNDEF_REAL )then
      call shr_orb_print( ClockInfo%orb_iyear_AD, ClockInfo%orb_eccen, &
                          ClockInfo%orb_obliq, ClockInfo%orb_mvelp )
      if ( ClockInfo%orb_obliqr /= SHR_ORB_UNDEF_REAL ) &
      write(6,F0F) 'orb_obliqr ', ClockInfo%orb_obliqr, " (radians)"
      if ( ClockInfo%orb_lambm0 /= SHR_ORB_UNDEF_REAL ) &
      write(6,F0F) 'orb_lambm0 ', ClockInfo%orb_lambm0, " (radians)"
      if ( ClockInfo%orb_mvelpp /= SHR_ORB_UNDEF_REAL ) &
      write(6,F0F) 'orb_mvelpp ', ClockInfo%orb_mvelpp, " (radians)"
   else if ( ClockInfo%orb_eccen  /= SHR_ORB_UNDEF_REAL &
      .and.  ClockInfo%orb_obliqr /= SHR_ORB_UNDEF_REAL &
      .and.  ClockInfo%orb_lambm0 /= SHR_ORB_UNDEF_REAL &
      .and.  ClockInfo%orb_mvelpp /= SHR_ORB_UNDEF_REAL )then
      write(6,F0F) 'orb_eccen  ', ClockInfo%orb_eccen,  " (unit-less)"
      write(6,F0F) 'orb_obliqr ', ClockInfo%orb_obliqr, " (radians)"
      write(6,F0F) 'orb_lambm0 ', ClockInfo%orb_lambm0, " (radians)"
      write(6,F0F) 'orb_mvelpp ', ClockInfo%orb_mvelpp, " (radians)"
   else
      write(6,FA) 'Orbit information not set'
   end if

END SUBROUTINE eshr_timemgr_clockInfoPrint

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoIsMaster -- Return true if this is a Master clock
!   
! !DESCRIPTION:
!   
!      Return true if the input clockInfo object is a Master Synchronization clock.
!      
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION eshr_timemgr_clockInfoIsMaster( ClockInfo )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(eshr_timemgr_clockInfoType), intent(in) :: ClockInfo  ! Input clock info

!EOP

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   eshr_timemgr_clockInfoIsMaster = ClockInfo%MasterSyncClock

END FUNCTION eshr_timemgr_clockInfoIsMaster

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoNoRest --- Return true if no restarts are done
!   
! !DESCRIPTION:
!   
!      Return true if restarts will NOT be done.
!      
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION eshr_timemgr_clockInfoNoRest( ClockInfo )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(eshr_timemgr_clockInfoType), intent(in) :: ClockInfo  ! Input clock info

!EOP

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   eshr_timemgr_clockInfoNoRest = ClockInfo%NoRestarts

END FUNCTION eshr_timemgr_clockInfoNoRest

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoIsPerpet --- Return true if in perpetual mode
!   
! !DESCRIPTION:
!   
!      Return true if in perpetual mode
!      
! !INTERFACE: ------------------------------------------------------------------

logical FUNCTION eshr_timemgr_clockInfoIsPerpet( ClockInfo )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(eshr_timemgr_clockInfoType), intent(in) :: ClockInfo  ! Input clock info

!EOP

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   eshr_timemgr_clockInfoIsPerpet = ClockInfo%perpetual_run

END FUNCTION eshr_timemgr_clockInfoIsPerpet

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoInitRest -- setup clockInfo restart information
!   
! !DESCRIPTION:
!   
!     Setup clockInfo restart information. Private to inside this module.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInfoInitRest( clockInfo )

! !USES:

   use shr_ncio_mod,   only: shr_ncio_descripInit, shr_ncio_descripSetDefault

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockInfoType), intent(INOUT) :: clockInfo    ! Input clockInfo object

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockInfoInitRest) '
    integer :: n    ! Index

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Clock info type restart variables
    !---------------------------------------------------------------------------
    call shr_ncio_descripSetDefault( nInfoVars, clockInfo%var )
    do n = 1, nInfoVars
       selectcase ( trim(ClockInfoSave(n)) )
          case( "desc" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n), &
                                   LongName="Description of this clock",     &
                                   StringData=.true. )
          case( "MasterSyncClock" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n), &
                                   LongName="Master synchronization clock"// &
                                            " or internal model clock",      &
                                   LogicalData=.true. )
          case( "calendar" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n), &
                                   LongName="Type of calendar used",         &
                                   units=eshr_timemgr_noLeap//" or "//        &
                                   eshr_timemgr_Gregorian, StringData=.true. )
          case( "perpetual_run" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n),   &
                                   LongName="Clock info flag for perpetual "// &
                                            "run mode",                        &
                                   LogicalData=.true. )
          case( "perpetual_ymd" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n),   &
                                   LongName="Clock info perpetual date",       &
                                   units = "date_[YYYYMMDD]",                  &
                                   IntegerData =.true.,                        &
                                   IntegerFill=0 )
          case( "orb_mode" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n), &
                                   LongName="Clock info orbital mode",       &
                                   units = "list", IntegerData =.true.,      &
                                   IntegerFill = SHR_ORB_UNDEF_INT,          &
                                   ListDescrips=eshr_timemgr_orbModeDescrips, &
                                   ListIntValues=eshr_timemgr_orbModeList )
          case( "orb_iyear_AD" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n), &
                                   LongName="Clock info orbital year",       &
                                   units = "Gregorian_year",                 &
                                   IntegerData =.true., &
                                   IntegerFill = SHR_ORB_UNDEF_INT )
          case( "orb_obliq" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n), &
                                   LongName="Clock info orbital obliquity",  &
                                   units = "degrees", RealR8Data=.true.,     &
                                   RealR8Fill = SHR_ORB_UNDEF_REAL )
          case( "orb_eccen" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n),   &
                                   LongName="Clock info orbital eccentricity", &
                                   units = "ratio[0-.1]", RealR8Data=.true.,   &
                                   RealR8Fill = SHR_ORB_UNDEF_REAL )
          case( "orb_mvelp" )
             call shr_ncio_descripInit( clockInfo%var(n), ClockInfoSave(n), &
                                   LongName="Clock info orbital moving "//   &
                                   "vernal equinox at perihelion",           &
                                   units = "degrees", RealR8Data=.true.,     &
                                   RealR8Fill = SHR_ORB_UNDEF_REAL )
          case default
             call shr_sys_abort( subname//': Unknown data type to save to '// &
                                 ' clock info restart file:'//&
                                 trim(ClockInfoSave(n)) )
       end select
    end do

END SUBROUTINE eshr_timemgr_clockInfoInitRest

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoRestRead -- Read in the clockInfo restart file
!   
! !DESCRIPTION:
!   
! Read in clockInfo restart information from given input netCDF file
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInfoRestRead( restart_file, MPICom,   &
                                           MasterTask, ClockInfo,  &
                                           clockNMLinfo )
! !USES:

   USE shr_ncio_mod,   ONLY : shr_ncio_descripRead, shr_ncio_descripName,     &
                              shr_ncio_descripGetInteger,                     &
                              shr_ncio_descripGetRealR8,                      &
                              shr_ncio_descripGetString,                      &
                              shr_ncio_descripGetLogical
   USE shr_string_mod, ONLY : shr_string_listGetIndexF

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(SHR_KIND_CL), intent(IN) :: restart_file  ! Restart local filename
   integer, optional,      intent(IN) :: MPICom        ! MPI communicator
   logical, optional,      intent(IN) :: MasterTask    ! If root MPI task or not
   type(eshr_timemgr_clockInfoType), intent(INOUT) :: clockInfo    ! Input Clock info type
   type(eshr_timemgr_NMLinfoType),   intent(INOUT) :: clockNMLinfo ! Input clock NML info

!EOP

   !----- local -----
   character(len=*), parameter :: subname = '(eshr_timemgr_clockInfoRestRead) '
   logical :: MasterTask2                   ! If root MPI task or not
   integer :: i                             ! Index
   integer :: rc                            ! Return code
   integer :: perpetual_ymd                 ! Perpetual date (YYYYMMDD)
   integer :: orb_mode                      ! Orbital mode
   integer :: nSet                          ! Number of values set
   integer :: list                          ! List index of item in overRide list
   integer :: ncId                          ! NetCDF file id
   logical :: exists                        ! If file exists or not

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( present(MasterTask) )then
      MasterTask2 = MasterTask
   else
      MasterTask2 = .true.
   end if
   ! --- Setup restart information --------------------------------------------
   call eshr_timemgr_clockInfoInitRest( ClockInfo )

   ! --- Read in restart file -------
   call shr_ncio_open( restart_file, MasterTask2, FileType=prefix//"restart_file", &
                       ncId=ncId, exists=exists )
   if ( present(MPICom) )then
      call shr_ncio_descripRead( ncId, NInfoVars,  prefix//"ClockInfo_",     &
                                 MPICom, MasterTask2, var=ClockInfo%var )
   else
      call shr_ncio_descripRead( ncId, NInfoVars,  prefix//"ClockInfo_",     &
                                 var=ClockInfo%var )
   end if
   call shr_ncio_close( ncId, MasterTask2, type=prefix//"restart_file", &
                        NCFileName=restart_file )
   !----------------------------------------------------------------------------
   ! Put information read into derived type for clock-info variables
   ! Check if item is on over-ride list (if list_index != 0)
   ! If on over-ride list use value from clockNMLinfo else use the value from the restart
   !----------------------------------------------------------------------------
   do i = 1, NInfoVars
      if (      trim(shr_ncio_descripName(ClockInfo%var(i))) == "desc"   )then
         ClockNMLinfo%desc = shr_ncio_descripGetString( ClockInfo%var(i) )
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "calendar"   )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                    trim(shr_ncio_descripName(ClockInfo%var(i)) ) )
         if ( list == 0 )then
            clockNMLinfo%calendar = &
                      trim(shr_ncio_descripGetString( ClockInfo%var(i) ))
         end if
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "perpetual_run" &
      )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                          shr_ncio_descripName(ClockInfo%var(i)) )
         if ( list == 0 )then
            clockNMLinfo%perpetual_run = &
                  shr_ncio_descripGetLogical( ClockInfo%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "MasterSyncClock" &
      )then
         ClockNMLinfo%MasterSyncClock = shr_ncio_descripGetLogical( &
                                             ClockInfo%var(i) )
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "perpetual_ymd" &
      )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                trim(shr_ncio_descripName(ClockInfo%var(i)) ) )
         if ( list == 0 )then
            clockNMLinfo%perpetual_ymd = &
                           shr_ncio_descripGetInteger( ClockInfo%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "orb_mode" &
      )then
         orb_mode = shr_ncio_descripGetInteger( ClockInfo%var(i) )
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "orb_iyear_AD"  &
      )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                 trim(shr_ncio_descripName(ClockInfo%var(i)) ) )
         if ( list == 0 )then
            clockNMLinfo%orb_iyear_AD = &
                   shr_ncio_descripGetInteger( ClockInfo%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "orb_eccen" &
      )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                          shr_ncio_descripName(ClockInfo%var(i) ) )
         if ( list == 0 )then
            clockNMLinfo%orb_eccen =  &
                     shr_ncio_descripGetRealR8( ClockInfo%var(i) )
         end if
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "orb_obliq"    &
      )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                          shr_ncio_descripName(ClockInfo%var(i) ) )
         if ( list == 0 )then
            clockNMLinfo%orb_obliq = &
                  shr_ncio_descripGetRealR8( ClockInfo%var(i) ) 
         end if
      else if ( trim(shr_ncio_descripName(ClockInfo%var(i))) == "orb_mvelp"    &
      )then
         list = shr_string_listGetIndexF( clockNMLinfo%restart_file_TMGoverRide, &
                                          shr_ncio_descripName(ClockInfo%var(i) ) )
         if ( list == 0 )then
            clockNMLinfo%orb_mvelp = &
                  shr_ncio_descripGetRealR8( ClockInfo%var(i) )
         end if
      else
         call shr_sys_abort( subname// &
                             'ERROR: unknown clock info variable name: '// &
                             trim(shr_ncio_descripName(ClockInfo%var(i) ) ) )
      end if
   end do
   call eshr_timemgr_NMLinfoOrbReconcl( clockNMLinfo )

END SUBROUTINE eshr_timemgr_clockInfoRestRead       

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoRestWrit -- Write out the clockInfo restart file
!   
! !DESCRIPTION:
!   
! Write time-manager clockInfo information out to input netCDF file
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInfoRestWrit( restart_file, MPICom, MasterTask, &
                                           clockInfo )
! !USES:

   use shr_ncio_mod,  only: shr_ncio_descripName, shr_ncio_descripPutData, &
                            shr_ncio_descripWrite

   implicit none

! !INPUT/OUTPUT PARAMETERS:
   character(SHR_KIND_CL),           intent(IN)    :: restart_file ! Restart filename
   integer, optional,                intent(IN)    :: MPICom       ! MPI communicator
   logical, optional,                intent(IN)    :: MasterTask   ! If root MPI task
   type(eshr_timemgr_clockInfoType), intent(INOUT) :: clockInfo    ! Input clock

!EOP

  !----- local -----
  character(len=*), parameter :: subname = '(eshr_timemgr_clockRestWrite) '
  logical :: MasterTask2                ! If root MPI task or not
  integer :: perpetual_ymd              ! Perpetual date (YYYYMMDD)
  integer :: perpetual_tod              ! Perpetual time of day (seconds)
  integer :: i                          ! Index
  integer :: ncId                       ! NetCDF file id
  logical :: exists                     ! If file exists or not

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  if ( present(MasterTask) )then
     MasterTask2 = MasterTask
  else
     MasterTask2 = .true.
  end if

  !-----------------------------------------------------------------------------
  ! Loop over the info variables, set the NetCDF description data
  !-----------------------------------------------------------------------------
  do i = 1, NInfoVars
      if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) == "perpetual_ymd" &
      )then
         if ( eshr_timemgr_clockInfoIsPerpet(  clockInfo ) )then
            perpetual_ymd = eshr_timemgr_ETimeGetDataYMD( &
                                                clockInfo%perpetual_time )
            call shr_ncio_descripPutData( clockInfo%var(i), "perpetual_ymd",    &
                                      IntegerData=perpetual_ymd )
         end if
      else if ( trim(shr_ncio_descripName( clockInfo%var(i) ) ) == "MasterSyncClock" &
      )then
         call shr_ncio_descripPutData( clockInfo%var(i), "MasterSyncClock",    &
                                LogicalData=clockInfo%MasterSyncClock )
      else if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) == "desc" )then
         call shr_ncio_descripPutData( clockInfo%var(i), "desc", &
                                       StringData=clockInfo%desc )
      else if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) == "calendar" &
      )then
         call shr_ncio_descripPutData( clockInfo%var(i), "calendar", &
                                       StringData=clockInfo%calendar )
      else if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) =="perpetual_run"&
      )then
         call shr_ncio_descripPutData( clockInfo%var(i), "perpetual_run", &
                                       LogicalData=clockInfo%perpetual_run )
      else if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) == "orb_mode" &
      )then
         call shr_ncio_descripPutData( clockInfo%var(i), "orb_mode", &
                                       IntegerData=clockInfo%orb_mode )
      else if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) =="orb_iyear_AD"&
      )then
         call shr_ncio_descripPutData( clockInfo%var(i), "orb_iyear_AD", &
                                       IntegerData=clockInfo%orb_iyear_AD )
      else if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) == "orb_obliq" &
      )then
         call shr_ncio_descripPutData( clockInfo%var(i), "orb_obliq", &
                                       RealR8Data=clockInfo%orb_obliq )
      else if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) == "orb_eccen" &
      )then
         call shr_ncio_descripPutData( clockInfo%var(i), "orb_eccen", &
                                       RealR8Data=clockInfo%orb_eccen )
      else if (      trim(shr_ncio_descripName( clockInfo%var(i) ) ) == "orb_mvelp" &
      )then
         call shr_ncio_descripPutData( clockInfo%var(i), "orb_mvelp", &
                                       RealR8Data=clockInfo%orb_mvelp )
      else
         call shr_sys_abort( subname//': unknown info variable name: '// &
                             trim(shr_ncio_descripName( clockInfo%var(i) ) ) )
      end if
  end do
  ! --- Write the info variables out to the restart file ------
  call shr_ncio_open( NCFileName=restart_file, MasterTask=masterTask2,  &
                      FileType=trim(prefix)//"restart_file",            &
                      ncId=ncId, exists=exists, writing=.true. )
  if ( present(MPICom) )then
     call shr_ncio_descripWrite( ncId, NInfoVars,  prefix//"ClockInfo_", &
                                 mpicom, MasterTask2, exists=exists,     &
                                 var=clockInfo%var )
  else
     call shr_ncio_descripWrite( ncId, NInfoVars,  prefix//"ClockInfo_", &
                                 exists=exists, var=clockInfo%var )
  end if
  call shr_ncio_close( ncId, MasterTask2, type=prefix//"restart_file",  &
                       NCFileName=restart_file )

END SUBROUTINE eshr_timemgr_clockInfoRestWrit      

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoSetOrb -- Set the orbital mode based on ClockInfo
!   
! !DESCRIPTION:
!   
! Private method:
!
! Set the orbital mode in the clock info type
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInfoSetOrb( ClockInfo )

! !USES:

    USE shr_string_mod, only : shr_string_listGetIndexF, shr_string_listGetNum
    USE shr_orb_mod,    only : shr_orb_params

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_clockInfoType), intent(INOUT) :: ClockInfo

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockInfoSetOrb) '
    integer            :: list               ! Index
    integer            :: nList              ! Number of elements in list
    logical, parameter :: LogPrint = .false. ! Do not print out log info

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- Orbit is determined by fixed orbit characteristics -----
    if ( ClockInfo%orb_iyear_AD == SHR_ORB_UNDEF_INT )then
       ClockInfo%orb_mode = shr_string_listGetIndexF(             &
                                    eshr_timemgr_orbModeDescrips, &
                                    eshr_timemgr_orb_mode_fixed )
    ! --- Orbit is determined by fixed year -----
    else
       ClockInfo%orb_mode = shr_string_listGetIndexF(                &
                                    eshr_timemgr_orbModeDescrips,    &
                                    eshr_timemgr_orb_mode_fixed_yr )
    end if
 
    !---------------------------------------------------------------------------
    ! Calculate index values of each orbit mode descriptions
    ! This is used for the output description of the orbit mode on the restart file
    !---------------------------------------------------------------------------
    nList = shr_string_listGetNum( eshr_timemgr_orbModeDescrips )
    if ( NList /= NOrbList )then
       call shr_sys_abort( subname//' Number of orbModeDescrip list '// &
                           'different than NOrbList' )
    end if
    !
    ! ---- Set the orbital parameters ----
    !
    if (  (ClockInfo%orb_obliqr == SHR_ORB_UNDEF_REAL) &
    .and. (ClockInfo%orb_lambm0 == SHR_ORB_UNDEF_REAL) &
    .and. (ClockInfo%orb_mvelpp == SHR_ORB_UNDEF_REAL) ) &
       call shr_orb_params( ClockInfo%orb_iyear_AD, ClockInfo%orb_eccen,  &
                            ClockInfo%orb_obliq,    ClockInfo%orb_mvelp,  &
                            ClockInfo%orb_obliqr,   ClockInfo%orb_lambm0, &
                            ClockInfo%orb_mvelpp,   LogPrint )

END SUBROUTINE eshr_timemgr_clockInfoSetOrb

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockInfoInit -- Initialize the clockInfo object
!   
! !DESCRIPTION:
!   
! Private method:
!
! Setup the clockInfo object
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockInfoInit( clockNMLinfo, desc, MasterSyncClock, &
                                       ClockInfoIn, ClockInfoOut )
! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_NMLinfoType), intent(IN), optional :: clockNMLinfo ! input info obj
    character(len=*), optional,     intent(IN) :: desc            ! Clock description
    logical,          optional,     intent(IN) :: MasterSyncClock ! If Master Sync. clock
    ! --- Input clockInfo to copy in -------------------------------------------
    type(eshr_timemgr_clockInfoType), intent(IN), optional :: ClockInfoIn 
    type(eshr_timemgr_clockInfoType), intent(OUT) :: ClockInfoOut   ! Output clockInfo

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_clockInfoInit) '
    logical                     :: orb_notSet           ! Orbit not set yet

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(ClockInfoIn) )then
       ! --- Description -----
       ClockInfoOut%calendar   = clockInfoIn%calendar
       ClockInfoOut%NoRestarts = clockInfoIn%NoRestarts
       ! --- Perpetual information -----
       ClockInfoOut%perpetual_run  = clockInfoIn%perpetual_run
       if ( clockInfoOut%perpetual_run ) &
           ClockInfoOut%perpetual_time = clockInfoIn%perpetual_time
       ! --- Orbital information -----
       ClockInfoOut%orb_mvelp       = clockInfoIn%orb_mvelp
       ClockInfoOut%orb_eccen       = clockInfoIn%orb_eccen
       ClockInfoOut%orb_obliq       = clockInfoIn%orb_obliq
       ClockInfoOut%orb_iyear_AD    = clockInfoIn%orb_iyear_AD
       orb_notSet                   = .false.
       ! ----- Description -----------------------------------------------------
       if ( present(desc) )then
          ClockInfoOut%desc = desc
       else
          call shr_sys_abort( subname//': desc NOT present and MUST be' )
       end if
       if ( present(MasterSyncClock) )then
          ClockInfoOut%MasterSyncClock = MasterSyncClock
       else
          call shr_sys_abort( subname//': MasterSyncClock NOT present and MUST be' )
       end if
    else
       ! --- Perpetual information ------------------------------------------------
       ClockInfoOut%perpetual_run = clockNMLinfo%perpetual_run
       if ( clockNMLinfo%perpetual_run ) &
          ClockInfoOut%perpetual_time = eshr_timemgr_ETimeInit( &
                                            clockNMLinfo%perpetual_ymd, 0, &
                                            "Perpetual date" )
       ! ----- Calendar --------------------------------------------------------
       ClockInfoOut%calendar        = clockNMLinfo%calendar
       ! --- Orbital information -----------------------------------------------
       orb_notSet                 = clockNMLinfo%orb_notSet
       ClockInfoOut%orb_mvelp     = clockNMLinfo%orb_mvelp
       ClockInfoOut%orb_eccen     = clockNMLinfo%orb_eccen
       ClockInfoOut%orb_obliq     = clockNMLinfo%orb_obliq
       ClockInfoOut%orb_iyear_AD  = clockNMLinfo%orb_iyear_AD
       ! ----- Description -----------------------------------------------------
       ClockInfoOut%desc            = clockNMLinfo%desc
       ClockInfoOut%MasterSyncClock = clockNMLinfo%MasterSyncClock
    end if
    ! ----- Orbit ----
    if ( .not. orb_notSet ) &
       call eshr_timemgr_clockInfoSetOrb( ClockInfoOut )
    ! --- Setup restart information --------------------------------------------
    call eshr_timemgr_clockInfoInitRest( ClockInfoOut )

END SUBROUTINE eshr_timemgr_clockInfoInit

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockAlarmRstInit -- Set restart alarm from ClockclockNMLinfo
!   
! !DESCRIPTION:
!   
!     Setup the restart alarm from the input ClockclockNMLinfo
!      
! !INTERFACE: ------------------------------------------------------------------

subroutine eshr_timemgr_clockAlarmRstInit( clockNMLinfo, ClockOut )

    implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(eshr_timemgr_NMLinfoType), intent(IN)    :: clockNMLinfo ! clock clockNMLinfo
    type(eshr_timemgr_clockType),   intent(INOUT) :: ClockOut     ! Output clock

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_setRestartFreqFromSetup) '
    integer :: FirstYMD                       ! First alarm time date (YMD)
    integer :: FirstTOD                       ! First alarm time of day (seconds)
    integer :: rc                             ! Return code
    type(ESMF_Time) :: FirstAlarm             ! First restart alarm time
    type(ESMF_Time) :: CurrentTime            ! Current time
    type(ESMF_TimeInterval) :: AlarmInterval  ! Alarm interval

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- If no restarts ------
    if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNone )then
       call eshr_timemgr_clockInfoPutData( clockOut%info, NoRestarts=.true. )
    ! --- If restarts every restart_n steps ------
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNSteps &
    )then
       call ESMF_TimeIntervalSet( AlarmInterval, s=clockNMLinfo%atm_cpl_dt, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalSet" )
       AlarmInterval = AlarmInterval * clockNMLinfo%restart_n
       FirstYMD = clockNMLinfo%start_ymd
       FirstTOD = clockNMLinfo%start_tod
    ! --- If restarts every restart_n days ------
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNDays &
    )then
       call ESMF_TimeIntervalSet( AlarmInterval, d=1, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalSet" )
       AlarmInterval = AlarmInterval * clockNMLinfo%restart_n
       FirstYMD = clockNMLinfo%start_ymd
       FirstTOD = clockNMLinfo%start_tod
    ! --- If restarts are monthly at the beginning of the month ----------------
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optMonthly &
    )then
       FirstYMD = int(clockNMLinfo%start_ymd / 100)*100 + 1
       FirstTOD = 0
       call ESMF_TimeIntervalSet( AlarmInterval, mm=1, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalSet" )
    ! --- If restarts every restart_n months ------
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNMonths &
    )then
       call ESMF_TimeIntervalSet( AlarmInterval, mm=1, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalSet" )
       AlarmInterval = AlarmInterval * clockNMLinfo%restart_n
       FirstYMD = clockNMLinfo%start_ymd
       FirstTOD = clockNMLinfo%start_tod
    ! --- If restarts are yearly at the beginning of the year ------------------
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optYearly &
    )then
       FirstYMD = int(clockNMLinfo%start_ymd / 10000)*10000 + 101
       FirstTOD = 0
       call ESMF_TimeIntervalSet( AlarmInterval, yy=1, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalSet" )
    ! --- If restarts every restart_n years ------
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optNYears &
    )then
       call ESMF_TimeIntervalSet( AlarmInterval, yy=1, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalSet" )
       AlarmInterval = AlarmInterval * clockNMLinfo%restart_n
       FirstYMD = clockNMLinfo%start_ymd
       FirstTOD = clockNMLinfo%start_tod
    ! --- If restarts only at end of simulation ------
    else if ( trim(clockNMLinfo%restart_option) == eshr_timemgr_rest_optEnd &
    )then
       call eshr_timemgr_clockGet( ClockOut, stop_ymd=FirstYMD, stop_tod=FirstTOD )
       call ESMF_TimeIntervalSet( AlarmInterval, yy=9999, rc=rc )
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_TimeIntervalSet" )
    end if
    ! --- If restarts are done at some frequency ------
    if ( .not. eshr_timemgr_clockInfoNoRest( ClockOut%info ) )then
       if ( clockNMLinfo%NextRestYMD /= 0 )then
          FirstYMD = clockNMLinfo%NextRestYMD
          FirstTOD = clockNMLinfo%NextRestTOD
       end if

       FirstAlarm = eshr_timemgr_ETimeInit( FirstYMD, FirstTOD, &
                                             "First restart alarm" )
       call eshr_timemgr_clockGet( ClockOut, ECurrTime=CurrentTime )
       if ( FirstAlarm <= CurrentTime ) FirstAlarm = FirstAlarm + AlarmInterval
       ClockOut%restart = ESMF_AlarmCreate( name=eshr_timemgr_restartAlarmName, &
                                            clock=ClockOut%EClock,             &
                                            ringTime=FirstAlarm,               &
                                            ringInterval=AlarmInterval, rc=rc)
       call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                      "Error from ESMF_AlarmCreate" )
    end if

end subroutine eshr_timemgr_clockAlarmRstInit

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockSetupRest -- setup clock restart information
!   
! !DESCRIPTION:
!   
!     Setup clock restart information. Private to inside this module.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockSetupRest( clock )

! !USES:

   use shr_ncio_mod,   only: shr_ncio_descripInit, shr_ncio_descripSetDefault

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(INOUT) :: clock    ! Input CCSM clock

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_NMLinfoRestart) '
    integer :: n    ! Index

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Main clock restart variables
    !---------------------------------------------------------------------------
    call shr_ncio_descripSetDefault( nClockVars, clock%var )
    do n = 1, nClockVars
       selectcase (trim(ClockSave(n)) )
          case( "start_ymd" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),           &
                                         LongName="Clock start date",          &
                                         units = "date [YYYYMMDD]",            &
                                         IntegerData =.true. )
          case( "start_tod" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),           &
                                         LongName="Clock start time of day",   &
                                         units = "seconds", IntegerData =.true. )
          case( "ref_ymd" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),           &
                                         LongName="Clock reference date",      &
                                         units = "date [YYYYMMDD]",            &
                                         IntegerData =.true. )
          case( "ref_tod" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),             &
                                         LongName="Clock reference time of day", &
                                         units = "seconds", IntegerData =.true. )
          case( "CurrentYMD" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),           &
                                         LongName="Clock current date",        &
                                         units = "date [YYYYMMDD]",            &
                                         IntegerData =.true. )
          case( "CurrentTOD" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),           &
                                         LongName="Clock current time of day", &
                                         units = "seconds", IntegerData =.true. )
          case( "DTime" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),           &
                                         LongName="Clock time-step",           &
                                         units = "seconds", IntegerData =.true. )
          case( "RestartIntervalSec" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),              &
                                         LongName="Clock restart alarm interval", &
                                         units = "seconds", IntegerData =.true.,  &
                                         IntegerFill=0 )
          case( "RestartIntervalMonths" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),              &
                                         LongName="Clock restart alarm interval", &
                                         units = "months", IntegerData =.true.,   &
                                         IntegerFill=0 )
          case( "RestartIntervalYears" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),              &
                                         LongName="Clock restart alarm interval", &
                                         units = "years", IntegerData =.true.,    &
                                         IntegerFill=0 )
          case( "RestartNextAlarmYMD" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),               &
                                         LongName="Restart alarm next alarm date", &
                                         units = "date [YYYYMMDD]",                &
                                         IntegerData =.true., &
                                         IntegerFill=0 )
          case( "RestartNextAlarmTOD" )
              call shr_ncio_descripInit( clock%var(n), ClockSave(n),               &
                                         LongName="Restart alarm next alarm "//    &
                                                  "time of day",                   &
                                         units = "seconds", IntegerData =.true.,   &
                                         IntegerFill=-1 )
          case default
             call shr_sys_abort( subname//': Unknown data type to save to '// &
                                 ' clock restart file:'//trim(ClockSave(n)) )
       end select
    end do

END SUBROUTINE eshr_timemgr_clockSetupRest

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockGetESMFClock -- Get ESMF_Clock from a clock
!   
! !DESCRIPTION:
!   
!     Get an ESMF clock from a CCSM share clock.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockGetESMFClock( clock, ESMFClockOut )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(IN) :: clock         ! Input CCSM clock
   type(ESMF_Clock),             intent(OUT)::  ESMFClockOut ! Output ESMF Clock object

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ESMFClockOut = clock%EClock

END SUBROUTINE eshr_timemgr_clockGetESMFClock

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_clockGetInfo  -- Get the ClockInfo from the shared clock
!   
! !DESCRIPTION:
!   
! Get the information part of the shared clock
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_clockGetInfo( clock, ClockInfo )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(eshr_timemgr_clockType), intent(IN) :: clock          ! Input clock
   type(eshr_timemgr_clockInfoType), intent(OUT) :: ClockInfo ! clock info

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   ClockInfo = clock%info

END SUBROUTINE eshr_timemgr_clockGetInfo

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_EClockInit -- Initialize the ESMF clock in the shared clock
!   
! !DESCRIPTION:
!   
! Private method:
!
! Setup the ESMF clock inside the wrapped CCSM clock
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_EClockInit( TimeStep, StartTime, RefTime,              &
                                    CurrentTime, desc, clockNMLinfo, StopTime, &
                                    ESMFClockOut )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_TimeInterval), intent(IN) :: TimeStep            ! Time-step of clock
    type(ESMF_Time), intent(IN) :: StartTime                   ! Start time
    type(ESMF_Time), intent(IN) :: RefTime                     ! Reference time
    type(ESMF_Time), intent(IN) :: CurrentTime                 ! Current time
    character(len=*), intent(IN) :: desc                       ! Description of this clock
    type(eshr_timemgr_NMLinfoType), intent(IN), optional :: clockNMLinfo ! Input info obj.
    type(ESMF_Time), intent(IN), optional :: StopTime          ! Stop time
    type(ESMF_Clock), intent(OUT) :: ESMFClockOut              ! Output ESMF clock

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(eshr_timemgr_EClockInit) '
    integer :: rc                             ! ESMF return code
    character(len=SHR_KIND_CL) :: description ! Description of this clock
    type(ESMF_Time) :: current                ! Current time
    type(ESMF_Time) :: StopDate               ! Stop time
    type(ESMF_Time) :: StopFinalDate          ! Final stop time
    type(ESMF_TimeInterval) :: StopInterval   ! Time-interval until stop-time

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Get stop -time
    !---------------------------------------------------------------------------
    if ( present(clockNMLinfo) )then
       if (      trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNSteps   &
       )then
          StopInterval = TimeStep * clockNMLinfo%stop_n
          StopDate = CurrentTime + StopInterval
       else if ( trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNDays    &
       )then
          call ESMF_TimeIntervalSet( StopInterval, d=clockNMLinfo%stop_n, rc=rc )
          call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                         "Error from ESMF_TimeIntervalSet" )
          StopDate = CurrentTime + StopInterval
       else if ( trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNMons  &
       )then
          call ESMF_TimeIntervalSet( StopInterval, mm=clockNMLinfo%stop_n, rc=rc )
          call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                         "Error from ESMF_TimeIntervalSet" )
          StopDate = CurrentTime + StopInterval
       else if ( trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionNYears   &
       )then
          call ESMF_TimeIntervalSet( StopInterval, yy=clockNMLinfo%stop_n, rc=rc )
          call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                         "Error from ESMF_TimeIntervalSet" )
          StopDate = CurrentTime + StopInterval
       else if ( trim(clockNMLinfo%stop_option) == eshr_timemgr_stop_optionDate   &
       )then
          StopDate = eshr_timemgr_ETimeInit( clockNMLinfo%stop_ymd, &
                     clockNMLinfo%stop_tod, "Stop date" )
       else
          call shr_sys_abort( subname//' stop time not set correctly' )
       end if
       StopFinalDate = eshr_timemgr_ETimeInit( clockNMLinfo%stop_final_ymd, 0, &
                                               "Stop final date" )
       if ( CurrentTime >= StopFinalDate )then
          call shr_sys_abort( subname//' current time or start time '// &
                              'less than or equal to stop_final_ymd date' )
       end if
       if ( StopFinalDate < StopDate ) StopDate = StopFinalDate
    else if ( .not. present(StopTime) )then
       call shr_sys_abort( subname//' clockNMLinfo not present NOR StopTime' )
    else
      StopDate = StopTime
    end if
    description = 'CCSM shared Time-manager clock:'//trim(desc)
    ! ------ Create ESMF Clobk with input characteristics -------------------
    ESMFClockOut = ESMF_ClockCreate(trim(description), &
                                    TimeStep=TimeStep, startTime=StartTime,&
                                    stopTime=StopDate, refTime=RefTime,    &
                                    rc=rc)
    call eshr_timemgr_ErCodeCheck( rc, mes=subname// &
                                   "Error from ESMF_ClockCreate" )
    ! ------ Advance clock to the current time (in case of a restart) -------
     call ESMF_ClockGet(ESMFClockOut, currTime=current, rc=rc )
     call eshr_timemgr_ErCodeCheck(rc, subname// &
                                   ': error return from ESMF_ClockGet' )
     do while( CurrentTime > current )
        call ESMF_ClockAdvance( ESMFClockOut, rc=rc )
        call eshr_timemgr_ErCodeCheck(rc, subname// &
                                      ': error return from ESMF_ClockAdvance' )
        call ESMF_ClockGet( ESMFClockOut, currTime=current )
        call eshr_timemgr_ErCodeCheck(rc, subname// &
                                      ': error return from ESMF_ClockGet' )
     end do

END SUBROUTINE eshr_timemgr_EClockInit

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_ErCodeCheck -- Check return-code from ESMF -- abort if not
!   
! !DESCRIPTION:
!   
!     Check ESMF return code and abort if not successful.
! NOTE: SHOULD THIS BE A eeshr_CheckESMFrc subroutine? <<<<<<<<<<<<<<<<<<<<
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_ErCodeCheck( rc, mes )

! !USES:
   use eshr_rc_mod,   only: eshr_rc_check

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer, intent(in)          :: rc   ! return code from ESMF
   character(len=*), intent(in) :: mes  ! error message

!EOP
!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   call eshr_rc_check( rc, mes )

END SUBROUTINE eshr_timemgr_ErCodeCheck

!===============================================================================
!===============================================================================
! !IROUTINE: eshr_timemgr_initCalendar -- Initialize the calendar that will be used
!   
! !DESCRIPTION:
!   
!     Private method to initialize the calendar that will be used for all time-instants.
!      
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE eshr_timemgr_initCalendar( calendar )

  use shr_string_mod, only: shr_string_toUpper

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(len=*), intent(IN) :: calendar         ! Name of calendar to use

!EOP

   !----- local -----
   character(len=*), parameter :: subname = '(eshr_timemgr_initCalendar) '
   type(ESMF_CalendarType) :: cal_type  ! calendar type
   character(SHR_KIND_CS)  :: caltmp    ! Uppercase of calendar string
   integer                 :: rc        ! return code
   type(ESMF_Calendar)     :: cal       ! Local calendar

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   ! --- Figure out calendar type ------
   caltmp = shr_string_toUpper(calendar)
   if ( trim(caltmp) == eshr_timemgr_noLeap ) then
      cal_type = ESMF_CAL_NOLEAP
   else if ( trim(caltmp) == eshr_timemgr_Gregorian ) then
      cal_type = ESMF_CAL_GREGORIAN
   else
      write(6,FA) subname//': unrecognized calendar specified: '// &
                  trim(calendar)
      call shr_sys_abort( subname//'ERROR:: bad calendar' )
   end if
   ! --- Create the new calendar if not already set ------
   if ( .not. eshr_timemgr_setCalendar )then
      eshr_timemgr_cal = ESMF_CalendarCreate( name=caltmp, &
                                             calendarType=cal_type, rc=rc )
      call eshr_timemgr_ErCodeCheck( rc, subname//': error return '//&
                                     'from ESMF_CalendarSet' )
   else
      ! --- If calendar already created still create it ------
      cal = ESMF_CalendarCreate( name=caltmp, calendarType=cal_type, rc=rc )
      call eshr_timemgr_ErCodeCheck( rc, subname//': error return from '//&
                                     'ESMF_CalendarSet' )
      !if ( .not. (cal == eshr_timemgr_cal) )then
      !   write(6,FA)subname//': Trying to create clock with calendar type = '//&
      !             caltmp
      !   call shr_sys_abort( subname//'ERROR:: Can NOT create a different '//&
      !                       'calendar type' )
      !end if
   end if
   eshr_timemgr_setCalendar = .true. 

END SUBROUTINE eshr_timemgr_initCalendar

!===============================================================================
!===============================================================================

end module eshr_timemgr_mod
