!===============================================================================
! SVN $Id: cpl_control_mod.F90 1777 2006-09-01 15:47:11Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_control_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_control_mod -- basic coupler control function logic.
!
! !DESCRIPTION:
!   This module represents a major subsystem of cpl6.  It contains 
!   data type definitions and associated methods used for controlling
!   a coupled integration.  Here ``controlling" refers to issues such as:
!   \begin{itemize}
!   \item selecting integration start and stop dates, 
!   \item determining when history and restart data files should be generated,
!   \item determining when diagnostic data should be generated,
!   \item verifying that all coupled system component models (eg. atm, ocn, etc)
!         are synchronized in time.
!   \item reads and parses namelist variables
!   \item makes various simulation control variables available to other modules
!   \end{itemize}
!
! !REVISION HISTORY:
!    2002-Sep-18 - B. Kauffman -- reworked using shr_alarm_mod
!    2001-May-27 - T. Bettge -- initial prototype
!
! !INTERFACE:  -----------------------------------------------------------------

module cpl_control_mod

! !USES:

   use shr_sys_mod     ! wrappers around system calls
   use shr_cal_mod     ! calendar module
   use shr_date_mod    ! date data-type & methods
   use shr_alarm_mod   ! alarm data-type & methods
   use cpl_kind_mod    ! access to F90 kind declarations

   implicit none

   private ! except

! !PUBLIC TYPES:

   ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_control_readNList  ! read & parse namelist input
   public :: cpl_control_init       ! initialize alarms, etc.
   public :: cpl_control_update     ! update control flags

! !PUBLIC DATA MEMBERS:

   !----- rest, stop, hist, diag control/alarm flags -----
   logical      ,public :: cpl_control_stopNow    ! T => stop model now
   logical      ,public :: cpl_control_stopEOD    ! T => stop model at end of day
   logical      ,public :: cpl_control_restNow    ! T => create restart data now
   logical      ,public :: cpl_control_restEOD    ! T => create restart data at EOD
   logical      ,public :: cpl_control_histNow    ! T => create history data now
   logical      ,public :: cpl_control_histEOD    ! T => create history data at EOD
!  logical      ,public :: cpl_control_histSave   ! T => archive history data now
   logical      ,public :: cpl_control_hist64bit  ! T => use 64 bit netCDFfiles 
   logical      ,public :: cpl_control_avhistNow  ! T => create history data now
   logical      ,public :: cpl_control_avhistEOD  ! T => create history data at EOD
   logical      ,public :: cpl_control_diagNow    ! T => print diagnostic data now
   logical      ,public :: cpl_control_diagEOD    ! T => print diagnostic data at EOD
   logical      ,public :: cpl_control_avDiagNow  ! T => print tavg diag data now
   logical      ,public :: cpl_control_avDiagEOD  ! T => print tavg diag data at EOD
   logical      ,public :: cpl_control_bfbflag    ! T => bfb with different pes

   !----- case name & descriptive string -----
   character(CL),public :: cpl_control_caseName   ! case name
   character(CL),public :: cpl_control_caseDesc   ! case description

   !----- restart control -----
   character(16),public :: cpl_control_restType   ! restart type: init,cont,branch
   integer(IN)  ,public :: cpl_control_restCDate  ! restart cDate from namelist
   integer(IN)  ,public :: cpl_control_restDate   ! restart date
   character(CL),public :: cpl_control_restPFn    ! restart pointer file name
   character(CL),public :: cpl_control_restBFn    ! restart branch  file name
   logical      ,public :: cpl_control_lagOcn     ! T => lag the ocn at startup
   logical      ,public :: cpl_control_sendAtmAlb ! T => send albedo ICs to atm
   logical      ,public :: cpl_control_sendLndDom ! T => send lnd domain to lnd
   logical      ,public :: cpl_control_icData_a   ! T => use IC data provided by atm
   logical      ,public :: cpl_control_icData_i   ! T => use IC data provided by ice
   logical      ,public :: cpl_control_icData_l   ! T => use IC data provided by lnd
   logical      ,public :: cpl_control_icData_o   ! T => use IC data provided by ocn
   logical      ,public :: cpl_control_icData_r   ! T => use IC data provided by roff
   character(16),public :: cpl_control_avhistType ! tavg history file type

   !----- mapping file names -----
   character(CL),public :: cpl_control_mapFn_a2oF ! map data file: a->o fluxes
   character(CL),public :: cpl_control_mapFn_a2oS ! map data file: a->o states
   character(CL),public :: cpl_control_mapFn_o2aF ! map data file: o->a fluxes
   character(CL),public :: cpl_control_mapFn_o2aS ! map data file: o->a states
   character(CL),public :: cpl_control_mapFn_r2o  ! map data file: r->o runoff

   !----- flux & orbital options -----
   logical      ,public :: cpl_control_fluxAlbAv  ! T => NO diurnal cycle in albedos
   character(16),public :: cpl_control_fluxEPbal  ! selects E,P,R adjustment technique
   real(R8)     ,public :: cpl_control_fluxEPfac  ! E,P,R adjust factor recv'd from ocn
   integer(IN)  ,public :: cpl_control_fluxAShift ! albedo calc time-shift (seconds)
   real(R8)     ,public :: cpl_control_orbEccen   ! eccen of earth orbit (unitless)
   real(R8)     ,public :: cpl_control_orbObliqr  ! earth's obliquity (rad)
   real(R8)     ,public :: cpl_control_orbLambm0  ! mean lon perihelion @ vernal eq (rad)
   real(R8)     ,public :: cpl_control_orbMvelpp  ! moving vernal equinox longitude
                                                  ! of perihelion plus pi (rad)

   !----- info about which specific component models are in use -----
   logical      ,public :: cpl_control_dead_a     ! T => atm component is dead comp
   logical      ,public :: cpl_control_dead_i     ! T => ice component is dead comp
   logical      ,public :: cpl_control_dead_l     ! T => lnd component is dead comp
   logical      ,public :: cpl_control_dead_o     ! T => ocn component is dead comp
   logical      ,public :: cpl_control_dead_ao    ! T => atm and/or ocn are dead comp

   !----- date/time & timestep info -----
   integer(IN)  ,public :: cpl_control_nCpl_a     ! atm/cpl communications per day
   integer(IN)  ,public :: cpl_control_nCpl_i     ! ice/cpl communications per day
   integer(IN)  ,public :: cpl_control_nCpl_l     ! lnd/cpl communications per day
   integer(IN)  ,public :: cpl_control_nCpl_o     ! ocn/cpl communications per day
   integer(IN)  ,public :: cpl_control_nCpl_r     ! rof/cpl communications per day
   integer(IN)  ,public :: cpl_control_cDate_a    ! atm coded date
   integer(IN)  ,public :: cpl_control_cDate_i    ! ice coded date
   integer(IN)  ,public :: cpl_control_cDate_l    ! lnd coded date
   integer(IN)  ,public :: cpl_control_cDate_o    ! ocn coded date
   integer(IN)  ,public :: cpl_control_sec_a      ! atm secs on date
   integer(IN)  ,public :: cpl_control_sec_i      ! ice secs on date
   integer(IN)  ,public :: cpl_control_sec_l      ! lnd secs on date
   integer(IN)  ,public :: cpl_control_sec_o      ! ocn secs on date

   !----- grid checking -----
   real(R8)     ,public :: cpl_control_eps_almask ! epsilon for masks, atm/lnd
   real(R8)     ,public :: cpl_control_eps_algrid ! epsilon for grid coords, atm/lnd
   real(R8)     ,public :: cpl_control_eps_alarea ! epsilon for areas, atm/lnd
   real(R8)     ,public :: cpl_control_eps_oimask ! epsilon for masks, ocn/ice
   real(R8)     ,public :: cpl_control_eps_oigrid ! epsilon for grid coords, ocn/ice
   real(R8)     ,public :: cpl_control_eps_oiarea ! epsilon for areas, ocn/ice

   !----- decomposition settings -----
   integer(IN)  ,public :: cpl_control_decomp_a   ! atm decomposition type
   integer(IN)  ,public :: cpl_control_decomp_l   ! lnd decomposition type
   integer(IN)  ,public :: cpl_control_decomp_o   ! ocn decomposition type
   integer(IN)  ,public :: cpl_control_decomp_i   ! ice decomposition type
   integer(IN)  ,public :: cpl_control_decomp_r   ! rof decomposition type

   !----- other flags -----
   integer(IN)  ,public :: cpl_control_infoDBug=1 ! user specified dbug level
   logical      ,public :: cpl_control_infoBcheck ! T => do bit-check now

!EOP

   !----- local -----
   type(shr_alarm) :: stopAlarm    ! goes on when model should stop
   type(shr_alarm) :: restAlarm    ! goes on when a restart file should be made
   type(shr_alarm) :: histAlarm    ! goes on when a history file should be made
   type(shr_alarm) :: avhistAlarm  ! goes on when tavg history file is written
   type(shr_alarm) :: diagAlarm    ! goes on when diagnostics should be done
   type(shr_alarm) :: avDiagAlarm  ! goes on when tavg diags should be done

   type(shr_date)  :: startDate    ! keeps a record of the model start date

   integer(IN),parameter :: unset = -999   ! flags unspecified namelist var

   !----------------------------------------------------------------------------
   ! Following are namelist variables needed by this control module.   They are
   ! declared/saved here because some of them must be read-in and used by
   ! two different routines.
   !----------------------------------------------------------------------------
   character(CL)  :: case_name     ! case name for history files, etc.
   character(CL)  :: case_desc     ! short text description of case
   character(16)  :: start_type    ! initial, continue, branch
   integer(IN)    :: start_date    ! required for some start types
   character(CL)  :: start_pfile   ! restart pointer file
   character(CL)  :: start_bfile   ! restart branch  file
   character(16)  :: rest_option   ! frequency
   integer(IN)    :: rest_n        ! required for some options
   integer(IN)    :: rest_date     ! required for some options
   character(16)  :: stop_option   ! frequency
   integer(IN)    :: stop_n        ! required for some options
   integer(IN)    :: stop_date     ! required for some options
   character(16)  :: hist_option   ! frequency
   integer(IN)    :: hist_n        ! required for some options
   integer(IN)    :: hist_date     ! required for some options
   logical        :: hist_64bit    ! make 64 bit netCDF files
   character(16)  :: avhist_option ! frequency
   integer(IN)    :: avhist_n      ! required for some options
   integer(IN)    :: avhist_date   ! required for some options
   character(16)  :: diag_option   ! frequency
   integer(IN)    :: diag_n        ! required for some options
   integer(IN)    :: diag_date     ! required for some options
   character(16)  :: avDiag_option ! frequency
   integer(IN)    :: avDiag_n      ! required for some options
   integer(IN)    :: avDiag_date   ! required for some options
   character(CL)  :: map_a2of_fn   ! map data file: a->o fluxes
   character(CL)  :: map_a2os_fn   ! map data file: a->o states
   character(CL)  :: map_o2af_fn   ! map data file: o->a fluxes
!  character(CL)  :: map_o2as_fn   ! map data file: o->a states ! not an option
   character(CL)  :: map_r2o_fn    ! map data file: r->o runoff flux
   integer(IN)    :: orb_year      ! year (AD) wrt orbital parameters
   real(R8)       :: orb_eccen     ! eccentricity of earths orbit (unitless)
   real(R8)       :: orb_mvelp     ! moving vernel equinox of orbit (degrees)
   real(R8)       :: orb_obliq     ! obliquity of orbit (degrees)
   logical        :: flx_albav     ! T => NO diurnal cycle in albedos
   character(16)  :: flx_epbal     ! flag useage of E,P,R adjustment factor
   integer(IN)    :: info_dbug     ! dbug level setting (0,1,2,...)
   integer(IN)    :: info_bcheck   ! bitCheck level setting (0,1,2,...)
   real(R8)       :: eps_almask    ! epsilon for masks, atm/lnd
   real(R8)       :: eps_algrid    ! epsilon for grid coords, atm/lnd
   real(R8)       :: eps_alarea    ! epsilon for areas, atm/lnd
   real(R8)       :: eps_oimask    ! epsilon for masks, ocn/ice
   real(R8)       :: eps_oigrid    ! epsilon for grid coords, ocn/ice
   real(R8)       :: eps_oiarea    ! epsilon for areas, ocn/ice
   integer(IN)    :: decomp_al     ! decomposition for atm/lnd
   integer(IN)    :: decomp_oi     ! decomposition for ocn/ice
   integer(IN)    :: decomp_r      ! decomposition for rof
   logical        :: bfbflag       ! bfb flag for different pes

   namelist / inparm /  &
   &     case_name   , case_desc   ,               & 
   &     start_type  , start_date  ,               & 
   &     start_pfile , start_bfile ,               &
   &     rest_option , rest_n      , rest_date   , &
   &     stop_option , stop_n      , stop_date   , &
   &     hist_option , hist_n      , hist_date   , hist_64bit , &
   &   avhist_option , avhist_n    , avhist_date , &
   &     diag_option , diag_n      , diag_date   , &
   &   avDiag_option , avDiag_n    , avDiag_date , &
   &     map_a2of_fn , map_a2os_fn , map_o2af_fn , map_r2o_fn , &
   &     orb_year    , orb_eccen   , orb_mvelp   , orb_obliq  , &
   &     flx_albav   , flx_epbal   , info_dbug   , info_bcheck, &
   &     eps_almask  , eps_algrid  , eps_alarea  , &
   &     eps_oimask  , eps_oigrid  , eps_oiarea  , &
   &     decomp_al   , decomp_oi   , decomp_r    , bfbflag

   character(*),parameter :: modName = 'cpl_control_mod'

   save

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_control_readNList - initialize and read namelist values.
!
! !DESCRIPTION:
!    Initialize and read namelist values.
! 
! !REVISION HISTORY:
!    2002-Sep-18 - B. Kauffman -- 1st version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_control_readNList()

! !USES:

   use shr_orb_mod

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   ! modifies private module data for later use (the namelist variables)

!EOP

   !----- local -----
   integer(IN) :: rcode       ! return code
   real(R8):: orb_mvelpp  ! moving vernal equinox lon of perihelion + pi (rad)
   real(R8):: orb_lambm0  ! Mean longitude of perihelion at the
   real(R8):: orb_obliqr  ! Earth's obliquity (rad)
   logical :: orb_log     ! flag to print-out orbital info

   !----- formats -----
   character(*),parameter :: subName = "('cpl_control_readNList') "
   character(*),parameter :: F00 = "('(cpl_control_readNList) ',8a)"
   character(*),parameter :: F01 = "('(cpl_control_readNList) ',a,i6)"
   character(*),parameter :: F12 = "('(cpl_control_readNList) ',a,f11.6)"
   character(*),parameter :: F90 = "('(cpl_control_readNList) ',60('-'))"

!-------------------------------------------------------------------------------
! NOTES:
!   Here we merely set default values and then read-in name-list values.  
!   They will be accessed and used elsewhere.
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! set default values for namelist vars
   !----------------------------------------------------------------------------
   case_name   = "unset"   ! case name (used in output file names)
   case_desc   = "unset"   ! case description text string
   start_type  = "initial" ! initial, continue, or branch
   start_date  = 010101    ! required for some types
   start_pfile = "unset"   ! restart pointer file
   start_bfile = "unset"   ! restart branch  file
   rest_option = "monthly" ! frequency
   rest_n      = 3         ! required for some options
   rest_date   = unset     ! required for some options
   stop_option = "monthly" ! frequency
   stop_n      = 3         ! required for some options
   stop_date   = unset     ! required for some options
   hist_option = "monthly" ! frequency
   hist_n      = 3         ! required for some options
   hist_date   = unset     ! required for some options
   hist_64bit  = .false.   ! false => 32 bit netCDF files
   avhist_option = "none"  ! frequency
   avhist_n      = 1       ! required for some options
   avhist_date   = unset   ! required for some options
   diag_option = "monthly" ! frequency
   diag_n      = 3         ! required for some options
   diag_date   = unset     ! required for some options
   avDiag_option = "yearly"! frequency
   avDiag_n      = 3       ! required for some options
   avDiag_date   = unset   ! required for some options
   map_a2of_fn = "unknown" ! map data file: a->o fluxes
   map_a2os_fn = "unknown" ! map data file: a->o states
   map_o2af_fn = "unknown" ! map data file: o->a fluxes
!  map_o2as_fn = "unknown" ! map data file: o->a states ! not an option
   map_r2o_fn  = "unknown" ! map data file: r->o runoff flux
   orb_year    = unset     ! orbital year wrt solar angle 
   orb_eccen   = unset     
   orb_mvelp   = unset
   orb_obliq   = unset
   flx_albav   = .false.   ! compute "time average" albedos
   flx_epbal   = 'off'     ! activate artificial evap/precip/runoff balance
   info_dbug   = 1         ! level of debug info to stdout
   info_bcheck = 0         ! level of debug info to stdout
   eps_almask  = 1.0e-6_R8 ! epsilon for masks
   eps_algrid  = 1.0e-2_R8 ! epsilon for grid coords
   eps_alarea  = 1.0e-1_R8 ! epsilon for areas
   eps_oimask  = 1.0e-6_R8 ! epsilon for masks
   eps_oigrid  = 1.0e-2_R8 ! epsilon for grid coords
   eps_oiarea  = 1.0e-1_R8 ! epsilon for areas
   decomp_al   = 1
   decomp_oi   = 1
   decomp_r    = 1
   bfbflag     = .false.

   !----------------------------------------------------------------------------
   ! document namelist values
   !----------------------------------------------------------------------------
   write(6,F90)
   write(6,F00) "Namelist values BEFORE reading file..."
   write(6,inparm)
   write(6,F90)

   !----------------------------------------------------------------------------
   ! read input namelist
   !----------------------------------------------------------------------------
   open(15,file="cpl.nml",status="old",action="read")
   read(15,nml=inparm,iostat=rcode)
   close(unit=15)
   if (rcode > 0) then
      write(6,F01) 'ERROR: reading input namelist, iostat=',rcode
      call shr_sys_abort(subName//": namelist read error")
   end if

   !----------------------------------------------------------------------------
   ! document namelist values
   !----------------------------------------------------------------------------
   write(6,F90)
   write(6,F00) "Namelist values AFTER reading file..."
   write(6,inparm)
   write(6,F90)

   !----------------------------------------------------------------------------
   ! must confirm an acceptable start-type
   !----------------------------------------------------------------------------
   if (  start_type /= "initial"  &
   .and. start_type /= "continue" &
   .and. start_type /= "branch"   ) then
      write(6,F00) 'ERROR: invalid start_type =',trim(start_type)
      call shr_sys_abort(subName//": namelist value error")
   end if
   if ( .not. shr_cal_validDate(start_date) ) then
      write(6,F00) 'ERROR: invalid start_date =',start_date
      call shr_sys_abort(subName//": namelist value error")
   end if

   !----------------------------------------------------------------------------
   ! orbital year and/or parameter selection
   !----------------------------------------------------------------------------
   if ( orb_year /= unset ) then
      write(6,F01) 'orbit based on orb_year = ',orb_year
      orb_eccen = SHR_ORB_UNDEF_REAL
      orb_obliq = SHR_ORB_UNDEF_REAL
      orb_mvelp = SHR_ORB_UNDEF_REAL
   else if ( orb_eccen /= unset &
   &   .and. orb_obliq /= unset &
   &   .and. orb_mvelp /= unset ) then
      write(6,F12) 'orbit based on orb_eccen = ',orb_eccen
      write(6,F12) 'orbit based on orb_obliq = ',orb_obliq
      write(6,F12) 'orbit based on orb_mvelp = ',orb_mvelp
      orb_year  = SHR_ORB_UNDEF_INT
   else
      write(6,F00) 'orbit ERROR: must set input parm orb_year or'
      write(6,F00) 'all three of orb_eccen ,orb_obliq, & orb_mvelp.'
      call shr_sys_abort(subName//": namelist value error")
   endif
   orb_log = .true.
   call shr_orb_params( orb_year  , orb_eccen , orb_obliq , orb_mvelp, &
   &                    orb_obliqr, orb_lambm0, orb_mvelpp, orb_log  )

   !----------------------------------------------------------------------------
   ! set default restart pointer file name, if necessary
   !----------------------------------------------------------------------------
   if ( trim(start_pfile) == 'unset') then
      if ( trim(case_name) == 'unset') then
         start_pfile = 'rpointer'
      else
         start_pfile = '$HOME/cpl.'//trim(case_name)//'.rpointer'
      endif
   else if ( trim(start_pfile) == './rpointer ') then
      start_pfile = 'rpointer'
   endif

   !----------------------------------------------------------------------------
   ! make local/private control variables available via public variables
   !----------------------------------------------------------------------------
   cpl_control_caseName   = case_name
   cpl_control_caseDesc   = case_desc
   cpl_control_restType   = start_type
   cpl_control_restCDate  = start_date
   cpl_control_restPFn    = start_pfile
   cpl_control_restBFn    = start_bfile
   cpl_control_avhistTYPE = avhist_option
   cpl_control_hist64bit  = hist_64bit
   cpl_control_orbEccen   = orb_eccen
   cpl_control_orbObliqr  = orb_obliqr
   cpl_control_orbLambm0  = orb_lambm0
   cpl_control_orbMvelpp  = orb_mvelpp
   cpl_control_fluxAlbav  = flx_albav 
   cpl_control_fluxEpbal  = flx_epbal
   cpl_control_mapFn_a2oF = map_a2of_fn 
   cpl_control_mapFn_a2oS = map_a2os_fn 
   cpl_control_mapFn_o2aF = map_o2af_fn 
   cpl_control_mapFn_o2aS = map_o2af_fn ! note: same as o2aF
   cpl_control_mapFn_r2o  = map_r2o_fn 
   cpl_control_fluxEpfac  = 0
   cpl_control_infoDBug   = info_dbug
   if (start_type == "initial") then
      cpl_control_lagOcn = .true.
   else
      cpl_control_lagOcn = .false.
   endif
   cpl_control_eps_almask = eps_almask
   cpl_control_eps_algrid = eps_algrid
   cpl_control_eps_alarea = eps_alarea
   cpl_control_eps_oimask = eps_oimask
   cpl_control_eps_oigrid = eps_oigrid
   cpl_control_eps_oiarea = eps_oiarea
   cpl_control_decomp_a   = decomp_al
   cpl_control_decomp_l   = decomp_al
   cpl_control_decomp_o   = decomp_oi
   cpl_control_decomp_i   = decomp_oi
   cpl_control_decomp_r   = decomp_r
   cpl_control_bfbflag    = bfbflag

!  ---------------------------------------------------
!  check for valid values of decomp, 
!  currently only 1 is valid for production
!  1002 is a special test value, convert it to 2 for testing
!  ---------------------------------------------------
   if (cpl_control_decomp_a /= 1 .and. cpl_control_decomp_a /= 1002) then
     write(6,*) trim(subName),' ERROR cpl decomp_a error ',cpl_control_decomp_a
     write(6,*) trim(subName),' ERROR decomp_al must be 1 in cpl namelist'
     call shr_sys_abort()
   endif
   if (cpl_control_decomp_l /= 1 .and. cpl_control_decomp_l /= 1002) then
     write(6,*) trim(subName),' ERROR cpl decomp_l error ',cpl_control_decomp_l
     write(6,*) trim(subName),' ERROR decomp_al must be 1 in cpl namelist'
     call shr_sys_abort()
   endif
   if (cpl_control_decomp_o /= 1 .and. cpl_control_decomp_o /= 1002) then
     write(6,*) trim(subName),' ERROR cpl decomp_o error ',cpl_control_decomp_o
     write(6,*) trim(subName),' ERROR decomp_oi must be 1 in cpl namelist'
     call shr_sys_abort()
   endif
   if (cpl_control_decomp_i /= 1 .and. cpl_control_decomp_i /= 1002) then
     write(6,*) trim(subName),' ERROR cpl decomp_i error ',cpl_control_decomp_i
     write(6,*) trim(subName),' ERROR decomp_oi must be 1 in cpl namelist'
     call shr_sys_abort()
   endif
   if (cpl_control_decomp_r /= 1 .and. cpl_control_decomp_r /= 1002) then
     write(6,*) trim(subName),' ERROR cpl decomp_r error ',cpl_control_decomp_r
     write(6,*) trim(subName),' ERROR decomp_r must be 1 in cpl namelist'
     call shr_sys_abort()
   endif
   if (cpl_control_decomp_a == 1002) cpl_control_decomp_a = 2
   if (cpl_control_decomp_l == 1002) cpl_control_decomp_l = 2
   if (cpl_control_decomp_o == 1002) cpl_control_decomp_o = 2
   if (cpl_control_decomp_i == 1002) cpl_control_decomp_i = 2
   if (cpl_control_decomp_r == 1002) cpl_control_decomp_r = 2

end subroutine cpl_control_readNList

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_control_init - initializes flags for stopping, restart, etc.
!
! !DESCRIPTION:
!    Set the module variable {\tt startDate} to the input argument {\tt date}.
!    Also set the module variables {\tt stopAlarm}, {\tt restAlarm},
!    {\tt histAlarm}, {\tt avhistAlarm}, {\tt diagAlarm} and {\tt avdiagAlarm}.
!    If any of {\tt rest\_date, stop\_date, etc.} are unset or negative,
!    set them to the input {\tt date}.
! 
! !REVISION HISTORY:
!    2002-Sep-18 - B. Kauffman -- 1st version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_control_init(date)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

  !output: modifies data declared in this module's header
  type(shr_date),intent(in) :: date ! model start date 

!EOP

   !----- local -----
   integer(IN) :: cDate     ! coded date (yymmdd)
   integer(IN) :: sec       ! seconds
   integer(IN) :: rcode     ! return code
   integer(IN) :: ns        ! number of steps per model day

   character(*),parameter :: subName = "('cpl_control_init') "

   !----- formats -----
   character(*),parameter :: F00 = "('(cpl_control_init) ',8a)"

!-------------------------------------------------------------------------------
! ASSUMPTIONS: 
!   alarms turn on on specified days, but alarms are not aware of the time of 
!   day. Thus these control flags will always turn on at the beginning of a day.
!   This functionality can be generalized when it is clear what alternative is
!   more desirable.
!-------------------------------------------------------------------------------

   !--- clobber model start date with date from namelist ?? ---
   startDate = date ! keep a record of the model start date
   call shr_date_getCDate(date,cDate,sec)

   !--- any unset or negative dates get set to the model start date ---
   if ((  rest_date == unset) .or. (   rest_date < 0))   rest_date = cDate
   if ((  stop_date == unset) .or. (   stop_date < 0))   stop_date = cDate
   if ((  hist_date == unset) .or. (   hist_date < 0))   hist_date = cDate
   if ((avhist_date == unset) .or. ( avhist_date < 0)) avhist_date = cDate
   if ((  diag_date == unset) .or. (   diag_date < 0))   diag_date = cDate
   if ((avDiag_date == unset) .or. ( avDiag_date < 0)) avDiag_date = cDate
   
   !--- initialize the alarms ---
   ns = shr_date_getStepsPerDay(date)
   stopAlarm   = alarmInit(  stop_option,  stop_n,shr_date_initCDate(  stop_date,ns))
   restAlarm   = alarmInit(  rest_option,  rest_n,shr_date_initCDate(  rest_date,ns))
   histAlarm   = alarmInit(  hist_option,  hist_n,shr_date_initCDate(  hist_date,ns))
   avhistAlarm = alarmInit(avhist_option,avhist_n,shr_date_initCDate(avhist_date,ns))
   diagAlarm   = alarmInit(  diag_option,  diag_n,shr_date_initCDate(  diag_date,ns))
   avDiagAlarm = alarmInit(avDiag_option,avDiag_n,shr_date_initCDate(avDiag_date,ns))

   !--- debugging info ---
   if (cpl_control_infoDBug > 2) then
      write(6,F00) "stopAlarm details..."
      call shr_alarm_dump(stopAlarm)
      write(6,F00) "restAlarm details..."
      call shr_alarm_dump(restAlarm)
      write(6,F00) "histAlarm details..."
      call shr_alarm_dump(histAlarm)
      write(6,F00) "avhistAlarm details..."
      call shr_alarm_dump(avhistAlarm)
      write(6,F00) "diagAlarm details..."
      call shr_alarm_dump(diagAlarm)
      write(6,F00) "avDiagAlarm details..."
      call shr_alarm_dump(avDiagAlarm)
   end if

end subroutine cpl_control_init

!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE: alarmInit - parse alarm namelist variables and initialize alarms.
!
! !DESCRIPTION:
!     Local/private routine to parse alarm namelist vars and initialize alarms.
!
! !REVISION HISTORY:
!    2002-Sep-18 - B. Kauffman -- 1st version
!
! !INTERFACE: ------------------------------------------------------------------

type(shr_alarm) function alarmInit(option,n,date)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

  !modifies public cpl_control_mod data
   character(*)   ,intent(inout):: option  ! "date","ndays","monthly", etc.
   integer(IN)    ,intent(in)   :: n       ! n wrt ndays or nmonths
   type(shr_date) ,intent(in)   :: date    ! alarm-on or offset date

! !EOP

   !----- formats ---
   character(*),parameter :: subName = '('//modName//':alarmInit) '
   character(*),parameter :: F00 = "('(cpl_control:alarmInit) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if      ( option(1:4) == "date"   ) then
      alarmInit = shr_alarm_initDate(date)
   else if ( option(1:4) == "nday" ) then
      alarmInit = shr_alarm_initNDays(n,date)
   else if ( option(1:6) == "nmonth" ) then
      alarmInit = shr_alarm_initNMonths(n,date)
   else if ( option(1:5) == "daily"  ) then
      alarmInit = shr_alarm_initNDays(1,date)
   else if ( option(1:5) == "month") then
      alarmInit = shr_alarm_initMonthly()
   else if ( option(1:4) == "year" ) then
      alarmInit = shr_alarm_initYearly()
   else if ( option(1:5) == "nstep" ) then
      alarmInit = shr_alarm_initNStep(n)
   else if ( option(1:5) == "ifsec" ) then
      alarmInit = shr_alarm_initifsec(n)
   else if ( option(1:7) == "ifdays0" ) then
      alarmInit = shr_alarm_initifdays0(n)
   else if ( option(1:5) == "ifday" ) then
      alarmInit = shr_alarm_initifday(n)
   else if ( option(1:5) == "ifmon" ) then
      alarmInit = shr_alarm_initifmon(n)
   else if ( option(1:5) == "ifyear") then
      alarmInit = shr_alarm_initifyear(n)
   else if ( option(1:5) == "none") then
      alarmInit = shr_alarm_initNone()
   else if ( option(1:5) == "never") then
      option = "none"
      alarmInit = shr_alarm_initNone()
   else 
      write(6,F00) "ERROR: unknown option = ",option
      write(6,F00) "Recognized option strings are: " &
      & // "none, never, date, ndays, nmonths, daily, monthly, yearly" &
      & // ", nstep, ifsec, ifdays0, ifday, ifmon, ifyear"
      write(6,F00) "String recognition is case sensitive."
      call shr_sys_abort(subName)
   end if

end function alarmInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_control_update - sets control flags for stopping,restart,etc.
!
! !DESCRIPTION:
!   Update all the module {\tt cpl\_control\_*Now} and {\tt cpl\_control\_*EOD} flags
!   using input argument {\tt currentDate}.
!
!   Also set {\tt cpl\_control\_infoBCheck}.
! 
! !REVISION HISTORY:
!    2002-Sep-18 - B. Kauffman -- 1st version using shr_alarm_mod
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_control_update(currentDate)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(in) :: currentDate  ! current model date

   !OUTPUT: modifies public cpl_control_mod data

!EOP

   !----- local -----
   type(shr_date) :: tomorrowsDate  ! day following current model date 
   logical        :: old_stopNow    ! value before update
   logical        :: old_stopEOD    ! value before update
   logical        :: old_restNow    ! value before update
   logical        :: old_histNow    ! value before update
   logical        :: old_avHistNow  ! value before update
   logical        :: old_diagNow    ! value before update
   logical        :: old_avDiagNow  ! value before update
   integer        :: cDate,sec      ! coded date (yymmdd) and seconds
   integer        :: y,m,d,s        ! year, month, day, seconds

   !----- formate -----
   character(*),parameter :: subName = "(cpl_control_update)"
   character(*),parameter :: F00 = "('(cpl_control_update) ',8a)"

!-------------------------------------------------------------------------------
! ASSUMPTIONS: 
!   alarms turn on on specified days, but alarms are not aware of the time of 
!   day. Thus these control flags will always turn on at the beginning of a day.
!   This functionality can be generalized when it is clear what alternative is
!   more desirable.
!-------------------------------------------------------------------------------

   !--- save old values so we can test if they've changed (below) ---
   old_stopNow   = cpl_control_stopNow 
   old_stopEOD   = cpl_control_stopEOD 
   old_restNow   = cpl_control_restNow 
   old_histNow   = cpl_control_histNow 
   old_avHistNow = cpl_control_avhistNow 
   old_diagNow   = cpl_control_diagNow 
   old_avDiagNow = cpl_control_avDiagNow 

   !--- determine tomorrow's date ---
   tomorrowsDate = currentDate
   call shr_date_advNextDay(tomorrowsDate)

   !----------------------------------------------------------------------------
   ! stop flags
   !----------------------------------------------------------------------------
   cpl_control_stopNow = shr_alarm_isOn(  currentDate,stopAlarm)
   cpl_control_stopEOD = shr_alarm_isOn(tomorrowsDate,stopAlarm)

   !--- handle special cases on startup ---
   if (currentDate == startDate) then

      !--- we don't want periodic alarms to stop model on start date ---
      if (shr_alarm_getType(stopAlarm) /= shr_alarm_date) &
      &   cpl_control_stopNow = .false.

      !--- if type is stop-on-date, then stop date must be after start date ---
      if (shr_alarm_getType(stopAlarm) == shr_alarm_date ) then
         call shr_date_getCDate(startDate,cDate,sec)
         if (stop_date < cDate) then
            write(6,*) "ERROR: user specified stop date prior to start date"
            call shr_sys_abort(subName//": stop date is prior to start date")
         end if
         if (stop_date == cDate) then
            write(6,*) "ERROR: user specified stop date equal to start date"
            call shr_sys_abort(subName//": stop date is equal to start date")
         end if
      end if

      !--- model must run at least two days so ocn can get stopEOD warning ---
      !--- note: this is only necessary if ocn is lagged by one whole day ---
      !--- note: could put more elaborate logic in here if requested ---
      if ( cpl_control_stopNow .or. cpl_control_stopEOD )  then
         write(6,*) "ERROR: model must run at least two days"
         call shr_sys_abort(subName//": model must run at least two days")
      end if

   end if

   !----------------------------------------------------------------------------
   ! restart flags
   !----------------------------------------------------------------------------
   cpl_control_restNow = shr_alarm_isOn(  currentDate,restAlarm)
   cpl_control_restEOD = shr_alarm_isOn(tomorrowsDate,restAlarm)
   if (cpl_control_restNow) call shr_alarm_set(restAlarm,.false.)

   !----------------------------------------------------------------------------
   ! history flag ---
   !----------------------------------------------------------------------------
   cpl_control_histNow = shr_alarm_isOn(  currentDate,histAlarm)
   cpl_control_histEOD = shr_alarm_isOn(tomorrowsDate,histAlarm)
   if (cpl_control_histNow) call shr_alarm_set(histAlarm,.false.)

   cpl_control_avhistNow = shr_alarm_isOn(  currentDate,avhistAlarm)
   cpl_control_avhistEOD = shr_alarm_isOn(tomorrowsDate,avhistAlarm)
   if (cpl_control_avhistNow) call shr_alarm_set(avhistAlarm,.false.)

   !----------------------------------------------------------------------------
   ! diagnostic flags
   !----------------------------------------------------------------------------
   cpl_control_diagNow = shr_alarm_isOn(  currentDate,diagAlarm)
   cpl_control_diagEOD = shr_alarm_isOn(tomorrowsDate,diagAlarm)
   if (cpl_control_diagNow) call shr_alarm_set(diagAlarm,.false.)

   cpl_control_avDiagNow = shr_alarm_isOn(  currentDate,avDiagAlarm)
   cpl_control_avDiagEOD = shr_alarm_isOn(tomorrowsDate,avDiagAlarm)
   if (cpl_control_avDiagNow) call shr_alarm_set(avDiagAlarm,.false.)

   !----------------------------------------------------------------------------
   ! bit-check flag
   !----------------------------------------------------------------------------
   call shr_date_getYMD(currentDate,y,m,d,s)
   cpl_control_infoBCheck = .false.
   if (     info_bCheck == 1) then
      if (d == 1 .and. s == 0) cpl_control_infoBCheck = .true.
   else if (info_bCheck == 2) then
      if (s == 0) cpl_control_infoBCheck = .true.
   else if (info_bCheck >= 3) then
                  cpl_control_infoBCheck = .true.
   end if

   !----------------------------------------------------------------------------
   ! log status of control flags
   !----------------------------------------------------------------------------
   if ( cpl_control_infoDBug < 3) then
      if (.not.old_stopNow   .and. cpl_control_stopNow   )  &
                   write(6,F00) "cpl_control_stopNow   = .true."
      if (.not.old_stopEOD   .and. cpl_control_stopEOD   )  &
                   write(6,F00) "cpl_control_stopEOD   = .true."
      if (.not.old_histNow   .and. cpl_control_histNow   )  &
                   write(6,F00) "cpl_control_histNow   = .true."
      if (.not.old_avHistNow   .and. cpl_control_avHistNow)  &
                   write(6,F00) "cpl_control_avHistNow = .true."
      if (.not.old_diagNow   .and. cpl_control_diagNow   )  &
                   write(6,F00) "cpl_control_diagNow   = .true."
      if (.not.old_avDiagNow   .and. cpl_control_avDiagNow)  &
                   write(6,F00) "cpl_control_avDiagNow = .true."
   else
      if (cpl_control_stopNow)   write(6,F00) "cpl_control_stopNow   = .true."
      if (cpl_control_stopEOD)   write(6,F00) "cpl_control_stopEOD   = .true."
      if (cpl_control_restNow)   write(6,F00) "cpl_control_restNow   = .true."
      if (cpl_control_restEOD)   write(6,F00) "cpl_control_restEOD   = .true."
      if (cpl_control_histNow)   write(6,F00) "cpl_control_histNow   = .true."
      if (cpl_control_histEOD)   write(6,F00) "cpl_control_histEOD   = .true."
      if (cpl_control_avhistNow) write(6,F00) "cpl_control_avhistNow = .true."
      if (cpl_control_avhistEOD) write(6,F00) "cpl_control_avhistEOD = .true."
      if (cpl_control_diagNow)   write(6,F00) "cpl_control_diagNow   = .true."
      if (cpl_control_diagEOD)   write(6,F00) "cpl_control_diagEOD   = .true."
      if (cpl_control_avDiagNow) write(6,F00) "cpl_control_avDiagNow = .true."
      if (cpl_control_avDiagEOD) write(6,F00) "cpl_control_avDiagEOD = .true."
   end if

end subroutine cpl_control_update

!===============================================================================
!===============================================================================
end module cpl_control_mod
