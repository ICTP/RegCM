#include <misc.h>
#include <preproc.h>

module clm_time_manager

   use shr_kind_mod, only: r8 => shr_kind_r8
   use spmdMod    , only: masterproc, iam
   use abortutils , only: endrun
   use ESMF_Mod

! Just use RegCM date library here.

   use mod_date
   use mod_intkinds, only : ik8
   use mod_mpmessage
   use mod_runparams, only : idate0, idate1, idate2, rcmtimer, &
                  syncro_srf, dtsec, dtsrf, doing_restart

   use mod_constants, only : secpd
   use clm_varsur  , only : r2coutfrq

   implicit none
   private
   save

   type (rcm_time_and_date) :: cordex_refdate
   integer :: ioffset

! Public methods

   public ::&
      timemgr_init,             &! time manager initialization
      timemgr_restart_io,       &! read/write time manager restart info and restart time manager
      timemgr_restart,          &! restart the time manager using info from timemgr_restart
      timemgr_datediff,         &! calculate difference between two time instants
      advance_timestep,         &! increment timestep number
      get_step_size,            &! return step size in seconds
      get_nstep,                &! return timestep number
      get_curr_date,            &! return date components at end of current timestep
      get_prev_date,            &! return date components at beginning of current timestep
      get_start_date,           &! return components of the start date
      get_ref_date,             &! return components of the reference date
      get_perp_date,            &! return components of the perpetual date, and current time of day
      get_curr_time,            &! return components of elapsed time since reference date at end of current timestep
      get_prev_time,            &! return components of elapsed time since reference date at beg of current timestep
      get_curr_calday,          &! return calendar day at end of current timestep
      get_calday,               &! return calendar day from input date
      is_first_step,            &! return true on first step of initial run
      is_first_restart_step,    &! return true on first step of restart or branch run
      is_end_curr_day,          &! return true on last timestep in current day
      is_end_curr_month,        &! return true on last timestep in current month
      is_last_step,             &! return true on last timestep
      is_perpetual               ! return true if perpetual calendar is in use

! Public data for namelist input

   character(len=ESMF_MAXSTR), public ::&
      calendar   = 'NO_LEAP'     ! Calendar to use in date calculations.
                                 ! 'NO_LEAP' or 'GREGORIAN'
   integer, parameter :: uninit_int = -999999999

! Namelist read in all modes
   integer, public ::&
      dtime         = uninit_int    ! timestep in seconds

! Namelist read in only in ccsm and offline modes
   integer, public ::&
      nelapse       = uninit_int,  &! number of timesteps (or days if negative) to extend a run
      start_ymd     = uninit_int,  &! starting date for run in yearmmdd format
      start_tod     = 0,           &! starting time of day for run in seconds
      stop_ymd      = uninit_int,  &! stopping date for run in yearmmdd format
      stop_tod      = 0,           &! stopping time of day for run in seconds
      ref_ymd       = uninit_int,  &! reference date for time coordinate in yearmmdd format
      ref_tod       = 0             ! reference time of day for time coordinate in seconds

! Private module data

   logical :: tm_first_restart_step = .false.  ! true for first step of a restart or branch run
   logical :: tm_perp_calendar = .false.       ! true when using perpetual calendar
   integer :: cal_type = uninit_int            ! calendar type

! Private module methods

!=========================================================================================
contains
!=========================================================================================

subroutine timemgr_init( calendar_in, start_ymd_in, start_tod_in, ref_ymd_in, &
                         ref_tod_in, stop_ymd_in, stop_tod_in,                &
                         perpetual_run_in, perpetual_ymd_in )

  !---------------------------------------------------------------------------------
  ! Initialize the ESMF time manager from the sync clock
  !
  ! Arguments
  character(len=*), optional, intent(IN) :: calendar_in       ! Calendar type
  integer        , optional, intent(IN) :: start_ymd_in      ! Start date (YYYYMMDD)
  integer        , optional, intent(IN) :: start_tod_in      ! Start time of day (sec)
  integer        , optional, intent(IN) :: ref_ymd_in        ! Reference date (YYYYMMDD)
  integer        , optional, intent(IN) :: ref_tod_in        ! Reference time of day (sec)
  integer        , optional, intent(IN) :: stop_ymd_in       ! Stop date (YYYYMMDD)
  integer        , optional, intent(IN) :: stop_tod_in       ! Stop time of day (sec)
  logical        , optional, intent(IN) :: perpetual_run_in  ! If in perpetual mode or not
  integer        , optional, intent(IN) :: perpetual_ymd_in  ! Perpetual date (YYYYMMDD)
  !
  calendar = calstr(rcmtimer%idate%calendar)
  cordex_refdate = 1949120100
  call setcal(cordex_refdate,rcmtimer%idate)
  ioffset = -(dtsrf-dtsec)

end subroutine timemgr_init

!=========================================================================================

subroutine timemgr_restart_io( ncid, flag )

  !---------------------------------------------------------------------------------
  ! Read/Write information needed on restart to a netcdf file.

  ! Arguments
  integer        , intent(in) :: ncid  ! netcdf id
  character(len=*), intent(in) :: flag  ! 'read' or 'write'
  !

  calendar = calstr(rcmtimer%idate%calendar)
  cordex_refdate = 1949120100
  call setcal(cordex_refdate,rcmtimer%idate)
  ioffset = -(dtsrf-dtsec)

end subroutine timemgr_restart_io

!=========================================================================================

subroutine timemgr_restart( stop_ymd_synclock, stop_tod_synclock )

  !---------------------------------------------------------------------------------
  ! Restart the ESMF time manager using the synclock for ending date.
  !
  integer, optional, intent(in) :: stop_ymd_synclock
  integer, optional, intent(in) :: stop_tod_synclock

  ! Will do nothing here to be safe
  ! Set flag that this is the first timestep of the restart run.

  tm_first_restart_step = .true.

  ! Calculate ending time step

  calendar = calstr(rcmtimer%idate%calendar)
  cordex_refdate = 1949120100
  call setcal(cordex_refdate,rcmtimer%idate)

  ! Print configuration summary to log file (stdout).

end subroutine timemgr_restart

!=========================================================================================

!=========================================================================================

subroutine init_calendar( )

  !---------------------------------------------------------------------------------
  ! Initialize calendar
  !
  ! Local variables
  !
  calendar = calstr(rcmtimer%idate%calendar)

end subroutine init_calendar

!=========================================================================================

subroutine advance_timestep()

  ! Increment the timestep number.

  character(len=*), parameter :: sub = 'advance_timestep'
  integer :: rc

  tm_first_restart_step = .false.

end subroutine advance_timestep

!=========================================================================================

integer function get_step_size()

  ! Return the step size in seconds.

  get_step_size = idnint(dtsrf)

end function get_step_size

!=========================================================================================

integer function get_nstep()

  ! Return the timestep number.

   get_nstep = syncro_srf%lcount + 1

end function get_nstep

!=========================================================================================

subroutine get_curr_date(yr, mon, day, tod, offset)

  !-----------------------------------------------------------------------------------------
  ! Return date components valid at end of current timestep with an optional
  ! offset (positive or negative) in seconds.

  integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative
                                            ! for previous times.
    type (rcm_time_and_date) :: id
    type (rcm_time_interval) :: tdif
    integer :: ih
    id = rcmtimer%idate
    ih = ioffset
    if ( present(offset) ) then
      ih = ih + offset
    end if
    tdif = ih
    id = id + tdif
    call split_idate(id,yr,mon,day,ih)
    tod = id%second_of_day

end subroutine get_curr_date

!=========================================================================================

subroutine get_perp_date(yr, mon, day, tod, offset)

  integer, intent(in) :: yr, mon, day, tod, offset

  call fatal(__FILE__,__LINE__,'NOT IMPLEMENTED get_perp_date')

end subroutine get_perp_date

!=========================================================================================

subroutine get_prev_date(yr, mon, day, tod)

! Return date components valid at beginning of current timestep.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

    type (rcm_time_and_date) :: id
    type (rcm_time_interval) :: tdif
    integer :: ih
    id = rcmtimer%idate
    tdif = (ioffset - r2coutfrq*60)
    id = id + tdif
    call split_idate(id,yr,mon,day,ih)
    tod = id%second_of_day

end subroutine get_prev_date

!=========================================================================================

subroutine get_start_date(yr, mon, day, tod)

! Return date components valid at beginning of initial run.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

    integer :: ih

    call split_idate(idate0,yr,mon,day,ih)
    tod = idate0%second_of_day

end subroutine get_start_date

!=========================================================================================

subroutine get_ref_date(yr, mon, day, tod)

! Return date components of the reference date.

! Arguments
   integer, intent(out) ::&
      yr,    &! year
      mon,   &! month
      day,   &! day of month
      tod     ! time of day (seconds past 0Z)

    integer :: ih

    yr = 1949
    mon = 12
    day = 1
    tod = 0

end subroutine get_ref_date

!=========================================================================================

subroutine get_curr_time(days, seconds)

! Return time components valid at end of current timestep.
! Current time is the time interval between the current date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

    type (rcm_time_interval) :: tdif
    real(r8) :: rh

    tdif = rcmtimer%idate - cordex_refdate
    rh = tohours(tdif) + dble(ioffset)/3600.0D0
    days = idint(rh/24.0D0)
    seconds = (rh-days*24.0D0)*3600.0D0

end subroutine get_curr_time

!=========================================================================================

subroutine get_prev_time(days, seconds)

! Return time components valid at beg of current timestep.
! prev time is the time interval between the prev date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

   type (rcm_time_interval) :: tdif
   real(r8) :: rh

   tdif = rcmtimer%idate - cordex_refdate
   rh = tohours(tdif)
   if ( doing_restart ) then
     rh = rh - r2coutfrq
   end if
   days = idint(rh/24.0D0)
   seconds = (rh-days*24.0D0)*3600.0D0

end subroutine get_prev_time

!=========================================================================================

real(r8) function get_curr_calday(offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative
                                            ! for previous times.
   type (rcm_time_and_date) :: id
   type (rcm_time_interval) :: tdif
   id = rcmtimer%idate
   if ( present(offset) ) then
     tdif = offset
     id = id + tdif
   end if
   get_curr_calday = yeardayfrac(id)

end function get_curr_calday

!=========================================================================================

real(r8) function get_calday(ymd, tod)

  integer, intent(in) :: ymd, tod

! Return calendar day corresponding to specified time instant.
! Calendar day 1.0 = 0Z on Jan 1.

   type (rcm_time_and_date) :: id

   id = ymd
   id%second_of_day = tod
   call setcal(id,rcmtimer%idate)
   get_calday = yeardayfrac(id)

end function get_calday

!=========================================================================================

logical function is_end_curr_day()

! Return true if current timestep is last timestep in current day.

   is_end_curr_day = (rcmtimer%idate%second_of_day == int(secpd-dtsec, ik8))

end function is_end_curr_day

!=========================================================================================

logical function is_end_curr_month()

! Return true if current timestep is last timestep in current month.

   type (rcm_time_and_date) :: id
   id = monlast(rcmtimer%idate)
   is_end_curr_month = ( &
      id%days_from_reference == rcmtimer%idate%days_from_reference .and. &
      is_end_curr_day() )
end function is_end_curr_month

!=========================================================================================

logical function is_first_step()

! Return true on first step of initial run only.

   is_first_step = rcmtimer%start( )

end function is_first_step

!=========================================================================================

logical function is_first_restart_step()

! Return true on first step of restart run only.

  is_first_restart_step = doing_restart

end function is_first_restart_step

!=========================================================================================

logical function is_last_step()

! Return true on last timestep.

   is_last_step = rcmtimer%reached_endtime

end function is_last_step

!=========================================================================================

logical function is_perpetual()

   is_perpetual = .false.

end function is_perpetual

!=========================================================================================

subroutine timemgr_datediff(ymd1, tod1, ymd2, tod2, days)

! Calculate the difference (ymd2,tod2) - (ymd1,tod1) and return the result in days.

! Arguments
   integer, intent(in) ::&
      ymd1,    &! date1 in yyyymmdd format
      tod1,    &! time of day relative to date1 (seconds past 0Z)
      ymd2,    &! date2 in yyyymmdd format
      tod2      ! time of day relative to date2 (seconds past 0Z)

   real(r8) :: days ! (ymd2,tod2)-(ymd1,tod1) in days

   type (rcm_time_and_date) :: id1, id2
   type (rcm_time_interval) :: tdif

   id1 = ymd1
   id2 = ymd2

   call setcal(id1,rcmtimer%idate)
   call setcal(id2,rcmtimer%idate)

   tdif = id2-id1

   days = tohours(tdif)/24.0D0 + dble(tod2-tod1)/86400.0D0

end subroutine timemgr_datediff

end module clm_time_manager
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
