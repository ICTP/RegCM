module mod_clm_time_manager

! Just use RegCM date library here.

   use mod_date
   use mod_realkinds
   use mod_mpmessage
   use mod_runparams , only : idate0 , idate1 , idate2 , idatex , &
                  ktau , mtau , dtsec , dtsrf , ntsrf , doing_restart , &
                  dayspy
   use mod_clm_varpar , only : outfrq

   implicit none
   private

  save
   save

   real(rk8), parameter :: uninit_r8  = -999999999.0
   integer, parameter :: uninit_int = -999999999

   type (rcm_time_and_date) :: cordex_refdate
   integer :: ioffset
   integer :: rst_nstep_rad_prev  ! nstep of previous radiation call

! Public methods

   public ::&
      timemgr_init,             &! time manager initialization
      timemgr_restart,          &! restart the time manager using info from timemgr_restart
      timemgr_datediff,         &! calculate difference between two time instants
      advance_timestep,         &! increment timestep number
      get_step_size,            &! return step size in seconds
      get_nstep,                &! return timestep number
      curr_date,            &! return date components at end of current timestep
      get_prev_date,            &! return date components at beginning of current timestep
      get_start_date,           &! return components of the start date
      get_driver_start_ymd,     &! return year/month/day (as integer in YYYYMMDD format) of driver start date
      get_prev_time,            &! return components of elapsed time since reference date at beg of current timestep
      get_curr_calday,          &! return calendar day at end of current timestep
      get_calday,               &! return calendar day from input date
      is_first_step,            &! return true on first step of initial run
      is_first_restart_step,    &! return true on first step of restart run
      is_end_curr_day,          &! return true on last timestep in current day
      is_end_curr_month,        &! return true on last timestep in current month
      is_last_step,             &! return true on last timestep
      is_restart,               &! return true if this is a restart run
      get_days_per_year,        &
      getdatetime,              &
      get_rad_step_size

! Public data for namelist input

   character(len=32), public ::&
      calendar   = 'NO_LEAP'     ! Calendar to use in date calculations.
                                 ! 'NO_LEAP' or 'GREGORIAN'

! Namelist read in all modes
   integer, public ::&
     dtime          = uninit_int,  &! timestep in seconds
     dtime_rad      = uninit_int,  &! radiation interval in seconds
     nstep_rad_prev = uninit_int    ! radiation interval in seconds

! Namelist read in only in ccsm and offline modes
   integer, public ::&
      nestep        = uninit_int,  &! final timestep (or day if negative) number
      nelapse       = uninit_int,  &! number of timesteps (or days if negative) to extend a run
      start_ymd     = uninit_int,  &! starting date for run in yearmmdd format
      start_tod     = 0,           &! starting time of day for run in seconds
      stop_ymd      = uninit_int,  &! stopping date for run in yearmmdd format
      stop_tod      = 0,           &! stopping time of day for run in seconds
      ref_ymd       = uninit_int,  &! reference date for time coordinate in yearmmdd format
      ref_tod       = 0             ! reference time of day for time coordinate in seconds

! Private module data

   logical :: tm_first_restart_step = .false.  ! true for first step of a restart run
   integer :: cal_type = uninit_int            ! calendar type

! Private module methods

   private :: calc_nestep

!=========================================================================================
contains
!=========================================================================================

subroutine timemgr_init( calendar_in, start_ymd_in, start_tod_in, ref_ymd_in, &
                         ref_tod_in, stop_ymd_in, stop_tod_in)

  !---------------------------------------------------------------------------------
  ! 
  ! Arguments
  character(len=*), optional, intent(IN) :: calendar_in       ! Calendar type
  integer         , optional, intent(IN) :: start_ymd_in      ! Start date (YYYYMMDD)
  integer         , optional, intent(IN) :: start_tod_in      ! Start time of day (sec)
  integer         , optional, intent(IN) :: ref_ymd_in        ! Reference date (YYYYMMDD)
  integer         , optional, intent(IN) :: ref_tod_in        ! Reference time of day (sec)
  integer         , optional, intent(IN) :: stop_ymd_in       ! Stop date (YYYYMMDD)
  integer         , optional, intent(IN) :: stop_tod_in       ! Stop time of day (sec)
  !
  calendar = calstr(idatex%calendar)
  cordex_refdate = 1949120100
  call setcal(cordex_refdate,idatex)
  ioffset = -idnint(dtsec)*(ntsrf-1)

end subroutine timemgr_init

!=========================================================================================

subroutine timemgr_restart( stop_ymd_synclock, stop_tod_synclock )
 
  !---------------------------------------------------------------------------------
  !
  integer, optional, intent(in) :: stop_ymd_synclock
  integer, optional, intent(in) :: stop_tod_synclock

  ! Will do nothing here to be safe
  ! Set flag that this is the first timestep of the restart run.

  tm_first_restart_step = .true.

  ! Calculate ending time step

  calendar = calstr(idatex%calendar)
  cordex_refdate = 1949120100
  call setcal(cordex_refdate,idatex)

  call calc_nestep( )

  ! Print configuration summary to log file (stdout).

end subroutine timemgr_restart

!=========================================================================================

subroutine calc_nestep()
  !---------------------------------------------------------------------------------
  !
  ! Calculate ending timestep number
  ! Calculation of ending timestep number (nestep) assumes a constant stepsize.
  !
  
  nestep = (mtau-ktau)/ntsrf

end subroutine calc_nestep

!=========================================================================================

subroutine init_calendar( )

  !---------------------------------------------------------------------------------
  ! Initialize calendar
  !
  ! Local variables
  !
  calendar = calstr(idatex%calendar)

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

   get_nstep = (ktau/ntsrf)

end function get_nstep

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
    id = idatex
    tdif = (ioffset - outfrq*60)
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

integer function get_driver_start_ymd( tod )

! Return date of start of simulation from driver (i.e. NOT from restart file)
! Note: get_start_date gets you the date from the beginning of the simulation
!       on the restart file.

! Arguments
   integer, intent(out) , optional::&
      tod     ! time of day (seconds past 0Z)

    integer :: ih , yr , mon , day

    call split_idate(idate0,yr,mon,day,ih)

    get_driver_start_ymd = yr*10000+mon*100+day

    if (present(tod)) tod = idate0%second_of_day

end function get_driver_start_ymd

!=========================================================================================

subroutine get_prev_time(days, seconds)

! Return time components valid at beg of current timestep.
! prev time is the time interval between the prev date and the reference date.

! Arguments
   integer, intent(out) ::&
      days,   &! number of whole days in time interval
      seconds  ! remaining seconds in time interval

   type (rcm_time_interval) :: tdif
   real(rk8) :: rh

   tdif = idatex - cordex_refdate
   rh = tohours(tdif)
   if ( doing_restart ) then
     rh = rh - outfrq
   end if
   days = idint(rh/24.0D0)
   seconds = (rh-days*24.0D0)*3600.0D0

end subroutine get_prev_time

!=========================================================================================

real(rk8) function get_curr_calday(offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
   integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
   type (rcm_time_and_date) :: id
   type (rcm_time_interval) :: tdif
   id = idatex
   if ( present(offset) ) then
     tdif = offset
     id = id + tdif
   end if
   get_curr_calday = yeardayfrac(id)

end function get_curr_calday

!=========================================================================================

real(rk8) function get_calday(ymd, tod)

  integer , intent(in) :: ymd , tod

! Return calendar day corresponding to specified time instant.
! Calendar day 1.0 = 0Z on Jan 1.

   type (rcm_time_and_date) :: id

   id = ymd
   id%second_of_day = tod
   call setcal(id,idatex)
   get_calday = yeardayfrac(id)

end function get_calday

!=========================================================================================
 
logical function is_end_curr_day()

! Return true if current timestep is last timestep in current day.

   is_end_curr_day = (idatex%second_of_day == 0)

end function is_end_curr_day

!=========================================================================================

logical function is_end_curr_month()

! Return true if current timestep is last timestep in current month.

   integer :: iy , im , id , ih
   call split_idate(idatex,iy,im,id,ih)
   is_end_curr_month = ( id == 1 .and. idatex%second_of_day == 0 )

end function is_end_curr_month

!=========================================================================================

logical function is_first_step()

! Return true on first step of initial run only.

   is_first_step = ( ktau == 0 )

end function is_first_step

!=========================================================================================

logical function is_first_restart_step()

! Return true on first step of restart run only.

  is_first_restart_step = doing_restart

end function is_first_restart_step

!=========================================================================================

logical function is_last_step()

! Return true on last timestep.

   is_last_step = ( ktau == mtau )

end function is_last_step

!=========================================================================================

subroutine timemgr_datediff(ymd1, tod1, ymd2, tod2, days)

! Calculate the difference (ymd2,tod2) - (ymd1,tod1) and return the result in days.

! Arguments
   integer, intent(in) ::&
      ymd1,    &! date1 in yyyymmdd format
      tod1,    &! time of day relative to date1 (seconds past 0Z)
      ymd2,    &! date2 in yyyymmdd format
      tod2      ! time of day relative to date2 (seconds past 0Z)

   real(rk8) :: days ! (ymd2,tod2)-(ymd1,tod1) in days

   type (rcm_time_and_date) :: id1 , id2
   type (rcm_time_interval) :: tdif

   id1 = ymd1
   id2 = ymd2

   call setcal(id1,idatex)
   call setcal(id2,idatex)

   tdif = id2-id1

   days = tohours(tdif)/24.0D0 + dble(tod2-tod1)/86400.0D0

end subroutine timemgr_datediff

  real(rk8) function get_days_per_year( )
    implicit none
    get_days_per_year = dayspy
  end function get_days_per_year

  integer function get_rad_step_size()
    implicit none
    if (nstep_rad_prev == uninit_int ) then
      get_rad_step_size=get_step_size()
    else
      get_rad_step_size=dtime_rad
    end if
  end function get_rad_step_size

  logical function is_restart( )
    ! Determine if restart run
    use mod_clm_varctl, only : nsrest, nsrContinue
    if (nsrest == nsrContinue) then
      is_restart = .true.
    else
      is_restart = .false.
    end if
  end function is_restart

  subroutine getdatetime(currdate,currtime)
    implicit none
    character(len=8) , intent(out) :: currdate
    character(len=8) , intent(out) :: currtime
    currdate = '0000000'
    currtime = '0000000'
  end subroutine getdatetime

end module mod_clm_time_manager
