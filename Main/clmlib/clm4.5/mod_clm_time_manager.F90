module mod_clm_time_manager

   use mod_date
   use mod_intkinds
   use mod_realkinds
   use mod_mpmessage
   use mod_constants , only : secpd
   use mod_runparams , only : idate0 , idatex , dtsec , dtsrf , dayspy

   implicit none

   private

   save

   ! return components of the start date
   public :: get_start_date
   ! return year/month/day (as integer(ik4) in YYYYMMDD format) of driver start
   public :: get_driver_start_ymd
   ! return calendar day at end of current timestep
   public :: get_curr_calday
   ! return calendar day from input date
   public :: get_calday
   ! return true on last timestep in current day
   public :: is_end_curr_day
   ! return true on last timestep in current month
   public :: is_end_curr_month
   ! return true if this is a restart run
   public :: is_restart

   public :: getdatetime

  contains

  subroutine get_start_date(yr, mon, day, tod)
    implicit none
    ! Return date components valid at beginning of initial run.
    integer(ik4), intent(out) ::&
       yr,    &! year
       mon,   &! month
       day,   &! day of month
       tod     ! time of day (seconds past 0Z)
    integer(ik4) :: ih

    call split_idate(idate0,yr,mon,day,ih)
    tod = idate0%second_of_day
  end subroutine get_start_date

  integer(ik4) function get_driver_start_ymd( tod )
    implicit none
    ! Return date of start of simulation from driver
    ! (i.e. NOT from restart file)
    ! Note: get_start_date gets you the date from the beginning of
    ! the simulation on the restart file.
    integer(ik4), intent(out) , optional::&
      tod     ! time of day (seconds past 0Z)
    integer(ik4) :: ih , yr , mon , day

    call split_idate(idate0,yr,mon,day,ih)

    get_driver_start_ymd = yr*10000+mon*100+day

    if ( present(tod) ) tod = idate0%second_of_day
  end function get_driver_start_ymd

  real(rk8) function get_curr_calday(offset)
    implicit none
    ! Return calendar day at end of current timestep with optional offset.
    ! Calendar day 1.0 = 0Z on Jan 1.
    ! Offset from current time in seconds.
    ! Positive for future times, negative for previous times.
    integer(ik4), optional, intent(in) :: offset
    type (rcm_time_and_date) :: id
    type (rcm_time_interval) :: tdif

    id = idatex
    if ( present(offset) ) then
      tdif = offset
      id = id + tdif
    end if
    get_curr_calday = yeardayfrac(id)
  end function get_curr_calday

  real(rk8) function get_calday(ymd, tod)
    implicit none
    integer(ik4) , intent(in) :: ymd , tod
    ! Return calendar day corresponding to specified time instant.
    ! Calendar day 1.0 = 0Z on Jan 1.
    type (rcm_time_and_date) :: id

    id = ymd
    id%second_of_day = tod
    call setcal(id,idatex)
    get_calday = yeardayfrac(id)
  end function get_calday

  logical function is_end_curr_day()
    implicit none
    ! Return true if current timestep is last timestep in current day.
    is_end_curr_day = (idatex%second_of_day == int(secpd-dtsec, ik8))
  end function is_end_curr_day

  logical function is_end_curr_month()
    implicit none
    ! Return true if current timestep is last timestep in current month.
    type (rcm_time_and_date) :: id
    id = monlast(idatex)
    is_end_curr_month = ( &
            id%days_from_reference == idatex%days_from_reference .and. &
            is_end_curr_day() )
  end function is_end_curr_month

  logical function is_restart( )
    use mod_clm_varctl, only : nsrest, nsrContinue
    implicit none
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
