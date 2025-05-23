module mod_clm_time_manager

   use mod_date
   use mod_intkinds
   use mod_realkinds
   use mod_mpmessage
   use mod_constants, only : secpd
   use mod_dynparam, only : dayspy
   use mod_runparams, only : idate0
   use mod_clm_varctl, only : nextdate

   implicit none

   private

   save

   ! return components of the start date
   public :: get_start_date
   ! return year/month/day (as integer(ik4) in YYYYMMDD format) of driver start
   public :: get_driver_start_ymd
   ! return calendar day at end of current timestep
   public :: get_curr_calday
   ! return position in the earth orbit in days
   public :: get_curr_yearpoint
   ! return calendar day from input date
   public :: get_calday
   ! return true on last timestep in current day
   public :: is_end_curr_day
   ! return true on last timestep in current month
   public :: is_end_curr_month
   ! return true if this is a restart run
   public :: is_restart
   ! return true if year end
   public :: is_end_curr_year
   ! return true if year middle
   public :: is_middle_curr_year

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
    integer(ik4), intent(out), optional::&
      tod     ! time of day (seconds past 0Z)
    integer(ik4) :: ih, yr, mon, day

    call split_idate(idate0,yr,mon,day,ih)

    get_driver_start_ymd = yr*10000+mon*100+day

    if ( present(tod) ) tod = idate0%second_of_day
  end function get_driver_start_ymd

  real(rk8) function get_curr_calday( )
    implicit none
    ! Return calendar day at end of current timestep with optional offset.
    ! Calendar day 1.0 = 0Z on Jan 1.
    ! Offset from current time in seconds.
    ! Positive for future times, negative for previous times.
    get_curr_calday = yeardayfrac(nextdate)
  end function get_curr_calday

  real(rk8) function get_curr_yearpoint(offset)
   implicit none
   ! Return calendar day at end of current timestep with optional offset.
   ! Calendar day 1.0 = 0Z on Jan 1.
    ! Offset from current time in seconds.
    ! Positive for future times, negative for previous times.
    integer(ik4), optional, intent(in) :: offset
    type (rcm_time_and_date) :: id
    type (rcm_time_interval) :: tdif

    id = nextdate
    if ( present(offset) ) then
      tdif = offset
      id = id + tdif
    end if
    get_curr_yearpoint = yearpoint(id)
  end function get_curr_yearpoint

  real(rk8) function get_calday(ymd, tod)
    implicit none
    integer(ik4), intent(in) :: ymd, tod
    ! Return calendar day corresponding to specified time instant.
    ! Calendar day 1.0 = 0Z on Jan 1.
    type (rcm_time_and_date) :: id

    id = i4wcal(ymd,nextdate%calendar)
    id%second_of_day = tod
    get_calday = yeardayfrac(id)
  end function get_calday

  logical function is_end_curr_day()
    implicit none
    ! Return true if current timestep is last timestep in current day.
    is_end_curr_day = (nextdate%second_of_day == int(0, ik8))
  end function is_end_curr_day

  logical function is_end_curr_month()
    implicit none
    ! Return true if current timestep is last timestep in current month.
    integer(ik4) :: iy, im, id, ih, imm, iss
    call split_idate(nextdate,iy,im,id,ih,imm,iss)
    is_end_curr_month = ( id == 1 .and. ih == 0 .and. &
                          imm == 0 .and. iss == 0 )
  end function is_end_curr_month

  logical function is_end_curr_year( )
    implicit none
    is_end_curr_year = date_is(nextdate,1,1) .and. time_is(nextdate,0)
  end function is_end_curr_year

  logical function is_middle_curr_year( )
    implicit none
    is_middle_curr_year = date_is(nextdate,7,1) .and. time_is(nextdate,0)
  end function is_middle_curr_year

  logical function is_restart( )
    use mod_clm_varctl, only : nsrest, nsrContinue, nsrStartup, DoForceRestart
    implicit none
    if (nsrest == nsrContinue .or. &
        (nsrest == nsrStartup .and. DoForceRestart) ) then
      is_restart = .true.
    else
      is_restart = .false.
    end if
  end function is_restart

  subroutine getdatetime(currdate,currtime)
    implicit none
    character(len=8), intent(out) :: currdate
    character(len=8), intent(out) :: currtime
    integer(ik4), dimension(8) :: tval
    call date_and_time(values=tval)
    write(currdate,'(i0.4,i0.2,i0.2)') tval(1), tval(2), tval(3)
    write(currtime,'(i0.2,i0.2,i0.2,i0.2)') tval(5), tval(6), tval(7), 0
  end subroutine getdatetime

end module mod_clm_time_manager
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
