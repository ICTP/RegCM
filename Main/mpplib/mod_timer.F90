!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_timer

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_date

  implicit none

  private

  public :: rcm_timer , rcm_alarm

  type rcm_timer
    integer(ik8) :: model_initial_time
    integer(ik8) :: model_start_time
    integer(ik8) :: model_stop_time
    integer(ik8) :: model_timestep
    integer(ik8) :: model_internal_time
    logical :: reached_endtime
    type(rcm_time_and_date) :: idate
    type(rcm_time_interval) :: intmdl
    integer(ik4) :: nowinday
    integer(ik4) :: year , month , day , hour , minute , second
  contains
    procedure :: step => step_timer
    procedure :: str => nowstring
    procedure :: from_start => time_from_start
    procedure :: ktau => step_from_start
  end type rcm_timer

  interface rcm_timer
    module procedure init_timer
  end interface rcm_timer

  ! Supported precision is 1 second (minimum dt is 1)

  type rcm_alarm
    real(rkx) :: dt
    integer(ik8) :: actint
    integer(ik8) :: lastact
    real(rkx) , dimension(2) :: wt
    integer(ik8) :: now
    logical :: triggered
    type(rcm_timer) , pointer :: timer
  contains
    procedure :: act => alarm_act
  end type rcm_alarm

  interface rcm_alarm
    module procedure init_alarm
  end interface rcm_alarm

  contains

  function init_timer(mdate0,mdate1,mdate2,mdt)
    implicit none
    type(rcm_time_and_date) , intent(in) :: mdate0 , mdate1 , mdate2
    type(rcm_timer) , pointer :: init_timer
    real(rkx) , intent(in) :: mdt
    type(rcm_time_interval) :: tdif
    allocate(init_timer)
    init_timer%model_initial_time = 0
    tdif = mdate1 - mdate0
    init_timer%model_start_time = int(tohours(tdif)*secph,ik8)
    tdif = mdate2 - mdate0
    init_timer%model_stop_time = int(tohours(tdif)*secph,ik8)
    init_timer%model_timestep = int(mdt,ik8)
    init_timer%model_internal_time = init_timer%model_start_time
    init_timer%reached_endtime = &
       init_timer%model_internal_time >= init_timer%model_stop_time
    init_timer%intmdl = rcm_time_interval(mdt,usec)
    init_timer%idate = mdate1
    call split_idate(init_timer%idate,init_timer%year, &
                     init_timer%month,init_timer%day,  &
                     init_timer%hour,init_timer%minute,init_timer%second)
    init_timer%nowinday = init_timer%idate%second_of_day
  end function init_timer

  subroutine step_timer(t)
    implicit none
    class(rcm_timer) , intent(inout) :: t
    integer(ik4) :: tmp1 , tmp2
    t%model_internal_time = t%model_internal_time + t%model_timestep
    t%nowinday = t%nowinday + t%model_timestep
    t%idate = t%idate + t%intmdl
    if ( t%nowinday > 86400 ) then
      call split_idate(t%idate,t%year,t%month,t%day,t%hour,t%minute,t%second)
      t%nowinday = t%idate%second_of_day
    else
      t%hour = t%nowinday/3600
      t%minute = mod(t%nowinday,3600)/60
      t%second = mod(t%nowinday,60)
    end if
    t%reached_endtime = t%model_internal_time >= t%model_stop_time
  end subroutine step_timer

  character (len=32) function nowstring(t) result(ns)
    implicit none
    class(rcm_timer) , intent(in) :: t
    ns = tochar(t%idate)
  end function nowstring

  integer(ik8) function step_from_start(t)
    implicit none
    class(rcm_timer) , intent(in) :: t
    step_from_start = t%model_internal_time/t%model_timestep
  end function step_from_start

  integer(ik8) function time_from_start(t)
    implicit none
    class(rcm_timer) , intent(in) :: t
    time_from_start = t%model_internal_time
  end function time_from_start

  function init_alarm(t,dt,act0)
    implicit none
    class(rcm_timer) , pointer , intent(in) :: t
    type(rcm_alarm) :: init_alarm
    real(rkx) , intent(in) :: dt
    logical , intent(in) , optional :: act0
    logical :: lact0
    init_alarm%timer => null( )
    if ( .not. associated(t) ) return
    init_alarm%timer => t
    lact0 = .false.
    if ( present(act0) ) lact0 = act0
    init_alarm%dt = dt
    init_alarm%actint = int(dt,ik8)
    init_alarm%lastact = t%model_internal_time
    init_alarm%triggered = lact0
  end function init_alarm

  logical function alarm_act(alarm) result(res)
    implicit none
    real(rkx) :: t1 , t2
    class(rcm_alarm) , intent(inout) :: alarm
    res = .false.
    t1 = real(alarm%timer%model_internal_time,rkx)
    t2 = real(alarm%lastact+alarm%actint,rkx)
    if ( t1 >= t2 ) then
      alarm%lastact = (alarm%timer%model_internal_time/alarm%actint) * &
                       alarm%actint
      alarm%triggered = .true.
    end if
    if ( alarm%triggered ) then
      res = .true.
      alarm%triggered = .false.
      alarm%wt(1) = (t1 - t2)/alarm%dt
      alarm%wt(2) = d_one - alarm%wt(1)
      alarm%now = alarm%timer%model_internal_time
    end if
  end function alarm_act

end module mod_timer

#ifdef TESTME

subroutine myabort
  implicit none
  call abort
end subroutine myabort

program test_timing
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use mod_timer
  implicit none

  integer(ik8) , dimension(3) :: idates = [ 1950010100_ik8, &
                                            1950010100_ik8, &
                                            1951010100_ik8 ]

  type(rcm_time_and_date) :: mdate0 , mdate1 , mdate2

  type(rcm_timer) , pointer :: timer
  type(rcm_alarm) :: srf_alarm , rad_alarm , cum_alarm
  type(rcm_alarm) :: srf_output

  mdate0 = idates(1)
  mdate1 = idates(2)
  mdate2 = idates(3)

  timer => rcm_timer(mdate0,mdate1,mdate2,213.0_rkx)

  print *, timer%str( ) , timer%ktau( )

  srf_alarm = rcm_alarm(timer,600.0_rkx,.true.)
  srf_output = rcm_alarm(timer,3600.0_rkx*3.0,.true.)
  rad_alarm = rcm_alarm(timer,1800.0_rkx,.true.)
  cum_alarm = rcm_alarm(timer,300.0_rkx)

  do while ( .not. timer%reached_endtime )
    call timer%step( )
    print *, timer%year,timer%month,timer%day, &
             timer%hour,timer%minute,timer%second
    if ( srf_alarm%act( ) ) then
      print *, 'SRF ', srf_alarm%now , srf_alarm%wt(1)
    end if
    if ( rad_alarm%act( ) ) then
      print *, 'RAD ', rad_alarm%now , rad_alarm%wt(1)
    end if
    if ( cum_alarm%act( ) ) then
      print *, 'CUM ', cum_alarm%now , cum_alarm%wt(1)
    end if
    if ( srf_output%act( ) ) then
      print *, 'OUT_SRF' , srf_output%now , srf_output%wt(1)
    end if
  end do

  print *, timer%str( ) , timer%ktau( )

  deallocate(timer)

end program test_timing
#endif

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
