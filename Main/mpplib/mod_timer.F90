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

  public :: init_timer
  public :: step_timer
  public :: reached_endtime
  public :: nowstring
  public :: ktau
  public :: time_from_start
  public :: rcm_alarm

  ! Supported precision is 1 second (minimum dt is 1)

  integer(ik8) :: model_initial_time
  integer(ik8) :: model_start_time
  integer(ik8) :: model_stop_time
  integer(ik8) :: model_timestep
  integer(ik8) :: model_internal_time

  logical :: reached_endtime

  type(rcm_time_and_date) :: idate
  type(rcm_time_interval) :: intmdl

  type rcm_alarm
    real(rkx) :: dt
    integer(ik8) :: actint
    integer(ik8) :: lastact
    real(rkx) , dimension(2) :: wt
    integer(ik8) :: now
    logical :: triggered
  contains
    procedure :: act => alarm_act
  end type rcm_alarm

  interface rcm_alarm
    module procedure init_alarm
  end interface rcm_alarm

  contains

  subroutine init_timer(mdate0,mdate1,mdate2,mdt)
    implicit none
    type(rcm_time_and_date) , intent(in) :: mdate0 , mdate1 , mdate2
    real(rkx) , intent(in) :: mdt
    type(rcm_time_interval) :: tdif
    model_initial_time = 0
    tdif = mdate1 - mdate0
    model_start_time = int(tohours(tdif)*secph,ik8)
    tdif = mdate2 - mdate0
    model_stop_time = int(tohours(tdif)*secph,ik8)
    model_timestep = int(mdt,ik8)
    model_internal_time = model_start_time
    reached_endtime = model_internal_time >= model_stop_time
    intmdl = rcm_time_interval(mdt,usec)
    idate = mdate1
  end subroutine init_timer

  subroutine step_timer
    model_internal_time = model_internal_time + model_timestep
    reached_endtime = model_internal_time >= model_stop_time
    idate = idate + intmdl
  end subroutine step_timer

  character (len=32) function nowstring result(ns)
    implicit none
    ns = tochar(idate)
  end function nowstring

  integer(ik8) function ktau
    implicit none
    ktau = model_internal_time/model_timestep
  end function ktau

  integer(ik8) function time_from_start
    implicit none
    time_from_start = model_internal_time
  end function time_from_start

  function init_alarm(dt,act0)
    implicit none
    type(rcm_alarm) :: init_alarm
    real(rkx) , intent(in) :: dt
    logical , intent(in) , optional :: act0
    logical :: lact0
    lact0 = .false.
    if ( present(act0) ) lact0 = act0
    init_alarm%dt = dt
    init_alarm%actint = int(dt,ik8)
    init_alarm%lastact = model_internal_time
    init_alarm%triggered = lact0
  end function init_alarm

  logical function alarm_act(alarm) result(res)
    implicit none
    real(rkx) :: t1 , t2
    class(rcm_alarm) , intent(inout) :: alarm
    res = .false.
    t1 = real(model_internal_time,rkx)
    t2 = real(alarm%lastact+alarm%actint,rkx)
    if ( t1 >= t2 ) then
      alarm%lastact = (model_internal_time/alarm%actint) * alarm%actint
      alarm%triggered = .true.
    end if
    if ( alarm%triggered ) then
      res = .true.
      alarm%triggered = .false.
      alarm%wt(1) = (t1 - t2)/alarm%dt
      alarm%wt(2) = d_one - alarm%wt(1)
      alarm%now = model_internal_time
    end if
  end function alarm_act

end module mod_timer

#ifdef TESTME

subroutine myabort
  implicit none
  call abort
end subroutine myabort

program test_timing
  use mod_realkinds
  use mod_date
  use mod_timer
  implicit none

  type(rcm_time_and_date) :: mdate0 , mdate1 , mdate2
  type(rcm_alarm) :: srf_alarm , rad_alarm , cum_alarm

  mdate0 = 1950010100
  mdate1 = 1950010100
  mdate2 = 2300010100

  call init_timer(mdate0,mdate1,mdate2,213.0_rkx)

  print *, nowstring( ) , ktau( )

  srf_alarm = rcm_alarm(600.0_rkx,.true.)
  rad_alarm = rcm_alarm(1800.0_rkx,.true.)
  cum_alarm = rcm_alarm(300.0_rkx)

  do while ( .not. reached_endtime )
    call step_timer
    if ( srf_alarm%act( ) ) then
      print *, 'SRF ', srf_alarm%now , srf_alarm%wt(1)
    end if
    if ( rad_alarm%act( ) ) then
      print *, 'RAD ', rad_alarm%now , rad_alarm%wt(1)
    end if
    if ( cum_alarm%act( ) ) then
      print *, 'CUM ', cum_alarm%now , cum_alarm%wt(1)
    end if
  end do

  print *, nowstring( ) , ktau( )

end program test_timing
#endif

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
