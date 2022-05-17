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

  ! Supported precision is 1 second (minimum dt is 1)

  public :: rcm_timer , rcm_alarm , rcm_syncro
  public :: operator(/)

  integer(ik4) , parameter :: maxalarms = 64
  integer(ik4) , parameter :: maxsyncro = 8

  type alarmp
    type(rcm_alarm) , pointer :: ap
  end type alarmp

  type syncrop
    type (rcm_syncro) , pointer :: sp
  end type syncrop

  type rcm_timer
    integer(ik8) :: model_initial_time
    integer(ik8) :: model_start_time
    integer(ik8) :: model_stop_time
    integer(ik8) :: model_timestep
    integer(ik8) :: model_internal_time
    integer(ik8) :: lcount
    logical :: reached_endtime
    logical :: next_is_endtime
    type(rcm_time_and_date) :: idate
    type(rcm_time_interval) :: intmdl
    integer(ik8) :: nowinday
    integer(ik4) :: year , month , day , hour , minute , second
    character(len=32) :: model_timestring
    integer(ik4) :: nalarm = 0
    integer(ik4) :: nsyncro = 0
    type(alarmp) , dimension(maxalarms) :: ap
    type(syncrop) , dimension(maxsyncro) :: sp
  contains
    procedure :: advance => step_timer
    procedure :: str => nowstring
    procedure :: from_start => time_from_start
    procedure :: step => step_from_start
    procedure :: start => is_start
    procedure :: integrating => is_integrating
    procedure :: dismiss => cleanup
  end type rcm_timer

  interface rcm_timer
    module procedure init_timer
  end interface rcm_timer

  type rcm_syncro
    integer(ik8) :: frq
    integer(ik8) :: lcount
    real(rkx) :: rw
    class(rcm_timer) , pointer :: timer
  contains
    procedure :: check => syncro_check
    procedure :: act => syncro_act
    procedure :: will_act => syncro_willact
  end type rcm_syncro

  interface rcm_syncro
    module procedure init_syncro
  end interface rcm_syncro

  interface operator(/)
    module procedure ratio_freq_syncro
    module procedure ratio_freq_alarm
    module procedure ratio_freq_syncro_alarm
    module procedure ratio_freq_alarm_syncro
  end interface operator(/)

  type rcm_alarm
    real(rkx) :: dt
    integer(ik8) :: now
    integer(ik8) :: actint
    integer(ik8) :: lastact
    integer(ik8) :: lcount
    real(rkx) :: rw
    real(rkx) , dimension(2) :: wt
    logical :: triggered
    class(rcm_timer) , pointer :: timer
    type(rcm_time_and_date) :: idate
    type(rcm_time_interval) :: intalm
  contains
    procedure :: check => alarm_check
    procedure :: act => alarm_act
    procedure :: will_act => alarm_willact
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
    type(rcm_timer) , pointer :: t
    allocate(t)
    t%model_initial_time = 0_ik8
    tdif = mdate1 - mdate0
    t%model_start_time = int(tohours(tdif)*secph,ik8)
    tdif = mdate2 - mdate0
    t%model_stop_time = int(tohours(tdif)*secph,ik8)
    t%model_timestep = int(mdt,ik8)
    t%model_internal_time = t%model_start_time
    t%reached_endtime = t%model_internal_time >= t%model_stop_time
    t%next_is_endtime = (t%model_internal_time + &
                         t%model_timestep) >= t%model_stop_time
    t%intmdl = rcm_time_interval(int(mdt,ik8),usec)
    t%lcount = t%model_internal_time/t%model_timestep
    t%idate = mdate1
    call split_idate(t%idate,t%year,t%month,t%day,t%hour,t%minute,t%second)
    t%nowinday = t%idate%second_of_day
    t%model_timestring = tochar(mdate1)
    init_timer => t
  end function init_timer

  subroutine step_timer(t)
    implicit none
    class(rcm_timer) , intent(inout) :: t
    integer(ik4) :: i
    t%model_internal_time = t%model_internal_time + t%model_timestep
    t%nowinday = t%nowinday + int(t%model_timestep,ik4)
    t%idate = t%idate + t%intmdl
    t%lcount = t%lcount + 1
    if ( t%nowinday >= 86400 ) then
      call split_idate(t%idate,t%year,t%month,t%day,t%hour,t%minute,t%second)
      t%nowinday = t%idate%second_of_day
    else
      t%hour = t%nowinday/3600
      t%minute = mod(t%nowinday,3600)/60
      t%second = mod(t%nowinday,60)
    end if
    t%reached_endtime = t%model_internal_time >= t%model_stop_time
    t%next_is_endtime = (t%model_internal_time + &
                         t%model_timestep) >= t%model_stop_time
    do i = 1 , t%nalarm
      call t%ap(i)%ap%check( )
    end do
    do i = 1 , t%nsyncro
      call t%sp(i)%sp%check( )
    end do
  end subroutine step_timer

  subroutine cleanup(t)
    implicit none
    class(rcm_timer) , intent(inout) :: t
    integer(ik4) :: i
    do i = 1 , t%nalarm
      deallocate(t%ap(i)%ap)
      nullify(t%ap(i)%ap)
    end do
    do i = 1 , t%nsyncro
      deallocate(t%sp(i)%sp)
      nullify(t%sp(i)%sp)
    end do
    t%nalarm = 0
    t%nsyncro = 0
  end subroutine cleanup

  logical function is_start(t)
    implicit none
    class(rcm_timer) , intent(in) :: t
    is_start = t%lcount == 0
  end function is_start

  logical function is_integrating(t)
    implicit none
    class(rcm_timer) , intent(in) :: t
    is_integrating = t%lcount > 0
  end function is_integrating

  character (len=32) function nowstring(t) result(ns)
    implicit none
    class(rcm_timer) , intent(inout) :: t
    integer(ik8) , save :: last
    if ( t%model_internal_time /= last ) then
      t%model_timestring = tochar(t%idate)
      last = t%model_internal_time
    end if
    ns = t%model_timestring
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

  function init_syncro(t,dt)
    implicit none
    type(rcm_syncro) , pointer :: init_syncro
    type(rcm_timer) , pointer , intent(inout) :: t
    real(rkx) , intent(in) :: dt
    type(rcm_syncro) , pointer :: syncro
    allocate(syncro)
    syncro%timer => null( )
    if ( associated(t) ) then
      syncro%timer => t
      syncro%frq = int(dt,ik8)
      syncro%rw = real(syncro%timer%model_timestep,rkx)/dt
      syncro%timer%nsyncro = syncro%timer%nsyncro + 1
      syncro%timer%sp(syncro%timer%nsyncro)%sp => syncro
      syncro%lcount = syncro%timer%model_internal_time/syncro%frq
    end if
    init_syncro => syncro
  end function init_syncro

  subroutine syncro_check(s)
    implicit none
    class(rcm_syncro) , intent(inout) :: s
    if ( mod(s%timer%model_internal_time,s%frq) == 0 ) then
      s%lcount = s%lcount + 1
    end if
  end subroutine syncro_check

  logical function syncro_act(s) result(res)
    implicit none
    class(rcm_syncro) , intent(in) :: s
    res = (mod(s%timer%model_internal_time,s%frq) == 0)
  end function syncro_act

  logical function syncro_willact(s,dt) result(res)
    implicit none
    class(rcm_syncro) , intent(in) :: s
    real(rkx) , optional , intent(in) :: dt
    integer(ik8) :: idt
    if ( present(dt) ) then
      idt = int(dt,ik8)
      res = (mod(s%timer%model_internal_time+ &
                 idt+s%timer%model_timestep,s%frq) == 0)
    else
      res = (mod(s%timer%model_internal_time+s%timer%model_timestep,s%frq) == 0)
    end if
  end function syncro_willact

  function init_alarm(t,dt,act0)
    implicit none
    type(rcm_alarm) , pointer :: init_alarm
    type(rcm_timer) , pointer , intent(inout) :: t
    real(rkx) , intent(in) :: dt
    logical , intent(in) , optional :: act0
    logical :: lact0
    type(rcm_alarm) , pointer :: alarm
    allocate(alarm)
    alarm%timer => null( )
    if ( associated(t) ) then
      alarm%timer => t
      lact0 = .false.
      if ( present(act0) ) lact0 = act0
      alarm%dt = dt
      alarm%idate = alarm%timer%idate
      alarm%intalm = rcm_time_interval(int(dt,ik8),usec)
      alarm%actint = int(dt,ik8)
      alarm%lastact = alarm%timer%model_internal_time
      alarm%triggered = lact0
      alarm%timer%nalarm = alarm%timer%nalarm + 1
      alarm%timer%ap(alarm%timer%nalarm)%ap => alarm
      alarm%lcount = alarm%timer%model_internal_time/alarm%actint
      alarm%rw = real(alarm%timer%model_timestep,rkx)/dt
      alarm%now = alarm%timer%model_internal_time
    end if
    init_alarm => alarm
  end function init_alarm

  subroutine alarm_check(alarm)
    class(rcm_alarm) , intent(inout) :: alarm
    real(rkx) :: t1 , t2
    if ( alarm%now == alarm%timer%model_internal_time ) then
      return
    end if
    t1 = real(alarm%timer%model_internal_time,rkx)
    t2 = real(alarm%lastact+alarm%actint,rkx)
    if ( t1 >= t2 ) then
      alarm%lastact = ((alarm%timer%model_internal_time) / &
      alarm%actint)*alarm%actint
      alarm%triggered = .true.
    end if
    if ( alarm%triggered ) then
      alarm%triggered = .false.
      alarm%idate = alarm%idate + alarm%intalm
      alarm%wt(1) = (t1 - t2)/alarm%dt
      alarm%wt(2) = d_one - alarm%wt(1)
      alarm%now = alarm%timer%model_internal_time
      alarm%lcount = alarm%lcount + 1
    end if
  end subroutine alarm_check

  logical function alarm_act(alarm) result(res)
    implicit none
    class(rcm_alarm) , intent(in) :: alarm
    res = .false.
    if ( alarm%now == alarm%timer%model_internal_time ) then
      res = .true.
      return
    end if
  end function alarm_act

  logical function alarm_willact(alarm,dt) result(res)
    implicit none
    class(rcm_alarm) , intent(in) :: alarm
    real(rkx) , intent(in) :: dt
    integer(ik8) :: t1 , t2 , idt
    res = .false.
    idt = int(dt,ik8)
    t1 = alarm%timer%model_internal_time+idt
    t2 = alarm%lastact+alarm%actint
    if ( t1 >= t2 ) then
      res = .true.
    end if
  end function alarm_willact

  real(rkx) function ratio_freq_syncro(s1,s2) result(res)
    implicit none
    type(rcm_syncro) , intent(in) :: s1 , s2
    res = real(s1%frq,rkx)/real(s2%frq,rkx)
  end function ratio_freq_syncro

  real(rkx) function ratio_freq_alarm(a1,a2) result(res)
    implicit none
    type(rcm_alarm) , intent(in) :: a1 , a2
    res = real(a1%actint,rkx)/real(a2%actint,rkx)
  end function ratio_freq_alarm

  real(rkx) function ratio_freq_syncro_alarm(s1,a2) result(res)
    implicit none
    type(rcm_syncro) , intent(in) :: s1
    type(rcm_alarm) , intent(in) :: a2
    res = real(s1%frq,rkx)/real(a2%actint,rkx)
  end function ratio_freq_syncro_alarm

  real(rkx) function ratio_freq_alarm_syncro(a1,s2) result(res)
    implicit none
    type(rcm_alarm) , intent(in) :: a1
    type(rcm_syncro) , intent(in) :: s2
    res = real(a1%actint,rkx)/real(s2%frq,rkx)
  end function ratio_freq_alarm_syncro

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
  type(rcm_alarm) , pointer :: srf_alarm , rad_alarm , cum_alarm
  type(rcm_alarm) , pointer :: srf_output

  mdate0 = idates(1)
  mdate1 = idates(2)
  mdate2 = idates(3)

  timer => rcm_timer(mdate0,mdate1,mdate2,213.0_rkx)

  print *, timer%str( ) , timer%step( )

  srf_alarm  => rcm_alarm(timer,600.0_rkx,.true.)
  srf_output => rcm_alarm(timer,3600.0_rkx*3.0,.true.)
  rad_alarm  => rcm_alarm(timer,1800.0_rkx,.true.)
  cum_alarm  => rcm_alarm(timer,300.0_rkx)

  do while ( .not. timer%reached_endtime )
    call timer%advance( )
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

  print *, timer%str( ) , timer%step( )

  deallocate(timer)

end program test_timing
#endif

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
