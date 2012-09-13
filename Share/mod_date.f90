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

module mod_date

  use mod_realkinds
  use mod_intkinds
  use mod_stdio
  use mod_message

  private

  integer(kind=ik8) , parameter :: i8spw = 604800_8
  integer(kind=ik8) , parameter :: i8spd = 86400_8
  integer(kind=ik8) , parameter :: i8mpd = 1440_8
  integer(kind=ik8) , parameter :: i8hpd = 24_8
  integer(kind=ik8) , parameter :: i8sph = 3600_8
  integer(kind=ik8) , parameter :: i8spm = 60_8
  integer(kind=ik8) , parameter :: i8mph = 60_8
  integer(kind=ik8) , parameter :: i8mpy = 12_8
  integer(kind=ik8) , parameter :: i8mpc = 1200_8
  integer(kind=ik8) , parameter :: i8ypc = 100_8

  integer(ik4) , public , parameter :: gregorian = 1
  integer(ik4) , public , parameter :: noleap    = 2
  integer(ik4) , public , parameter :: y360      = 3

  integer(ik4) , public , parameter :: usec = 1
  integer(ik4) , public , parameter :: umin = 2
  integer(ik4) , public , parameter :: uhrs = 3
  integer(ik4) , public , parameter :: uday = 4
  integer(ik4) , public , parameter :: umnt = 5
  integer(ik4) , public , parameter :: uyrs = 6
  integer(ik4) , public , parameter :: ucnt = 7

  integer(ik4) , parameter :: reference_year = 1900

  character (len=16) , public , dimension(7) :: cintstr
  character (len=12) , public , dimension(3) :: calstr

  integer(ik4) , dimension(12) :: mlen

  type rcm_time_and_date
    integer(ik4) :: calendar = gregorian
    integer(ik4) :: days_from_reference = 0
    integer(ik4) :: second_of_day = 0
  end type rcm_time_and_date

  type rcm_time_interval
    integer(kind=ik8) :: ival = 0
    integer(ik4) :: iunit = usec
  end type rcm_time_interval

  type iadate
    integer(ik4) :: calendar = gregorian
    integer(ik4) :: year = reference_year
    integer(ik4) :: month = 1
    integer(ik4) :: day = 1
  end type iadate

  type iatime
    integer(ik4) :: hour = 0
    integer(ik4) :: minute = 0
    integer(ik4) :: second = 0
  end type iatime

  interface assignment(=) 
    module procedure initfromintdt , initfromtypedt
    module procedure initfromintit , initfromdbleit , initfromtypeit
  end interface assignment(=)

  interface operator(==)
    module procedure isequaldt , isequalidt , isequalit
  end interface operator(==)

  interface operator(+)
    module procedure add_interval
  end interface operator(+)

  interface operator(-)
    module procedure diffdate , sub_interval
  end interface operator(-)

  interface operator(>)
    module procedure date_greater , interval_greater , idate_greater
  end interface

  interface operator(<)
    module procedure date_less , interval_less , idate_less
  end interface

  interface operator(>=)
    module procedure date_ge , idate_ge
  end interface

  interface operator(<=)
    module procedure date_le , idate_le
  end interface

  interface operator(/=)
    module procedure date_ne , idate_ne
  end interface

  interface setcal
    module procedure set_calint , set_calstring , set_caltype
  end interface

  interface split_idate
    module procedure split_i10 , split_rcm_time_and_date
  end interface

  public :: timeval2ym
  public :: rcm_time_and_date , assignment(=) , operator(==)
  public :: rcm_time_interval , operator(+) , operator(-)
  public :: operator(>) , operator(<) , operator(>=) , &
            operator(<=) , operator(/=)
  public :: print_rcm_time_and_date , print_rcm_time_interval
  public :: setcal
  public :: tochar , toint10 , tohours
  public :: lsamemonth , imondiff , lfhomonth , monfirst , monlast , monmiddle
  public :: nextmon , prevmon , yrfirst , nextwk , prevwk
  public :: lsameweek , iwkdiff , idayofweek , ifdoweek , ildoweek
  public :: timeval2date , lfdomonth , lfdoyear , lmidnight , yeardayfrac
  public :: split_idate , julianday
  public :: getyear , getmonth , getday
  public :: gethour , getminute , getsecond

  data mlen /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data calstr /'gregorian','noleap','360_day'/
  data cintstr /'seconds', 'minutes', 'hours', 'days', &
                'months', 'years', 'centuries'/

  contains

  function lleap(iyear)
    implicit none
    logical :: lleap
    integer(ik4) , intent(in) :: iyear
    if ( mod(iyear,400) == 0 .or.  &
        ( mod(iyear,4) == 0 .and. mod(iyear,100) /= 0 ) ) then
      lleap = .true.
    else
      lleap = .false.
    end if
  end function lleap

  integer function mdays_leap(iyear, imon) result(mdays)
    implicit none
    integer(ik4) , intent(in) :: iyear , imon
    if (imon /= 2) then
      mdays = mlen(imon)
    else
      mdays = mlen(2)
      if (lleap(iyear)) then
        mdays = mdays + 1
      end if
    end if
  end function mdays_leap

  integer function yeardays(y,c)
    implicit none
    integer(ik4) , intent(in) :: y , c
    select case (c)
      case (noleap)
        yeardays = 365
      case (y360)
        yeardays = 360
      case default
        yeardays = 365
        if ( lleap(y) ) yeardays = 366
    end select
  end function yeardays

  subroutine idayofyear_to_monthdate(j,y,c,m,d)
    implicit none
    integer(ik4) , intent(in) :: j , y , c
    integer(ik4) , intent(out) :: m , d
    integer(ik4) :: id , md
    id = j
    m = 1
    d = 1
    select case (c)
      case (noleap)
        do while (id > mlen(m))
          id = id - mlen(m)
          m = m + 1
        end do
      case (y360)
        do while (id > 30)
          id = id - 30
          m  = m + 1
        end do
      case default
        md = mdays_leap(y,1)
        do while (id > md)
          id = id - md
          m = m + 1
          md = mdays_leap(y,m)
        end do
    end select
    d = id
  end subroutine idayofyear_to_monthdate

  integer function idayofyear(x) result(id)
    implicit none
    type (iadate) , intent(in) :: x
    integer(ik4) :: i
    id = x%day
    select case (x%calendar)
      case (gregorian)
        do i = 1 , x%month-1
          id = id + mdays_leap(x%year, i)
        end do
      case (noleap)
        id = id + sum(mlen(1:x%month-1))
      case (y360)
        id = id + 30*(x%month-1)
    end select
  end function idayofyear

  subroutine date_to_days_from_reference(d,x)
    type(iadate) , intent(in) :: d
    type(rcm_time_and_date) , intent(inout) :: x
    integer(ik4) :: ny , ipm , id , iy
    x%calendar = d%calendar
    ny = d%year - reference_year
    ipm = isign(1,ny)
    ny = ny*ipm
    if ( ipm < 0 ) ny = ny-1
    iy = reference_year
    id = 0
    do while (ny > 0)
      id = id + ipm*yeardays(iy,d%calendar)
      iy = iy + ipm
      ny = ny - 1
    end do
    if ( ipm > 0 ) then
      x%days_from_reference = id + idayofyear(d) - 1
    else
      x%days_from_reference = id - (yeardays(iy,d%calendar)-idayofyear(d)) - 1
    end if
  end subroutine date_to_days_from_reference

  subroutine days_from_reference_to_date(x,d)
    type(rcm_time_and_date) , intent(in) :: x
    type(iadate) , intent(out) :: d
    integer(ik4) :: id , ipm , iy
    d%calendar = x%calendar
    id = x%days_from_reference
    ipm = isign(1,id)
    d%year = reference_year
    iy = yeardays(d%year,d%calendar)
    do while ((id*ipm) >= iy)
      d%year = d%year + ipm
      id = id - ipm*iy
      iy = yeardays(d%year,d%calendar)
    end do
    if (id >= 0) then
      call idayofyear_to_monthdate(id+1,d%year,d%calendar,d%month,d%day)
    else
      d%year = d%year - 1
      id = yeardays(d%year,d%calendar)+id
      call idayofyear_to_monthdate(id,d%year,d%calendar,d%month,d%day)
    end if
  end subroutine days_from_reference_to_date

  subroutine time_to_second_of_day(t,x)
    type(iatime) , intent(in) :: t
    type(rcm_time_and_date) , intent(inout) :: x
    x%second_of_day = t%hour*3600+t%minute*60+t%second
  end subroutine time_to_second_of_day

  subroutine second_of_day_to_time(x,t)
    type(rcm_time_and_date) , intent(in) :: x
    type(iatime) , intent(out) :: t
    integer(ik4) :: i1 , i2
    i1 = x%second_of_day
    i2 = i1/3600
    t%hour = i2
    i1 = i1-i2*3600
    i2 = i1/60
    t%minute = i2
    i1 = i1-i2*60
    t%second = i1
  end subroutine second_of_day_to_time

  subroutine date_time_to_internal(d,t,x)
    type(iadate) , intent(in) :: d
    type(iatime) , intent(in) :: t
    type(rcm_time_and_date) , intent(out) :: x
    call time_to_second_of_day(t,x)
    call date_to_days_from_reference(d,x)
  end subroutine date_time_to_internal

  subroutine internal_to_date_time(x,d,t)
    type(rcm_time_and_date) , intent(in) :: x
    type(iadate) , intent(out) :: d
    type(iatime) , intent(out) :: t
    call second_of_day_to_time(x,t)
    call days_from_reference_to_date(x,d)
  end subroutine internal_to_date_time

  subroutine adjustpm(a,b,i)
    implicit none
    integer(ik4) , intent(inout) :: a , b
    integer(ik4) , intent(in) :: i
    if ( a > i ) then
      a = a - i
      b = b + 1
    end if
  end subroutine adjustpm

  subroutine adjustmp(a,b,i)
    implicit none
    integer(ik4) , intent(inout) :: a , b
    integer(ik4) , intent(in) :: i
    if ( a < 1 ) then
      a = i - a
      b = b - 1
    end if
  end subroutine adjustmp

  subroutine normidate(idate)
    implicit none
    integer(ik4) , intent(inout) :: idate
    if (idate < 10000) idate = idate*1000000+10100
    if (idate < 1000000) idate = idate*10000+100
    if (idate < 100000000) idate = idate*100
  end subroutine normidate

  subroutine split_i10(idate, iy, im, id, ih)
    implicit none
    integer(ik4) , intent(in) :: idate
    integer(ik4) , intent(out) :: iy , im , id , ih
    integer(ik4) :: base , iidate
    iidate = idate
    call normidate(iidate)
    base = iidate
    iy = base/1000000
    base = base-iy*1000000
    im = base/10000
    base = base-im*10000
    id = base/100
    base = base-id*100
    ih = base
  end subroutine split_i10

  subroutine initfromintdt(x, i)
    implicit none
    integer(ik4) , intent(in) :: i
    type (rcm_time_and_date) , intent(out) :: x
    type (iadate) :: d
    type (iatime) :: t
    call split_i10(i, d%year, d%month, d%day, t%hour)
    d%calendar = gregorian
    call date_time_to_internal(d,t,x)
  end subroutine initfromintdt

  subroutine initfromintit(x, i)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    integer(ik4) , intent(in) :: i
    x%ival = i
    x%iunit = usec
  end subroutine initfromintit

  subroutine initfromdbleit(x, d)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    real(rk8) , intent(in) :: d
    x%ival = d
    x%iunit = usec
  end subroutine initfromdbleit

  subroutine set_caltype(x, y)
    implicit none
    type (rcm_time_and_date) , intent(inout) :: x
    type (rcm_time_and_date) , intent(in) :: y
    x%calendar = y%calendar
  end subroutine set_caltype

  subroutine set_calint(x, c)
    implicit none
    type (rcm_time_and_date) , intent(inout) :: x
    integer(ik4) , intent(in) :: c
    type (iadate) :: d
    if ( c /= x%calendar) then
      call days_from_reference_to_date(x,d)
      select case (c)
        case (gregorian)
          d%calendar = c
        case (noleap)
          d%calendar = c
        case (y360)
          d%calendar = c
        case default
          write (stderr,*) 'Unknown calendar, using Julian/Gregorian'
          d%calendar = gregorian
      end select
      call date_to_days_from_reference(d,x)
    end if
  end subroutine set_calint

  subroutine set_calstring(x, c)
    implicit none
    type (rcm_time_and_date) , intent(inout) :: x
    character (len=*) , intent(in) :: c
    integer(ik4) :: ic
    if ( c == 'gregorian' ) then
      ic = gregorian
    else if ( c(1:6) == 'noleap' .or.   &
              c(1:8) == 'days_365' .or. &
              c(1:7) == '365_day' ) then
      ic = noleap
    else if ( c(1:7) == 'days_360' ) then
      ic = y360
    else
      write (stderr,*) 'Unknown calendar, using Julian/Gregorian'
      ic = gregorian
    end if
    call set_calint(x,ic)
  end subroutine set_calstring

  subroutine set_timeunit(x, u)
    implicit none
    type (rcm_time_interval) , intent(inout) :: x
    integer(ik4) , intent(in) :: u
    select case (u)
      case (usec)
        x%iunit = u
      case (umin)
        x%iunit = u
      case (uhrs)
        x%iunit = u
      case (uday)
        x%iunit = u
      case (umnt)
        x%iunit = u
      case (uyrs)
        x%iunit = u
      case (ucnt)
        x%iunit = u
      case default
        write (stderr,*) 'Unknown time unit, assuming hours'
        x%iunit = uhrs
    end select
  end subroutine set_timeunit

  subroutine initfromtypedt(x, y)
    implicit none
    type (rcm_time_and_date) , intent(out) :: x
    type (rcm_time_and_date) , intent(in) :: y
    x%calendar            = y%calendar
    x%days_from_reference = y%days_from_reference
    x%second_of_day       = y%second_of_day
  end subroutine initfromtypedt

  subroutine initfromtypeit(x, y)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    type (rcm_time_interval) , intent(in) :: y
    x%iunit = y%iunit
    x%ival  = y%ival
  end subroutine initfromtypeit

  function tochar(x) result(cdat)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    character (len=32) :: cdat
    type (iadate) :: d
    type (iatime) :: t
    call internal_to_date_time(x,d,t)
    write (cdat,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2," UTC")') &
       d%year, d%month, d%day, t%hour, t%minute, t%second
  end function tochar

  subroutine print_rcm_time_and_date(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    write (stdout,*) tochar(x)
  end subroutine print_rcm_time_and_date

  subroutine print_rcm_time_interval(x)
    implicit none
    type (rcm_time_interval) , intent(in) :: x
    write (stdout,'(i16," ",a)') x%ival, trim(cintstr(x%iunit))
  end subroutine print_rcm_time_interval

  logical function isequaldt(x, y)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) , intent(in) :: y
    call check_cal(x,y)
    isequaldt = ( x%days_from_reference == y%days_from_reference) .and. &
                ( x%second_of_day       == y%second_of_day)
  end function isequaldt

  logical function isequalidt(x, y)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call set_calint(yy,x%calendar)
    isequalidt = isequaldt(x,yy)
  end function isequalidt

  logical function isequalit(x, y)
    implicit none
    type (rcm_time_interval) , intent(in) :: x
    type (rcm_time_interval) , intent(in) :: y
    select case (x%iunit)
      case (usec)
        select case (y%iunit)
          case (usec)
            isequalit = ( x%ival == y%ival )
          case (umin)
            isequalit = ( x%ival == y%ival*i8spm )
          case (uhrs)
            isequalit = ( x%ival == y%ival*i8sph )
          case (uday)
            isequalit = ( x%ival == y%ival*i8spd )
          case default
            isequalit = .false.
        end select
      case (umin)
        select case (y%iunit)
          case (usec)
            isequalit = ( x%ival == y%ival/i8spm )
          case (umin)
            isequalit = ( x%ival == y%ival )
          case (uhrs)
            isequalit = ( x%ival == y%ival*i8mph )
          case (uday)
            isequalit = ( x%ival == y%ival*i8mpd )
          case default
            isequalit = .false.
        end select
      case (uhrs)
        select case (y%iunit)
          case (usec)
            isequalit = ( x%ival == y%ival/i8sph )
          case (umin)
            isequalit = ( x%ival == y%ival/i8mph )
          case (uhrs)
            isequalit = ( x%ival == y%ival )
          case (uday)
            isequalit = ( x%ival == y%ival*i8hpd )
          case default
            isequalit = .false.
        end select
      case (uday)
        select case (y%iunit)
          case (usec)
            isequalit = ( x%ival == y%ival/i8spd )
          case (umin)
            isequalit = ( x%ival == y%ival/i8mpd )
          case (uhrs)
            isequalit = ( x%ival == y%ival*i8hpd )
          case (uday)
            isequalit = ( x%ival == y%ival )
          case default
            isequalit = .false.
        end select
      case (umnt)
        select case (y%iunit)
          case (umnt)
            isequalit = ( x%ival == y%ival )
          case (uyrs)
            isequalit = ( x%ival == y%ival*i8mpy )
          case (ucnt)
            isequalit = ( x%ival == y%ival*i8mpc )
          case default
            isequalit = .false.
        end select
      case (uyrs)
        select case (y%iunit)
          case (umnt)
            isequalit = ( x%ival == y%ival*i8mpy )
          case (uyrs)
            isequalit = ( x%ival == y%ival )
          case (ucnt)
            isequalit = ( x%ival == y%ival*i8ypc )
          case default
            isequalit = .false.
        end select
      case (ucnt)
        select case (y%iunit)
          case (umnt)
            isequalit = ( x%ival == y%ival*i8mpc )
          case (uyrs)
            isequalit = ( x%ival == y%ival*i8ypc )
          case (ucnt)
            isequalit = ( x%ival == y%ival )
          case default
            isequalit = .false.
        end select
      case default
        isequalit = ( x%iunit == y%iunit ) .and. ( x%ival == y%ival )
    end select
  end function isequalit

  recursive subroutine add_days_leap(d,m,y,a)
    integer(ik4) , intent(inout) :: d , m , y
    integer(ik4) , intent(in) :: a
    integer(ik4) :: tmp
    tmp = a
    if ( tmp > (mdays_leap(y,m)-d) ) then
      tmp = tmp - (mdays_leap(y,m)-d+1)
      m = m+1
      call adjustpm(m,y,12)
      d = 1
      call add_days_leap(d,m,y,tmp)
    else
      d = d + tmp
    end if
  end subroutine add_days_leap

  recursive subroutine sub_days_leap(d,m,y,a)
    integer(ik4) , intent(inout) :: d , m , y
    integer(ik4) , intent(in) :: a
    integer(ik4) :: tmp
    tmp = a
    if ( tmp >= d ) then
      tmp = tmp-d
      m = m-1
      call adjustmp(m,y,12)
      d = mdays_leap(y,m)
      call sub_days_leap(d,m,y,tmp)
    else
      d = d - tmp
    end if
  end subroutine sub_days_leap

  recursive subroutine add_days_noleap(d,m,y,a)
    integer(ik4) , intent(inout) :: d , m , y
    integer(ik4) , intent(in) :: a
    integer(ik4) :: tmp
    tmp = a
    if ( tmp > mlen(m)-d ) then
      tmp = tmp - mlen(m)-d+1
      m = m + 1
      call adjustpm(m,y,12)
      d = 1
      call add_days_noleap(d,m,y,tmp)
    else
      d = d + tmp
    end if
  end subroutine add_days_noleap

  recursive subroutine sub_days_noleap(d,m,y,a)
    integer(ik4) , intent(inout) :: d , m , y
    integer(ik4) , intent(in) :: a
    integer(ik4) :: tmp
    tmp = a
    if ( tmp > d ) then
      tmp = tmp-d
      m = m-1
      call adjustmp(m,y,12)
      d = mlen(m)
      call sub_days_noleap(d,m,y,tmp)
    else
      d = d - tmp
    end if
  end subroutine sub_days_noleap

  function add_interval(x, y) result (z)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_interval) , intent(in) :: y
    type (rcm_time_and_date) :: z
    type (iadate) :: d
    real(rk8) :: dm
    integer(ik8) :: tmp
    z = x
    tmp = y%ival
    select case (y%iunit)
      case (usec)
        tmp = tmp + z%second_of_day
        z%second_of_day = mod(tmp, i8spd)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = i8spd+z%second_of_day
          tmp = tmp-i8spd
        end if
        tmp = tmp/i8spd
        z%days_from_reference = z%days_from_reference + tmp
      case (umin)
        tmp = tmp*i8spm+z%second_of_day
        z%second_of_day = mod(tmp, i8spd)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = i8spd+z%second_of_day
          tmp = tmp-i8spd
        end if
        tmp = tmp/i8spd
        z%days_from_reference = z%days_from_reference + tmp
      case (uhrs)
        tmp = tmp*i8sph+z%second_of_day
        z%second_of_day = mod(tmp, i8spd)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = i8spd+z%second_of_day
          tmp = tmp-i8spd
        end if
        tmp = tmp/i8spd
        z%days_from_reference = z%days_from_reference + tmp
      case (uday)
        z%days_from_reference = z%days_from_reference + tmp
      case (umnt)
        call days_from_reference_to_date(x,d)
        ! Adjust date of the month: This is really a trick...
        select case (z%calendar)
          case (gregorian)
            dm = dble(d%day)/dble(mdays_leap(d%year, d%month))
          case (noleap)
            dm = dble(d%day)/dble(mlen(d%month))
        end select
        d%month = d%month+mod(tmp,i8mpy)
        call adjustpm(d%month,d%year,12)
        tmp = tmp/i8mpy
        d%year = d%year+tmp
        ! Adjust date of the month: This is really a trick...
        select case (z%calendar)
          case (gregorian)
            d%day = idnint(dble(mdays_leap(d%year, d%month)) * dm)
          case (noleap)
            d%day = idnint(dble(mlen(d%month)) * dm)
        end select
        call date_to_days_from_reference(d,z)
      case (uyrs)
        call days_from_reference_to_date(x,d)
        d%year = d%year+tmp
        call date_to_days_from_reference(d,z)
      case (ucnt)
        call days_from_reference_to_date(x,d)
        d%year = d%year+i8ypc*tmp
        call date_to_days_from_reference(d,z)
    end select
  end function add_interval

  function julianday(iy, im, id)
    implicit none
    real(rk8) :: julianday
    integer(ik4) , intent(in) :: iy , im , id
    integer(ik4) :: ia , ib , iiy , iim
    iiy = iy
    iim = im
    if (iim <= 2) then
      iiy = iiy - 1
      iim = iim + 12
    end if
    if (lcaltype(iy,im,id)) then
      ia = iiy/100
      ib = 2 - ia + ia / 4
    else
      ib = 0
    end if
    julianday = dint(365.25D+00*dble(iiy+4716)) + &
                dint(30.6001D+00*dble(iim+1))   + &
                dble(id + ib) - 1524.5D+00
    contains
      function lcaltype(iy, im, id)
      implicit none
      logical :: lcaltype
      integer(ik4) :: icaltype
      integer(ik4) , intent(in) :: iy , im , id
      ! Standard Julian/Gregorian switch
      ! Return true  if before 1582-10-04
      !        false if after  1582-10-15
      icaltype = 0
      if (iy < 1582) then
        icaltype = 1
      else if (iy == 1582) then
        if (im < 10) then
          icaltype = 1
        else if (im == 10) then
          if (id <= 4) then
            icaltype = 1
          end if
        end if
      end if
      if (iy > 1582) then
        icaltype = 2
      else if (iy == 1582) then
        if (im > 10) then
          icaltype = 2
        else if (im == 10) then
          if (id >= 15) then
            icaltype = 2
          end if
        end if
      end if
      if (icaltype == 1) then
        lcaltype = .false.
      else if (icaltype == 2) then
        lcaltype = .true.
      else
        write (stderr, *) 'year  = ', iy
        write (stderr, *) 'month = ', im
        write (stderr, *) 'day   = ', id
        call die('mod_date','Day non existent, inside Julian/Gregorian jump',1)
      end if
    end function lcaltype
  end function julianday

  subroutine check_cal(x,y)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) , intent(in) :: y
    if ( x%calendar /= y%calendar ) then
      write (stderr,*) 'X calendar = ', calstr(x%calendar)
      write (stderr,*) 'Y calendar = ', calstr(y%calendar)
      call die('mod_date','Dates not comparable: not using same calendar.',1)
    end if
  end subroutine check_cal

  function diffdate(x,y) result(z)
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) , intent(in) :: y
    type (rcm_time_interval) :: z
    call check_cal(x,y)
    z%iunit = usec
    z%ival = (x%second_of_day-y%second_of_day) + &
             (x%days_from_reference-y%days_from_reference)*i8spd
  end function diffdate

  function sub_interval(x,y) result (z)
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_interval) , intent(in) :: y
    type (rcm_time_and_date) :: z
    type (iadate) :: d
    integer(kind=ik8) :: tmp
    z = x
    tmp = y%ival
    select case (y%iunit)
      case (usec)
        tmp = tmp - z%second_of_day
        z%second_of_day = mod(tmp,i8spd)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = i8spd+z%second_of_day
          tmp = tmp+i8spd
        end if
        tmp = tmp/i8spd
        z%days_from_reference = z%days_from_reference - tmp
      case (umin)
        tmp = tmp*i8spm - z%second_of_day
        z%second_of_day = mod(tmp, i8spd)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = i8spd+z%second_of_day
          tmp = tmp+i8spd
        end if
        tmp = tmp/i8spd
        z%days_from_reference = z%days_from_reference - tmp
      case (uhrs)
        tmp = tmp*i8sph - z%second_of_day
        z%second_of_day = mod(tmp, i8spd)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = i8spd+z%second_of_day
          tmp = tmp+i8spd
        end if
        tmp = tmp/i8spd
        z%days_from_reference = z%days_from_reference - tmp
      case (uday)
        z%days_from_reference = z%days_from_reference - tmp
      case (umnt)
        call days_from_reference_to_date(x,d)
        d%month = d%month-mod(tmp,i8mpy)
        call adjustmp(d%month,d%year,12)
        d%year = d%year-tmp/i8mpy
        call date_to_days_from_reference(d,z)
      case (uyrs)
        call days_from_reference_to_date(x,d)
        d%year = d%year-tmp
        call date_to_days_from_reference(d,z)
      case (ucnt)
        call days_from_reference_to_date(x,d)
        d%year = d%year-i8ypc*tmp
        call date_to_days_from_reference(d,z)
    end select
  end function sub_interval

  integer function toint10(x) result(z)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    type (iatime) :: t
    call internal_to_date_time(x,d,t)
    z = d%year*1000000+d%month*10000+d%day*100+t%hour;
  end function toint10

  logical function lfdoyear(x) result(l11)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    l11 = (d%month == 1 .and. d%day == 1)
  end function

  logical function lfdomonth(x) result(l1)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    l1 = (d%day == 1)
  end function

  logical function lmidnight(x) result(lm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    lm = (x%second_of_day == 0)
  end function

  logical function lsamemonth(x,y) result(lm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    type (iadate) :: d1 , d2
    call check_cal(x,y)
    call days_from_reference_to_date(x,d1)
    call days_from_reference_to_date(y,d2)
    lm = (d1%month == d2%month)
  end function

  integer function imondiff(x,y) result(z)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    type (iadate) :: d1 , d2
    call check_cal(x,y)
    call days_from_reference_to_date(x,d1)
    call days_from_reference_to_date(y,d2)
    z = (d1%year-d2%year)*12+(d1%month-d2%month)
  end function imondiff

  logical function lfhomonth(x) result(lf)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mm
    mm = monmiddle(x)
    lf = (x < mm)
  end function lfhomonth

  integer function idayofweek(x) result(iday)
    ! Sun Mon Tue Wed Thu Fri Sat
    ! 1   2   3   4   5   6   7
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    real(rk8) :: jd
    if (x%calendar /= gregorian) then
      call die('mod_date', &
               'Error: week concept works only for gregorian calendar.',1)
    end if
    call days_from_reference_to_date(x,d)
    jd = julianday(d%year, d%month, d%day)
    iday = idint(dmod(jd+1.5D+00, 7.0D+00))+1
  end function idayofweek

  function setmidnight(x) result(mx)
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mx
    mx = x
    mx%second_of_day = 0
  end function setmidnight

  logical function lsameweek(x,y) result(ls)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    type (rcm_time_and_date) :: xx , x1 , x2
    integer(ik4) :: idwk
    xx = setmidnight(x)
    idwk = idayofweek(xx)
    x1 = xx - rcm_time_interval(idwk,uday)
    x2 = xx + rcm_time_interval(7-idwk,uday)
    ls = (y > x1) .and. (y < x2)
  end function lsameweek

  function ifdoweek(x) result(ifd)
    implicit none
    integer(ik4) :: iwkday
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_interval) :: z
    type (rcm_time_and_date) :: ifd
    iwkday = idayofweek(x) - 1
    z = rcm_time_interval(iwkday,uday)
    ifd = setmidnight(x - z)
  end function ifdoweek

  function ildoweek(x) result(ild)
    implicit none
    integer(ik4) :: iwkday
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_interval) :: z
    type (rcm_time_and_date) :: ild
    iwkday = 7 - idayofweek(x)
    z = rcm_time_interval(iwkday,uday)
    ild = setmidnight(x + z)
  end function ildoweek

  integer function iwkdiff(x,y) result(iwk)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    type (rcm_time_interval) :: z
    call check_cal(x,y)
    z = x - y
    iwk = int(z%ival/i8spw) + 1
  end function iwkdiff

  function monfirst(x) result(mf)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mf
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    d%day = 1
    call date_to_days_from_reference(d,mf)
  end function monfirst

  function yrfirst(x) result(yf)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: yf
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    d%month = 1
    d%day = 1
    call date_to_days_from_reference(d,yf)
  end function yrfirst

  function monlast(x) result(ml)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: ml
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    select case (d%calendar)
      case (gregorian)
        d%day = mdays_leap(d%year, d%month)
      case (noleap)
        d%day = mlen(d%month)
      case (y360)
        d%day = 30
    end select
    call date_to_days_from_reference(d,ml)
  end function monlast

  function monmiddle(x) result(mm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mm
    integer(ik4) :: dm
    type (iadate) :: d
    type (iatime) :: t
    call days_from_reference_to_date(x,d)
    select case (x%calendar)
      case (gregorian)
        dm = mdays_leap(d%year, d%month)
        d%day = dm/2
        t%hour = (dm-d%day*2)*12
      case (noleap)
        dm = mlen(d%month)
        d%day = dm/2
        t%hour = (dm-d%day*2)*12
      case (y360)
        d%day = 15
    end select
    call date_time_to_internal(d,t,mm)
  end function monmiddle

  function nextmon(x) result(nm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: nm
    type (rcm_time_interval) :: z
    z = rcm_time_interval(1,umnt)
    nm = x + z
  end function nextmon

  function prevmon(x) result(pm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: pm
    type (rcm_time_interval) :: z
    z = rcm_time_interval(1,umnt)
    pm = x - z
  end function prevmon

  function nextwk(x) result(nm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: nm
    type (rcm_time_interval) :: z
    z = rcm_time_interval(7,uday)
    nm = x + z
  end function nextwk

  function prevwk(x) result(pm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: pm
    type (rcm_time_interval) :: z
    z = rcm_time_interval(7,uday)
    pm = x - z
  end function prevwk

  subroutine timeval2ym(xval,cunit,year,month,julday)
    implicit none
    real(rk8) , intent(in) :: xval
    character(*) , intent(in) :: cunit
    integer(ik4) , intent(out) :: year , month, julday

    type (iadate) , save :: d
    type (iatime) :: t
    character(16) :: cdum

    character(64) , save :: csave,dsave
    data csave /'months since XXXX-XX-XX XX:XX:XX XXX'/
    data dsave /'days since XXXX-XX-XX XX:XX:XX XXX'/

    if (cunit(1:6) == 'months') then

    
    if (csave == cunit) then
      year = d%year+idint(xval/12.0D0)
      month = d%month+idint(mod(xval,12.0D0))
      if ( month > 12 ) then
        month = month - 12
        year = year + 1
      else if ( month < 1 ) then
        month = month + 12
        year = year - 1
      end if
      return
    end if

      ! Unit is months since reference
      if (len_trim(cunit) >= 31) then
        read(cunit,'(a13,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') &
          cdum, d%year, cdum, d%month, cdum, d%day, &
          cdum, t%hour, cdum, t%minute, cdum, t%second
      else if (len_trim(cunit) >= 28) then
        read(cunit,'(a13,i4,a1,i2,a1,i2,a1,i2,a1,i2)') &
          cdum, d%year, cdum, d%month, cdum, d%day, &
          cdum, t%hour, cdum, t%minute
      else if (len_trim(cunit) >= 25) then
        read(cunit,'(a13,i4,a1,i2,a1,i2,a1,i2)') &
          cdum, d%year, cdum, d%month, cdum, d%day, cdum, t%hour
      else if (len_trim(cunit) >= 22) then
        read(cunit,'(a13,i4,a1,i2,a1,i2)') &
          cdum, d%year, cdum, d%month, cdum, d%day
      else if (len_trim(cunit) >= 19) then
        read(cunit,'(a13,i4,a1,i2)') &
          cdum, d%year, cdum, d%month
      else
        call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2YM')
      end if

    csave = cunit
    year = d%year+idint(xval/12.0D0)
    month = d%month+idint(mod(xval,12.0D0))
    if ( month > 12 ) then
      month = month - 12
      year = year + 1
    else if ( month < 1 ) then
      month = month + 12
      year = year - 1
    end if


    elseif (cunit(1:4) == 'months') then

      if (dsave == cunit) then
       year = d%year+idint(xval/366)
       julday = d%day  + xval - idint(xval/366)* 366     
       month = 0 
       return
       end if

      if (len_trim(cunit) >= 31) then
        read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') &
          cdum, d%year, cdum, d%month, cdum, d%day, &
          cdum, t%hour, cdum, t%minute, cdum, t%second
      else if (len_trim(cunit) >= 28) then
        read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2,a1,i2)') &
          cdum, d%year, cdum, d%month, cdum, d%day, &
          cdum, t%hour, cdum, t%minute
      else if (len_trim(cunit) >= 25) then
        read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2)') &
          cdum, d%year, cdum, d%month, cdum, d%day, cdum, t%hour
      else if (len_trim(cunit) >= 22) then
        read(cunit,'(a11,i4,a1,i2,a1,i2)') &
          cdum, d%year, cdum, d%month, cdum, d%day
      else if (len_trim(cunit) >= 19) then
        read(cunit,'(a11,i4,a1,i2)') &
          cdum, d%year, cdum, d%month
      else
        call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2YM')
      end if

      dsave = cunit
      year = d%year+ idint(xval/366)
      julday = d%day  + xval - idint(xval/366) * 366     
      month=0
    else 
       call die('mod_date','TIME UNIT IN TIMEVAL2YM MUST BE MONTHS or DAYS')
    end if

  end subroutine timeval2ym

  function timeval2date(xval,cunit,ccal) result(dd)
    implicit none
    real(rk8) , intent(in) :: xval
    character(*) , intent(in) :: cunit
    character(*) , intent(in) :: ccal

    type (rcm_time_and_date) :: dd
    type (rcm_time_interval) :: z , zz
    character(64) , save :: csave
    character(16) :: safeccal
    integer(ik4) , save :: iunit
    type (rcm_time_and_date) , save :: dref
    character(16) :: cdum
    type (iadate) :: d
    type (iatime) :: t

    data csave/'none'/

    t%hour = 0
    t%minute = 0
    t%second = 0

    safeccal = ccal

    if (csave == cunit) then
      z = xval
      z%iunit = iunit
      dd = dref + z
      if (iunit == uday) then
        zz%ival = idint((xval-dble(z%ival))*24.0D0)
        zz%iunit = uhrs
        if ( zz%ival /= 0 ) then
          dd = dd + zz
        end if
      end if
    else
      csave = cunit
      z = xval
      if (cunit(1:5) == 'hours') then
        ! Unit is hours since reference
        if (len_trim(cunit) >= 30) then
          read(cunit,'(a12,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, &
            cdum, t%hour, cdum, t%minute, cdum, t%second
        else if (len_trim(cunit) >= 27) then
          read(cunit,'(a12,i4,a1,i2,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, &
            cdum, t%hour, cdum, t%minute
        else if (len_trim(cunit) >= 24) then
          read(cunit,'(a12,i4,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, cdum, t%hour
        else
          call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
        end if
        iunit = uhrs
      else if (cunit(1:7) == 'minutes') then
        ! Unit is minutes since reference
        if (len_trim(cunit) >= 32) then
          read(cunit,'(a14,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, &
            cdum, t%hour, cdum, t%minute, cdum, t%second
        else if (len_trim(cunit) >= 29) then
          read(cunit,'(a14,i4,a1,i2,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, &
            cdum, t%hour, cdum, t%minute
        else if (len_trim(cunit) >= 26) then
          read(cunit,'(a14,i4,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, cdum, t%hour
        else
          call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
        end if
        iunit = umin
      else if (cunit(1:7) == 'seconds') then
        ! Unit is seconds since reference
        if (len_trim(cunit) >= 32) then
          read(cunit,'(a14,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, &
            cdum, t%hour, cdum, t%minute, cdum, t%second
        else if (len_trim(cunit) >= 29) then
          read(cunit,'(a14,i4,a1,i2,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, &
            cdum, t%hour, cdum, t%minute
        else if (len_trim(cunit) >= 26) then
          read(cunit,'(a14,i4,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, cdum, t%hour
        else
          call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
        end if
        iunit = usec
      else if (cunit(1:4) == 'days') then
        ! Unit is days since reference
        if (len_trim(cunit) >= 30) then
          read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, &
            cdum, t%hour, cdum, t%minute, cdum, t%second
        else if (len_trim(cunit) >= 27) then
          read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, &
            cdum, t%hour, cdum, t%minute
        else if (len_trim(cunit) >= 24) then
          read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day, cdum, t%hour
        else if (len_trim(cunit) >= 21) then
          read(cunit,'(a11,i4,a1,i2,a1,i2)') &
            cdum, d%year, cdum, d%month, cdum, d%day
        else if (len_trim(cunit) >= 19) then
          read(cunit,'(a11,i4,a1,i1,a1,i1)') &
            cdum, d%year, cdum, d%month, cdum, d%day
        else
          call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
        end if
        iunit = uday
      else
        call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
      end if
      if ( safeccal(1:6) == 'noleap' .or.   &
           safeccal(1:8) == 'days_365' .or. &
           safeccal(1:7) == '365_day' ) then
        d%calendar = noleap
      else if (safeccal(1:7) == '360_day') then
        d%calendar = y360
      else
        d%calendar = gregorian
      end if
      call date_time_to_internal(d,t,dref)
      z%iunit = iunit
      dd = dref + z
      if (iunit == uday) then
        zz%ival = idint((xval-dble(z%ival))*24.0D0)
        zz%iunit = uhrs
        if (zz%ival /= 0) then
          dd = dd + zz
        end if
      end if
    end if
  end function timeval2date
!
  logical function interval_greater(x, y)
    implicit none
    type (rcm_time_interval) , intent(in) :: x
    type (rcm_time_interval) , intent(in) :: y
    select case (x%iunit)
      case (usec)
        select case (y%iunit)
          case (usec)
            interval_greater = ( x%ival > y%ival )
          case (umin)
            interval_greater = ( x%ival > y%ival*i8spm )
          case (uhrs)
            interval_greater = ( x%ival > y%ival*i8sph )
          case (uday)
            interval_greater = ( x%ival > y%ival*i8spd )
          case default
            interval_greater = .false.
        end select
      case (umin)
        select case (y%iunit)
          case (usec)
            interval_greater = (y%ival > i8sph) .and. ( x%ival > y%ival/i8sph )
          case (umin)
            interval_greater = ( x%ival > y%ival )
          case (uhrs)
            interval_greater = ( x%ival > y%ival*i8mph )
          case (uday)
            interval_greater = ( x%ival > y%ival*i8mpd )
          case default
            interval_greater = .false.
        end select
      case (uhrs)
        select case (y%iunit)
          case (usec)
            interval_greater = (y%ival > i8sph) .and. ( x%ival > y%ival/i8sph )
          case (umin)
            interval_greater = (y%ival > i8mph) .and. ( x%ival > y%ival/i8mph )
          case (uhrs)
            interval_greater = ( x%ival > y%ival )
          case (uday)
            interval_greater = ( x%ival > y%ival*i8hpd )
          case default
            interval_greater = .false.
        end select
      case (uday)
        select case (y%iunit)
          case (usec)
            interval_greater = (y%ival > i8spd) .and. ( x%ival > y%ival/i8spd )
          case (umin)
            interval_greater = (y%ival > i8mpd) .and. ( x%ival > y%ival/i8mpd )
          case (uhrs)
            interval_greater = (y%ival > i8hpd) .and. ( x%ival > y%ival/i8hpd )
          case (uday)
            interval_greater = ( x%ival > y%ival )
          case default
            interval_greater = .false.
        end select
      case (umnt)
        select case (y%iunit)
          case (umnt)
            interval_greater = ( x%ival > y%ival )
          case (uyrs)
            interval_greater = ( x%ival > y%ival*i8mpy )
          case (ucnt)
            interval_greater = ( x%ival > y%ival*i8mpc )
          case default
            interval_greater = .false.
        end select
      case (uyrs)
        select case (y%iunit)
          case (umnt)
            interval_greater = ( x%ival > y%ival*i8mpy )
          case (uyrs)
            interval_greater = ( x%ival > y%ival )
          case (ucnt)
            interval_greater = ( x%ival > y%ival*i8ypc )
          case default
            interval_greater = .false.
        end select
      case (ucnt)
        select case (y%iunit)
          case (umnt)
            interval_greater = ( x%ival > y%ival*i8mpc )
          case (uyrs)
            interval_greater = ( x%ival > y%ival*i8ypc )
          case (ucnt)
            interval_greater = ( x%ival == y%ival )
          case default
            interval_greater = .false.
        end select
      case default
        interval_greater = ( x%iunit > y%iunit ) .and. ( x%ival > y%ival )
    end select
  end function interval_greater

  logical function interval_less(x, y)
    implicit none
    type (rcm_time_interval) , intent(in) :: x
    type (rcm_time_interval) , intent(in) :: y
    select case (x%iunit)
      case (usec)
        select case (y%iunit)
          case (usec)
            interval_less = ( x%ival < y%ival )
          case (umin)
            interval_less = ( x%ival < y%ival*i8spm )
          case (uhrs)
            interval_less = ( x%ival < y%ival*i8sph )
          case (uday)
            interval_less = ( x%ival < y%ival*i8spd )
          case default
            interval_less = .false.
        end select
      case (umin)
        select case (y%iunit)
          case (usec)
            interval_less = (y%ival > i8spm) .and. ( x%ival < y%ival/i8spm )
          case (umin)
            interval_less = ( x%ival < y%ival )
          case (uhrs)
            interval_less = ( x%ival < y%ival*i8mph )
          case (uday)
            interval_less = ( x%ival < y%ival*i8mpd )
          case default
            interval_less = .false.
        end select
      case (uhrs)
        select case (y%iunit)
          case (usec)
            interval_less = (y%ival > i8sph) .and. ( x%ival < y%ival/i8sph )
          case (umin)
            interval_less = (y%ival > i8mph) .and. ( x%ival < y%ival/i8mph )
          case (uhrs)
            interval_less = ( x%ival < y%ival )
          case (uday)
            interval_less = ( x%ival < y%ival*i8hpd )
          case default
            interval_less = .false.
        end select
      case (uday)
        select case (y%iunit)
          case (usec)
            interval_less = (y%ival > i8spd) .and. ( x%ival < y%ival/i8spd )
          case (umin)
            interval_less = (y%ival > i8mpd) .and. ( x%ival < y%ival/i8mpd )
          case (uhrs)
            interval_less = (y%ival > i8hpd) .and. ( x%ival < y%ival/i8hpd )
          case (uday)
            interval_less = ( x%ival < y%ival )
          case default
            interval_less = .false.
        end select
      case (umnt)
        select case (y%iunit)
          case (umnt)
            interval_less = ( x%ival < y%ival )
          case (uyrs)
            interval_less = ( x%ival < y%ival*i8mpy )
          case (ucnt)
            interval_less = ( x%ival < y%ival*i8mpc )
          case default
            interval_less = .false.
        end select
      case (uyrs)
        select case (y%iunit)
          case (umnt)
            interval_less = ( x%ival < y%ival*i8mpy )
          case (uyrs)
            interval_less = ( x%ival < y%ival )
          case (ucnt)
            interval_less = ( x%ival < y%ival*i8ypc )
          case default
            interval_less = .false.
        end select
      case (ucnt)
        select case (y%iunit)
          case (umnt)
            interval_less = ( x%ival < y%ival*i8mpc )
          case (uyrs)
            interval_less = ( x%ival < y%ival*i8ypc )
          case (ucnt)
            interval_less = ( x%ival == y%ival )
          case default
            interval_less = .false.
        end select
      case default
        interval_less = ( x%iunit < y%iunit ) .and. ( x%ival < y%ival )
    end select
  end function interval_less

  logical function date_greater(x,y) result(gt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    gt = (x%days_from_reference > y%days_from_reference)
    if (gt) return
    if (x%days_from_reference == y%days_from_reference) then
      gt = (x%second_of_day > y%second_of_day)
    end if
  end function date_greater

  logical function idate_greater(x,y) result(gt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call set_calint(yy,x%calendar)
    gt = date_greater(x,yy)
  end function idate_greater

  logical function date_less(x,y) result(lt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    lt = (x%days_from_reference < y%days_from_reference)
    if (lt) return
    if (x%days_from_reference == y%days_from_reference) then
      lt = (x%second_of_day < y%second_of_day)
    end if
  end function date_less

  logical function idate_less(x,y) result(lt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call set_calint(yy,x%calendar)
    lt = date_less(x,yy)
  end function idate_less

  logical function date_ge(x,y) result(gt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    gt = (x%days_from_reference == y%days_from_reference .and. &
          x%second_of_day       == y%second_of_day)
    if (gt) return
    gt = date_greater(x,y)
  end function date_ge

  logical function idate_ge(x,y) result(gt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call set_calint(yy,x%calendar)
    gt = date_ge(x,yy)
  end function idate_ge

  logical function date_le(x,y) result(lt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    lt = (x%days_from_reference == y%days_from_reference .and. &
          x%second_of_day       == y%second_of_day)
    if (lt) return
    lt = date_less(x,y)
  end function date_le

  logical function idate_le(x,y) result(lt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call set_calint(yy,x%calendar)
    lt = date_le(x,yy)
  end function idate_le

  logical function date_ne(x,y) result(ln)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    ln = (x%days_from_reference /= y%days_from_reference .or. &
          x%second_of_day       /= y%second_of_day)
  end function date_ne

  logical function idate_ne(x,y) result(ln)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call set_calint(yy,x%calendar)
    ln = date_ne(x,yy)
  end function idate_ne

  real(rk8) function tohours(x) result(hs)
    implicit none
    type (rcm_time_interval) , intent(in) :: x
    select case (x%iunit)
      case (usec)
        hs = dble(x%ival)/3600.0D0
      case (umin)
        hs = dble(x%ival)/60.0D0
      case (uhrs)
        hs = dble(x%ival)
      case (uday)
        hs = dble(x%ival)*24.0D0
      case default
        call die('mod_date','Interval unit conversion depend on calendar',1)
    end select
  end function tohours

  real(rk8) function yeardayfrac(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    type (iatime) :: t
    call internal_to_date_time(x,d,t)
    yeardayfrac = dble(idayofyear(d)) + dble(t%hour)/24.0D+00 + &
                  dble(t%minute/1440.0D0) + dble(t%second/86400.0D0)
  end function yeardayfrac

  integer function getyear(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    getyear = d%year
  end function getyear

  integer function getmonth(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    getmonth = d%month
  end function getmonth

  integer function getday(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    getday = d%day
  end function getday

  integer function gethour(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iatime) :: t
    call second_of_day_to_time(x,t)
    gethour = t%hour
  end function gethour

  integer function getminute(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iatime) :: t
    call second_of_day_to_time(x,t)
    getminute = t%minute
  end function getminute

  integer function getsecond(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iatime) :: t
    call second_of_day_to_time(x,t)
    getsecond = t%second
  end function getsecond

  subroutine split_rcm_time_and_date(x,iy,im,id,ih)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(out) :: iy , im , id , ih
    type(iadate) :: d
    type(iatime) :: t
    call internal_to_date_time(x,d,t)
    iy = d%year
    im = d%month
    id = d%day
    ih = t%hour
  end subroutine split_rcm_time_and_date

end module mod_date
