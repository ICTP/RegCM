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

  implicit none

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
  integer(ik4) , parameter :: cordex_refdate = 1949120100

  character (len=16) , public , dimension(7) :: cintstr
  character (len=12) , public , dimension(3) :: calstr

  integer(ik4) , dimension(12) :: mlen

  type i4wcal
    integer(ik4) :: i
    integer(ik4) :: cal
  end type i4wcal

  type i8wcal
    integer(ik8) :: i
    integer(ik4) :: cal
  end type i8wcal

  public :: i4wcal , i8wcal

  type rcm_time_and_date
    integer(ik4) :: calendar = gregorian
    integer(ik8) :: days_from_reference = 0
    integer(ik8) :: second_of_day = 0
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
    module procedure initfromintdt ,    &
                     initfromintdtwc ,  &
                     initfromint8dt ,   &
                     initfromint8dtwc , &
                     initfromtypedt
    module procedure initfromintit ,    &
                     initfromsingleit , &
                     initfromdoubleit , &
                     initfromtypeit
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
    module procedure split_i10 , split_i10_8 , split_rcm_time_and_date , &
                    split_rcm_time_and_date_complete , split_rcm_date
  end interface

  interface date_is
    module procedure full_date_is , daymon_is
  end interface date_is

  interface time_is
    module procedure time_complete_is , time_of_day_is
  end interface time_is

  interface timeval2date
    module procedure timeval2date_single , timeval2date_double
  end interface

  interface julianday
    module procedure ymd_julianday
    module procedure idate_julianday
  end interface julianday

  public :: timeval2ym
  public :: rcm_time_and_date , assignment(=) , operator(==)
  public :: rcm_time_interval , operator(+) , operator(-)
  public :: operator(>) , operator(<) , operator(>=) , &
            operator(<=) , operator(/=)
  public :: print_rcm_time_and_date , print_rcm_time_interval
  public :: setcal , set_timeunit
  public :: tochar , toint10 , tochar10 , tohours , toiso8601
  public :: lsamemonth , imondiff , iyeardiff
  public :: lfhomonth , monfirst , monlast , monmiddle
  public :: hourdiff , nextmon , prevmon , yrfirst , nextwk , prevwk
  public :: lsameweek , iwkdiff , idayofweek , ifdoweek , ildoweek
  public :: timeval2date , lfdomonth , lfdoyear , lmidnight , yeardayfrac
  public :: split_idate , julianday , dayofyear
  public :: getyear , getmonth , getday
  public :: gethour , getminute , getsecond
  public :: date_is , time_is
  public :: curr_date , ref_date , curr_time , calendar_str , calendar_int
  public :: date_in_scenario
  public :: ndaypm , yeardays , yearpoint

  data mlen /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data calstr /'gregorian','noleap','360_day'/
  data cintstr /'seconds', 'minutes', 'hours', 'days', &
                'months', 'years', 'centuries'/

  contains

  pure integer(ik4) function ndaypm(y,m,cal)
    implicit none
    integer(ik4) , intent(in) :: y , m , cal
    select case (cal)
      case (gregorian)
        ndaypm = mdays_leap(y, m)
      case (noleap)
        ndaypm = mlen(m)
      case (y360)
        ndaypm = 30
    end select
  end function ndaypm

  pure logical function lleap(iyear)
    implicit none
    integer(ik4) , intent(in) :: iyear
    if ( mod(iyear,400) == 0 .or.  &
        ( mod(iyear,4) == 0 .and. mod(iyear,100) /= 0 ) ) then
      lleap = .true.
    else
      lleap = .false.
    end if
  end function lleap

  pure integer(ik4) function mdays_leap(iyear, imon) result(mdays)
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

  pure integer(ik4) function yeardays(y,c) result(yd)
    implicit none
    integer(ik4) , intent(in) :: y , c
    select case (c)
      case (noleap)
        yd = 365
      case (y360)
        yd = 360
      case default
        yd = 365
        if ( lleap(y) ) yd = 366
    end select
  end function yeardays

  pure subroutine idayofyear_to_monthdate(j,y,c,m,d)
    implicit none
    integer(ik8) , intent(in) :: j
    integer(ik4) , intent(in) :: y , c
    integer(ik4) , intent(out) :: m , d
    integer(ik8) :: id
    integer(ik4) :: md
    id = j + 1
    m = 1
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
    d = int(id,ik4)
  end subroutine idayofyear_to_monthdate

  pure integer(ik4) function idayofyear(x) result(id)
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

  integer(ik4) function dayofyear(x) result(id)
    implicit none
    type(rcm_time_and_date) , intent(in) :: x
    type(iadate) :: d
    type(iatime) :: t
    call internal_to_date_time(x,d,t)
    id = idayofyear(d)
  end function dayofyear

  pure subroutine date_to_days_from_reference(d,x)
    type(iadate) , intent(in) :: d
    type(rcm_time_and_date) , intent(inout) :: x
    integer(ik4) :: ny , iy
    integer(ik8) :: id
    x%calendar = d%calendar
    ny = d%year - reference_year
    iy = reference_year
    id = 0
    if ( ny > 0 ) then
      do while ( ny > 0 )
        id = id + yeardays(iy,d%calendar)
        iy = iy + 1
        ny = ny - 1
      end do
    else
      do while ( ny < 0 )
        id = id - yeardays(iy-1,d%calendar)
        iy = iy - 1
        ny = ny + 1
      end do
    end if
    x%days_from_reference = id + idayofyear(d) - 1
  end subroutine date_to_days_from_reference

  pure subroutine days_from_reference_to_date(x,d)
    type(rcm_time_and_date) , intent(in) :: x
    type(iadate) , intent(out) :: d
    integer(ik4) :: iy
    integer(ik8) :: id
    d%calendar = x%calendar
    id = x%days_from_reference
    d%year = reference_year
    if ( id > 0 ) then
      iy = yeardays(d%year,d%calendar)
      do while ( id >= iy )
        d%year = d%year + 1
        id = id - iy
        iy = yeardays(d%year,d%calendar)
      end do
    else if ( id < 0 ) then
      d%year = d%year - 1
      iy = yeardays(d%year,d%calendar)
      do while ( id < -iy )
        d%year = d%year - 1
        id = id + iy
        iy = yeardays(d%year,d%calendar)
      end do
      if ( id <= 0 ) id = id + iy
    end if
    call idayofyear_to_monthdate(id,d%year,d%calendar,d%month,d%day)
  end subroutine days_from_reference_to_date

  subroutine time_to_second_of_day(t,x)
    type(iatime) , intent(in) :: t
    type(rcm_time_and_date) , intent(inout) :: x
    x%second_of_day = t%hour*3600+t%minute*60+t%second
  end subroutine time_to_second_of_day

  subroutine second_of_day_to_time(x,t)
    type(rcm_time_and_date) , intent(in) :: x
    type(iatime) , intent(out) :: t
    t%hour = int(x%second_of_day/3600,ik4)
    t%minute = int(mod(x%second_of_day,3600_ik8)/60,ik4)
    t%second = int(mod(x%second_of_day,60_ik8),ik4)
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

  subroutine split_i10(idate, iy, im, id, ih)
    implicit none
    integer(ik4) , intent(in) :: idate
    integer(ik4) , intent(out) :: iy , im , id , ih
    integer(ik8) :: base , iidate
    integer(ik8) :: iy8 , im8 , id8 , ih8
    iidate = int(idate,ik8)
    base = iidate
    iy8 = base/1000000_ik8
    base = abs(base-iy8*1000000_ik8)
    im8 = base/10000_ik8
    base = base-im8*10000_ik8
    id8 = base/100_ik8
    base = base-id8*100_ik8
    ih8 = base
    iy = int(iy8,ik4)
    im = int(im8,ik4)
    id = int(id8,ik4)
    ih = int(ih8,ik4)
  end subroutine split_i10

  subroutine split_i10_8(idate, iy, im, id, ih)
    implicit none
    integer(ik8) , intent(in) :: idate
    integer(ik4) , intent(out) :: iy , im , id , ih
    integer(ik8) :: base , iidate
    integer(ik8) :: iy8 , im8 , id8 , ih8
    iidate = idate
    base = iidate
    iy8 = base/1000000_ik8
    base = abs(base-iy8*1000000_ik8)
    im8 = base/10000_ik8
    base = base-im8*10000_ik8
    id8 = base/100_ik8
    base = base-id8*100_ik8
    ih8 = base
    iy = int(iy8,ik4)
    im = int(im8,ik4)
    id = int(id8,ik4)
    ih = int(ih8,ik4)
  end subroutine split_i10_8

  subroutine initfromintdtwc(x, iwc)
    implicit none
    type(i4wcal) , intent(in) :: iwc
    type (rcm_time_and_date) , intent(out) :: x
    type (iadate) :: d
    type (iatime) :: t
    call split_i10(iwc%i, d%year, d%month, d%day, t%hour)
    if ( d%year <= -10000 .or. d%year >= 10000 ) then
      write(stderr,*) 'YEAR  = ', d%year
      write(stderr,*) 'MONTH = ', d%month
      write(stderr,*) 'DAY   = ', d%day
      write(stderr,*) 'HOUR  = ', t%hour
      call die('mod_date','Inconsistent year in date')
    end if
    if ( d%month < 1 .or. d%month > 12 ) then
      write(stderr,*) 'YEAR  = ', d%year
      write(stderr,*) 'MONTH = ', d%month
      write(stderr,*) 'DAY   = ', d%day
      write(stderr,*) 'HOUR  = ', t%hour
      call die('mod_date','Inconsistent month in date')
    end if
    if ( d%day < 1 .or. d%day > 31 ) then
      write(stderr,*) 'YEAR  = ', d%year
      write(stderr,*) 'MONTH = ', d%month
      write(stderr,*) 'DAY   = ', d%day
      write(stderr,*) 'HOUR  = ', t%hour
      call die('mod_date','Inconsistent day in date')
    end if
    if ( t%hour+1 < 1 .or. t%hour+1 > 24 ) then
      write(stderr,*) 'YEAR  = ', d%year
      write(stderr,*) 'MONTH = ', d%month
      write(stderr,*) 'DAY   = ', d%day
      write(stderr,*) 'HOUR  = ', t%hour
      call die('mod_date','Inconsistent hour in date')
    end if
    d%calendar = iwc%cal
    call date_time_to_internal(d,t,x)
  end subroutine initfromintdtwc

  subroutine initfromintdt(x, i)
    implicit none
    integer(ik4) , intent(in) :: i
    type (rcm_time_and_date) , intent(out) :: x
    type (iadate) :: d
    type (iatime) :: t
    call split_i10(i, d%year, d%month, d%day, t%hour)
    if ( d%year <= -10000 .or. d%year >= 10000 ) then
      write(stderr,*) 'YEAR  = ', d%year
      write(stderr,*) 'MONTH = ', d%month
      write(stderr,*) 'DAY   = ', d%day
      write(stderr,*) 'HOUR  = ', t%hour
      call die('mod_date','Inconsistent year in date')
    end if
    if ( d%month < 1 .or. d%month > 12 ) then
      write(stderr,*) 'YEAR  = ', d%year
      write(stderr,*) 'MONTH = ', d%month
      write(stderr,*) 'DAY   = ', d%day
      write(stderr,*) 'HOUR  = ', t%hour
      call die('mod_date','Inconsistent month in date')
    end if
    if ( d%day < 1 .or. d%day > 31 ) then
      write(stderr,*) 'YEAR  = ', d%year
      write(stderr,*) 'MONTH = ', d%month
      write(stderr,*) 'DAY   = ', d%day
      write(stderr,*) 'HOUR  = ', t%hour
      call die('mod_date','Inconsistent day in date')
    end if
    if ( t%hour+1 < 1 .or. t%hour+1 > 24 ) then
      write(stderr,*) 'YEAR  = ', d%year
      write(stderr,*) 'MONTH = ', d%month
      write(stderr,*) 'DAY   = ', d%day
      write(stderr,*) 'HOUR  = ', t%hour
      call die('mod_date','Inconsistent hour in date')
    end if
    d%calendar = gregorian
    call date_time_to_internal(d,t,x)
  end subroutine initfromintdt

  subroutine initfromint8dt(x, i)
    implicit none
    integer(ik8) , intent(in) :: i
    type (rcm_time_and_date) , intent(out) :: x
    type (iadate) :: d
    type (iatime) :: t
    call split_i10_8(i, d%year, d%month, d%day, t%hour)
    if ( d%year <= -10000 .or. d%year >= 10000 ) then
      call die('mod_date','Inconsistent year in date')
    end if
    if ( d%month < 1 .or. d%month > 12 ) then
      call die('mod_date','Inconsistent month in date')
    end if
    if ( d%day < 1 .or. d%day > 31 ) then
      call die('mod_date','Inconsistent day in date')
    end if
    if ( t%hour+1 < 1 .or. t%hour+1 > 24 ) then
      call die('mod_date','Inconsistent hour in date')
    end if
    d%calendar = gregorian
    call date_time_to_internal(d,t,x)
  end subroutine initfromint8dt

  subroutine initfromint8dtwc(x, iwc)
    implicit none
    type(i8wcal) , intent(in) :: iwc
    type (rcm_time_and_date) , intent(out) :: x
    type (iadate) :: d
    type (iatime) :: t
    call split_i10_8(iwc%i, d%year, d%month, d%day, t%hour)
    if ( d%year <= -10000 .or. d%year >= 10000 ) then
      call die('mod_date','Inconsistent year in date')
    end if
    if ( d%month < 1 .or. d%month > 12 ) then
      call die('mod_date','Inconsistent month in date')
    end if
    if ( d%day < 1 .or. d%day > 31 ) then
      call die('mod_date','Inconsistent day in date')
    end if
    if ( t%hour+1 < 1 .or. t%hour+1 > 24 ) then
      call die('mod_date','Inconsistent hour in date')
    end if
    d%calendar = iwc%cal
    call date_time_to_internal(d,t,x)
  end subroutine initfromint8dtwc

  subroutine initfromintit(x, i)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    integer(ik4) , intent(in) :: i
    x%ival = i
    x%iunit = usec
  end subroutine initfromintit

  subroutine initfromsingleit(x, d)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    real(rk4) , intent(in) :: d
    x%ival = int(d,ik8)
    x%iunit = usec
  end subroutine initfromsingleit

  subroutine initfromdoubleit(x, d)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    real(rk8) , intent(in) :: d
    x%ival = int(d,ik8)
    x%iunit = usec
  end subroutine initfromdoubleit

  subroutine set_caltype(x, y)
    implicit none
    type (rcm_time_and_date) , intent(inout) :: x
    type (rcm_time_and_date) , intent(in) :: y
    if ( y%calendar /= x%calendar) then
      call set_calint(x,y%calendar)
    end if
  end subroutine set_caltype

  pure subroutine set_calint(x, c)
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
    if ( c == 'gregorian' .or. &
         c == 'standard' ) then
      ic = gregorian
    else if ( c == 'noleap' .or.   &
              c == 'days_365' .or. &
              c == 'day_365' .or.  &
              c == '365_days' .or. &
              c == '365_day' ) then
      ic = noleap
    else if ( c == 'days_360' .or. &
              c == 'day_360' .or.  &
              c == '360_days' .or. &
              c == '360_day' ) then
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
    if ( d%year >= 0 ) then
      write (cdat, &
        '(" ",i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2," UTC")') &
         d%year, d%month, d%day, t%hour, t%minute, t%second
    else
      write (cdat, &
        '("-",i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2," UTC")') &
         -d%year, d%month, d%day, t%hour, t%minute, t%second
    end if
  end function tochar

  function toiso8601(x) result(cdat)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    character (len=21) :: cdat ! Accomodate for - sign
    type (iadate) :: d
    type (iatime) :: t
    call internal_to_date_time(x,d,t)
    if ( d%year >= 0 ) then
      write (cdat, &
         '(" ",i0.4,"-",i0.2,"-",i0.2,"T",i0.2,":",i0.2,":",i0.2,"Z")') &
         d%year, d%month, d%day, t%hour, t%minute, t%second
    else
      write (cdat, &
         '("-",i0.4,"-",i0.2,"-",i0.2,"T",i0.2,":",i0.2,":",i0.2,"Z")') &
         -d%year, d%month, d%day, t%hour, t%minute, t%second
    end if
  end function toiso8601

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

  pure logical function isequalit(x, y)
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
        ! It is to have 14 as February middle month for Gregorian.
        select case (z%calendar)
          case (gregorian)
            dm = real(d%day-1,rk8)/real(mdays_leap(d%year, d%month),rk8)
          case (noleap)
            dm = real(d%day-1,rk8)/real(mlen(d%month),rk8)
          case default
            dm = real(d%day-1,rk8)/30.0_rk8
        end select
        d%month = d%month+int(mod(tmp,i8mpy),ik4)
        call adjustpm(d%month,d%year,12)
        tmp = tmp/i8mpy
        d%year = d%year+int(tmp,ik4)
        ! Adjust date of the month: This is really a trick...
        ! It is to have 14 as February middle month for Gregorian.
        select case (z%calendar)
          case (gregorian)
            d%day = nint(real(mdays_leap(d%year, d%month),rk8) * dm)+1
          case (noleap)
            d%day = nint(real(mlen(d%month),rk8) * dm)+1
          case default
            d%day = nint(30.0_rk8*dm)+1
        end select
        call date_to_days_from_reference(d,z)
      case (uyrs)
        call days_from_reference_to_date(x,d)
        d%year = d%year+int(tmp,ik4)
        call date_to_days_from_reference(d,z)
      case (ucnt)
        call days_from_reference_to_date(x,d)
        d%year = d%year+int(i8ypc*tmp,ik4)
        call date_to_days_from_reference(d,z)
    end select
  end function add_interval

  real(rkx) function idate_julianday(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) :: iy , im , id
    call split_rcm_date(x,iy,im,id)
    idate_julianday = real(ymd_julianday(iy, im, id),rkx)
  end function idate_julianday

  real(rk8) function ymd_julianday(iy, im, id)
    implicit none
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
    ymd_julianday = int(365.25_rk8*real(iiy+4716,rk8)) + &
                int(30.6001_rk8*real(iim+1,rk8))   + &
                real(id + ib,rk8) - 1524.5_rk8
    contains
      logical function lcaltype(iy, im, id)
      implicit none
      integer(ik4) , intent(in) :: iy , im , id
      integer(ik4) :: icaltype
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
      lcaltype = .false.
      if (icaltype == 2) then
        lcaltype = .true.
      else
        write (stderr, *) 'year  = ', iy
        write (stderr, *) 'month = ', im
        write (stderr, *) 'day   = ', id
        call die('mod_date','Day non existent, inside Julian/Gregorian jump')
      end if
    end function lcaltype
  end function ymd_julianday

  subroutine check_cal(x,y)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) , intent(in) :: y
    if ( x%calendar /= y%calendar ) then
      write (stderr,*) 'X calendar = ', calstr(x%calendar)
      write (stderr,*) 'Y calendar = ', calstr(y%calendar)
      call die('mod_date','Dates not comparable: not using same calendar.')
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
        if ( z%second_of_day > 0 ) then
          z%second_of_day = i8spd-z%second_of_day
          tmp = tmp+i8spd
        end if
        tmp = tmp/i8spd
        z%days_from_reference = z%days_from_reference - tmp
      case (uhrs)
        tmp = tmp*i8sph - z%second_of_day
        z%second_of_day = mod(tmp, i8spd)
        if ( z%second_of_day > 0 ) then
          z%second_of_day = i8spd-z%second_of_day
          tmp = tmp+i8spd
        end if
        tmp = tmp/i8spd
        z%days_from_reference = z%days_from_reference - tmp
      case (uday)
        z%days_from_reference = z%days_from_reference - tmp
      case (umnt)
        call days_from_reference_to_date(x,d)
        d%month = d%month-int(mod(tmp,i8mpy),ik4)
        call adjustmp(d%month,d%year,12)
        d%year = d%year-int(tmp/i8mpy,ik4)
        call date_to_days_from_reference(d,z)
      case (uyrs)
        call days_from_reference_to_date(x,d)
        d%year = d%year-int(tmp,ik4)
        call date_to_days_from_reference(d,z)
      case (ucnt)
        call days_from_reference_to_date(x,d)
        d%year = d%year-int(i8ypc*tmp,ik4)
        call date_to_days_from_reference(d,z)
    end select
  end function sub_interval

  function tochar10(x) result(cdat)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    character (len=11) :: cdat
    integer(ik8) :: ival
    ival = toint10(x)
    if ( ival > 0 ) then
      write(cdat,'(i0.10)') ival
    else
      write(cdat,'(i0.10)') ival
    end if
    cdat = adjustl(cdat)
  end function tochar10

  integer(ik8) function toint10(x) result(z)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    type (iatime) :: t
    call internal_to_date_time(x,d,t)
    if ( d%year >= 0 ) then
      z = d%year*1000000_ik8+d%month*10000+d%day*100+t%hour
    else
      z = -((-d%year)*1000000_ik8+d%month*10000+d%day*100+t%hour)
    end if
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

  integer(ik4) function iyeardiff(x,y) result(z)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    type (iadate) :: d1 , d2
    call check_cal(x,y)
    call days_from_reference_to_date(x,d1)
    call days_from_reference_to_date(y,d2)
    z = abs((d1%year-d2%year))
  end function iyeardiff

  integer(ik4) function imondiff(x,y) result(z)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    type (iadate) :: d1 , d2
    call check_cal(x,y)
    call days_from_reference_to_date(x,d1)
    call days_from_reference_to_date(y,d2)
    z = abs((d1%year-d2%year))*12+(d1%month-d2%month)
  end function imondiff

  real(rk8) function hourdiff(x,y) result(z)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    z = real(x%days_from_reference-y%days_from_reference,rk8)*24.0_rk8 + &
        real(x%second_of_day-y%second_of_day,rk8)/3600.0_rk8
  end function hourdiff

  logical function lfhomonth(x) result(lf)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mm
    mm = monmiddle(x)
    lf = (x < mm)
  end function lfhomonth

  integer(ik4) function idayofweek(x) result(iday)
    ! Sun Mon Tue Wed Thu Fri Sat
    ! 1   2   3   4   5   6   7
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    real(rk8) :: jd
    if (x%calendar /= gregorian) then
      call die('mod_date', &
               'Error: week concept works only for gregorian calendar.')
    end if
    call days_from_reference_to_date(x,d)
    jd = julianday(d%year, d%month, d%day)
    iday = int(mod(jd+1.5_rk8, 7.0_rk8))+1
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

  integer(ik4) function iwkdiff(x,y) result(iwk)
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

  subroutine timeval2ym(xval,cunit,year,month)
    implicit none
    real(rk8) , intent(in) :: xval
    character(*) , intent(in) :: cunit
    integer(ik4) , intent(out) :: year , month
    integer(ik4) :: istat , aunit
    type (iadate) , save :: d
    type (iatime) :: t

    character(len=64) , save :: csave
    data csave /'months since XXXX-XX-XX XX:XX:XX XXX'/

    if (csave == cunit) then
      year = d%year+int(xval/12.0_rk8)
      month = d%month+int(mod(xval,12.0_rk8))
      if ( month > 12 ) then
        month = month - 12
        year = year + 1
      else if ( month < 1 ) then
        month = month + 12
        year = year - 1
      end if
      return
    end if
    istat = parse_timestring(cunit,aunit,d%year,d%month,d%day, &
                             t%hour,t%minute,t%second)
    if ( istat /= 0 .or. aunit /= umnt ) then
      call die('mod_date','TIME UNIT ERROR IN timeval2ym')
    end if
    csave = cunit
    year = d%year+int(xval/12.0_rk8)
    month = d%month+int(mod(xval,12.0_rk8))
    if ( month > 12 ) then
      month = month - 12
      year = year + 1
    else if ( month < 1 ) then
      month = month + 12
      year = year - 1
    end if
  end subroutine timeval2ym

  function timeval2date_single(xval,cunit,ccal) result(dd)
    implicit none
    real(rk4) , intent(in) :: xval
    character(*) , intent(in) :: cunit
    character(*) , intent(in) :: ccal

    type (rcm_time_and_date) :: dd
    type (rcm_time_interval) :: z , zz
    character(len=64) , save :: csave
    character(len=16) , save :: csavecal
    integer(ik4) , save :: iunit
    type (rcm_time_and_date) , save :: dref
    type (iadate) :: d
    type (iatime) :: t
    integer(ik4) :: istat

    data csave/'none'/
    data csavecal/'none'/

    t%hour = 0
    t%minute = 0
    t%second = 0

    if ( csave == cunit .and. csavecal == ccal ) then
      z = abs(xval)
      z%iunit = iunit
      if ( xval >= 0.0_rk8 ) then
        dd = dref + z
      else
        dd = dref - z
      end if
      ! Some datasets have fraction of days as hours.
      ! I.e. 15.5 days from a date means 15 days + 12 hours.
      ! Never seen fraction of minutes or hours.
      ! Months fraction are almost nonsense.
      if (iunit == uday) then
        zz%ival = (int((abs(xval)-abs(real(z%ival,rk8)))*24.0_rk8))
        zz%iunit = uhrs
        if ( zz%ival /= 0 ) then
          if ( xval >= 0.0_rk8 ) then
            dd = dd + zz
          else
            dd = dd - zz
          end if
        end if
      end if
    else
      csave = cunit
      csavecal = ccal
      z = abs(xval)
      istat = parse_timestring(cunit,iunit,d%year,d%month,d%day, &
                               t%hour,t%minute,t%second)
      if ( istat /= 0 ) then
        call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
      end if
      if ( csavecal(1:6) == 'noleap' .or.   &
           csavecal(1:8) == 'days_365' .or. &
           csavecal(1:7) == '365_day' ) then
        d%calendar = noleap
      else if (csavecal(1:7) == '360_day') then
        d%calendar = y360
      else
        d%calendar = gregorian
      end if
      call date_time_to_internal(d,t,dref)
      z%iunit = iunit
      if ( xval >= 0.0_rk8 ) then
        dd = dref + z
      else
        dd = dref - z
      end if
      ! Some datasets have fraction of days as hours.
      ! I.e. 15.5 days from a date means 15 days + 12 hours.
      ! Never seen fraction of minutes or hours.
      ! Months fraction are almost nonsense.
      if (iunit == uday) then
        zz%ival = int((abs(xval)-abs(real(z%ival,rk8)))*24.0_rk8)
        zz%iunit = uhrs
        if ( zz%ival > 0 ) then
          if ( xval > 0 ) then
            dd = dd + zz
          else if ( xval < 0 ) then
            dd = dd -zz
          end if
        end if
      end if
    end if
  end function timeval2date_single

  function timeval2date_double(xval,cunit,ccal) result(dd)
    implicit none
    real(rk8) , intent(in) :: xval
    character(*) , intent(in) :: cunit
    character(*) , intent(in) :: ccal

    type (rcm_time_and_date) :: dd
    type (rcm_time_interval) :: z , zz
    character(len=64) , save :: csave
    character(len=16) , save :: csavecal
    integer(ik4) , save :: iunit
    type (rcm_time_and_date) , save :: dref
    type (iadate) :: d
    type (iatime) :: t
    integer(ik4) :: istat

    data csave/'none'/
    data csavecal/'none'/

    t%hour = 0
    t%minute = 0
    t%second = 0

    if ( csave == cunit .and. csavecal == ccal ) then
      z = abs(xval)
      z%iunit = iunit
      if ( xval >= 0.0_rk8 ) then
        dd = dref + z
      else
        dd = dref - z
      end if
      ! Some datasets have fraction of days as hours.
      ! I.e. 15.5 days from a date means 15 days + 12 hours.
      ! Never seen fraction of minutes or hours.
      ! Months fraction are almost nonsense.
      if (iunit == uday) then
        zz%ival = (int((abs(xval)-abs(real(z%ival,rk8)))*24.0_rk8))
        zz%iunit = uhrs
        if ( zz%ival /= 0 ) then
          if ( xval >= 0.0_rk8 ) then
            dd = dd + zz
          else
            dd = dd - zz
          end if
        end if
      end if
    else
      csave = cunit
      csavecal = ccal
      z = abs(xval)
      istat = parse_timestring(cunit,iunit,d%year,d%month,d%day, &
                               t%hour,t%minute,t%second)
      if ( istat /= 0 ) then
        call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
      end if
      if ( csavecal(1:6) == 'noleap' .or.   &
           csavecal(1:8) == 'days_365' .or. &
           csavecal(1:7) == '365_day' ) then
        d%calendar = noleap
      else if (csavecal(1:7) == '360_day') then
        d%calendar = y360
      else
        d%calendar = gregorian
      end if
      call date_time_to_internal(d,t,dref)
      z%iunit = iunit
      if ( xval >= 0.0_rk8 ) then
        dd = dref + z
      else
        dd = dref - z
      end if
      ! Some datasets have fraction of days as hours.
      ! I.e. 15.5 days from a date means 15 days + 12 hours.
      ! Never seen fraction of minutes or hours.
      ! Months fraction are almost nonsense.
      if (iunit == uday) then
        zz%ival = int((abs(xval)-abs(real(z%ival,rk8)))*24.0_rk8)
        zz%iunit = uhrs
        if ( zz%ival > 0 ) then
          if ( xval > 0 ) then
            dd = dd + zz
          else if ( xval < 0 ) then
            dd = dd -zz
          end if
        end if
      end if
    end if
  end function timeval2date_double
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
        hs = real(x%ival,rk8)/3600.0_rk8
      case (umin)
        hs = real(x%ival,rk8)/60.0_rk8
      case (uhrs)
        hs = real(x%ival,rk8)
      case (uday)
        hs = real(x%ival,rk8)*24.0_rk8
      case default
        hs = 0.0_rk8
        call die('mod_date','Interval unit conversion depend on calendar')
    end select
  end function tohours

  real(rk8) function yearpoint(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    type (iatime) :: t
    integer(ik4) :: id , i
    real(rk8) :: lc
    call internal_to_date_time(x,d,t)
    id = d%day
    select case (x%calendar)
      case (gregorian)
        do i = 1 , d%month-1
          id = id + mdays_leap(d%year, i)
        end do
        yearpoint = real(id,rk8)
        if ( lleap(d%year) ) then
          lc = -(yearpoint+1095_rk8)/1461_rk8
        else if ( lleap(d%year+1) ) then
          lc = -(yearpoint+730_rk8)/1461_rk8
        else if ( lleap(d%year+2) ) then
          lc = -(yearpoint+365_rk8)/1461_rk8
        else
          lc = -yearpoint/1461_rk8
        end if
        yearpoint = yearpoint+lc
      case (noleap)
        id = id + sum(mlen(1:d%month-1))
        yearpoint = real(id,rk8)
      case (y360)
        id = id + 30*(d%month-1)
        yearpoint = real(id,rk8)
    end select
  end function yearpoint

  real(rk8) function yeardayfrac(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    type (iatime) :: t
    call internal_to_date_time(x,d,t)
    yeardayfrac = real(idayofyear(d),rk8) + &
                  real(t%hour,rk8)/24.0_rk8 + &
                  real(t%minute/1440.0_rk8,rk8) + &
                  real(t%second/86400.0_rk8,rk8) - 1.0_rk8
  end function yeardayfrac

  integer(ik4) function getyear(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    getyear = d%year
  end function getyear

  integer(ik4) function getmonth(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    getmonth = d%month
  end function getmonth

  integer(ik4) function getday(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iadate) :: d
    call days_from_reference_to_date(x,d)
    getday = d%day
  end function getday

  integer(ik4) function gethour(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iatime) :: t
    call second_of_day_to_time(x,t)
    gethour = t%hour
  end function gethour

  integer(ik4) function getminute(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iatime) :: t
    call second_of_day_to_time(x,t)
    getminute = t%minute
  end function getminute

  integer(ik4) function getsecond(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (iatime) :: t
    call second_of_day_to_time(x,t)
    getsecond = t%second
  end function getsecond

  subroutine ref_date(iy,im,id,isec)
    implicit none
    integer(ik4) , intent(out) :: iy , im , id , isec
    iy = 1949
    im = 12
    id = 1
    isec = 0
  end subroutine ref_date

  subroutine curr_date(x,iy,im,id,isec,offset)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(out) :: iy , im , id , isec
    integer(ik4) , intent(in) , optional :: offset
    type (rcm_time_interval) :: tdif
    type (rcm_time_and_date) :: y
    type(iadate) :: d
    if ( present(offset) ) then
      tdif = offset
      y = x + tdif
    else
      y = x
    end if
    call days_from_reference_to_date(y,d)
    iy = d%year
    im = d%month
    id = d%day
    isec = int(y%second_of_day,ik4)
  end subroutine curr_date

  ! Returns days and seconds from cordex reference date
  subroutine curr_time(x,iday,isec)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(out) :: iday , isec
    type (rcm_time_and_date) , save :: y
    logical , save :: isinit = .false.
    if ( .not. isinit ) then
      y = cordex_refdate
      call set_caltype(y,x)
      isinit = .true.
    else
      if ( x%calendar /= y%calendar ) then
        y = cordex_refdate
        call set_caltype(y,x)
      end if
    end if
    iday = int(x%days_from_reference - y%days_from_reference,ik4)
    isec = int(x%second_of_day,ik4)
  end subroutine curr_time

  subroutine split_rcm_date(x,iy,im,id)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(out) :: iy , im , id
    type(iadate) :: d
    call days_from_reference_to_date(x,d)
    iy = d%year
    im = d%month
    id = d%day
  end subroutine split_rcm_date

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

  subroutine split_rcm_time_and_date_complete(x,iy,im,id,ih,imm,iss)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(out) :: iy , im , id , ih , imm , iss
    type(iadate) :: d
    type(iatime) :: t
    call internal_to_date_time(x,d,t)
    iy = d%year
    im = d%month
    id = d%day
    ih = t%hour
    imm = t%minute
    iss = t%second
  end subroutine split_rcm_time_and_date_complete

  logical function full_date_is(x,y,m,d)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: y , m , d
    integer(ik4) :: iy , im , id
    call split_rcm_date(x,iy,im,id)
    full_date_is = ( y == iy .and. m == im .and. d == id )
  end function full_date_is

  logical function daymon_is(x,m,d)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: m , d
    integer(ik4) :: iy , im , id
    call split_rcm_date(x,iy,im,id)
    daymon_is = ( m == im .and. d == id )
  end function daymon_is

  logical function time_complete_is(x,h,m,s,delta)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: h , m , s
    real(rkx) , optional , intent(in) :: delta
    type(iatime) :: t
    type (rcm_time_and_date) :: y
    integer(ik8) :: id
    t%hour = h
    t%minute = m
    t%second = s
    call time_to_second_of_day(t,y)
    if ( present(delta) ) then
      id = int(delta,ik8)
      time_complete_is = ( abs(y%second_of_day - x%second_of_day) < id )
    else
      time_complete_is = ( y%second_of_day == x%second_of_day )
    end if
  end function time_complete_is

  logical function time_of_day_is(x,s,delta)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: s
    real(rkx) , optional , intent(in) :: delta
    integer(ik8) :: id
    if ( present(delta) ) then
      id = int(delta,ik8)
      time_of_day_is = ( abs(s - x%second_of_day) < id )
    else
      time_of_day_is = ( s == x%second_of_day )
    end if
  end function time_of_day_is

  character(len=16) function calendar_str(x) result(calstr)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    select case (x%calendar)
      case(gregorian)
        calstr = 'gregorian'
      case(noleap)
        calstr = 'noleap'
      case(y360)
        calstr = '360_day'
      case default
        call die('mod_date','Unknown Calendar!!!')
    end select
  end function calendar_str

  pure integer function calendar_int(calstr) result(calint)
    implicit none
    character(len=*) , intent(in) :: calstr
    select case (calstr)
      case('gregorian', 'standard')
        calint = gregorian
      case('noleap', '365_days', '365_day', 'days_365', 'day_365')
        calint = noleap
      case('days_360', 'day_360', '360_days', '360_day')
        calint = y360
      case default
        calint = gregorian
    end select
  end function calendar_int

  logical function date_in_scenario(x,cmip_cycle,l2006)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    integer(ik4) , intent(in) :: cmip_cycle
    logical , optional , intent(in) :: l2006
    integer(ik4) :: iy , im , id
    call split_rcm_date(x,iy,im,id)
    if ( cmip_cycle == 3 ) then
      date_in_scenario = ( iy > 2000 )
    else ! Assume CMIP5
      if ( present(l2006) ) then
        if ( l2006 ) then
          date_in_scenario = ( iy > 2005 )
        else
          date_in_scenario = ( iy > 2005 .or. ( iy == 2005 .and. im == 12 ) )
        end if
      else
        date_in_scenario = ( iy > 2005 .or. ( iy == 2005 .and. im == 12 ) )
      end if
    end if
  end function date_in_scenario

  integer(ik4) function parse_timestring(string,aunit, &
                                         y,m,d,h,minu,seco) result(ires)
    implicit none
    character(len=*) :: string
    integer(ik4) , intent(out) :: aunit
    integer(ik4) , intent(out) :: y , m , d , h , minu , seco
    integer(ik4) :: ip , istat
    character(len=128) :: cdum , cdum1
    character(len=1) :: cdiv

    ires = -1
    aunit = -1

    cdum(:) = ' '
    cdum1(:) = ' '

    y = 1900
    m = 1
    d = 1
    h = 0
    minu = 0
    seco = 0

    if ( string(1:5) == 'hours' ) then
      aunit = uhrs
      ip = 13
    else if ( string(1:6) == 'months' ) then
      aunit = umnt
      ip = 14
    else if ( string(1:4) == 'days' ) then
      aunit = uday
      ip = 12
    else if ( string(1:7) == 'minutes' ) then
      aunit = umin
      ip = 15
    else if ( string(1:7) == 'seconds' ) then
      aunit = usec
      ip = 15
    else if ( string(1:5) == 'years' ) then
      aunit = uyrs
      ip = 13
    else
      write(stderr,*) 'Unrecognized unit of measure for time string.'
      write(stderr,*) 'Offending timeunit string :',trim(string)
    end if

    if ( len_trim(string(ip:)) >= 6 ) then
      read(string(ip:),'(i4,a1,a)',iostat=istat) y , cdiv , cdum
      if ( istat /= 0 ) then
        write(stderr,*) 'Unrecognized year read from time string.'
        write(stderr,*) 'Offending timeunit string :',trim(string)
        return
      end if
      if ( cdiv /= '-' ) then
        write(stderr,*) 'Unrecognized year read from time string.'
        write(stderr,*) 'Offending timeunit string :',trim(string)
        return
      end if
    else
      read(string(ip:),'(i4)',iostat=istat) y
      if ( istat /= 0 ) then
        write(stderr,*) 'Unrecognized year read from time string.'
        write(stderr,*) 'Offending timeunit string :',trim(string)
        return
      end if
      ires = 0
      return
    end if

    if ( len_trim(cdum) > 1 ) then
      if ( cdum(2:2) == '-' ) then
        read(cdum,'(i1,a1,a)',iostat=istat) m , cdiv , cdum1
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized month read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      else if ( len_trim(cdum) > 2 .and. cdum(3:3) == '-' ) then
        read(cdum,'(i2,a1,a)',iostat=istat) m , cdiv , cdum1
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized month read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      else
        read(cdum,'(i2)',iostat=istat) m
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized month read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
        ires = 0
        return
      end if
    else
      read(cdum,'(i1)',iostat=istat) m
      if ( istat /= 0 ) then
        write(stderr,*) 'Unrecognized month read from time string.'
        write(stderr,*) 'Offending timeunit string :',trim(string)
        return
      end if
      ires = 0
      return
    end if

    cdum = cdum1
    cdum1(:) = ' '

    if ( len_trim(cdum) > 1 ) then
      if ( cdum(2:2) == ' ' .or. cdum(2:2) == 'T' ) then
        read(cdum,'(i1,a1,a)',iostat=istat) d , cdiv , cdum1
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized day read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      else if ( len_trim(cdum) > 2 .and. &
                cdum(3:3) == ' ' .or. cdum(3:3) == 'T' ) then
        read(cdum,'(i2,a1,a)',iostat=istat) d , cdiv , cdum1
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized day read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      else
        read(cdum,'(i2)',iostat=istat) d
        if ( istat /= 0 ) then
          read(cdum,'(i1)',iostat=istat) d
          if ( istat /= 0 ) then
            write(stderr,*) 'Unrecognized day read from time string.'
            write(stderr,*) 'Offending timeunit string :',trim(string)
            return
          end if
        end if
        ires = 0
        return
      end if
    else
      read(cdum,'(i1)',iostat=istat) d
      if ( istat /= 0 ) then
        write(stderr,*) 'Unrecognized day read from time string.'
        write(stderr,*) 'Offending timeunit string :',trim(string)
        return
      end if
      ires = 0
      return
    end if

    cdum = cdum1
    cdum1(:) = ' '

    if ( len_trim(cdum) > 1 ) then
      if ( cdum(2:2) == ':' ) then
        read(cdum,'(i1,a1,a)',iostat=istat) h , cdiv , cdum1
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized hour read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      else if ( len_trim(cdum) > 2 .and. cdum(3:3) == '-' ) then
        read(cdum,'(i2,a1,a)',iostat=istat) h , cdiv , cdum1
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized hour read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      else
        read(cdum,'(i2)',iostat=istat) h
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized hour read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
        ires = 0
        return
      end if
    else
      read(cdum,'(i1)',iostat=istat) h
      if ( istat /= 0 ) then
        write(stderr,*) 'Unrecognized hour read from time string.'
        write(stderr,*) 'Offending timeunit string :',trim(string)
        return
      end if
      ires = 0
      return
    end if

    cdum = cdum1
    cdum1(:) = ' '

    if ( len_trim(cdum) > 1 ) then
      if ( cdum(2:2) == ':' ) then
        read(cdum,'(i1,a1,a)',iostat=istat) minu , cdiv , cdum1
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized minute read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      else if ( len_trim(cdum) > 2 .and. cdum(3:3) == '-' ) then
        read(cdum,'(i2,a1,a)',iostat=istat) minu , cdiv , cdum1
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized minute read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      else
        read(cdum,'(i2)',iostat=istat) minu
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized minute read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
        ires = 0
        return
      end if
    else
      read(cdum,'(i1)',iostat=istat) minu
      if ( istat /= 0 ) then
        write(stderr,*) 'Unrecognized minute read from time string.'
        write(stderr,*) 'Offending timeunit string :',trim(string)
        return
      end if
      ires = 0
      return
    end if

    if ( len_trim(cdum1) > 1 ) then
      read(cdum1,'(i2)',iostat=istat) seco
      if ( istat /= 0 ) then
        read(cdum1,'(i1)', iostat=istat) seco
        if ( istat /= 0 ) then
          write(stderr,*) 'Unrecognized second read from time string.'
          write(stderr,*) 'Offending timeunit string :',trim(string)
          return
        end if
      end if
    end if
    ires = 0
  end function parse_timestring

end module mod_date
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
