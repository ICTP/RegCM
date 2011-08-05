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

  use m_realkinds
  use m_stdio
  use m_die
  use mpi

  private

  integer , public , parameter :: gregorian = 1
  integer , public , parameter :: noleap    = 2
  integer , public , parameter :: y360      = 3

  integer , public , parameter :: usec = 1
  integer , public , parameter :: umin = 2
  integer , public , parameter :: uhrs = 3
  integer , public , parameter :: uday = 4
  integer , public , parameter :: umnt = 5
  integer , public , parameter :: uyrs = 6
  integer , public , parameter :: ucnt = 7

  integer , parameter :: reference_year = 1900

  character (len=16) , public , dimension(7) :: cintstr
  character (len=12) , public , dimension(3) :: calstr

  integer , dimension(12) :: mlen

  type rcm_time_and_date
    integer :: calendar = gregorian
    integer :: year = reference_year
    integer :: month = 1
    integer :: day = 1
    integer :: hour = 0
    integer :: minute = 0
    integer :: second = 0
    integer :: days_from_reference = 0
    integer :: second_of_day = 0
    contains
      procedure , pass(x) :: printdate => print_rcm_time_and_date
      procedure , pass(x) :: toidate => toint
      procedure , pass(x) :: tostring => tochar
      procedure , pass(x) :: setup => internal_setup
      procedure , pass(x) :: recalculate => recalculate_from_internal
      procedure , pass(x) :: setcal => set_calint
      procedure , pass(x) :: setccal => set_calstring
      procedure , pass(x) :: broadcast => date_bcast
  end type rcm_time_and_date

  type rcm_time_interval
    integer :: ival = 0
    integer :: iunit = usec
    contains
      procedure , pass(x) :: printint => print_rcm_time_interval
      procedure , pass(x) :: setunit => set_timeunit
      procedure , pass(x) :: hours => tohours
  end type rcm_time_interval

  interface assignment(=) 
    module procedure initfromintdt , initfromtypedt
    module procedure initfromintit , initfromtypeit
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

  public :: rcm_time_and_date , assignment(=) , operator(==)
  public :: rcm_time_interval , operator(+) , operator(-)
  public :: operator(>) , operator(<) , operator(>=) , &
            operator(<=) , operator(/=)
  public :: lsamemonth , imondiff , lfhomonth , monfirst , monlast , monmiddle
  public :: nextmon , prevmon , yrfirst , nextwk , prevwk
  public :: lsameweek , iwkdiff , idayofweek , ifdoweek , ildoweek , idayofyear
  public :: timeval2date , lfdomonth , lfdoyear , lmidnight , yeardayfrac
  public :: julianday

  data mlen /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data calstr /'gregorian','noleap','360_day'/
  data cintstr /'seconds', 'minutes', 'hours', 'days', &
                'months', 'years', 'centuries'/

  contains

  function lleap(iyear)
    implicit none
    logical :: lleap
    integer , intent(in) :: iyear
    if ( mod(iyear,400) == 0 .or.  &
        ( mod(iyear,4) == 0 .and. mod(iyear,100) /= 0 ) ) then
      lleap = .true.
    else
      lleap = .false.
    end if
  end function lleap

  integer function mdays_leap(iyear, imon) result(mdays)
    implicit none
    integer , intent(in) :: iyear , imon
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
    integer , intent(in) :: y , c
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
    integer , intent(in) :: j , y , c
    integer , intent(out) :: m , d
    integer :: id , md
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
    class (rcm_time_and_date) , intent(in) :: x
    integer :: i
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

  subroutine days_from_reference(x)
    class(rcm_time_and_date) , intent(inout) :: x
    integer :: ny , ipm , id , iy
    ny = x%year - reference_year
    ipm = isign(1,ny)
    ny = ny*ipm
    iy = reference_year
    id = 0
    do while (ny > 0)
      id = id + ipm*yeardays(iy,x%calendar)
      iy = iy + ipm
      ny = ny - 1
    end do
    x%days_from_reference = id + ipm*idayofyear(x) - ipm
  end subroutine days_from_reference

  subroutine days_from_reference_to_date(x)
    class(rcm_time_and_date) , intent(inout) :: x
    integer :: id , ipm , iy
    id = x%days_from_reference
    ipm = isign(1,id)
    x%year = reference_year
    iy = yeardays(x%year,x%calendar)
    do while ((id*ipm) >= iy)
      x%year = x%year + ipm
      id = id - ipm*iy
      iy = yeardays(x%year,x%calendar)
    end do
    if (id >= 0) then
      call idayofyear_to_monthdate(id+1,x%year,x%calendar,x%month,x%day)
    else
      x%year = x%year - 1
      id = yeardays(x%year,x%calendar)+id
      call idayofyear_to_monthdate(id,x%year,x%calendar,x%month,x%day)
    end if
  end subroutine days_from_reference_to_date

  subroutine seconds_from_midnight(x)
    class(rcm_time_and_date) , intent(inout) :: x
    x%second_of_day = x%hour*3600+x%minute*60+x%second
  end subroutine seconds_from_midnight

  subroutine seconds_from_midnight_to_time(x)
    class(rcm_time_and_date) , intent(inout) :: x
    integer :: i1 , i2
    i1 = x%second_of_day
    i2 = i1/3600
    x%hour = i2
    i1 = i1-i2*3600
    i2 = i1/60
    x%minute = i2
    i1 = i1-i2*60
    x%second = i1
  end subroutine seconds_from_midnight_to_time

  subroutine internal_setup(x)
    class(rcm_time_and_date) , intent(inout) :: x
    call days_from_reference(x)
    call seconds_from_midnight(x)
  end subroutine internal_setup

  subroutine recalculate_from_internal(x)
    class(rcm_time_and_date) , intent(inout) :: x
    call days_from_reference_to_date(x)
    call seconds_from_midnight_to_time(x)
  end subroutine recalculate_from_internal

  subroutine adjustpm(a,b,i)
    implicit none
    integer , intent(inout) :: a , b
    integer , intent(in) :: i
    if ( a > i ) then
      a = a - i
      b = b + 1
    end if
  end subroutine adjustpm

  subroutine adjustmp(a,b,i)
    implicit none
    integer , intent(inout) :: a , b
    integer , intent(in) :: i
    if ( a < 1 ) then
      a = i - a
      b = b - 1
    end if
  end subroutine adjustmp

  subroutine normidate(idate)
    implicit none
    integer , intent(inout) :: idate
    if (idate < 10000) idate = idate*1000000+10100
    if (idate < 1000000) idate = idate*10000+100
    if (idate < 100000000) idate = idate*100
  end subroutine normidate

  subroutine split_idate(idate, iy, im, id, ih)
    implicit none
    integer , intent(in) :: idate
    integer , intent(out) :: iy , im , id , ih
    integer :: base , iidate
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
  end subroutine split_idate

  subroutine initfromintdt(x, i)
    implicit none
    type (rcm_time_and_date) , intent(out) :: x
    integer , intent(in) :: i
    call split_idate(i, x%year, x%month, x%day, x%hour)
    x%minute = 0
    x%second = 0
    x%calendar = gregorian
    call x%setup( )
  end subroutine initfromintdt

  subroutine initfromintit(x, i)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    integer , intent(in) :: i
    x%ival = i
    x%iunit = uhrs
  end subroutine initfromintit

  subroutine set_calint(x, c)
    implicit none
    class (rcm_time_and_date) , intent(inout) :: x
    integer , intent(in) :: c
    select case (c)
      case (gregorian)
        x%calendar = c
      case (noleap)
        x%calendar = c
      case (y360)
        x%calendar = c
      case default
        write (stderr,*) 'Unknown calendar, using Julian/Gregorian'
        x%calendar = gregorian
    end select
    call x%setup()
  end subroutine set_calint

  subroutine set_calstring(x, c)
    implicit none
    class (rcm_time_and_date) , intent(inout) :: x
    character (len=*) , intent(in) :: c

    if ( c == 'gregorian' ) then
      x%calendar = gregorian
    else if ( c(1:6) == 'noleap' .or.   &
              c(1:8) == 'days_365' .or. &
              c(1:7) == '365_day' ) then
      x%calendar = noleap
    else if ( c(1:7) == 'days_360' ) then
      x%calendar = y360
    else
      write (stderr,*) 'Unknown calendar, using Julian/Gregorian'
      x%calendar = gregorian
    end if
    call x%setup( )
  end subroutine set_calstring

  subroutine set_timeunit(x, u)
    implicit none
    class (rcm_time_interval) , intent(inout) :: x
    integer , intent(in) :: u
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
    x%calendar = y%calendar
    x%year = y%year
    x%month = y%month
    x%day = y%day
    x%hour = y%hour
    x%minute = y%minute
    x%second = y%second
    call x%setup()
  end subroutine initfromtypedt

  subroutine initfromtypeit(x, y)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    type (rcm_time_interval) , intent(in) :: y
    x%ival = y%ival
    x%iunit = y%iunit
  end subroutine initfromtypeit

  subroutine print_rcm_time_and_date(x)
    implicit none
    class (rcm_time_and_date) , intent(in) :: x
    write (stdout,'(i0.4,"-",i0.2,"-",i0.2,"T",i0.2,":",i0.2,":",i0.2)') &
       x%year, x%month, x%day, x%hour, x%minute, x%second
  end subroutine print_rcm_time_and_date

  function tochar(x) result(cdat)
    implicit none
    class (rcm_time_and_date) , intent(in) :: x
    character (len=32) :: cdat
    write (cdat,'(i0.4,"-",i0.2,"-",i0.2," ",i0.2,":",i0.2,":",i0.2," UTC")') &
       x%year, x%month, x%day, x%hour, x%minute, x%second
  end function tochar

  subroutine print_rcm_time_interval(x)
    implicit none
    class (rcm_time_interval) , intent(in) :: x
    write (stdout,'(i10," ",a)') x%ival, trim(cintstr(x%iunit))
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
    integer , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call yy%setcal(x%calendar)
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
            isequalit = ( x%ival == y%ival*60 )
          case (uhrs)
            isequalit = ( x%ival == y%ival*3600 )
          case (uday)
            isequalit = ( x%ival == y%ival*86400 )
          case default
            isequalit = .false.
        end select
      case (umin)
        select case (y%iunit)
          case (usec)
            isequalit = ( x%ival == y%ival/60 )
          case (umin)
            isequalit = ( x%ival == y%ival )
          case (uhrs)
            isequalit = ( x%ival == y%ival*60 )
          case (uday)
            isequalit = ( x%ival == y%ival*1440 )
          case default
            isequalit = .false.
        end select
      case (uhrs)
        select case (y%iunit)
          case (usec)
            isequalit = ( x%ival == y%ival/3600 )
          case (umin)
            isequalit = ( x%ival == y%ival/60 )
          case (uhrs)
            isequalit = ( x%ival == y%ival )
          case (uday)
            isequalit = ( x%ival == y%ival*24 )
          case default
            isequalit = .false.
        end select
      case (uday)
        select case (y%iunit)
          case (usec)
            isequalit = ( x%ival == y%ival/86400 )
          case (umin)
            isequalit = ( x%ival == y%ival/1440 )
          case (uhrs)
            isequalit = ( x%ival == y%ival*24 )
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
            isequalit = ( x%ival == y%ival*12 )
          case (ucnt)
            isequalit = ( x%ival == y%ival*1200 )
          case default
            isequalit = .false.
        end select
      case (uyrs)
        select case (y%iunit)
          case (umnt)
            isequalit = ( x%ival == y%ival*12 )
          case (uyrs)
            isequalit = ( x%ival == y%ival )
          case (ucnt)
            isequalit = ( x%ival == y%ival*100 )
          case default
            isequalit = .false.
        end select
      case (ucnt)
        select case (y%iunit)
          case (umnt)
            isequalit = ( x%ival == y%ival*1200 )
          case (uyrs)
            isequalit = ( x%ival == y%ival*100 )
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
    integer , intent(inout) :: d , m , y
    integer , intent(in) :: a
    integer :: tmp
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
    integer , intent(inout) :: d , m , y
    integer , intent(in) :: a
    integer :: tmp
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
    integer , intent(inout) :: d , m , y
    integer , intent(in) :: a
    integer :: tmp
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
    integer , intent(inout) :: d , m , y
    integer , intent(in) :: a
    integer :: tmp
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
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_interval) , intent(in) :: y
    type (rcm_time_and_date) :: z
    integer :: tmp
    z = x
    tmp = y%ival
    select case (y%iunit)
      case (usec)
        z%second_of_day = z%second_of_day+mod(tmp,86400)
        if ( z%second_of_day >= 86400 ) then
          z%second_of_day = z%second_of_day - 86400
          z%days_from_reference = z%days_from_reference + 1
        end if
        z%days_from_reference = z%days_from_reference + (tmp/86400)
        call z%recalculate()
      case (umin)
        tmp = tmp*60
        z%second_of_day = z%second_of_day+mod(tmp,86400)
        if ( z%second_of_day >= 86400 ) then
          z%second_of_day = z%second_of_day - 86400
          z%days_from_reference = z%days_from_reference + 1
        end if
        z%days_from_reference = z%days_from_reference + (tmp/86400)
        call z%recalculate()
      case (uhrs)
        tmp = tmp*3600
        z%second_of_day = z%second_of_day+mod(tmp,86400)
        if ( z%second_of_day >= 86400 ) then
          z%second_of_day = z%second_of_day - 86400
          z%days_from_reference = z%days_from_reference + 1
        end if
        z%days_from_reference = z%days_from_reference + (tmp/86400)
        call z%recalculate()
      case (uday)
        z%days_from_reference = z%days_from_reference + tmp
        call z%recalculate()
      case (umnt)
        z%month = z%month+mod(tmp,12)
        call adjustpm(z%month,z%year,12)
        tmp = tmp/12
        z%year = z%year+tmp
        call z%setup()
      case (uyrs)
        z%year = z%year+tmp
        call z%setup()
      case (ucnt)
        z%year = z%year+100*tmp
        call z%setup()
    end select
  end function add_interval

  function julianday(iy, im, id)
    implicit none
    real(dp) :: julianday
    integer , intent(in) :: iy , im , id
    integer :: ia , ib , iiy , iim
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
      integer :: icaltype
      integer , intent(in) :: iy , im , id
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
    if (abs(x%year - y%year) < 64) then
      z%iunit = usec
      z%ival = (x%second_of_day-y%second_of_day) +  &
               (x%days_from_reference-y%days_from_reference)*86400
    else
      z%iunit = uday
      z%ival = x%days_from_reference-y%days_from_reference
    end if
  end function diffdate

  function sub_interval(x,y) result (z)
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_interval) , intent(in) :: y
    type (rcm_time_and_date) :: z
    integer :: tmp
    z = x
    tmp = y%ival
    select case (y%iunit)
      case (usec)
        z%second_of_day = z%second_of_day-mod(tmp,86400)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = 86400 - z%second_of_day
          z%days_from_reference = z%days_from_reference - 1
        end if
        z%days_from_reference = z%days_from_reference - (tmp/86400)
        call z%recalculate()
      case (umin)
        tmp = tmp * 60
        z%second_of_day = z%second_of_day-mod(tmp,86400)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = 86400 - z%second_of_day
          z%days_from_reference = z%days_from_reference - 1
        end if
        z%days_from_reference = z%days_from_reference - (tmp/86400)
        call z%recalculate()
      case (uhrs)
        tmp = tmp * 3600
        z%second_of_day = z%second_of_day-mod(tmp,86400)
        if ( z%second_of_day < 0 ) then
          z%second_of_day = 86400 - z%second_of_day
          z%days_from_reference = z%days_from_reference - 1
        end if
        z%days_from_reference = z%days_from_reference - (tmp/86400)
        call z%recalculate()
      case (uday)
        z%days_from_reference = z%days_from_reference - tmp
        call z%recalculate()
      case (umnt)
        z%month = z%month-mod(tmp,12)
        call adjustmp(z%month,z%year,12)
        z%year = z%year-tmp/12
        call z%setup()
      case (uyrs)
        z%year = z%year-tmp
        call z%setup()
      case (ucnt)
        z%year = z%year-100*tmp
        call z%setup()
    end select
  end function sub_interval

  integer function toint(x) result(z)
    implicit none
    class (rcm_time_and_date) , intent(in) :: x
    z = x%year*1000000+x%month*10000+x%day*100+x%hour;
  end function toint

  logical function lfdoyear(x) result(l11)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    l11 = (x%month == 1 .and. x%day == 1)
  end function

  logical function lfdomonth(x) result(l1)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    l1 = (x%day == 1)
  end function

  logical function lmidnight(x) result(lm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    lm = (x%hour == 0 .and. x%minute == 0 .and. x%second == 0)
  end function

  logical function lsamemonth(x,y) result(lm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    lm = (x%month == y%month)
  end function

  integer function imondiff(x,y) result(z)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    z = (x%year-y%year)*12+(x%month-y%month)
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
    real(dp) :: jd
    if (x%calendar /= gregorian) then
      call die('mod_date','Error: week concept works only for gregorian calendar.',1)
    end if
    jd = julianday(x%year, x%month, x%day)
    iday = idint(dmod(jd+1.5D+00, 7.0D+00))+1
  end function idayofweek

  function setmidnight(x) result(mx)
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mx
    mx = x
    mx%hour = 0
    mx%minute = 0
    mx%second = 0
    call mx%setup()
  end function setmidnight

  logical function lsameweek(x,y) result(ls)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    type (rcm_time_and_date) :: xx , x1 , x2
    integer :: idwk
    xx = setmidnight(x)
    idwk = idayofweek(xx)
    x1 = xx - rcm_time_interval(idwk,uday)
    x2 = xx + rcm_time_interval(7-idwk,uday)
    ls = (y > x1) .and. (y < x2)
  end function lsameweek

  function ifdoweek(x) result(ifd)
    implicit none
    integer :: iwkday
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_interval) :: z
    type (rcm_time_and_date) :: ifd
    iwkday = idayofweek(x) - 1
    z = rcm_time_interval(iwkday,uday)
    ifd = setmidnight(x - z)
  end function ifdoweek

  function ildoweek(x) result(ild)
    implicit none
    integer :: iwkday
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
    iwk = z%ival/604800 + 1
  end function iwkdiff

  function monfirst(x) result(mf)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mf
    mf = rcm_time_and_date(x%calendar,x%year,x%month,1,0,0,0,0,0)
    call mf%setup( )
  end function monfirst

  function yrfirst(x) result(yf)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: yf
    yf = rcm_time_and_date(x%calendar,x%year,1,1,0,0,0,0,0)
    call yf%setup( )
  end function yrfirst

  function monlast(x) result(ml)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: ml
    select case (x%calendar)
      case (gregorian)
        ml = rcm_time_and_date(x%calendar,x%year,x%month, &
                               mdays_leap(x%year, x%month),0,0,0,0,0)
      case (noleap)
        ml = rcm_time_and_date(x%calendar,x%year,x%month,mlen(x%month),0,0,0,0,0)
      case (y360)
        ml = rcm_time_and_date(x%calendar,x%year,x%month,30,0,0,0,0,0)
    end select
    call ml%setup( )
  end function monlast

  function monmiddle(x) result(mm)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mm
    real(dp) :: rmom
    integer :: imom
    select case (x%calendar)
      case (gregorian)
        imom = mlen(x%month)/2
        rmom = dble(mdays_leap(x%year, x%month)) * 0.5D0
        mm = rcm_time_and_date(x%calendar,x%year,x%month, imom, &
                               idint((rmom-dble(imom))*24.0D0),0,0,0,0)
      case (noleap)
        imom = mlen(x%month)/2
        mm = rcm_time_and_date(x%calendar,x%year,x%month,imom, &
                               12*(mlen(x%month)-imom*2),0,0,0,0)
      case (y360)
        mm = rcm_time_and_date(x%calendar,x%year,x%month,15,0,0,0,0,0)
    end select
    call mm%setup( )
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

  function timeval2date(xval,cunit,ccal) result(dd)
    implicit none
    real(dp) , intent(in) :: xval
    character(*) , intent(in) :: cunit
    character(*) , intent(in) :: ccal

    type (rcm_time_and_date) :: dd
    type (rcm_time_interval) :: z , zz
    character(64) , save :: csave
    integer , save :: iunit
    integer :: year , month , day , hour , minute , second
    type (rcm_time_and_date) , save :: dref
    character(12) :: cdum

    data csave/'none'/

    if (csave == cunit) then
      z = idint(xval)
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
      z = idint(xval)
      if (cunit(1:5) == 'hours') then
        ! Unit is hours since reference
        if (len_trim(cunit) >= 30) then
          read(cunit,'(a12,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') cdum, year, &
            cdum, month, cdum, day, cdum, hour, cdum, minute, cdum, second
        else if (len_trim(cunit) >= 27) then
          read(cunit,'(a12,i4,a1,i2,a1,i2,a1,i2,a1,i2)') cdum, year, &
            cdum, month, cdum, day, cdum, hour, cdum, minute
          second = 0
        else if (len_trim(cunit) >= 24) then
          read(cunit,'(a12,i4,a1,i2,a1,i2,a1,i2)') cdum, year, &
            cdum, month, cdum, day, cdum, hour
          minute = 0
          second = 0
        else
          call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
        end if
        iunit = uhrs
      else if (cunit(1:4) == 'days') then
        ! Unit is days since reference
        if (len_trim(cunit) >= 30) then
          read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') cdum, year, &
            cdum, month, cdum, day, cdum, hour, cdum, minute, cdum, second
        else if (len_trim(cunit) >= 27) then
          read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2,a1,i2)') cdum, year, &
            cdum, month, cdum, day, cdum, hour, cdum, minute
          second = 0
        else if (len_trim(cunit) >= 24) then
          read(cunit,'(a11,i4,a1,i2,a1,i2,a1,i2)') cdum, year, &
            cdum, month, cdum, day, cdum, hour
          minute = 0
          second = 0
        else if (len_trim(cunit) >= 21) then
          read(cunit,'(a11,i4,a1,i2,a1,i2)') cdum, year, cdum, month, cdum, day
          hour = 0
          minute = 0
          second = 0
        else if (len_trim(cunit) >= 19) then
          read(cunit,'(a11,i4,a1,i1,a1,i1)') cdum, year, cdum, month, cdum, day
          hour = 0
          minute = 0
          second = 0
        else
          call die('mod_date','CANNOT PARSE TIME UNIT IN TIMEVAL2DATE')
        end if
        iunit = uday
      end if
      if ( ccal(1:6) == 'noleap' .or.   &
           ccal(1:8) == 'days_365' .or. &
           ccal(1:7) == '365_day' ) then
        dref = rcm_time_and_date(noleap,year,month,day,hour,minute,second,0,0)
      else if (ccal(1:7) == '360_day') then
        dref = rcm_time_and_date(y360,year,month,day,hour,minute,second,0,0)
      else
        dref = rcm_time_and_date(gregorian,year,month,day,hour,minute,second,0,0)
      end if
      call dref%setup( )
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
            interval_greater = ( x%ival > y%ival*60 )
          case (uhrs)
            interval_greater = ( x%ival > y%ival*3600 )
          case (uday)
            interval_greater = ( x%ival > y%ival*86400 )
          case default
            interval_greater = .false.
        end select
      case (umin)
        select case (y%iunit)
          case (usec)
            interval_greater = (y%ival > 60) .and. ( x%ival > y%ival/60 )
          case (umin)
            interval_greater = ( x%ival > y%ival )
          case (uhrs)
            interval_greater = ( x%ival > y%ival*60 )
          case (uday)
            interval_greater = ( x%ival > y%ival*1440 )
          case default
            interval_greater = .false.
        end select
      case (uhrs)
        select case (y%iunit)
          case (usec)
            interval_greater = (y%ival > 3600) .and. ( x%ival > y%ival/3600 )
          case (umin)
            interval_greater = (y%ival > 60) .and. ( x%ival > y%ival/60 )
          case (uhrs)
            interval_greater = ( x%ival > y%ival )
          case (uday)
            interval_greater = ( x%ival > y%ival*24 )
          case default
            interval_greater = .false.
        end select
      case (uday)
        select case (y%iunit)
          case (usec)
            interval_greater = (y%ival > 86400) .and. ( x%ival > y%ival/86400 )
          case (umin)
            interval_greater = (y%ival > 1440) .and. ( x%ival > y%ival/1440 )
          case (uhrs)
            interval_greater = (y%ival > 24) .and. ( x%ival > y%ival/24 )
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
            interval_greater = ( x%ival > y%ival*12 )
          case (ucnt)
            interval_greater = ( x%ival > y%ival*1200 )
          case default
            interval_greater = .false.
        end select
      case (uyrs)
        select case (y%iunit)
          case (umnt)
            interval_greater = ( x%ival > y%ival*12 )
          case (uyrs)
            interval_greater = ( x%ival > y%ival )
          case (ucnt)
            interval_greater = ( x%ival > y%ival*100 )
          case default
            interval_greater = .false.
        end select
      case (ucnt)
        select case (y%iunit)
          case (umnt)
            interval_greater = ( x%ival > y%ival*1200 )
          case (uyrs)
            interval_greater = ( x%ival > y%ival*100 )
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
            interval_less = ( x%ival < y%ival*60 )
          case (uhrs)
            interval_less = ( x%ival < y%ival*3600 )
          case (uday)
            interval_less = ( x%ival < y%ival*86400 )
          case default
            interval_less = .false.
        end select
      case (umin)
        select case (y%iunit)
          case (usec)
            interval_less = (y%ival > 60) .and. ( x%ival < y%ival/60 )
          case (umin)
            interval_less = ( x%ival < y%ival )
          case (uhrs)
            interval_less = ( x%ival < y%ival*60 )
          case (uday)
            interval_less = ( x%ival < y%ival*1440 )
          case default
            interval_less = .false.
        end select
      case (uhrs)
        select case (y%iunit)
          case (usec)
            interval_less = (y%ival > 3600) .and. ( x%ival < y%ival/3600 )
          case (umin)
            interval_less = (y%ival > 60) .and. ( x%ival < y%ival/60 )
          case (uhrs)
            interval_less = ( x%ival < y%ival )
          case (uday)
            interval_less = ( x%ival < y%ival*24 )
          case default
            interval_less = .false.
        end select
      case (uday)
        select case (y%iunit)
          case (usec)
            interval_less = (y%ival > 86400) .and. ( x%ival < y%ival/86400 )
          case (umin)
            interval_less = (y%ival > 1440) .and. ( x%ival < y%ival/1440 )
          case (uhrs)
            interval_less = (y%ival > 24) .and. ( x%ival < y%ival/24 )
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
            interval_less = ( x%ival < y%ival*12 )
          case (ucnt)
            interval_less = ( x%ival < y%ival*1200 )
          case default
            interval_less = .false.
        end select
      case (uyrs)
        select case (y%iunit)
          case (umnt)
            interval_less = ( x%ival < y%ival*12 )
          case (uyrs)
            interval_less = ( x%ival < y%ival )
          case (ucnt)
            interval_less = ( x%ival < y%ival*100 )
          case default
            interval_less = .false.
        end select
      case (ucnt)
        select case (y%iunit)
          case (umnt)
            interval_less = ( x%ival < y%ival*1200 )
          case (uyrs)
            interval_less = ( x%ival < y%ival*100 )
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
    integer , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call yy%setcal(x%calendar)
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
    integer , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call yy%setcal(x%calendar)
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
    integer , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call yy%setcal(x%calendar)
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
    integer , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call yy%setcal(x%calendar)
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
    integer , intent(in) :: y
    type (rcm_time_and_date) :: yy
    yy = y
    call yy%setcal(x%calendar)
    ln = date_ne(x,yy)
  end function idate_ne

  real(dp) function tohours(x) result(hs)
    implicit none
    class (rcm_time_interval) , intent(in) :: x
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

  real(dp) function yeardayfrac(x)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    yeardayfrac = dble(idayofyear(x)) + dble(x%hour)/24.0D+00 + &
                  dble(x%minute/1440.0D0) + dble(x%second/86400.0D0)
  end function yeardayfrac

  subroutine date_bcast(x,from,comm,ierr)
    class (rcm_time_and_date) , intent(inout) :: x
    integer , intent(in) :: from , comm
    integer , intent(out) :: ierr
    integer :: lerr
    ierr = 0
    call mpi_bcast(x%calendar,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%year,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%month,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%day,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%hour,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%minute,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%second,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%days_from_reference,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%second_of_day,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
  end subroutine date_bcast

end module mod_date
