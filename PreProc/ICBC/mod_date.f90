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

  character (len=16) , dimension(7) :: cintstr

  integer , dimension(12) :: mlen
  data mlen /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

  type rcm_time_and_date
    integer :: calendar
    integer :: year
    integer :: month
    integer :: day
    integer :: hour
    integer :: minute
    integer :: second
    contains
      procedure , pass :: printdate => print_rcm_time_and_date
      procedure , pass :: toidate => toint
      procedure , pass :: tostring => tochar
      procedure , pass(x) :: setcal => set_calendar
  end type rcm_time_and_date

  type rcm_time_interval
    integer :: ival
    integer :: iunit
    contains
      procedure , pass :: printint => print_rcm_time_interval
      procedure , pass(x) :: setunit => set_timeunit
      procedure , pass(x) :: hours => tohours
  end type rcm_time_interval

  interface assignment(=) 
    module procedure initfromintdt , initfromtypedt
    module procedure initfromintit , initfromtypeit
  end interface assignment(=)

  interface operator(==)
    module procedure isequaldt , isequalit
  end interface operator(==)

  interface operator(+)
    module procedure add_interval
  end interface operator(+)

  interface operator(-)
    module procedure diffdate , sub_interval
  end interface operator(-)

  interface operator(>)
    module procedure date_greater , interval_greater
  end interface

  interface operator(<)
    module procedure date_less , interval_less
  end interface

  public :: rcm_time_and_date , assignment(=) , operator(==)
  public :: rcm_time_interval , operator(+) , operator(-)
  public :: operator(>) , operator(<)
  public :: lsamemonth , imondiff , lfhomonth , monfirst , monlast , monmiddle
  public :: nextmon , prevmon , yrfirst
  public :: lsameweek , iwkdiff , idayofweek , ifdoweek , idayofyear
  public :: timeval2date

  data cintstr /'seconds', 'minutes', 'hours', 'days', &
                'months', 'years', 'centuries'/

  contains

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
  end subroutine initfromintdt

  subroutine initfromintit(x, i)
    implicit none
    type (rcm_time_interval) , intent(out) :: x
    integer , intent(in) :: i
    x%ival = i
    x%iunit = uhrs
  end subroutine initfromintit

  subroutine set_calendar(x, c)
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
  end subroutine set_calendar

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
    character (len=23) :: cdat
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
    isequaldt = ( x%calendar == y%calendar ) .and. &
                ( x%year == y%year )         .and. &
                ( x%month == y%month )       .and. &
                ( x%day == y%day )           .and. &
                ( x%hour == y%hour )         .and. &
                ( x%minute == y%minute )     .and. &
                ( x%second == y%second )
  end function isequaldt

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
  end function mdays_leap

  recursive subroutine add_days_leap(d,m,y,a)
    integer , intent(inout) :: d , m , y
    integer , intent(in) :: a
    integer :: tmp
    tmp = a
    if ( tmp > (mdays_leap(y,m)-d) ) then
      tmp = tmp - (mdays_leap(y,m)-d+1)
      d = 1
      m = mod(m,12)+1
      y = y+m/12
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
      if ( m == 0 ) then
        m = 12
        y = y-1
      end if
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
      d = 1
      m = mod(m,12)+1
      y = y+m/12
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
      if ( m == 0 ) then
        m = 12
        y = y-1
      end if
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
    select case (x%calendar)
      case (gregorian)
        select case (y%iunit)
          case (usec)
            z%second = z%second+mod(tmp,60)
            tmp = tmp/60
            z%minute = z%minute+mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            call add_days_leap(z%day, z%month, z%year, tmp)
          case (umin)
            z%minute = z%minute+mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            call add_days_leap(z%day, z%month, z%year, tmp)
          case (uhrs)
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            call add_days_leap(z%day, z%month, z%year, tmp)
          case (uday)
            call add_days_leap(z%day, z%month, z%year, tmp)
          case (umnt)
            z%month = z%month+mod(tmp,12)
            tmp = tmp/12
            z%year = z%year+tmp
          case (uyrs)
            z%year = z%year+tmp
          case (ucnt)
            z%year = z%year+100*tmp
        end select
      case (noleap)
        select case (y%iunit)
          case (usec)
            z%second = z%second+mod(tmp,60)
            tmp = tmp/60
            z%minute = z%minute+mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            call add_days_noleap(z%day, z%month, z%year, tmp)
          case (umin)
            z%minute = z%minute+mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            call add_days_noleap(z%day, z%month, z%year, tmp)
          case (uhrs)
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            call add_days_noleap(z%day, z%month, z%year, tmp)
          case (uday)
            call add_days_noleap(z%day, z%month, z%year, tmp)
          case (umnt)
            z%month = z%month+mod(tmp,12)
            tmp = tmp/12
            z%year = z%year+tmp
          case (uyrs)
            z%year = z%year+tmp
          case (ucnt)
            z%year = z%year+100*tmp
        end select
      case (y360)
        select case (y%iunit)
          case (usec)
            z%second = z%second+mod(tmp,60)
            tmp = tmp/60
            z%minute = z%minute+mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            z%day = z%day+mod(tmp,30)
            tmp = tmp/30
            z%month = z%month+mod(tmp,12)
            tmp = tmp/12
            z%year = z%year+tmp
          case (umin)
            z%minute = z%minute+mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            z%day = z%day+mod(tmp,30)
            tmp = tmp/30
            z%month = z%month+mod(tmp,12)
            tmp = tmp/12
            z%year = z%year+tmp
          case (uhrs)
            z%hour = z%hour+mod(tmp,24)
            tmp = tmp/24
            z%day = z%day+mod(tmp,30)
            tmp = tmp/30
            z%month = z%month+mod(tmp,12)
            tmp = tmp/12
            z%year = z%year+tmp
          case (uday)
            z%day = z%day+mod(tmp,30)
            tmp = tmp/30
            z%month = z%month+mod(tmp,12)
            tmp = tmp/12
            z%year = z%year+tmp
          case (umnt)
            z%month = z%month+mod(tmp,12)
            tmp = tmp/12
            z%year = z%year+tmp
          case (uyrs)
            z%year = z%year+tmp
          case (ucnt)
            z%year = z%year+100*tmp
        end select
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
      call die('mod_date','Dates not comparable: not using same calendar.',1)
    end if
  end subroutine check_cal

  function diffdate(x,y) result(z)
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) , intent(in) :: y
    type (rcm_time_interval) :: z
    real(dp) :: jd1 , jd2
    integer :: it1 , it2
    call check_cal(x,y)
    z%ival = 0
    z%iunit = usec
    select case (x%calendar)
      case (gregorian)
        jd2 = julianday(x%year, x%month, x%day)
        jd1 = julianday(y%year, y%month, y%day)
        z%ival = x%second-y%second +           &
                 (x%minute-y%minute) * 60 +    &
                 (x%hour-y%hour) * 3600 +      &
                 idnint(jd2-jd1)*86400
      case (noleap)
        it2 = (x%year-2000)*31536000+sum(mlen(1:x%month-1))*86400+(x%day-1)*86400
        it1 = (y%year-2000)*31536000+sum(mlen(1:y%month-1))*86400+(y%day-1)*86400
        z%ival = it2-it1 +                  &
                 x%second-y%second +        &
                 (x%minute-y%minute) * 60 + &
                 (x%hour-y%hour) * 3600
      case (y360)
        z%ival = x%second-y%second +           &
                 (x%minute-y%minute) * 60 +    &
                 (x%hour-y%hour) * 3600 +      &
                 (x%day-y%day) * 86400 +       &
                 (x%month-y%month) * 2592000 + &
                 (x%year-y%year) * 31104000
    end select

  end function diffdate

  function sub_interval(x,y) result (z)
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_interval) , intent(in) :: y
    type (rcm_time_and_date) :: z
    integer :: tmp
    z = x
    tmp = y%ival
    select case (x%calendar)
      case (gregorian)
        select case (y%iunit)
          case (usec)
            z%second = z%second-mod(tmp,60)
            tmp = tmp/60
            z%minute = z%minute-mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            call sub_days_leap(z%day, z%month, z%year, tmp)
          case (umin)
            z%minute = z%minute-mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            call sub_days_leap(z%day, z%month, z%year, tmp)
          case (uhrs)
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            call sub_days_leap(z%day, z%month, z%year, tmp)
          case (uday)
            call sub_days_leap(z%day, z%month, z%year, tmp)
          case (umnt)
            z%month = z%month-mod(tmp,12)
            tmp = tmp/12
            z%year = z%year-tmp
          case (uyrs)
            z%year = z%year-tmp
          case (ucnt)
            z%year = z%year-100*tmp
        end select
      case (noleap)
        select case (y%iunit)
          case (usec)
            z%second = z%second-mod(tmp,60)
            tmp = tmp/60
            z%minute = z%minute-mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            call sub_days_noleap(z%day, z%month, z%year, tmp)
          case (umin)
            z%minute = z%minute-mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            call sub_days_noleap(z%day, z%month, z%year, tmp)
          case (uhrs)
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            call sub_days_noleap(z%day, z%month, z%year, tmp)
          case (uday)
            call sub_days_noleap(z%day, z%month, z%year, tmp)
          case (umnt)
            z%month = z%month-mod(tmp,12)
            tmp = tmp/12
            z%year = z%year-tmp
          case (uyrs)
            z%year = z%year-tmp
          case (ucnt)
            z%year = z%year-100*tmp
        end select
      case (y360)
        select case (y%iunit)
          case (usec)
            z%second = z%second-mod(tmp,60)
            tmp = tmp/60
            z%minute = z%minute-mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            z%day = z%day-mod(tmp,30)
            tmp = tmp/30
            z%month = z%month-mod(tmp,12)
            tmp = tmp/12
            z%year = z%year-tmp
          case (umin)
            z%minute = z%minute-mod(tmp,60)
            tmp = tmp/60
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            z%day = z%day-mod(tmp,30)
            tmp = tmp/30
            z%month = z%month-mod(tmp,12)
            tmp = tmp/12
            z%year = z%year-tmp
          case (uhrs)
            z%hour = z%hour-mod(tmp,24)
            tmp = tmp/24
            z%day = z%day-mod(tmp,30)
            tmp = tmp/30
            z%month = z%month-mod(tmp,12)
            tmp = tmp/12
            z%year = z%year-tmp
          case (uday)
            z%day = z%day-mod(tmp,30)
            tmp = tmp/30
            z%month = z%month-mod(tmp,12)
            tmp = tmp/12
            z%year = z%year-tmp
          case (umnt)
            z%month = z%month-mod(tmp,12)
            tmp = tmp/12
            z%year = z%year-tmp
          case (uyrs)
            z%year = z%year-tmp
          case (ucnt)
            z%year = z%year-100*tmp
        end select
    end select
  end function sub_interval

  integer function toint(x) result(z)
    implicit none
    class (rcm_time_and_date) , intent(in) :: x
    z = x%year*1000000+x%month*10000+x%day*100+x%hour;
  end function toint

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

  integer function iwkdiff(x,y) result(iwk)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    type (rcm_time_interval) :: z
    call check_cal(x,y)
    z = x - y
    iwk = z%ival/604800
  end function iwkdiff

  function monfirst(x) result(mf)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mf
    mf = rcm_time_and_date(x%calendar,x%year,x%month,1,0,0,0)
  end function monfirst

  function yrfirst(x) result(mf)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: mf
    mf = rcm_time_and_date(x%calendar,x%year,1,1,0,0,0)
  end function yrfirst

  function monlast(x) result(ml)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
    type (rcm_time_and_date) :: ml
    select case (x%calendar)
      case (gregorian)
        ml = rcm_time_and_date(x%calendar,x%year,x%month, &
                               mdays_leap(x%year, x%month), 0, 0 ,0)
      case (noleap)
        ml = rcm_time_and_date(x%calendar,x%year,x%month,mlen(x%month), 0, 0 ,0)
      case (y360)
        ml = rcm_time_and_date(x%calendar,x%year,x%month,30, 0, 0 ,0)
    end select
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
                               idint((rmom-dble(imom))*24.0D0),0,0)
      case (noleap)
        imom = mlen(x%month)/2
        mm = rcm_time_and_date(x%calendar,x%year,x%month,imom, &
                               12*(mlen(x%month)-imom*2),0,0)
      case (y360)
        mm = rcm_time_and_date(x%calendar,x%year,x%month,15,0,0,0)
    end select
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

  integer function idayofyear(x) result(id)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x
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
  end function

  function timeval2date(xval,cunit,ccal) result(dd)
    implicit none
    real(dp) , intent(in) :: xval
    character(*) , intent(in) :: cunit
    character(*) , intent(in) :: ccal

    type (rcm_time_and_date) :: dd
    type (rcm_time_interval) :: z
    character(35) , save :: csave
    integer :: year , month , day , hour , minute , second
    type (rcm_time_and_date) , save :: dref
    character(12) :: cdum

    data csave/'none'/

    dd = dref
    if (csave == cunit) then
      z = idnint(xval)
      dd = dref + z
    else
      if (len_trim(cunit) < 35) then
        dd = 0
      else
        csave = cunit
        read(cunit,'(a12,i4,a1,i2,a1,i2,a1,i2,a1,i2,a1,i2)') cdum, year, &
          cdum, month, cdum, day, cdum, hour, cdum, minute, cdum, second
        if (ccal == 'noleap' .or. ccal == '365_day') then
          dref = rcm_time_and_date(noleap,year,month,day,hour,minute,second)
        else if (ccal == '360_day') then
          dref = rcm_time_and_date(y360,year,month,day,hour,minute,second)
        else
          dref = rcm_time_and_date(gregorian,year,month,day,hour,minute,second)
        end if
        z = idnint(xval)
        dd = dref + z
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
            interval_greater = ( x%ival > y%ival/60 )
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
            interval_greater = ( x%ival > y%ival/3600 )
          case (umin)
            interval_greater = ( x%ival > y%ival/60 )
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
            interval_greater = ( x%ival > y%ival/86400 )
          case (umin)
            interval_greater = ( x%ival > y%ival/1440 )
          case (uhrs)
            interval_greater = ( x%ival > y%ival*24 )
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
            interval_less = ( x%ival < y%ival/60 )
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
            interval_less = ( x%ival < y%ival/3600 )
          case (umin)
            interval_less = ( x%ival < y%ival/60 )
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
            interval_less = ( x%ival < y%ival/86400 )
          case (umin)
            interval_less = ( x%ival < y%ival/1440 )
          case (uhrs)
            interval_less = ( x%ival < y%ival*24 )
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
    gt = ((x - y) > rcm_time_interval(1,usec))
  end function date_greater

  logical function date_less(x,y) result(lt)
    implicit none
    type (rcm_time_and_date) , intent(in) :: x , y
    call check_cal(x,y)
    lt = ((x - y) < rcm_time_interval(1,usec))
  end function date_less

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

end module mod_date
