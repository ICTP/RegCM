!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Unit test of Share/mod_date.F90: intervals in different units, month
! arithmetic, subtraction as the inverse of addition, time units with one
! or two digit hours, the proleptic gregorian julian day, and day by day
! stepping on the three calendars.

program test_date
  use mod_intkinds
  use mod_realkinds
  use mod_date
  use share_test_utils
  implicit none
  type(rcm_time_interval) :: a
  type(rcm_time_and_date) :: d, e, d0
  integer(ik4), dimension(12), parameter :: mlen = &
    [31,28,31,30,31,30,31,31,30,31,30,31]
  integer(ik4), dimension(3), parameter :: units = [usec, umin, uhrs]
  character(len=10), dimension(3), parameter :: cals = &
    ['gregorian ','noleap    ','360_day   ']
  integer(ik4) :: y, m, k, ml, nbad, ndy, ic, iy0
  real(rk8) :: prev, jd

  a = rcm_time_interval(1,uday) + rcm_time_interval(1,uhrs)
  call check(a%ival == 90000_ik8,'1 day + 1 hour is 90000 s')
  call check(rcm_time_interval(1,uday) == rcm_time_interval(24,uhrs), &
    '1 day == 24 hours')
  call check(.not. (rcm_time_interval(1,uday) == &
                    rcm_time_interval(1,uhrs)),'1 day /= 1 hour')
  call check(.not. (rcm_time_interval(1,umin) == &
                    rcm_time_interval(61,usec)),'1 minute /= 61 s')
  call check(rcm_time_interval(2,umin) > rcm_time_interval(90,usec), &
    '2 minutes > 90 s')
  call check(rcm_time_interval(1,umin) < rcm_time_interval(90,usec), &
    '1 minute < 90 s')
  call check(rcm_time_interval(2,ucnt) > rcm_time_interval(1,ucnt), &
    '2 centuries > 1 century')
  call check(.not. (rcm_time_interval(1,uhrs) > &
                    rcm_time_interval(1,uday)),'not 1 hour > 1 day')

  d = 2001011500
  e = d - rcm_time_interval(1,umnt)
  call check(toint10(e) == 2000121500_ik8,'2001-01-15 - 1 month')
  e = d - rcm_time_interval(2,umnt)
  call check(toint10(e) == 2000111500_ik8,'2001-01-15 - 2 months')
  e = d - rcm_time_interval(14,umnt)
  call check(toint10(e) == 1999111500_ik8,'2001-01-15 - 14 months')
  ! Adding months keeps the fraction of the month: the 15th of January
  ! becomes the 14th of February, the middle of the month.
  e = d + rcm_time_interval(13,umnt)
  call check(toint10(e) == 2002021400_ik8,'2001-01-15 + 13 months')
  e = d + rcm_time_interval(-1,umnt)
  call check(toint10(e) == 2000121500_ik8,'2001-01-15 + (-1) month')
  e = d + rcm_time_interval(-25,umnt)
  call check(toint10(e) == 1998121500_ik8,'2001-01-15 + (-25) months')
  ! Subtracting months keeps the day, or the last day of a shorter month.
  d = 2001033100
  call check(toint10(prevmon(d)) == 2001022800_ik8, &
    'previous month of 2001-03-31')
  d = 2000022900
  e = d - rcm_time_interval(12,umnt)
  call check(toint10(e) == 1999022800_ik8,'2000-02-29 - 12 months')
  d = 2001033000
  call setcal(d,'360_day')
  e = d - rcm_time_interval(1,umnt)
  call check(toint10(e) == 2001023000_ik8,'360_day: 03-30 - 1 month')

  ! Subtraction is the addition of the opposite interval.
  nbad = 0
  do k = 0, 23
    d = 2001011500 + k
    do m = 0, 200000, 37
      do ic = 1, 3
        a = rcm_time_interval(m,units(ic))
        e = d - a
        d0 = d + rcm_time_interval(-m,units(ic))
        if ( .not. (e == d0) ) nbad = nbad + 1
        if ( .not. (e + a == d) ) nbad = nbad + 1
      end do
    end do
  end do
  call check(nbad == 0,'subtraction of seconds, minutes and hours')

  d = timeval2date(0.0_rk8,'hours since 1979-01-01 12:30:45','gregorian')
  call check(index(tochar(d),'1979-01-01 12:30:45') > 0, &
    'time unit with a two digit hour')
  d = timeval2date(2.0_rk8,'hours since 1979-1-1 6:00:00','365_day')
  call check(index(tochar(d),'1979-01-01 08:00:00') > 0, &
    'time unit with a one digit hour')
  d = timeval2date(1.0_rk8,'days since 1949-12-01 00:00:00','gregorian')
  call check(toint10(d) == 1949120200_ik8,'days since')

  call check(abs(julianday(2000,1,1)-2451544.5_rk8) < 1.0e-9_rk8, &
    'julian day of 2000-01-01')
  call check(abs(julianday(1582,10,15)-2299160.5_rk8) < 1.0e-9_rk8, &
    'julian day of 1582-10-15')
  call check(abs(julianday(1,1,1)-1721425.5_rk8) < 1.0e-9_rk8, &
    'julian day of 0001-01-01 (proleptic gregorian)')
  nbad = 0
  prev = julianday(-400,1,1) - 1.0_rk8
  do y = -400, 2800
    do m = 1, 12
      ml = mlen(m)
      if ( m == 2 .and. leap(y) ) ml = 29
      do k = 1, ml
        jd = julianday(y,m,k)
        if ( abs(jd-prev-1.0_rk8) > 1.0e-9_rk8 ) nbad = nbad + 1
        prev = jd
      end do
    end do
  end do
  call check(nbad == 0,'julian day is continuous from -400 to 2800')

  ! Day by day: every year must have the length of its calendar.
  do ic = 1, 3
    d0 = 1950010100
    call setcal(d0,trim(cals(ic)))
    d = d0
    nbad = 0
    iy0 = 1950
    ndy = 0
    do k = 1, 150*366
      d = d + rcm_time_interval(1,uday)
      ndy = ndy + 1
      if ( getyear(d) /= iy0 ) then
        if ( ndy /= yearlen(ic,iy0) ) nbad = nbad + 1
        iy0 = getyear(d)
        ndy = 0
      end if
      if ( abs(hourdiff(d,d0)-24.0_rk8*real(k,rk8)) > 1.0e-6_rk8 ) then
        nbad = nbad + 1
      end if
    end do
    call check(nbad == 0,'day stepping on the '//trim(cals(ic))// &
      ' calendar')
  end do

  call test_end('test_date')

  contains

  logical function leap(y)
    implicit none
    integer(ik4), intent(in) :: y
    leap = (modulo(y,4) == 0 .and. modulo(y,100) /= 0) .or. &
           modulo(y,400) == 0
  end function leap

  integer(ik4) function yearlen(ic,y)
    implicit none
    integer(ik4), intent(in) :: ic, y
    select case (ic)
      case (1)
        yearlen = 365
        if ( leap(y) ) yearlen = 366
      case (2)
        yearlen = 365
      case default
        yearlen = 360
    end select
  end function yearlen

end program test_date
