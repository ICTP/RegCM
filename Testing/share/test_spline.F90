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

! Unit test of Share/mod_spline.F90. A cubic spline with the exact end
! conditions reproduces a cubic polynomial: spline and splint with the end
! first derivatives, spline (Janjic) and spline1d with the end second
! derivatives, and splini integrates it exactly. A natural spline is exact
! for linear data. Points outside the input range are also checked.

program test_spline
  use mod_intkinds
  use mod_realkinds
  use mod_spline
  use share_test_utils
  implicit none
  integer(ik4), parameter :: n = 9, m = 7
  real(rkx), dimension(n) :: xa, ya, y2a, p, q
  real(rkx), dimension(m) :: xo, yo
  real(rkx), dimension(n) :: xn, yn
  real(rkx), dimension(3) :: x3, y3, y23
  real(rkx) :: y, yi, err
  integer(ik4) :: k

  xa = [0.0_rkx, 0.5_rkx, 1.2_rkx, 2.0_rkx, 3.1_rkx, 4.0_rkx, &
        5.5_rkx, 6.0_rkx, 8.0_rkx]
  xo = [0.1_rkx, 0.9_rkx, 2.5_rkx, 3.1_rkx, 4.7_rkx, 5.9_rkx, 7.3_rkx]
  ya = f(xa)

  call spline(xa,ya,df(xa(1)),df(xa(n)),y2a)
  err = 0.0_rkx
  do k = 1, m
    call splint(xa,ya,y2a,xo(k),y)
    err = max(err,abs(y-f(xo(k))))
  end do
  call check(err < 1.0e-10_rkx,'spline and splint reproduce a cubic')
  call check(maxval(abs(y2a-d2f(xa))) < 1.0e-10_rkx, &
    'spline gives the second derivatives of a cubic')

  err = 0.0_rkx
  do k = 1, m
    call splini(xa,ya,y2a,xo(k),yi)
    err = max(err,abs(yi-(intf(xo(k))-intf(xa(1)))))
  end do
  call check(err < 1.0e-10_rkx,'splini integrates a cubic')

  ya = 3.0_rkx*xa - 2.0_rkx
  call spline(xa,ya,1.0e30_rkx,1.0e30_rkx,y2a)
  err = 0.0_rkx
  do k = 1, m
    call splint(xa,ya,y2a,xo(k),y)
    err = max(err,abs(y-(3.0_rkx*xo(k)-2.0_rkx)))
  end do
  call check(err < 1.0e-10_rkx,'natural spline is exact for linear data')

  ! Janjic spline (mod_capecin): the end second derivatives are given.
  ya = f(xa)
  y2a = 0.0_rkx
  y2a(1) = d2f(xa(1))
  y2a(n) = d2f(xa(n))
  ! All its arrays have the same dimension jtb.
  xn = 0.0_rkx
  xn(1:m) = xo
  call spline(n,n,xa,ya,y2a,m,xn,yn,p,q)
  call check(maxval(abs(yn(1:m)-f(xo))) < 1.0e-10_rkx, &
    'Janjic spline reproduces a cubic')
  ! A point above the last one, then one in the first interval: the
  ! coefficients of the interval must be computed again.
  ya = 2.0_rkx*xa + 1.0_rkx
  y2a = 0.0_rkx
  xn = 0.0_rkx
  xn(1:2) = [9.0_rkx, 0.25_rkx]
  call spline(n,n,xa,ya,y2a,2,xn,yn,p,q)
  call check(abs(yn(1)-ya(n)) < 1.0e-12_rkx .and. &
             abs(yn(2)-1.5_rkx) < 1.0e-12_rkx, &
    'Janjic spline after a point above the last one')
  ya = f(xa)

  y2a = 0.0_rkx
  y2a(1) = d2f(xa(1))
  y2a(n) = d2f(xa(n))
  call spline1d(n,xa,ya,y2a,m,xo,yo)
  call check(maxval(abs(yo-f(xo))) < 1.0e-10_rkx, &
    'spline1d reproduces a cubic')
  call spline1d(n,xa,ya,y2a,2,[-1.0_rkx, 9.0_rkx],yo(1:2))
  call check(yo(1) == ya(1) .and. yo(2) == ya(n), &
    'spline1d is constant outside the points')

  ! Three points, as used for the boundary nudging profile.
  x3 = [0.0_rkx, 0.4_rkx, 1.0_rkx]
  y3 = f(x3)
  y23 = [d2f(x3(1)), 0.0_rkx, d2f(x3(3))]
  call spline1d(3,x3,y3,y23,m,xo/8.0_rkx,yo)
  call check(maxval(abs(yo-f(xo/8.0_rkx))) < 1.0e-10_rkx, &
    'spline1d with three points reproduces a cubic')

  call test_end('test_spline')

  contains

  elemental real(rkx) function f(x)
    implicit none
    real(rkx), intent(in) :: x
    f = 0.5_rkx*x**3 - 2.0_rkx*x**2 + x + 3.0_rkx
  end function f

  elemental real(rkx) function df(x)
    implicit none
    real(rkx), intent(in) :: x
    df = 1.5_rkx*x**2 - 4.0_rkx*x + 1.0_rkx
  end function df

  elemental real(rkx) function d2f(x)
    implicit none
    real(rkx), intent(in) :: x
    d2f = 3.0_rkx*x - 4.0_rkx
  end function d2f

  elemental real(rkx) function intf(x)
    implicit none
    real(rkx), intent(in) :: x
    intf = 0.125_rkx*x**4 - 2.0_rkx/3.0_rkx*x**3 + 0.5_rkx*x**2 + 3.0_rkx*x
  end function intf

end program test_spline
