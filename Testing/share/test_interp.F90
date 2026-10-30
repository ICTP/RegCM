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

! Unit test of Share/mod_interp.F90 (interp1d): exact for linear data with
! any tension, the same result for decreasing coordinates and for the
! scalar, single precision and OpenACC variants, constant and linear
! extrapolation, coincident or unordered input coordinates, and inputs of
! one or two points.

program test_interp
  use mod_intkinds
  use mod_realkinds
  use mod_interp
  use share_test_utils
  implicit none
  integer(ik4), parameter :: n = 9, m = 7
  real(rk8), dimension(n) :: xi, g, zi, zg
  real(rk8), dimension(m) :: xo, f, f2, ex
  real(rk4), dimension(n) :: xi4, g4
  real(rk4), dimension(m) :: xo4, f4
  real(rk8), dimension(3) :: f3
  real(rk8) :: f1
  real(rk8), dimension(3), parameter :: alfas = [0.0_rk8, 0.5_rk8, 1.0_rk8]
  integer(ik4) :: k
  character(len=80) :: msg

  xi = [0.0_rk8, 0.5_rk8, 1.2_rk8, 2.0_rk8, 3.1_rk8, 4.0_rk8, &
        5.5_rk8, 6.0_rk8, 8.0_rk8]
  xo = [0.1_rk8, 0.9_rk8, 2.5_rk8, 3.1_rk8, 4.7_rk8, 5.9_rk8, 7.3_rk8]

  g = 3.0_rk8*xi - 2.0_rk8
  ex = 3.0_rk8*xo - 2.0_rk8
  do k = 1, 3
    call interp1d(xi,g,xo,f,alfas(k),1.0_rk8,1.0_rk8)
    write(msg,'(a,f3.1)') 'linear data, tension ', alfas(k)
    call check(maxval(abs(f-ex)) < 1.0e-12_rk8,trim(msg))
  end do

  ! A smooth function: the cubic (no tension) is close in the interior.
  g = sin(xi)
  call interp1d(xi,g,xo,f,0.0_rk8,1.0_rk8,1.0_rk8)
  call check(maxval(abs(f(1:6)-sin(xo(1:6)))) < 0.05_rk8, &
    'sin(x) without tension')
  call check(abs(f(4)-sin(xo(4))) < 1.0e-12_rk8, &
    'value at an input coordinate')

  call interp1d(xi(n:1:-1),g(n:1:-1),xo,f2,0.0_rk8,1.0_rk8,1.0_rk8)
  call check(maxval(abs(f-f2)) < 1.0e-12_rk8,'decreasing coordinates')

  call interp1d(xi,g,xo(3),f1,0.0_rk8,1.0_rk8,1.0_rk8)
  call check(abs(f1-f(3)) < 1.0e-12_rk8,'scalar variant')

  call interp1d(xi,g,xo,f2,0.0_rk8,1.0_rk8,1.0_rk8,zi,zg,n,m)
  call check(maxval(abs(f-f2)) < 1.0e-12_rk8,'OpenACC variant')

  xi4 = real(xi,rk4)
  g4 = real(g,rk4)
  xo4 = real(xo,rk4)
  call interp1d(xi4,g4,xo4,f4,0.0_rk4,1.0_rk4,1.0_rk4)
  call check(maxval(abs(real(f4,rk8)-f)) < 1.0e-5_rk8, &
    'single precision variant')

  ! Below the first point constant (ex1 = 0), above the last linear.
  g = 3.0_rk8*xi - 2.0_rk8
  call interp1d(xi,g,[-1.0_rk8, 9.0_rk8],f(1:2),0.0_rk8,0.0_rk8,1.0_rk8)
  call check(abs(f(1)+2.0_rk8) < 1.0e-12_rk8 .and. &
             abs(f(2)-25.0_rk8) < 1.0e-12_rk8,'extrapolation')

  ! Coincident coordinates are merged.
  xi(4) = xi(3)
  g = 3.0_rk8*xi - 2.0_rk8
  call interp1d(xi,g,xo,f,0.0_rk8,1.0_rk8,1.0_rk8)
  call check(maxval(abs(f-ex)) < 1.0e-12_rk8,'coincident coordinates')

  call interp1d([0.0_rk8, 1.0_rk8],[0.0_rk8, 2.0_rk8],[0.5_rk8], &
                f(1:1),0.0_rk8,1.0_rk8,1.0_rk8)
  call check(abs(f(1)-1.0_rk8) < 1.0e-12_rk8,'two input points')
  call interp1d([2.0_rk8],[5.0_rk8],[1.0_rk8, 2.0_rk8, 3.0_rk8], &
                f3,0.0_rk8,1.0_rk8,1.0_rk8)
  call check(all(abs(f3-5.0_rk8) < 1.0e-12_rk8),'one input point')

  ! Coordinates out of order are merged until they increase.
  call interp1d([0.0_rk8, 2.0_rk8, 1.0_rk8, 3.0_rk8], &
                [0.0_rk8, 4.0_rk8, 2.0_rk8, 6.0_rk8], &
                [0.5_rk8, 1.5_rk8, 2.5_rk8],f3,1.0_rk8,1.0_rk8,1.0_rk8)
  call check(all(abs(f3-[1.0_rk8, 3.0_rk8, 5.0_rk8]) < 1.0e-12_rk8), &
    'one coordinate out of order')
  call interp1d([0.0_rk8, 3.0_rk8, 2.0_rk8, 1.0_rk8, 4.0_rk8], &
                [0.0_rk8, 6.0_rk8, 4.0_rk8, 2.0_rk8, 8.0_rk8], &
                [0.5_rk8, 2.0_rk8, 3.5_rk8],f3,1.0_rk8,1.0_rk8,1.0_rk8)
  call check(all(abs(f3-[1.0_rk8, 4.0_rk8, 7.0_rk8]) < 1.0e-12_rk8), &
    'several coordinates out of order')

  call test_end('test_interp')
end program test_interp
