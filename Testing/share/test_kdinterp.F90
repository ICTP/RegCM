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

! Unit test of Share/mod_kdinterp.F90: from a global 1 degree grid to a
! curvilinear grid, and between point lists.

program test_kdinterp
  use mod_intkinds
  use mod_realkinds
  use mod_kdinterp
  use share_test_utils
  implicit none
  integer, parameter :: nlon = 360, nlat = 181, nj = 40, ni = 30, nk = 3
  real(rk8), parameter :: d2r = acos(-1.0_rk8)/180.0_rk8
  real(rkx), dimension(nlat) :: slat
  real(rkx), dimension(nlon) :: slon
  real(rkx), dimension(nj,ni) :: tlat, tlon
  real(rk8), dimension(nlon,nlat) :: g
  real(rk8), dimension(nj,ni) :: f, ex
  real(rk8), dimension(nlon,nlat,nk) :: g3
  real(rk8), dimension(nj,ni,nk) :: f3, f3n, ftmp
  real(rkx), dimension(nlon,nlat,1) :: gs
  real(rkx), dimension(nj,ni,1) :: fs
  real(rkx), dimension(50) :: s1lat, s1lon
  real(rkx), dimension(5) :: t1lat, t1lon
  real(rk8), dimension(50) :: g1
  real(rk8), dimension(5) :: f1
  real(rkx), dimension(nj,3) :: sl, so
  type(h_interpolator) :: hi, hs
  integer :: i, j, k

  do j = 1, nlat
    slat(j) = -90.0_rkx + real(j-1,rkx)
  end do
  do i = 1, nlon
    slon(i) = -180.0_rkx + real(i-1,rkx)
  end do
  do i = 1, ni
    do j = 1, nj
      tlat(j,i) = -35.0_rkx + 0.5_rkx*real(i-1,rkx) + 0.05_rkx*real(j-1,rkx)
      tlon(j,i) = -70.0_rkx + 0.5_rkx*real(j-1,rkx) - 0.05_rkx*real(i-1,rkx)
    end do
  end do
  call h_interpolator_create(hi,slat,slon,tlat,tlon)

  g = 7.25_rk8
  call h_interpolate_cont(hi,g,f)
  call check(maxval(abs(f-7.25_rk8)) < 1.0e-12_rk8,'constant field')

  do j = 1, nlat
    do i = 1, nlon
      g(i,j) = 280.0_rk8 + 20.0_rk8*cos(slat(j)*d2r)*cos(2.0_rk8*slon(i)*d2r)
    end do
  end do
  ex = 280.0_rk8 + 20.0_rk8*cos(tlat*d2r)*cos(2.0_rk8*tlon*d2r)
  call h_interpolate_cont(hi,g,f)
  call check(maxval(abs(f-ex)) < 0.01_rk8,'smooth field within 0.01 K')

  do k = 1, nk
    g3(:,:,k) = g*real(k,rk8)
  end do
  call h_interpolate_cont(hi,g3,f3)
  call h_interpolate_cont(hi,g3,f3n,ftmp,nlon,nlat,nj,ni,nk)
  call check(maxval(abs(f3-f3n)) < 1.0e-10_rk8, &
    '3-D level by level and OpenACC variants agree')
  call h_interpolator_destroy(hi)

  ! Point lists, with missing values in the source.
  do i = 1, 50
    s1lat(i) = -30.0_rkx + 0.2_rkx*real(i,rkx)
    s1lon(i) = -60.0_rkx + 0.1_rkx*real(i,rkx)
    g1(i) = 10.0_rk8 + real(i,rk8)
  end do
  g1(20:50) = h_missing_value
  t1lat = [-29.0_rkx, -27.0_rkx, -25.0_rkx, -23.0_rkx, -21.0_rkx]
  t1lon = s1lon(5:45:10)
  call h_interpolator_create(hs,50,s1lat,s1lon,5,t1lat,t1lon)
  call h_interpolate_cont(hs,g1,f1)
  call check(all(f1(1:2) > 10.0_rk8 .and. f1(1:2) < 30.0_rk8), &
    'point list interpolation')
  call check(all(f1(3:5) == h_missing_value), &
    'point list keeps missing values')
  call h_interpolator_destroy(hs)

  ! The class fraction mean ignores the missing source points.
  gs = 50.0_rkx
  gs(100:130,50:70,1) = h_missing_value
  call h_interpolator_create(hs,slat,slon,tlat,tlon,300.0_rkx)
  call h_interpolate_sum(hs,gs,fs)
  call check(all(abs(fs-50.0_rkx) < 1.0e-4_rkx .or. &
                 fs == h_missing_value), &
    'sum interpolation ignores missing source points')
  call h_interpolator_destroy(hs)

  ! Fewer than 10 target rows: divided by zero in the progress dots.
  do i = 1, 3
    do j = 1, nj
      sl(j,i) = -20.0_rkx + real(i,rkx)
      so(j,i) = -60.0_rkx + 0.5_rkx*real(j,rkx)
    end do
  end do
  call h_interpolator_create(hs,slat,slon,sl,so)
  call check(.true.,'target grid with 3 rows')
  call h_interpolator_destroy(hs)

  call test_end('test_kdinterp')
end program test_kdinterp
