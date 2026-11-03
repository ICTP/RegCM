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
! curvilinear grid, and between point lists; the nearest neighbour against
! a brute force search, and the real class interpolation equal to it;
! from a coarser grid, close to bilinear also across the date line.

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
  real(rkx), dimension(nlon,nlat) :: gc
  real(rkx), dimension(nj,ni) :: fc
  real(rkx), dimension(nj,ni,1) :: fs
  real(rkx), dimension(50) :: s1lat, s1lon
  real(rkx), dimension(5) :: t1lat, t1lon
  real(rk8), dimension(50) :: g1
  real(rk8), dimension(5) :: f1
  real(rkx), dimension(nj,3) :: sl, so
  type(h_interpolator) :: hi, hs
  integer :: i, j, k, ii, jj, ib, jb, jp
  real(rk8) :: c, cmax, wx, wy, xl
  logical :: lnear

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

  ! Nearest neighbour: the source index of the closest point.
  do j = 1, nlat
    do i = 1, nlon
      g(i,j) = real(i + 1000*j,rk8)
    end do
  end do
  call h_interpolate_nn(hi,g,f)
  lnear = .true.
  do i = 1, ni
    do j = 1, nj
      cmax = -2.0_rk8
      do jj = 1, nlat
        do ii = 1, nlon
          c = sin(tlat(j,i)*d2r)*sin(slat(jj)*d2r) + &
              cos(tlat(j,i)*d2r)*cos(slat(jj)*d2r) * &
              cos((tlon(j,i)-slon(ii))*d2r)
          if ( c > cmax ) then
            cmax = c
            ib = ii
            jb = jj
          end if
        end do
      end do
      ! Accept ties: the chosen point must be as close as the closest.
      ii = nint(f(j,i)) - 1000*(nint(f(j,i))/1000)
      jj = nint(f(j,i))/1000
      c = sin(tlat(j,i)*d2r)*sin(slat(jj)*d2r) + &
          cos(tlat(j,i)*d2r)*cos(slat(jj)*d2r) * &
          cos((tlon(j,i)-slon(ii))*d2r)
      lnear = lnear .and. c > cmax - 1.0e-10_rk8
    end do
  end do
  call check(lnear,'nearest neighbour is the closest source point')
  gc = real(g,rkx)
  call h_interpolate_class(hi,gc,fc)
  call check(all(real(fc,rk8) == f), &
    'real class interpolation takes the nearest source point')
  g(ib,jb) = h_missing_value
  call h_interpolate_nn(hi,g,f)
  call check(f(nj,ni) /= h_missing_value .and. &
             f(nj,ni) /= real(ib + 1000*jb,rk8), &
    'nearest neighbour skips a missing source point')
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

  ! From a coarser grid the source points are the corners of the cell
  ! around the target: close to bilinear, also across the date line.
  do j = 1, nlat
    do i = 1, nlon
      g(i,j) = 280.0_rk8 + 0.1_rk8*slat(j) + &
               20.0_rk8*cos(slat(j)*d2r)*cos(2.0_rk8*slon(i)*d2r)
    end do
  end do
  do i = 1, ni
    do j = 1, nj
      tlat(j,i) = 10.0_rkx + 0.37_rkx*real(i-1,rkx)
      tlon(j,i) = 178.03_rkx + 0.11_rkx*real(j-1,rkx)
      if ( tlon(j,i) >= 180.0_rkx ) tlon(j,i) = tlon(j,i) - 360.0_rkx
    end do
  end do
  call h_interpolator_create(hi,slat,slon,tlat,tlon)
  call h_interpolate_cont(hi,g,f)
  do i = 1, ni
    do j = 1, nj
      ii = int(tlat(j,i)+90.0_rkx) + 1
      xl = real(tlon(j,i),rk8) + 180.0_rk8
      jj = int(xl) + 1
      jp = jj + 1
      if ( jp > nlon ) jp = 1
      wy = real(tlat(j,i),rk8) + 90.0_rk8 - real(ii-1,rk8)
      wx = xl - real(jj-1,rk8)
      ex(j,i) = (1.0_rk8-wx)*(1.0_rk8-wy)*g(jj,ii) + &
                wx*(1.0_rk8-wy)*g(jp,ii) + &
                (1.0_rk8-wx)*wy*g(jj,ii+1) + wx*wy*g(jp,ii+1)
    end do
  end do
  call check(maxval(abs(f-ex)) < 2.0e-3_rk8, &
    'cell corners: bilinear within 2e-3 K across the date line')
  call h_interpolator_destroy(hi)

  call test_end('test_kdinterp')
end program test_kdinterp
