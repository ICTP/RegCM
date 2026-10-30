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

! Unit test of Share/mod_intldtr.F90 (interpolation of the terrain and
! surface data on regular grids): bilinear and bicubic exact for low
! order polynomials, distance weighting of a constant, the most frequent
! class with the water threshold for land use and soil texture, class
! percentages, a box across the longitude wrap of a global grid, the
! removal of small water puddles and the domain extrema of mxmnll.

program test_intldtr
  use mod_intkinds
  use mod_realkinds
  use mod_intldtr
  use share_test_utils
  implicit none
  integer(ik4), parameter :: m = 30, n = 20, jx = 7, iy = 5
  real(rkx), dimension(m,n) :: imt
  real(rkx), dimension(jx,iy) :: xlat, xlon, omt, ex
  real(rkx), dimension(36,18) :: gmt
  real(rkx), dimension(3,3) :: glat, glon, gout
  real(rkx), dimension(10,8) :: lu
  real(rkx), dimension(12,10) :: dlat, dlon
  integer(ik4) :: i, j, k
  character(len=40) :: msg

  ! A regular source grid of 1 degree, indices 1..m and 1..n.
  call set_grid(m,n,0.0_rk8,real(m-1,rk8),0.0_rk8,real(n-1,rk8),.false.)
  do j = 1, iy
    do i = 1, jx
      xlon(i,j) = 5.3_rkx + 2.1_rkx*real(i,rkx)
      xlat(i,j) = 3.7_rkx + 2.3_rkx*real(j,rkx)
    end do
  end do
  do k = 1, 4
    do j = 1, n
      do i = 1, m
        imt(i,j) = f(real(i-1,rkx),real(j-1,rkx),k)
      end do
    end do
    ex = f(xlon,xlat,k)
    if ( k <= 2 ) then
      call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,1)
      write(msg,'(a,i0)') 'bilinear exact, field ', k
      call check(maxval(abs(omt-ex)) < 1.0e-10_rkx,trim(msg))
    end if
    call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,2)
    write(msg,'(a,i0)') 'bicubic exact, field ', k
    call check(maxval(abs(omt-ex)) < 1.0e-10_rkx,trim(msg))
  end do
  imt = 7.5_rkx
  call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,6)
  call check(maxval(abs(omt-7.5_rkx)) < 1.0e-10_rkx, &
    'distance weighting of a constant')
  call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,3)
  call check(all(omt == 7.5_rkx),'nearest point')

  ! Classes: a box of 3 ocean (15) and 1 forest (6) points (75 % water) is
  ! ocean with h2opct = 50, and land with h2opct = 80 (land below 80 %).
  imt = 6.0_rkx
  do j = 1, n
    do i = 1, m
      if ( mod(i,2) == 0 .or. mod(j,2) == 0 ) imt(i,j) = 15.0_rkx
    end do
  end do
  call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,4,ibnty=1, &
              h2opct=50.0_rkx)
  call check(all(omt(2:jx-1,2:iy-1) == 15.0_rkx), &
    'land use: 75 % water is water at 50 %')
  call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,4,ibnty=1, &
              h2opct=80.0_rkx)
  call check(all(omt(2:jx-1,2:iy-1) == 6.0_rkx), &
    'land use: 75 % water is land at 80 %')
  ! Soil texture: water is class 14, bedrock 15 stays.
  where ( imt == 15.0_rkx ) imt = 14.0_rkx
  call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,4,ibnty=2, &
              h2opct=80.0_rkx)
  call check(all(omt(2:jx-1,2:iy-1) == 6.0_rkx), &
    'texture: water dropped below the threshold')
  where ( imt == 14.0_rkx ) imt = 15.0_rkx
  call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,4,ibnty=2, &
              h2opct=80.0_rkx)
  call check(all(omt(2:jx-1,2:iy-1) == 15.0_rkx), &
    'texture: bedrock is not water')
  call interp(50.0_rkx,jx,iy,xlat,xlon,omt,imt,5,ival=15)
  call check(all(abs(omt(2:jx-1,2:iy-1)-75.0_rkx) < 1.0e-10_rkx), &
    'percentage of a class')

  ! Global grid of 10 degrees with longitude wrap; the box around
  ! longitude 0 (8 columns) takes points on both sides: 5 water columns
  ! (330-360 and 0 degrees) and 3 land columns (10-30 degrees).
  call set_grid(36,18,0.0_rk8,350.0_rk8,-85.0_rk8,85.0_rk8,.true.)
  gmt = 6.0_rkx
  gmt(33:36,:) = 15.0_rkx
  gmt(1,:) = 15.0_rkx
  do j = 1, 3
    do i = 1, 3
      glon(i,j) = -2.0_rkx + 2.0_rkx*real(i-1,rkx)
      glat(i,j) = -10.0_rkx + 10.0_rkx*real(j-1,rkx)
    end do
  end do
  call interp(7000.0_rkx,3,3,glat,glon,gout,gmt,4,ibnty=1, &
              h2opct=50.0_rkx)
  call check(gout(2,2) == 15.0_rkx,'class box across the longitude wrap')
  do j = 1, 18
    do i = 1, 36
      gmt(i,j) = real(cos(real(i-1,rk8)*10.0_rk8*acos(-1.0_rk8)/180.0_rk8), &
                      rkx)
    end do
  end do
  glon = 355.0_rkx
  call interp(50.0_rkx,3,3,glat,glon,gout,gmt,1)
  call check(maxval(abs(gout-real(0.5_rk8*(cos(350.0_rk8*acos(-1.0_rk8)/ &
             180.0_rk8)+1.0_rk8),rkx))) < 1.0e-10_rkx, &
    'bilinear across the longitude wrap')

  ! Puddles: an isolated water point becomes the most frequent neighbour.
  lu = 6.0_rkx
  lu(3,3) = 15.0_rkx
  lu(7:10,4:8) = 15.0_rkx
  lu(2,2:4) = 2.0_rkx
  lu(3,2) = 2.0_rkx
  lu(4,3) = 2.0_rkx
  call puddlefilter(10,8,lu)
  call check(lu(3,3) == 2.0_rkx .and. all(lu(7:10,4:8) == 15.0_rkx), &
    'small water puddles removed, large ones kept')

  ! Domain extrema: a domain across the dateline.
  do j = 1, 10
    do i = 1, 12
      dlat(i,j) = 30.0_rkx + real(j,rkx)
      dlon(i,j) = 170.0_rkx + 2.0_rkx*real(i,rkx)
      if ( dlon(i,j) > 180.0_rkx ) dlon(i,j) = dlon(i,j) - 360.0_rkx
    end do
  end do
  call mxmnll(12,10,dlon,dlat,0)
  call check(lcrosstime .and. .not. lonwrap .and. xminlat == 31.0_rk8 &
             .and. xmaxlat == 40.0_rk8,'domain across the dateline')

  call test_end('test_intldtr')

  contains

  subroutine set_grid(nx,ny,lon1,lon2,lat1,lat2,lwrap)
    implicit none
    integer(ik4), intent(in) :: nx, ny
    real(rk8), intent(in) :: lon1, lon2, lat1, lat2
    logical, intent(in) :: lwrap
    nlonin = nx
    nlatin = ny
    grdlnmn = lon1
    grdlnma = lon2
    grdltmn = lat1
    grdltma = lat2
    lonwrap = lwrap
    lcrosstime = .false.
  end subroutine set_grid

  elemental real(rkx) function f(x,y,k)
    implicit none
    real(rkx), intent(in) :: x, y
    integer(ik4), intent(in) :: k
    select case (k)
      case (1)
        f = 3.0_rkx*x + 2.0_rkx*y + 1.0_rkx
      case (2)
        f = x*y
      case (3)
        f = x*x
      case default
        f = y*y
    end select
  end function f

end program test_intldtr
