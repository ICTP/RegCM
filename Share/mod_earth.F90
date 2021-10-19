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

module mod_earth

  use mod_realkinds
  use mod_intkinds
  use mod_constants

  private

  real(rk8) , parameter :: mindist = 1.0e-6_rk8

  public :: gcdist_simple , gcdist
  public :: ll2xyz , longitude_circle
  public :: global_domain , get_window

  interface get_window
    module procedure get_window_r4
    module procedure get_window_r8
  end interface get_window

  interface ll2xyz
    module procedure ll2xyz_values
    module procedure ll2xyz_array
    module procedure ll2xyz_1d
    module procedure ll2xyz_arrays
    module procedure ll2xyz_grid
  end interface

  type global_domain
    integer(ik4) :: global_ni
    integer(ik4) :: global_nj
    integer(ik4) :: ntiles
    integer(ik4) , dimension(2) :: ni
    integer(ik4) , dimension(2) :: igstart
    integer(ik4) , dimension(2) :: igstop
    integer(ik4) :: nj
    integer(ik4) :: jgstart
    integer(ik4) :: jgstop
  end type global_domain

  contains

  real(rkx) function gcdist_simple(lat1,lon1,lat2,lon2)
    implicit none
    real(rkx) , intent(in) :: lat1 , lon1 , lat2 , lon2
    real(rk8) :: clat1 , slat1 , clat2 , slat2 , cdlon , crd
    clat1 = cos(lat1*degrad)
    slat1 = sin(lat1*degrad)
    clat2 = cos(lat2*degrad)
    slat2 = sin(lat2*degrad)
    cdlon = cos((lon1-lon2)*degrad)
    crd   = slat1*slat2+clat1*clat2*cdlon
    ! Have it in km to avoid numerical problems :)
    gcdist_simple = sngl(erkm*acos(crd))
  end function gcdist_simple

  real(rkx) function gcdist(ds,lat1,lon1,lat2,lon2)
    implicit none
    real(rkx) , intent(in) :: ds , lat1 , lon1 , lat2 , lon2
    real(rk8) :: clat1 , slat1 , clat2 , slat2 , cdlon , sdlon
    real(rk8) :: y , x
    clat1 = cos(lat1*degrad)
    slat1 = sin(lat1*degrad)
    clat2 = cos(lat2*degrad)
    slat2 = sin(lat2*degrad)
    cdlon = cos((lon1-lon2)*degrad)
    sdlon = sin((lon1-lon2)*degrad)
    y = sqrt((clat2*sdlon)**2+(clat1*slat2-slat1*clat2*cdlon)**2)
    x = slat1*slat2+clat1*clat2*cdlon
    gcdist = sngl(max(erkm*atan2(y,x)/ds,mindist))
  end function gcdist

  subroutine ll2xyz_values(lat,lon,x,y,z)
    implicit none
    real(rkx) , intent(in) :: lat , lon
    real(rk8) , intent(out) :: x , y , z
    real(rk8) :: rlat , rlon
    rlat = max(min(dble(lat),89.999_rk8),-89.999_rk8)*degrad
    rlon = dble(lon)*degrad
    x = cos(rlat) * sin(rlon)
    y = sin(rlat)
    z = cos(rlat) * cos(rlon)
  end subroutine ll2xyz_values

  subroutine ll2xyz_array(lat,lon,x)
    implicit none
    real(rkx) , intent(in) :: lat
    real(rkx) , intent(in) :: lon
    real(rk8) , intent(out) , dimension(3) :: x
    real(rk8) :: rlat , rlon
    rlat = max(min(dble(lat),89.999_rk8),-89.999_rk8)*degrad
    rlon = dble(lon)*degrad
    x(1) = cos(rlat) * sin(rlon)
    x(2) = sin(rlat)
    x(3) = cos(rlat) * cos(rlon)
  end subroutine ll2xyz_array

  subroutine ll2xyz_1d(ni,lat,lon,x)
    implicit none
    integer(ik4) , intent(in) :: ni
    real(rkx) , intent(in) , dimension(ni) :: lat
    real(rkx) , intent(in) , dimension(ni) :: lon
    real(rk8) , intent(out) , dimension(3,ni) :: x
    real(rk8) :: rlat , rlon
    integer(ik4) :: i
    do i = 1 , ni
      rlat = max(min(dble(lat(i)),89.999_rk8),-89.999_rk8)*degrad
      rlon = dble(lon(i))*degrad
      x(1,i) = cos(rlat) * sin(rlon)
      x(2,i) = sin(rlat)
      x(3,i) = cos(rlat) * cos(rlon)
    end do
  end subroutine ll2xyz_1d

  subroutine ll2xyz_arrays(lat,lon,x)
    implicit none
    real(rkx) , intent(in) , dimension(:) :: lat
    real(rkx) , intent(in) , dimension(:) :: lon
    real(rk8) , intent(out) , dimension(:,:) :: x
    real(rk8) :: rlat , rlon
    integer(ik4) :: i , j , n
    if ( size(x,1) /= 3 .or.  size(x,2) /= size(lat) * size(lon) ) then
      return
    end if
    n = 1
    do j = 1 , size(lat)
      do i = 1 , size(lon)
        rlat = max(min(dble(lat(j)),89.999_rk8),-89.999_rk8)*degrad
        rlon = dble(lon(i))*degrad
        x(1,n) = cos(rlat) * sin(rlon)
        x(2,n) = sin(rlat)
        x(3,n) = cos(rlat) * cos(rlon)
        n = n + 1
      end do
    end do
  end subroutine ll2xyz_arrays

  subroutine ll2xyz_grid(lat,lon,x)
    implicit none
    real(rkx) , intent(in) , dimension(:,:) :: lat
    real(rkx) , intent(in) , dimension(:,:) :: lon
    real(rk8) , intent(out) , dimension(:,:) :: x
    real(rk8) :: rlat , rlon
    integer(ik4) :: i , j , n
    if ( size(x,1) /= 3 .or.  size(x,2) /= product(shape(lat)) .or. &
         product(shape(lat)) /= product(shape(lon))  ) then
      return
    end if
    n = 1
    do j = 1 , size(lat,2)
      do i = 1 , size(lat,1)
        rlat = max(min(dble(lat(i,j)),89.999_rk8),-89.999_rk8)*degrad
        rlon = dble(lon(i,j))*degrad
        x(1,n) = cos(rlat) * sin(rlon)
        x(2,n) = sin(rlat)
        x(3,n) = cos(rlat) * cos(rlon)
        n = n + 1
      end do
    end do
  end subroutine ll2xyz_grid

  real(rk8) function longitude_circle(lat) result(er)
    implicit none
    real(rk8) , intent(in) :: lat
    er = d_two * mathpi * erkm * cos(lat*degrad)
  end function longitude_circle

  ! Assumptions :
  !  GLAT is regular ordered array of latitudes in a GLOBAL grid
  !  GLON is regular ordered array of longitudes in a GLOBAL grid
  !  XLAT is local grid latitudes going South to North
  !  XLON is local grid longitudes  going West to East
  ! ALL LONGITUDES ARE in the range -180.0 <-> 180.0
  subroutine get_window_r8(glat,glon,xlat,xlon,i_band,domain)
    implicit none
    real(rk8) , dimension(:) , intent(in) :: glat
    real(rk8) , dimension(:) , intent(in) :: glon
    real(rkx) , dimension(:,:) , intent(in) :: xlat
    real(rkx) , dimension(:,:) , intent(in) :: xlon
    integer(ik4) , intent(in) :: i_band
    type(global_domain) , intent(out) :: domain

    real(rk8) :: dlat , dlon
    real(rk8) , allocatable , dimension(:,:) :: xlon360
    real(rk8) :: maxlat , minlat , maxlon , minlon , d1 , d2
    integer :: gi , gj , xi , xj , l1 , l2 , i , j , itmp

    xi = size(xlon,1)
    xj = size(xlat,2)
    maxlat = maxval(xlat)
    minlat = minval(xlat)
    maxlon = maxval(xlon)
    minlon = minval(xlon)

    gi = size(glon)
    gj = size(glat)
    dlat = abs(glat(2) - glat(1))
    dlon = abs(glon(2) - glon(1))
    domain%global_ni = gi
    domain%global_nj = gj

    if ( i_band == 1 ) then
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else if ( glon(gi) > 350.0_rk8 ) then
      ! Input data is     0 : 360 , xlon is -180 : 180
      if ( (any(xlon(1,:) < 0.0_rk8) .and. any(xlon(1,:) > 0.0_rk8)) .or. &
           (any(xlon(xi,:) < 0.0_rk8) .and. any(xlon(xi,:) > 0.0_rk8)) .or. &
           (all(xlon(1,:) < 0.0_rk8) .and. all(xlon(xi,:) > 0.0_rk8)) ) then
        ! Cross Greenwich line
        minlon = minval(xlon(1,:))
        maxlon = maxval(xlon(xi,:))
        if ( minlon < 0.0_rk8 ) minlon = minlon + 360.0_rk8
        if ( maxlon < 0.0_rk8 ) maxlon = maxlon + 360.0_rk8
        domain%ntiles = 2
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      else
        allocate(xlon360(xi,xj))
        xlon360 = xlon
        where ( xlon < 0.0_rk8 )
          xlon360 = 360.0_rk8 + xlon
        end where
        domain%ntiles = 1
        minlon = minval(xlon360)
        maxlon = maxval(xlon360)
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
        deallocate(xlon360)
      end if
    else
      ! Input Data is -180 : 180 , xlon is -180 : 180
      d1 = 0.0_rk8
      d2 = 0.0_rk8
      do j = 1 , xj-1
        d1 = max(dble(abs(xlon(1,j+1)-xlon(1,j))),d1)
      end do
      do j = 1 , xj-1
        d2 = max(dble(abs(xlon(xi,j+1)-xlon(xi,j))),d2)
      end do
      if ( d1 < 350.0_rk8 .and. d2 < 350.0_rk8 .and. &
          .not. (all(xlon(1,:) > 0.0_rk8) .and. all(xlon(xi,:) < 0.0)) ) then
        domain%ntiles = 1
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
      else
        ! it is crossing timeline
        domain%ntiles = 2
        if ( d1 > 350.0_rk8 ) then
          minlon = minval(xlon(1 ,:),(xlon(1 ,:)>0.0_rk8))
          maxlon = maxval(xlon(xi,:))
        else if ( d2 > 350.0_rk8 ) then
          minlon = minval(xlon(1,:))
          maxlon = maxval(xlon(xi,:),(xlon(xi,:)<0.0_rk8))
        else
          maxlon = maxval(xlon(xi,:))
          minlon = minval(xlon(1,:))
        end if
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      end if
    end if
    if ( has_north_pole(xlat,xi/2) ) then
      ! North pole inside
      if ( glat(1) < glat(gj) ) then
        l1 = int((minval(xlat(:,1))-glat(1))/dlat) - 2
        l2 = int((minval(xlat(:,xj))-glat(1))/dlat) - 2
        l1 = min(l1,l2)
        l2 = int((minval(xlat(:,xj/2))-glat(1))/dlat) - 2
        domain%jgstart = min(l1,l2)
        domain%jgstop = gj
      else
        l1 = int((maxval(glat(1)-xlat(:,1)))/dlat) + 3
        l2 = int((maxval(glat(1)-xlat(:,xj)))/dlat) + 3
        l1 = max(l1,l2)
        l2 = int((maxval(glat(1)-xlat(:,xj/2)))/dlat) + 3
        domain%jgstart = 1
        domain%jgstop = max(l1,l2)
      end if
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else if ( has_south_pole(xlat,xi/2) ) then
      ! South Pole inside
      if ( glat(1) < glat(gj) ) then
        l1 = int((maxval(xlat(:,1))-glat(1))/dlat) + 3
        l2 = int((maxval(xlat(:,xj))-glat(1))/dlat) + 3
        l1 = max(l1,l2)
        l2 = int((maxval(xlat(:,xj/2))-glat(1))/dlat) + 3
        domain%jgstart = 1
        domain%jgstop = max(l1,l2)
      else
        l1 = int((minval(glat(1)-xlat(:,1)))/dlat) - 2
        l2 = int((minval(glat(1)-xlat(:,xj)))/dlat) - 2
        l1 = min(l1,l2)
        l2 = int((minval(glat(1)-xlat(:,xj/2)))/dlat) - 2
        domain%jgstart = min(l1,l2)
        domain%jgstop = gj
      end if
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else
      if ( glat(1) < glat(gj) ) then
        domain%jgstart = int((minlat-glat(1))/dlat) - 2
        domain%jgstop  = int((maxlat-glat(1))/dlat) + 3
      else
        domain%jgstart = int((minlat-glat(gj))/dlat) - 2
        domain%jgstop  = int((maxlat-glat(gj))/dlat) + 3
        domain%nj =  domain%jgstop - domain%jgstart + 1
        domain%jgstart = (gj+1)-(domain%jgstart+domain%nj-1)
        domain%jgstop = domain%jgstart + domain%nj - 1
      end if
    end if

    if ( domain%igstart(1) < 1 .and. domain%ntiles == 1 ) then
      domain%ntiles = 2
      itmp = domain%igstop(1)
      domain%igstart(1) = gi + domain%igstart(1) - 1
      domain%igstop(1) = gi
      domain%igstart(2) = 1
      domain%igstop(2) = itmp
    end if
    if ( domain%igstop(1) > gi .and. domain%ntiles == 1 ) then
      domain%ntiles = 2
      itmp = domain%igstop(1) - gi
      domain%igstart(1) = domain%igstart(1)
      domain%igstop(1) = gi
      domain%igstart(2) = 1
      domain%igstop(2) = itmp
    end if
    domain%ni = 0
    do i = 1 , domain%ntiles
      domain%ni(i) =  domain%igstop(i) - domain%igstart(i) + 1
    end do
    domain%jgstart = min(max(1,domain%jgstart),gj)
    domain%jgstop = max(min(domain%jgstop,gj),1)
    domain%nj =  domain%jgstop - domain%jgstart + 1

    contains

      logical function has_north_pole(l,i)
        real(rkx) , intent(in) , dimension(:,:) :: l
        integer(ik4) , intent(in) :: i
        integer(ik4) :: j
        has_north_pole = .false.
        if ( all(l(i,:) < 0.0) ) return
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_north_pole = .true.
            exit
          end if
        end do
      end function has_north_pole

      logical function has_south_pole(l,i)
        real(rkx) , intent(in) , dimension(:,:) :: l
        integer(ik4) , intent(in) :: i
        integer(ik4) :: j
        has_south_pole = .false.
        if ( all(l(i,:) > 0.0) ) return
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_south_pole = .true.
            exit
          end if
        end do
      end function has_south_pole

  end subroutine get_window_r8

  subroutine get_window_r4(glat,glon,xlat,xlon,i_band,domain)
    implicit none
    real(rk4) , dimension(:) , intent(in) :: glat
    real(rk4) , dimension(:) , intent(in) :: glon
    real(rkx) , dimension(:,:) , intent(in) :: xlat
    real(rkx) , dimension(:,:) , intent(in) :: xlon
    integer(ik4) , intent(in) :: i_band
    type(global_domain) , intent(out) :: domain

    real(rk8) :: dlat , dlon
    real(rk8) , allocatable , dimension(:,:) :: xlon360
    real(rk8) :: maxlat , minlat , maxlon , minlon , d1 , d2
    integer :: gi , gj , xi , xj , l1 , l2 , i , j , itmp

    xi = size(xlon,1)
    xj = size(xlat,2)
    maxlat = maxval(xlat)
    minlat = minval(xlat)
    maxlon = maxval(xlon)
    minlon = minval(xlon)

    gi = size(glon)
    gj = size(glat)
    dlat = abs(glat(2) - glat(1))
    dlon = abs(glon(2) - glon(1))
    domain%global_ni = gi
    domain%global_nj = gj

    if ( i_band == 1 ) then
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else if ( glon(gi) > 350.0_rk8 ) then
      ! Input data is     0 : 360 , xlon is -180 : 180
      if ( (any(xlon(1,:) < 0.0_rk8) .and. any(xlon(1,:) > 0.0_rk8)) .or. &
           (any(xlon(xi,:) < 0.0_rk8) .and. any(xlon(xi,:) > 0.0_rk8)) .or. &
           (all(xlon(1,:) < 0.0_rk8) .and. all(xlon(xi,:) > 0.0_rk8)) ) then
        ! Cross Greenwich line
        minlon = minval(xlon(1,:))
        maxlon = maxval(xlon(xi,:))
        if ( minlon < 0.0_rk8 ) minlon = minlon + 360.0_rk8
        if ( maxlon < 0.0_rk8 ) maxlon = maxlon + 360.0_rk8
        domain%ntiles = 2
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      else
        allocate(xlon360(xi,xj))
        xlon360 = xlon
        where ( xlon < 0.0_rk8 )
          xlon360 = 360.0_rk8 + xlon
        end where
        domain%ntiles = 1
        minlon = minval(xlon360)
        maxlon = maxval(xlon360)
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
        deallocate(xlon360)
      end if
    else
      ! Input Data is -180 : 180 , xlon is -180 : 180
      d1 = 0.0_rk8
      d2 = 0.0_rk8
      do j = 1 , xj-1
        d1 = max(dble(abs(xlon(1,j+1)-xlon(1,j))),d1)
      end do
      do j = 1 , xj-1
        d2 = max(dble(abs(xlon(xi,j+1)-xlon(xi,j))),d2)
      end do
      if ( d1 < 350.0_rk8 .and. d2 < 350.0_rk8 .and. &
          .not. (all(xlon(1,:) > 0.0_rk8) .and. all(xlon(xi,:) < 0.0)) ) then
        domain%ntiles = 1
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
      else
        ! it is crossing timeline
        domain%ntiles = 2
        if ( d1 > 350.0_rk8 ) then
          minlon = minval(xlon(1 ,:),(xlon(1 ,:)>0.0_rk8))
          maxlon = maxval(xlon(xi,:))
        else if ( d2 > 350.0_rk8 ) then
          minlon = minval(xlon(1,:))
          maxlon = maxval(xlon(xi,:),(xlon(xi,:)<0.0_rk8))
        else
          maxlon = maxval(xlon(xi,:))
          minlon = minval(xlon(1,:))
        end if
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      end if
    end if
    if ( has_north_pole(xlat,xi/2) ) then
      ! North pole inside
      if ( glat(1) < glat(gj) ) then
        l1 = int((minval(xlat(:,1))-glat(1))/dlat) - 2
        l2 = int((minval(xlat(:,xj))-glat(1))/dlat) - 2
        l1 = min(l1,l2)
        l2 = int((minval(xlat(:,xj/2))-glat(1))/dlat) - 2
        domain%jgstart = min(l1,l2)
        domain%jgstop = gj
      else
        l1 = int((maxval(glat(1)-xlat(:,1)))/dlat) + 3
        l2 = int((maxval(glat(1)-xlat(:,xj)))/dlat) + 3
        l1 = max(l1,l2)
        l2 = int((maxval(glat(1)-xlat(:,xj/2)))/dlat) + 3
        domain%jgstart = 1
        domain%jgstop = max(l1,l2)
      end if
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else if ( has_south_pole(xlat,xi/2) ) then
      ! South Pole inside
      if ( glat(1) < glat(gj) ) then
        l1 = int((maxval(xlat(:,1))-glat(1))/dlat) + 3
        l2 = int((maxval(xlat(:,xj))-glat(1))/dlat) + 3
        l1 = max(l1,l2)
        l2 = int((maxval(xlat(:,xj/2))-glat(1))/dlat) + 3
        domain%jgstart = 1
        domain%jgstop = max(l1,l2)
      else
        l1 = int((minval(glat(1)-xlat(:,1)))/dlat) - 2
        l2 = int((minval(glat(1)-xlat(:,xj)))/dlat) - 2
        l1 = min(l1,l2)
        l2 = int((minval(glat(1)-xlat(:,xj/2)))/dlat) - 2
        domain%jgstart = min(l1,l2)
        domain%jgstop = gj
      end if
      domain%ntiles = 1
      domain%igstart(1) = 1
      domain%igstop(1) = gi
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else
      if ( glat(1) < glat(gj) ) then
        domain%jgstart = int((minlat-glat(1))/dlat) - 2
        domain%jgstop  = int((maxlat-glat(1))/dlat) + 3
      else
        domain%jgstart = int((minlat-glat(gj))/dlat) - 2
        domain%jgstop  = int((maxlat-glat(gj))/dlat) + 3
        domain%nj =  domain%jgstop - domain%jgstart + 1
        domain%jgstart = (gj+1)-(domain%jgstart+domain%nj-1)
        domain%jgstop = domain%jgstart + domain%nj - 1
      end if
    end if

    if ( domain%igstart(1) < 1 .and. domain%ntiles == 1 ) then
      domain%ntiles = 2
      itmp = domain%igstop(1)
      domain%igstart(1) = gi + domain%igstart(1) - 1
      domain%igstop(1) = gi
      domain%igstart(2) = 1
      domain%igstop(2) = itmp
    end if
    if ( domain%igstop(1) > gi .and. domain%ntiles == 1 ) then
      domain%ntiles = 2
      itmp = domain%igstop(1) - gi
      domain%igstart(1) = domain%igstart(1)
      domain%igstop(1) = gi
      domain%igstart(2) = 1
      domain%igstop(2) = itmp
    end if
    domain%ni = 0
    do i = 1 , domain%ntiles
      domain%ni(i) =  domain%igstop(i) - domain%igstart(i) + 1
    end do
    domain%jgstart = min(max(1,domain%jgstart),gj)
    domain%jgstop = max(min(domain%jgstop,gj),1)
    domain%nj =  domain%jgstop - domain%jgstart + 1

    contains

      logical function has_north_pole(l,i)
        real(rkx) , intent(in) , dimension(:,:) :: l
        integer(ik4) , intent(in) :: i
        integer(ik4) :: j
        has_north_pole = .false.
        if ( all(l(i,:) < 0.0) ) return
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_north_pole = .true.
            exit
          end if
        end do
      end function has_north_pole

      logical function has_south_pole(l,i)
        real(rkx) , intent(in) , dimension(:,:) :: l
        integer(ik4) , intent(in) :: i
        integer(ik4) :: j
        has_south_pole = .false.
        if ( all(l(i,:) > 0.0) ) return
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_south_pole = .true.
            exit
          end if
        end do
      end function has_south_pole

  end subroutine get_window_r4

end module mod_earth

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
