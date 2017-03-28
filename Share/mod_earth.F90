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

  real(rkx) , parameter :: mindist = 1.0e-6_rkx

  public :: gcdist_simple , gcdist
  public :: ll2xyz
  public :: global_domain , get_window

  interface ll2xyz
    module procedure ll2xyz_values
    module procedure ll2xyz_array
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
    real(rkx) , intent(in) :: lat1 , lon1 , lat2, lon2
    real(rkx) :: clat1 , slat1 , clat2 , slat2 , cdlon , crd
    clat1 = cos(lat1*degrad)
    slat1 = sin(lat1*degrad)
    clat2 = cos(lat2*degrad)
    slat2 = sin(lat2*degrad)
    cdlon = cos((lon1-lon2)*degrad)
    crd   = slat1*slat2+clat1*clat2*cdlon
    ! Have it in km to avoid numerical problems :)
    gcdist_simple = erkm*acos(crd)
  end function gcdist_simple

  real(rkx) function gcdist(ds,lat1,lon1,lat2,lon2)
    implicit none
    real(rkx) , intent(in) :: ds , lat1 , lon1 , lat2, lon2
    real(rkx) :: clat1 , slat1 , clat2 , slat2 , cdlon , sdlon
    real(rkx) :: y , x
    clat1 = cos(lat1*degrad)
    slat1 = sin(lat1*degrad)
    clat2 = cos(lat2*degrad)
    slat2 = sin(lat2*degrad)
    cdlon = cos((lon1-lon2)*degrad)
    sdlon = sin((lon1-lon2)*degrad)
    y = sqrt((clat2*sdlon)**2+(clat1*slat2-slat1*clat2*cdlon)**2)
    x = slat1*slat2+clat1*clat2*cdlon
    gcdist = max(erkm*atan2(y,x)/ds,mindist)
  end function gcdist

  subroutine ll2xyz_values(lat,lon,x,y,z)
    implicit none
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: x , y , z
    real(rkx) :: rlat , rlon
    rlat = halfpi - lat*degrad
    rlon = lon*degrad
    x = sin(rlat) * sin(rlon)
    y = cos(rlat)
    z = sin(rlat) * cos(rlon)
  end subroutine ll2xyz_values

  subroutine ll2xyz_array(lat,lon,x)
    implicit none
    real(rkx) , intent(in) :: lat
    real(rkx) , intent(in) :: lon
    real(rkx) , intent(out) , dimension(3) :: x
    real(rkx) :: rlat , rlon
    rlat = halfpi - lat*degrad
    rlon = lon*degrad
    x(1) = sin(rlat) * sin(rlon)
    x(2) = cos(rlat)
    x(3) = sin(rlat) * cos(rlon)
  end subroutine ll2xyz_array

  subroutine ll2xyz_arrays(lat,lon,x)
    implicit none
    real(rkx) , intent(in) , dimension(:) :: lat
    real(rkx) , intent(in) , dimension(:) :: lon
    real(rkx) , intent(out) , dimension(:,:) :: x
    real(rkx) :: rlat , rlon
    integer(ik4) :: i , j , n
    if ( size(x,1) /= 3 .or.  size(x,2) /= size(lat) * size(lon) ) then
      return
    end if
    n = 1
    do j = 1 , size(lat)
      do i = 1 , size(lon)
        rlat = halfpi - lat(j)*degrad
        rlon = lon(i)*degrad
        x(1,n) = sin(rlat) * sin(rlon)
        x(2,n) = cos(rlat)
        x(3,n) = sin(rlat) * cos(rlon)
        n = n + 1
      end do
    end do
  end subroutine ll2xyz_arrays

  subroutine ll2xyz_grid(lat,lon,x)
    implicit none
    real(rkx) , intent(in) , dimension(:,:) :: lat
    real(rkx) , intent(in) , dimension(:,:) :: lon
    real(rkx) , intent(out) , dimension(:,:) :: x
    real(rkx) :: rlat , rlon
    integer(ik4) :: i , j , n
    if ( size(x,1) /= 3 .or.  size(x,2) /= product(shape(lat)) .or. &
         product(shape(lat)) /= product(shape(lon))  ) then
      return
    end if
    n = 1
    do j = 1 , size(lat,2)
      do i = 1 , size(lat,1)
        rlat = halfpi - lat(i,j)*degrad
        rlon = lon(i,j)*degrad
        x(1,n) = sin(rlat) * sin(rlon)
        x(2,n) = cos(rlat)
        x(3,n) = sin(rlat) * cos(rlon)
        n = n + 1
      end do
    end do
  end subroutine ll2xyz_grid

  real(rkx) function longitude_circle(lat) result(er)
    implicit none
    real(rkx) , intent(in) :: lat
    er = d_two * mathpi * erkm * cos(lat*degrad)
  end function longitude_circle

  ! Assumptions :
  !  GLAT is regular ordered array of latitudes in a GLOBAL grid
  !  GLON is regular ordered array of longitudes in a GLOBAL grid
  !  XLAT is local grid latitudes going South to North
  !  XLON is local grid longitudes  going West to East
  ! ALL LONGITUDES ARE in the range -180.0 <-> 180.0
  subroutine get_window(glat,glon,xlat,xlon,i_band,domain)
    implicit none
    real(rkx) , dimension(:) , intent(in) :: glat
    real(rkx) , dimension(:) , intent(in) :: glon
    real(rkx) , dimension(:,:) , intent(in) :: xlat
    real(rkx) , dimension(:,:) , intent(in) :: xlon
    integer(ik4) , intent(in) :: i_band
    type(global_domain) , intent(out) :: domain

    real(rkx) :: dlat , dlon
    real(rkx) , allocatable , dimension(:,:) :: xlon360
    real(rkx) :: maxlat
    real(rkx) :: minlat
    real(rkx) :: maxlon
    real(rkx) :: minlon
    integer :: gi , gj , xi , xj , l1 , l2 , i , j

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
    else if ( glon(gi) > 350.0_rkx ) then
      ! Input data is     0 : 360 , xlon is -180 : 180
      allocate(xlon360(xi,xj))
      xlon360 = xlon
      where ( xlon < 0.0_rkx )
        xlon360 = 360.0_rkx + xlon
      end where
      if ( minval(xlon360(1,:)) < maxval(xlon360(xi,:)) ) then
        domain%ntiles = 1
        minlon = minval(xlon360)
        maxlon = maxval(xlon360)
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
      else
        ! Cross Greenwich line
        minlon = minval(xlon360(1,:))
        maxlon = maxval(xlon(xi,:))
        domain%ntiles = 2
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = gi
        domain%igstart(2) = 1
        domain%igstop(2) = int((maxlon-glon(1))/dlon) + 3
      end if
      deallocate(xlon360)
    else
      ! Input Data is -180 : 180 , xlon is -180 : 180
      if ( xlon(1,xj)   <= xlon(xi,xj)   .and. &
           xlon(1,xj/2) <= xlon(xi,xj/2) .and. &
           xlon(1,1)    <= xlon(xi,1) ) then
        ! it is not crossing timeline
        domain%ntiles = 1
        domain%igstart(1) = int((minlon-glon(1))/dlon) - 2
        domain%igstop(1) = int((maxlon-glon(1))/dlon) + 3
        domain%igstart(2) = 0
        domain%igstop(2) = 0
      else
        domain%ntiles = 2
        minlon = 180.0_rkx
        do j = 1 , xj
          if ( xlon(1,j) > 0.0_rkx ) minlon = min(minlon,xlon(1,j))
        end do
        maxlon = -180.0_rkx
        do j = 1 , xj
          if ( xlon(xi,j) < 0.0_rkx ) maxlon = max(maxlon,xlon(xi,j))
        end do
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
        domain%jgstart = min(l1,l2)
        domain%jgstop = gj
      else
        l1 = int((maxval(glat(1)-xlat(:,1)))/dlat) + 3
        l2 = int((maxval(glat(1)-xlat(:,xj)))/dlat) + 3
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
        domain%jgstart = 1
        domain%jgstop = max(l1,l2)
      else
        l1 = int((maxval(glat(1)-xlat(:,1)))/dlat) - 2
        l2 = int((maxval(glat(1)-xlat(:,xj)))/dlat) - 2
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

    domain%ni = 0
    do i = 1 , domain%ntiles
      domain%igstart(i) = min(max(1,domain%igstart(i)),gi)
      domain%igstop(i) = max(min(domain%igstop(i),gi),1)
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

  end subroutine get_window

end module mod_earth

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
