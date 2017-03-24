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

  interface ll2xyz
    module procedure ll2xyz_values
    module procedure ll2xyz_array
    module procedure ll2xyz_arrays
  end interface

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

  real(rkx) function longitude_circle(lat) result(er)
    implicit none
    real(rkx) , intent(in) :: lat
    er = d_two * mathpi * erkm * cos(lat*degrad)
  end function longitude_circle

end module mod_earth

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
