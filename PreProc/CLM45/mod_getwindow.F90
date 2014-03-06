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
module mod_getwindow

  use mod_realkinds
  use mod_intkinds
  use mod_grid

  implicit none

  private

  type global_domain
    integer(ik4) :: ntiles
    integer(ik4) , dimension(2) :: ni
    integer(ik4) , dimension(2) :: igstart
    integer(ik4) , dimension(2) :: igstop
    integer(ik4) :: nj
    integer(ik4) :: jgstart
    integer(ik4) :: jgstop
  end type global_domain

  public :: global_domain , get_window

  contains

  ! Assumptions :
  !  GLAT is regular ordered array of latitudes in a GLOBAL grid
  !  GLON is regular ordered array of longitudes in a GLOBAL grid
  !  XLAT is local grid latitudes going South to North
  !  XLON is local grid longitudes  going West to East
  ! ALL LONGITUDES ARE in the range -180.0 <-> 180.0
  subroutine get_window(glat,glon,domain)
    implicit none
    real(rk4) , dimension(:) , intent(in) :: glat
    real(rk4) , dimension(:) , intent(in) :: glon
    type(global_domain) , intent(out) :: domain

    real(rk4) :: dlat , dlon
    real(rk4) :: maxlat
    real(rk4) :: minlat
    real(rk4) :: maxlon
    real(rk4) :: minlon
    integer :: gi , gj , xi , xj , l1 , l2 , i

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

    if ( xlon(1,xj)   < xlon(xi,xj)   .and. &
         xlon(1,xj/2) < xlon(xi,xj/2) .and. &
         xlon(1,1)    < xlon(xi,1) ) then
      ! it is not crossing timeline
      domain%ntiles = 1
      domain%igstart(1) = int((minlon-glon(1))/dlon) - 1
      domain%igstop(1) = int((maxlon-glon(1))/dlon) + 2
      domain%igstart(2) = 0
      domain%igstop(2) = 0
    else
      domain%ntiles = 2
      domain%igstart(1) = int((maxlon - glon(1))/dlon) - 1
      domain%igstop(1) = gj
      domain%igstart(2) = 1
      domain%igstop(2) = int((glon(1)-minlon)/dlon) + 2
    end if
    if ( has_north_pole(xlat,xi/2) ) then
      ! North pole inside
      l1 = int((minval(xlat(:,1))-glat(1))/dlat) - 1
      l2 = int((minval(xlat(:,xj))-glat(1))/dlat) - 1
      domain%jgstart = min(l1,l2)
      domain%jgstop = gj
    else if ( has_south_pole(xlat,xi/2) ) then
      ! South Pole inside
      l1 = int((maxval(xlat(:,1))-glat(1))/dlat) + 2
      l2 = int((maxval(xlat(:,xj))-glat(1))/dlat) + 2
      domain%jgstart = 1
      domain%jgstop = max(l1,l2)
    else
      domain%jgstart = int((minlat-glat(1))/dlat) - 1
      domain%jgstop  = int((maxlat-glat(1))/dlat) + 2
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
        real(rk4) , intent(in) , dimension(:,:) :: l
        integer(ik4) , intent(in) :: i
        integer(ik4) :: j
        has_north_pole = .false.
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_north_pole = .true.
            exit
          end if
        end do
      end function has_north_pole

      logical function has_south_pole(l,i)
        real(rk4) , intent(in) , dimension(:,:) :: l
        integer(ik4) , intent(in) :: i
        integer(ik4) :: j
        has_south_pole = .false.
        do j = 2 , size(l,2)
          if ( l(i,j) < l(i,j-1) ) then
            has_south_pole = .true.
            exit
          end if
        end do
      end function has_south_pole

  end subroutine get_window

end module mod_getwindow
