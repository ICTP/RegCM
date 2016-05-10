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

module mod_block

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_message
  use mod_constants

  real(rkx) :: grdlnmn , grdltmn , grdlnma , grdltma
  real(rkx) :: xmaxlat , xmaxlon , xminlat , xminlon
  integer(ik4) :: nlatin , nlonin
  logical :: lonwrap , lcrosstime

  contains

  subroutine mxmnll(jx,iy,xlon,xlat,iband)
    implicit none
    integer(ik4) , intent(in) :: iy , jx , iband
    real(rkx) , dimension(jx,iy) , intent(in) :: xlat , xlon

    real(rkx) :: xtstlon1 , xtstlon2
    real(rkx) , parameter :: coord_eps = 0.001_rkx
    !
    ! PURPOSE : FINDS THE MAXIMUM AND MINIMUM LATITUDE AND LONGITUDE
    !
    xminlat = floor(minval(xlat))
    xmaxlat = ceiling(maxval(xlat))

    if ( iband.eq.1 ) then
      xminlon  = -deg180
      xmaxlon  =  deg180
      xtstlon1 = xminlon
      xtstlon2 = xmaxlon
    else if ( abs(xminlat+deg90) < coord_eps .or. &
              abs(xmaxlat-deg90) < coord_eps ) then
      xminlon  = -deg180
      xmaxlon  =  deg180
      xtstlon1 = xminlon
      xtstlon2 = xmaxlon
    else
      xminlon  = floor(minval(xlon(1,:)))
      xmaxlon  = ceiling(maxval(xlon(jx,:)))
      xtstlon1 = floor(maxval(xlon(1,:)))
      xtstlon2 = ceiling(minval(xlon(jx,:)))
    end if

    if ( abs(xminlon-xmaxlon) < coord_eps ) then
      xminlon  = -deg180
      xmaxlon  =  deg180
      xtstlon1 = xminlon
      xtstlon2 = xmaxlon
    end if

    lonwrap = .false.
    lcrosstime = .false.
    if ( (xmaxlon-xminlon) > (deg360-coord_eps) ) then
      lonwrap = .true.
      write(stdout,*) 'Special case for longitude wrapping'
    end if
    if ( abs(xminlon - xtstlon1) > deg180 .or.   &
         abs(xmaxlon - xtstlon2) > deg180 .or.   &
         xminlon > deg00 .and. xmaxlon < deg00 ) then
      lcrosstime = .true.
      if ( xminlon < deg00 .and. xtstlon1 > deg00 ) xminlon = xtstlon1
      if ( xmaxlon > deg00 .and. xtstlon2 < deg00 ) xmaxlon = xtstlon2
      write(stdout,*) 'Special case for timeline crossing'
    end if

    write(stdout,*) 'Calculated large extrema:'
    write(stdout,*) '         MINLAT = ', xminlat
    write(stdout,*) '         MAXLAT = ', xmaxlat
    write(stdout,*) '         MINLON = ', xminlon
    write(stdout,*) '         MAXLON = ', xmaxlon

  end subroutine mxmnll

end module mod_block
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
