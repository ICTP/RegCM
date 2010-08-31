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
      implicit none
      integer :: nobs
      integer :: nnc
      integer :: iblk , iter, jter
      real(4) , allocatable , dimension(:) :: ht , ht2 , htsd
      real(8) , allocatable , dimension(:) :: xobs , yobs
      real(4) , dimension(50) :: stores
      real(4) :: xmaxlat , xmaxlon , xminlat , xminlon
      real(4) :: dsinm , rin , xnc
      real(8) :: grdlnmn , grdltmn , grdlnma , grdltma
      logical :: lonwrap , lcrosstime

      contains

      subroutine allocate_block(ni, nj)
        implicit none
        integer , intent(in) :: ni , nj
        iter = ni
        jter = nj
        iblk = (iter*jter)/2
        allocate(ht(iblk))
        allocate(ht2(iblk))
        allocate(htsd(iblk))
        allocate(xobs(iblk))
        allocate(yobs(iblk))
      end subroutine

      subroutine free_block
        iter = 0
        jter = 0
        iblk = 0
        if (allocated(ht)) deallocate(ht)
        if (allocated(ht2)) deallocate(ht2)
        if (allocated(htsd)) deallocate(htsd)
        if (allocated(xobs)) deallocate(xobs)
        if (allocated(yobs)) deallocate(yobs)
      end subroutine

      subroutine mxmnll(iy,jx,xlon,xlat)
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx
      real(4) , dimension(iy,jx) :: xlat , xlon
      intent (in) iy , jx , xlat , xlon
!
!     PURPOSE : FINDS THE MAXIMUM AND MINIMUM LATITUDE AND LONGITUDE
!
      xminlat = floor(minval(xlat))
      xmaxlat = ceiling(maxval(xlat))
      if (abs(xminlat+90.0)<0.0001 .or. abs(xmaxlat-90.0)<0.001) then
        xminlon = -180.0
        xmaxlon =  180.0
      else
        xminlon = floor(minval(xlon(:,1)))
        xmaxlon = ceiling(maxval(xlon(:,jx)))
      end if
 
      print *, 'Calculated large extrema:'
      print *, '         MINLAT = ', xminlat
      print *, '         MAXLAT = ', xmaxlat
      print *, '         MINLON = ', xminlon
      print *, '         MAXLON = ', xmaxlon

      lonwrap = .false.
      lcrosstime = .false.
      if ((xmaxlon-xminlon) > 359.99) then
        lonwrap = .true.
        print *, 'Special case for longitude wrapping'
      end if
      if (xminlon > 0.0 .and. xmaxlon < 0.0) then
        lcrosstime = .true.
        print *, 'Special case for timeline crossing'
      end if

!--------initialize minimum lat and lon of data from tape

      grdltmn = xminlat + 5.0
      grdltma = xmaxlat - 5.0
      grdlnmn = xminlon + 5.0
      grdlnma = xmaxlon - 5.0

      end subroutine mxmnll

      end module mod_block
