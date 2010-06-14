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
      real(4) , allocatable , dimension(:) :: xobs , yobs
      real(8) , allocatable , dimension(:,:) :: lnd8
      real(4) , dimension(50) :: stores
      real(4) :: grdlnmn , grdltmn
      real(4) :: xmaxlat , xmaxlon , xminlat , xminlon
      real(4) :: dsinm , rin , xn , xnc
      real(4) :: dxcen , dycen

      contains

      subroutine allocate_block(ni, nj)
        implicit none
        integer , intent(in) :: ni , nj
        iter = ni
        jter = nj
        iblk = (iter*jter)/2
        allocate(lnd8(iter,jter))
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
        if (allocated(lnd8)) deallocate(lnd8)
        if (allocated(ht)) deallocate(ht)
        if (allocated(ht2)) deallocate(ht2)
        if (allocated(htsd)) deallocate(htsd)
        if (allocated(xobs)) deallocate(xobs)
        if (allocated(yobs)) deallocate(yobs)
      end subroutine

      subroutine mxmnll(iy,jx,clon,xlon,xlat,ntypec)
      implicit none
!
! Dummy arguments
!
      real(4) :: clon
      integer :: iy , jx , ntypec
      real(4) , dimension(iy,jx) :: xlat , xlon
      intent (in) clon , iy , jx , ntypec , xlat , xlon
!
! Local variables
!
      integer :: i , j
!
!     PURPOSE : FINDS THE MAXIMUM AND MINIMUM LATITUDE AND LONGITUDE
!
      xmaxlat = -90
      xminlat = 90
      xminlon = 999999.
      xmaxlon = -999999.
!
      do i = 1 , iy
        do j = 1 , jx
          xminlat = amin1(xminlat,xlat(i,j))
          xmaxlat = amax1(xmaxlat,xlat(i,j))
        end do
      end do
      do i = 1 , iy
        do j = 1 , jx
          if ( clon>=0.0 ) then
            if ( xlon(i,j)>=0.0 ) then
              xminlon = amin1(xminlon,xlon(i,j))
              xmaxlon = amax1(xmaxlon,xlon(i,j))
            else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)+360.)) )  &
                    & then
              xminlon = amin1(xminlon,xlon(i,j))
              xmaxlon = amax1(xmaxlon,xlon(i,j))
            else
              xminlon = amin1(xminlon,xlon(i,j)+360.)
              xmaxlon = amax1(xmaxlon,xlon(i,j)+360.)
            end if
          else if ( xlon(i,j)<0.0 ) then
            xminlon = amin1(xminlon,xlon(i,j))
            xmaxlon = amax1(xmaxlon,xlon(i,j))
          else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)-360.)) )    &
                  & then
            xminlon = amin1(xminlon,xlon(i,j))
            xmaxlon = amax1(xmaxlon,xlon(i,j))
          else
            xminlon = amin1(xminlon,xlon(i,j)-360.)
            xmaxlon = amax1(xmaxlon,xlon(i,j)-360.)
          end if
        end do
      end do
 
      print 99001 , xminlat , xmaxlat , xminlon , xmaxlon , ntypec
!--------initialize minimum lat and lon of data from tape
      grdltmn = xminlat + 5.
      grdlnmn = xminlon + 5.
99001 format (1x,'xminlat,xmaxlat,xminlon,xmaxlon,ntypec= ',4F10.2,i10)
!
      end subroutine mxmnll

      end module mod_block
