!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine param(nx,ny,kz,np,xlat,xlon,vvarmin,vvarmax,           &
                &      xlat1d,xlon1d,iadim,ndim,plv)
 
      implicit none
!
! Dummy arguments
!
      integer :: kz , ndim , np , nx , ny
      logical :: plv
      integer , dimension(ndim) :: iadim
      real(4) , dimension(ndim) :: vvarmax , vvarmin
      real(4) , dimension(nx,ny) :: xlat , xlon
      real(4) , dimension(ny) :: xlat1d
      real(4) , dimension(nx) :: xlon1d
      intent (in) kz , ndim , np , nx , ny , plv , xlat , xlon
      intent (out) iadim , vvarmax , vvarmin , xlat1d , xlon1d
!
! Local variables
!
      integer :: i , j
!
      vvarmin(1) = xlon(1,ny/2)
      vvarmin(2) = xlat(nx/2,1)
      vvarmin(3) = 1050.
      vvarmax(1) = xlon(nx,ny/2)
      vvarmax(2) = xlat(nx/2,ny)
      vvarmax(3) = 0.
      iadim(1) = nx
      iadim(2) = ny
      if ( .not.plv ) then
        iadim(3) = kz
      else
        iadim(3) = np
      end if
      do i = 1 , nx
        xlon1d(i) = xlon(i,ny/2)
      end do
      do j = 1 , ny
        xlat1d(j) = xlat(nx/2,j)
      end do
 
      end subroutine param
