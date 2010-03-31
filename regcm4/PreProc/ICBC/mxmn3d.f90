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

      subroutine mxmn3d(var,cvar,jx,iy,np)
      implicit none
!
! Dummy arguments
!
      character(2) :: cvar
      integer :: iy , jx , np
      real , dimension(jx,iy,np) :: var
      intent (in) cvar , iy , jx , np , var
!
! Local variables
!
      integer :: i , j , k
      real :: smax , smin
!
      do k = 1 , np
        smax = -1.E8
        smin = 1.E8
        do j = 1 , iy
          do i = 1 , jx
            if ( smax<var(i,j,k) ) smax = var(i,j,k)
            if ( smin>var(i,j,k) ) smin = var(i,j,k)
          end do
        end do
        write (*,*) cvar , k , smax , smin
      end do
      end subroutine mxmn3d
