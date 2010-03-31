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

      subroutine xlon360(nx,ny,xlon)
      implicit none
!
! Dummy arguments
!
      integer :: nx , ny
      real(4) , dimension(ny,nx) :: xlon
      intent (in) nx , ny
      intent (inout) xlon
!
! Local variables
!
      integer :: i , j
      real(4) , dimension(ny,nx) :: xlon1
!
      do i = 1 , nx
        do j = 1 , ny
          xlon1(j,i) = xlon(j,i)
        end do
      end do
 
      do i = 1 , nx
        do j = 1 , ny
          if ( xlon1(j,i)<0. ) then
            xlon(j,i) = xlon1(j,i) + 360.
          else
            xlon(j,i) = xlon1(j,i)
          end if
        end do
      end do
 
      end subroutine xlon360
