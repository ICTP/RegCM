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

      subroutine two2three(i,nx,ny,kx,fin,fout)
 
      implicit none
!
! Dummy arguments
!
      integer :: i , nx , ny , kx
      real(4) , dimension(ny,kx) :: fin
      real(4) , dimension(nx,ny,kx) :: fout
      intent (in) i , fin , nx , ny , kx
      intent (out) fout
!
! Local variables
!
      integer :: j , k
!
      do k = 1 , kx
        do j = 1 , ny
          fout(i,j,k) = fin(j,k)
        end do
      end do
 
      end subroutine two2three
