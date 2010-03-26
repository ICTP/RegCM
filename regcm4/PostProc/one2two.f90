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

      subroutine one2two(i,nx,ny,fin,fout)
 
      implicit none
!
! Dummy arguments
!
      integer :: i , nx , ny
      real(4) , dimension(ny) :: fin
      real(4) , dimension(nx,ny) :: fout
      intent (in) i , fin , nx , ny
      intent (out) fout
!
! Local variables
!
      integer :: j
!
      do j = 1 , ny
        fout(i,j) = fin(j)
      end do
 
      end subroutine one2two
