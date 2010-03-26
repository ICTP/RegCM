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

      subroutine avgrain(favg,f,f1,nx,ny,nhr,ihr)
 
      implicit none
!
! Dummy arguments
!
      integer :: ihr , nhr , nx , ny
      real(4) , dimension(nx,ny) :: f , f1
      real(4) , dimension(nx,ny,nhr) :: favg
      intent (in) f , ihr , nhr , nx , ny
      intent (inout) f1 , favg
!
! Local variables
!
      integer :: i , j
!
      do j = 1 , ny
        do i = 1 , nx
          favg(i,j,ihr) = favg(i,j,ihr) + (f(i,j)-f1(i,j))
          f1(i,j) = f(i,j)
        end do
      end do
 
      end subroutine avgrain
