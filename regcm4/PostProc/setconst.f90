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

      subroutine setconst(f,const,n1,n2,n3,n4,n5,i1,ni1,i2,ni2)
 
      implicit none
!
! Dummy arguments
!
      real(4) :: const
      integer :: i1 , i2 , n1 , n2 , n3 , n4 , n5 , ni1 , ni2
      real(4) , dimension(n1,n2,n3,n4,n5) :: f
      intent (in) const , i1 , i2 , n1 , n2 , n3 , n4 , n5 , ni1 , ni2
      intent (out) f
!
! Local variables
!
      integer :: i , j , k , l , m
!
      do m = 1 , n5
        do l = 1 , n4
          do k = 1 , n3
            do j = i2 , ni2
              do i = i1 , ni1
                f(i,j,k,l,m) = const
              end do
            end do
          end do
        end do
      end do
 
      end subroutine setconst
