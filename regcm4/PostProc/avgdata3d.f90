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

      subroutine avgdata3d(favg,f,n1,n2,n3,n4,n5,ihr,vmisdat)
 
      implicit none
!
! Dummy arguments
!
      integer :: ihr , n1 , n2 , n3 , n4 , n5
      real(4) :: vmisdat
      real(4) , dimension(n1,n2,n3,n4) :: f
      real(4) , dimension(n1,n2,n3,n4,n5) :: favg
      intent (in) f , ihr , n1 , n2 , n3 , n4 , n5 , vmisdat
      intent (inout) favg
!
! Local variables
!
      integer :: i , j , k , l
!
      do l = 1 , n4
        do k = 1 , n3
          do j = 1 , n2
            do i = 1 , n1
              if ( f(i,j,k,l)>vmisdat ) then
                favg(i,j,k,l,ihr) = favg(i,j,k,l,ihr) + f(i,j,k,l)
              else
                favg(i,j,k,l,ihr) = vmisdat
              end if
            end do
          end do
        end do
      end do
 
      end subroutine avgdata3d
