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
 
      subroutine minimi(array,ix,kx,ks,kend,kt,istart,iend)
!
      implicit none
!
! Dummy arguments
!
      integer :: iend , istart , ix , kend , kx
      real(8) , dimension(ix,kx) :: array
      integer , dimension(ix) :: ks , kt
      intent (in) array , iend , istart , ix , kend , ks , kx
      intent (out) kt
!
! Local variables
!
      integer :: i , k
      real(8) :: x
!
      do i = istart , iend
        kt(i) = ks(i)
        x = array(i,ks(i))
        do k = ks(i) + 1 , kend
          if ( array(i,k).lt.x ) then
            x = array(i,k)
            kt(i) = k
          end if
        end do
      end do
!
      end subroutine minimi
