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
 
      function intmax(n,iy,inc)
!
      implicit none
!
! Dummy arguments
!
      integer :: inc , n
      integer :: intmax
      integer , dimension(*) :: iy
      intent (in) inc , iy , n
!
! Local variables
!
      integer :: i , mx
!
      mx = iy(1)
      intmax = 1
      do i = 1 + inc , inc*n , inc
        if ( iy(i).gt.mx ) then
          mx = iy(i)
          intmax = i
        end if
      end do
      end function intmax
