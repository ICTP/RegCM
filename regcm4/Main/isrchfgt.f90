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
 
      function isrchfgt(n,array,inc,rtarg)
      implicit none
!
! Dummy arguments
!
      integer :: inc , n
      real(8) :: rtarg
      real(8) , dimension(*) :: array
      integer :: isrchfgt
      intent (in) array , inc , n , rtarg
!
! Local variables
!
      integer :: i , ind
!
      if ( n.le.0 ) then
        isrchfgt = 0
        return
      end if
      ind = 1
      do i = 1 , n
        if ( array(ind).gt.rtarg ) then
          isrchfgt = i
          return
        else
          ind = ind + inc
        end if
      end do
      isrchfgt = n + 1
      end function isrchfgt
