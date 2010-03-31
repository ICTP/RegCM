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
 
      subroutine whenne(n,array,inc,itarg,indx,nval)
 
      implicit none
!
! Dummy arguments
!
      integer :: inc , itarg , n , nval
      integer , dimension(*) :: array , indx
      intent (in) array , inc , itarg , n
      intent (out) indx
      intent (inout) nval
!
! Local variables
!
      integer :: i , ina
!
      ina = 1
      nval = 0
      if ( inc.lt.0 ) ina = (-inc)*(n-1) + 1
      do i = 1 , n
        if ( array(ina).ne.itarg ) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do
      end subroutine whenne
