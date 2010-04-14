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
 
      module mod_when

      contains
!
! SUBROUTINE wheneq
!
      subroutine wheneq(n,array,inc,itarg,indx,nval)
 
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
        if ( array(ina).eq.itarg ) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do
 
      end subroutine wheneq
!
! SUBROUTINE whenne
!
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
!
! SUBROUTINE whenfgt
!
      subroutine whenfgt(n,array,inc,rtarg,indx,nval)
 
      implicit none
!
! Dummy arguments
!
      integer :: inc , n , nval
      real(8) :: rtarg
      real(8) , dimension(*) :: array
      integer , dimension(*) :: indx
      intent (in) array , inc , n , rtarg
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
        if ( array(ina).gt.rtarg ) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do
 
      end subroutine whenfgt
!
! SUBROUTINE whenflt
!
      subroutine whenflt(n,array,inc,rtarg,indx,nval)
 
      implicit none
!
! Dummy arguments
!
      integer :: inc , n , nval
      real(8) :: rtarg
      real(8) , dimension(*) :: array
      integer , dimension(*) :: indx
      intent (in) array , inc , n , rtarg
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
        if ( array(ina).lt.rtarg ) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do
      end subroutine whenflt

      end module mod_when
