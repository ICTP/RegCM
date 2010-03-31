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
 
      subroutine vcheki(ier,numerr,aname)
      implicit none
!
! Dummy arguments
!
      character(8) :: aname
      integer :: ier , numerr
      intent (in) aname , ier
      intent (inout) numerr
!
!  flag of detected errors in linear algebra routines
!
      if ( ier.ne.0 ) then
        numerr = numerr + 1
        print 99001 , aname , ier
99001   format ('0 error in determination of ',a8,                      &
               &' using library routine     ier=',i4)
      end if
!
      end subroutine vcheki
