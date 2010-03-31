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
 
      subroutine vcheke(er,ei,nk,numerr,aname)
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: tol = 1.E-9
!
! Dummy arguments
!
      character(8) :: aname
      integer :: nk , numerr
      real(8) , dimension(nk) :: ei , er
      intent (in) aname , ei , er , nk
      intent (inout) numerr
!
! Local variables
!
      real(8) :: emax
      integer :: n , nimag , numneg
!
!  check that eigenvalues are real and positive valued.
!
      numneg = 0
      emax = 0.
      do n = 1 , nk
        if ( er(n).le.0. ) numneg = numneg + 1
        if ( er(n).gt.emax ) emax = er(n)
      end do
!
      nimag = 0
      do n = 1 , nk
        if ( ei(n)/emax.gt.tol ) nimag = nimag + 1
      end do
!
      if ( numneg+nimag.eq.0 ) then
        return
      end if
!
      numerr = numerr + 1
      print 99001 , aname , numneg , nimag
99001 format ('0 problem with equivalent depths determined from ',a8,/, &
            & 10x,i3,' depths are nonpositive valued',10x,i3,           &
             &' are not real')
!
      end subroutine vcheke
