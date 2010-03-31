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
 
      subroutine vnorml(z,s,nk,nk1)
      implicit none
!
! Dummy arguments
!
      integer :: nk , nk1
      real(8) , dimension(nk1) :: s
      real(8) , dimension(nk,nk) :: z
      intent (in) nk , nk1 , s
      intent (inout) z
!
! Local variables
!
      real(8) :: a , v , zmax
      integer :: k , kmax , l
!
!  this routine normalizes the columns of z such that the component
!  with the largest absolute value is positive, and the sum of the
!  mass-weighted squares equals one.
!
      do l = 1 , nk
        zmax = -1.
        v = 0.
!
        do k = 1 , nk
          a = dabs(z(k,l))
          if ( a.gt.zmax ) then
            zmax = a
            kmax = k
          end if
          v = (s(k+1)-s(k))*a*a + v
        end do
!
        a = (z(kmax,l)/zmax)/dsqrt(v)
        do k = 1 , nk
          z(k,l) = a*z(k,l)
        end do
      end do
!
      end subroutine vnorml
