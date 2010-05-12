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
 
      subroutine vorder(z,hbar,wz,wh,nk)
      implicit none
!
! Dummy arguments
!
      integer :: nk
      real(8) , dimension(nk) :: hbar
      real(8) , dimension(nk,2) :: wh
      real(8) , dimension(nk,nk) :: wz , z
      intent (in) nk
      intent (inout) hbar , wh , wz , z
!
! Local variables
!
      real(8) :: hmax
      integer :: k , kmax , l
!
!  this routine orders the components of hbar so they are largest to
!  smallest valued.  the columns of z are reorded so that they
!  correspond to the same (but reordered) components of hbar.
!
      kmax = 1
      do k = 1 , nk
        wh(k,1) = hbar(k)
        wh(k,2) = 0.
        do l = 1 , nk
          wz(k,l) = z(k,l)
        end do
      end do
!
      do l = 1 , nk
        hmax = -1.D100
        do k = 1 , nk
          if ( (wh(k,2).eq.0.) .and. (wh(k,1).gt.hmax) ) then
            hmax = wh(k,1)
            kmax = k
          end if
        end do
!
        hbar(l) = hmax
        wh(kmax,2) = 1.
        do k = 1 , nk
          z(k,l) = wz(k,kmax)
        end do
      end do
!
      end subroutine vorder
