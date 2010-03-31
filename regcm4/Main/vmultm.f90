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
 
      subroutine vmultm(a1,a2,a3,nk)
      implicit none
!
! Dummy arguments
!
      real(8) :: c
      integer :: nk
      real(8) , dimension(nk,nk) :: a1 , a2 , a3
      intent (in) a2 , a3 , c , nk
      intent (inout) a1
!
! Local variables
!
      integer :: k , l , mm
!
      do l = 1 , nk
        do k = 1 , nk
          a1(k,l) = 0.
          do mm = 1 , nk
            a1(k,l) = a1(k,l) + a2(k,mm)*a3(mm,l)
          end do
        end do
      end do
      return
!
      entry vaddms(a1,a2,a3,nk)
      do l = 1 , nk
        do k = 1 , nk
          a1(k,l) = a2(k,l) + a3(k,l)
        end do
      end do
      return
!
      entry vsubtm(a1,a2,a3,nk)
      do l = 1 , nk
        do k = 1 , nk
          a1(k,l) = a2(k,l) - a3(k,l)
        end do
      end do
      return
!
      entry vmultc(a1,a2,nk,c)
      do l = 1 , nk
        do k = 1 , nk
          a1(k,l) = a2(k,l)*c
        end do
      end do
      end subroutine vmultm
