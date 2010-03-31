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

      subroutine xtdot(px,pd,ni,nj,nk,ni1,nj1)
 
      implicit none
!
! Dummy arguments
!
      integer :: ni , ni1 , nj , nj1 , nk
      real(4) , dimension(ni,nj,nk) :: pd , px
      intent (in) ni , ni1 , nj , nj1 , nk , px
      intent (out) pd
!
! Local variables
!
      integer :: i , j , k
!
!     this routine determines p(.) from p(x) by a 4-point interpolation.
!     on the x-grid, a p(x) point outside the grid domain is assumed to
!     satisfy p(0,j)=p(1,j); p(ni,j)=p(ni-1,j); and similarly for the
!     i's.
!
      do k = 1 , nk
        do j = 2 , nj1
          do i = 2 , ni1
            pd(i,j,k) = 0.25*(px(i,j,k)+px(i-1,j,k)+px(i,j-1,k)         &
                      & +px(i-1,j-1,k))
          end do
        end do
 
        do i = 2 , ni1
          pd(i,1,k) = 0.5*(px(i,1,k)+px(i-1,1,k))
          pd(i,nj,k) = 0.5*(px(i,nj1,k)+px(i-1,nj1,k))
        end do
 
        do j = 2 , nj1
          pd(1,j,k) = 0.5*(px(1,j,k)+px(1,j-1,k))
          pd(ni,j,k) = 0.5*(px(ni1,j,k)+px(ni1,j-1,k))
        end do
 
        pd(1,1,k) = px(1,1,k)
        pd(1,nj,k) = px(1,nj1,k)
        pd(ni,1,k) = px(ni1,1,k)
        pd(ni,nj,k) = px(ni1,nj1,k)
 
      end do
 
      end subroutine xtdot
