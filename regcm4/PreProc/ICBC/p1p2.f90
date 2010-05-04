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

      subroutine p1p2(pd,px,ni,nj)
      implicit none
!
! Dummy arguments
!
      integer :: ni , nj
      real(4) , dimension(ni,nj) :: pd , px
      intent (in) ni , nj , px
      intent (out) pd
!
! Local variables
!
      integer :: i , j , ni1 , nj1
!
!     THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
!     ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
!     SATISFY P(0,J)=P(1,J); P(NI,J)=P(NI-1,J); AND SIMILARLY FOR THE
!     I'S.
!
      ni1 = ni - 1
      nj1 = nj - 1
!
      do j = 2 , nj1
        do i = 2 , ni1
          pd(i,j) = 0.25*(px(i,j)+px(i-1,j)+px(i,j-1)+px(i-1,j-1))
        end do
      end do
!
      do i = 2 , ni1
        pd(i,1) = 0.5*(px(i,1)+px(i-1,1))
        pd(i,nj) = 0.5*(px(i,nj1)+px(i-1,nj1))
      end do
!
      do j = 2 , nj1
        pd(1,j) = 0.5*(px(1,j)+px(1,j-1))
        pd(ni,j) = 0.5*(px(ni1,j)+px(ni1,j-1))
      end do
!
      pd(1,1) = px(1,1)
      pd(1,nj) = px(1,nj1)
      pd(ni,1) = px(ni1,1)
      pd(ni,nj) = px(ni1,nj1)
!
      end subroutine p1p2
