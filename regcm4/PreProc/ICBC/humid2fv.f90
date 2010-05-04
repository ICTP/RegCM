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

      subroutine humid2fv(t,q,ps,pt,sigma,ni,nj,nk)
      use mod_constants , only : tzero , lh0 , lh1 , lsvp1 , lsvp2 , ep2
      implicit none
!
! PARAMETER definitions
!
      real(4) , parameter :: tr = 1./tzero
!
! Dummy arguments
!
      integer :: ni , nj , nk
      real(4) :: pt
      real(4) , dimension(ni,nj) :: ps
      real(4) , dimension(ni,nj,nk) :: q , t
      real(4) , dimension(nk) :: sigma
      intent (in) ni , nj , nk , ps , pt , sigma , t
      intent (inout) q
!
! Local variables
!
      real(4) :: hl , p , qs , satvp
      integer :: i , j , k
!
!     THIS ROUTINE REPLACES RELATIVE HUMIDITY BY SPECIFIC HUMIDITY
!
      do i = 1 , ni
        do j = 1 , nj
          do k = 1 , nk
            p = (pt+sigma(k)*ps(i,j))*10.
            hl = lh0 - lh1*(t(i,j,k)-tzero)
            satvp = lsvp1*exp(lsvp2*hl*(tr-1./t(i,j,k)))
            qs = ep2*satvp/(p-satvp)
            q(i,j,k) = amax1(q(i,j,k)*qs,0.0)
          end do
        end do
      end do
!
      end subroutine humid2fv
