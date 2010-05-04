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

      subroutine humid1fv(t,q,p3d,ni,nj,nk)
      use mod_constants , only : tzero , ep2 , lh0 , lh1 , lsvp1 , lsvp2
      implicit none
!
! PARAMETER definitions
!
      real(4) , parameter :: tr = 1./tzero
!
! Dummy arguments
!
      integer :: ni , nj , nk
      real(4) , dimension(ni,nj,nk) :: p3d , q , t
      intent (in) ni , nj , nk , p3d , t
      intent (inout) q
!
! Local variables
!
      real(4) :: hl , qs , satvp
      integer :: i , j , k
!
!     THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!
      do i = 1 , ni
        do j = 1 , nj
          do k = 1 , nk
            if ( p3d(i,j,k)>-9990. ) then
              hl = lh0 - lh1*(t(i,j,k)-tzero)  ! LATENT HEAT OF EVAP.
              satvp = lsvp1*exp(lsvp2*hl*(tr-1./t(i,j,k)))
                                                   ! SATURATION VAP PRESS.
              qs = ep2*satvp/(p3d(i,j,k)-satvp)   ! SAT. MIXING RATIO
              q(i,j,k) = amax1(q(i,j,k)/qs,0.0)    !ALREADY MIXING RATIO
            else
              q(i,j,k) = -9999.
            end if
          end do
        end do
      end do
      end subroutine humid1fv
