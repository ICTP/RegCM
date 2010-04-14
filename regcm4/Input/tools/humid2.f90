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

      subroutine humid2(t,q,preslv,im,jm,kp)
      use mod_constants , only : tzero , rtzero , lh0 , lh1 , lsvp1 ,   &
               &                 lsvp2 , ep2
      use mod_preproc_param , only : ptop
      implicit none
!
! PARAMETER definitions
!
      real(4) , parameter :: qmin = 0.0
!
! Dummy arguments
!
      integer :: im , jm , kp
      real(4) , dimension(kp) :: preslv
      real(4) , dimension(im,jm,kp) :: q , t
      intent (in) im , jm , kp , preslv , t
      intent (inout) q
!
! Local variables
!
      real(4) :: hl , qs , satvp
      integer :: i , j , k
!
!     THIS ROUTINE REPLACES SPECIFIC HUMIDITY BY RELATIVE HUMIDITY
!     DATA ON SIGMA LEVELS
!
      do k = 1 , kp          ! MINIMUM VALUE OF SPECIFIC HUMIDITY
        do j = 1 , jm
          do i = 1 , im
            hl = lh0 - lh1*(t(i,j,k)-tzero)
            satvp = lsvp1*exp(lsvp2*hl*(rtzero-1./t(i,j,k)))
            qs = ep2*satvp/(preslv(k)-satvp)          ! preslv (hPa)
            if ( q(i,j,k)<qmin ) q(i,j,k) = qmin      ! SPECIFIED MINIMUM
            q(i,j,k) = q(i,j,k)*qs
          end do
        end do
      end do
      end subroutine humid2
