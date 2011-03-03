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
 
      module mod_cumtran
!
! Tracer convective transport
!
      use mod_runparams
      use mod_trachem
      use mod_mainchem
!
      private
!
      public :: cumtran
!
      contains
!
      subroutine cumtran

      implicit none
!
! Local variables
!
      real(8) :: chiabar , chibbar , cumfrc , deltas
      integer :: i , j , k , kcumtop , n
!
!hy   the well-mixing only over cumulus cloud fraction (assuming 0.3)
!     add the well-mixing fraction in dt step (5%) - only updraft
!hy   cumfrc = 0.3
!
      cumfrc = 0.015D0
!
      do n = 1 , ntr
#ifdef MPP1
        do j = jbegin , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm1
#endif
#endif
          do i = 2 , iym1
            if ( icumtop(i,j) > 0 ) then
              deltas = d_zero
              chiabar = d_zero
              chibbar = d_zero
              kcumtop = max0(icumtop(i,j),4)
              do k = kcumtop , kz
                deltas = deltas + dsigma(k)
                chiabar = chiabar + chia(i,k,j,n)*dsigma(k)
                chibbar = chibbar + chib(i,k,j,n)*dsigma(k)
              end do
!?            do 95 k=icumtop(i,j),kz      ! yhuang, 12/98
!qhy
              do k = kcumtop , kz
                chia(i,k,j,n) = chia(i,k,j,n)*(d_one-cumfrc)  &
                              & + cumfrc*chiabar/deltas
                chib(i,k,j,n) = chib(i,k,j,n)*(d_one-cumfrc)  &
                              & + cumfrc*chibbar/deltas
              end do
            end if
          end do
        end do
      end do
!
      end subroutine cumtran
!
      end module mod_cumtran
