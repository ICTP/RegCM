!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine cumtran

      use mod_regcm_param
      use mod_trachem
      use mod_mainchem
      use mod_param3
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
      cumfrc = 0.015
 
      do n = 1 , ntr
#ifdef MPP1
        do j = jbegin , jendx
#else
        do j = 2 , jxm1
#endif
          do i = 2 , ixm1
            if ( icumtop(i,j).gt.0 ) then
              deltas = 0.
              chiabar = 0.
              chibbar = 0.
              kcumtop = max0(icumtop(i,j),4)
              do k = kcumtop , kx
                deltas = deltas + dsigma(k)
                chiabar = chiabar + chia(i,k,j,n)*dsigma(k)
                chibbar = chibbar + chib(i,k,j,n)*dsigma(k)
              end do
!?            do 95 k=icumtop(i,j),kx      ! yhuang, 12/98
!qhy
              do k = kcumtop , kx
                chia(i,k,j,n) = chia(i,k,j,n)*(1.-cumfrc)               &
                              & + cumfrc*chiabar/deltas
                chib(i,k,j,n) = chib(i,k,j,n)*(1.-cumfrc)               &
                              & + cumfrc*chibbar/deltas
              end do
            end if
          end do
        end do
      end do
      end subroutine cumtran
