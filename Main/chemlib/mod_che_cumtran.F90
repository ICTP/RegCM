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
 
module mod_che_cumtran
!
! Tracer convective transport
!
  use mod_dynparam
  use mod_realkinds
  use mod_che_common
!
  private
!
  public :: cumtran
!
  contains
!
  subroutine cumtran
    implicit none
    real(dp) :: chiabar , chibbar , deltas , cumfrc
    integer :: i , j , k , kctop , n

    if ( ichdiag == 1 ) chiten0 = chib
    do n = 1 , ntr
      do j = jci1 , jci2
        do i = ici1 , ici2
          if ( kcumtop(j,i) > 0 ) then
            deltas = d_zero
            chiabar = d_zero
            chibbar = d_zero
            kctop = max0(kcumtop(j,i),4)
            do k = kctop , kz
              deltas = deltas + cdsigma(k)
              chiabar = chiabar + chia(j,i,k,n)*cdsigma(k)
              chibbar = chibbar + chib(j,i,k,n)*cdsigma(k)
            end do
            do k = kctop , kz
 !             cumfrc =  ccldfra(j,i,k) - cfcc(j,i,k)
              cumfrc =  convcldfra (j,i,k)
              chia(j,i,k,n) = chia(j,i,k,n)*(d_one-cumfrc)+cumfrc*chiabar/deltas
              chib(j,i,k,n) = chib(j,i,k,n)*(d_one-cumfrc)+cumfrc*chibbar/deltas
            end do
          end if
        end do
      end do
    end do
    ! here calculate a pseudo tendency.
    ! factor 2 is added since we are out of leap frog 
    if ( ichdiag == 1 ) then
      do j = jci1 , jci2
        do i = ici1 , ici2
        cconvdiag(j,i,:,:)  = cconvdiag(j,i,:,:) + ( (chiten0(j,i,:,:) - chib(j,i,:,:)) / dtche ) * d_two * cdiagf 
       end do
     end do
    end if
  end subroutine cumtran
!
end module mod_che_cumtran
