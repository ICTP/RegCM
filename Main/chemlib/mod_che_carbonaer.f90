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

module mod_che_carbonaer
  
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_che_common
  use mod_che_species
  use mod_che_indices

  private

  ! Parameter usefull for wet and dry deposition of carbon aerosol 
  ! densities in kg/m3

  real(rk8) , public , parameter :: rhobc   = 2000.0D0
  real(rk8) , public , parameter :: rhooc   = 1200.0D0
  real(rk8) , public , parameter :: rhobchl = 1600.0D0
  real(rk8) , public , parameter :: rhoochl = 1200.0D0

  ! effctive dimaters ( and not radius!)  in micrometer
  ! ( should they be defined intercatively in the future ? ) 
  real(rk8) , public , parameter :: reffbc   = 0.05D0
  real(rk8) , public , parameter :: reffbchl = 0.3D0
  real(rk8) , public , parameter :: reffoc   = 0.2D0
  real(rk8) , public , parameter :: reffochl = 0.3D0

  ! aging efolding time (s), from hydrophobic to hydrophilic
  ! Cooke et al.
  real(rk8) , parameter :: chagct = 1.15D0 * 86400.0D0
  !
  ! solubility of carbon aer for rain out param of giorgi and chameides
  !
  real(rk8) , parameter :: solbc = 0.05D0
  real(rk8) , parameter :: solbchl = 0.8D0
  real(rk8) , parameter :: soloc = 0.05D0
  real(rk8) , parameter :: solochl = 0.8D0
  public :: aging_carb , solbc , solbchl , soloc , solochl

  ! bin size for carboneaceous aerosols
  ! ps add one dimension for sulfate too.
  real(rk8) , public , dimension(5) :: carbed

  contains

    subroutine aging_carb(j)
      implicit none
      integer, intent(in) :: j
      integer(ik4) :: i , k
      real(rk8) :: agingtend1 , agingtend2
      !
      ! aging o carbon species : Conversion from hydrophobic to 
      ! hydrophilic: Carbonaceopus species time constant
      ! ( 1.15 day cooke et al.,1999)
      !
      if ( ibchb > 0 .and. ibchl > 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            agingtend1 = -chib(j,i,k,ibchb)*(d_one-dexp(-dt/chagct))/dt
            agingtend2 = -agingtend1
            chiten(j,i,k,ibchb) = chiten(j,i,k,ibchb) + agingtend1
            chiten(j,i,k,ibchl) = chiten(j,i,k,ibchl) + agingtend2
          end do
        end do
      end if
      if ( iochb > 0  .and. iochl > 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            agingtend1 = -chib(j,i,k,iochb)*(d_one-dexp(-dt/chagct))/dt
            agingtend2 = -agingtend1
            chiten(j,i,k,iochb) = chiten(j,i,k,iochb) + agingtend1
            chiten(j,i,k,iochl) = chiten(j,i,k,iochl) + agingtend2
          end do
        end do
      end if
    end subroutine aging_carb
!
end module mod_che_carbonaer
