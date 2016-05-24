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

  implicit none

  private

  ! Parameter usefull for wet and dry deposition of carbon aerosol
  ! densities in kg/m3

  real(rkx) , public , parameter :: rhobc   = 2000.0_rkx
  real(rkx) , public , parameter :: rhooc   = 1200.0_rkx
  real(rkx) , public , parameter :: rhobchl = 1600.0_rkx
  real(rkx) , public , parameter :: rhoochl = 1200.0_rkx
  real(rkx) , public , parameter :: rhosm1 = 1200.0_rkx
  real(rkx) , public , parameter :: rhosm2 = 1200.0_rkx


  ! effctive dimaters ( and not radius!)  in micrometer
  ! ( should they be defined intercatively in the future ? )
  real(rkx) , public , parameter :: reffbc   = 0.05_rkx
  real(rkx) , public , parameter :: reffbchl = 0.3_rkx
  real(rkx) , public , parameter :: reffoc   = 0.2_rkx
  real(rkx) , public , parameter :: reffochl = 0.3_rkx
  real(rkx) , public , parameter :: reffsm1 = 0.3_rkx
  real(rkx) , public , parameter :: reffsm2 = 0.3_rkx


  ! aging efolding time (s), from hydrophobic to hydrophilic
  ! Cooke et al.
  real(rkx) , parameter :: chagct = 1.15_rkx * 86400.0_rkx
  real(rkx) , parameter :: chsmct = 1.15_rkx * 86400.0_rkx

  !
  ! solubility of carbon aer for rain out param of giorgi and chameides
  !
  real(rkx) , parameter :: solbc = 0.05_rkx
  real(rkx) , parameter :: solbchl = 0.8_rkx
  real(rkx) , parameter :: soloc = 0.05_rkx
  real(rkx) , parameter :: solochl = 0.8_rkx
  real(rkx) , parameter :: solsm1 = 0.05_rkx
  real(rkx) , parameter :: solsm2 = 0.8_rkx


  public :: aging_carb , solbc , solbchl , soloc , solochl, solsm1,solsm2

  ! bin size for carboneaceous aerosols
  ! ps add one dimension for sulfate too.
  real(rkx) , public , dimension(9) :: carbed

  contains

    subroutine aging_carb(j)
      implicit none
      integer, intent(in) :: j
      integer(ik4) :: i , k
      real(rkx) :: agingtend1 , agingtend2 , kav
      !
      ! aging o carbon species : Conversion from hydrophobic to
      ! hydrophilic: Carbonaceopus species time constant
      ! ( 1.15 day cooke et al.,1999)
      !
      if ( ibchb > 0 .and. ibchl > 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            kav = max(chib(j,i,k,ibchb)-mintr,d_zero)
            agingtend1 = -kav*(d_one-exp(-dt/chagct))/dt
            agingtend2 = -agingtend1
            chiten(j,i,k,ibchb) = chiten(j,i,k,ibchb) + agingtend1
            chiten(j,i,k,ibchl) = chiten(j,i,k,ibchl) + agingtend2
          end do
        end do
      end if
      if ( iochb > 0  .and. iochl > 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            kav = max(chib(j,i,k,iochb)-mintr,d_zero)
            agingtend1 = -kav*(d_one-exp(-dt/chagct))/dt
            agingtend2 = -agingtend1
            chiten(j,i,k,iochb) = chiten(j,i,k,iochb) + agingtend1
            chiten(j,i,k,iochl) = chiten(j,i,k,iochl) + agingtend2
          end do
        end do
      end if

      if ( ism1 > 0  .and. ism2 > 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            kav = max(chib(j,i,k,iochb)-mintr,d_zero)
            agingtend1 = -kav*(d_one-exp(-dt/chsmct))/dt
            agingtend2 = -agingtend1
            chiten(j,i,k,ism1) = chiten(j,i,k,ism1) + agingtend1
            chiten(j,i,k,ism2) = chiten(j,i,k,ism2) + agingtend2
          end do
        end do
      end if

    end subroutine aging_carb
!
end module mod_che_carbonaer
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
