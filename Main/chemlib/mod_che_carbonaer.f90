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
  
  use m_realkinds
  use mod_dynparam
  use mod_constants
  use mod_che_common
  use mod_che_species
  use mod_che_indices

  ! Parameter usefull for wet and dry deposition of carbon aerosol 
  ! densities in kg/m3
  real(dp) , parameter :: chrhobc   = 2000.0D0
  real(dp) , parameter :: chrhobchl = 1600.0D0
  real(dp) , parameter :: chrhooc   = 1200.0D0
  real(dp) , parameter :: chrhoochl = 1200.0D0

  ! effctive dimaters ( and not radius!)  in micrometer
  ! ( should they be defined intercatively in the future ? ) 
  real(dp) , parameter :: chreffbc   = 0.05D0
  real(dp) , parameter :: chreffbchl = 0.3D0
  real(dp) , parameter :: chreffoc   = 0.2D0
  real(dp) , parameter :: chreffochl = 0.3D0

  ! aging efolding time (s), from hydrophobic to hydrophilic
  ! Cooke et al.
  real(dp) , parameter :: chagct = 1.15D0 * 86400.0D0
  !
  ! bin size for carboneaceous aerosols
  ! ps add one dimension for sulfate too.
  real(dp) , dimension(5,2) :: carbsiz

  private

  public :: aging_carb
  public :: chrhobc , chrhobchl , chrhooc , chrhoochl
  public :: chreffbc , chreffbchl , chreffoc , chreffochl
  public :: carbsiz

  contains

  subroutine aging_carb(j)

    implicit none
!
    integer, intent(in) :: j
    integer :: i , k
    real(dp) :: agingtend1 , agingtend2

    ! aging o carbon species : Conversion from hydrophobic to 
    ! hydrophilic: Carbonaceopus species time constant
    ! ( 1.15 day cooke et al.,1999)

    if ( ibchb > 0 .and. ibchl > 0 ) then
      do k = 1 , kz
        do i = 2 , iym2
          agingtend1 = -chib(i,k,j,ibchb)*(d_one-dexp(-mdfrq/chagct))/mdfrq
          agingtend2 = -agingtend1
          chiten(i,k,j,ibchb) = chiten(i,k,j,ibchb) + agingtend1
          chiten(i,k,j,ibchl) = chiten(i,k,j,ibchl) + agingtend2
        end do
      end do
    end if
 
    if ( iochb > 0  .and. iochl > 0 ) then
      do k = 1 , kz
        do i = 2 , iym2
          agingtend1 = -chib(i,k,j,iochb)*(d_one-dexp(-mdfrq/chagct))/mdfrq
          agingtend2 = -agingtend1
          chiten(i,k,j,iochb) = chiten(i,k,j,iochb) + agingtend1
          chiten(i,k,j,iochl) = chiten(i,k,j,iochl) + agingtend2
        end do
      end do
    end if
  end subroutine aging_carb

end module mod_che_carbonaer
