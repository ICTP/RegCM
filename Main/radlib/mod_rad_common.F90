!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_rad_common

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam, only : nspi

  implicit none

  public

  ! absnxt  - Nearest layer absorptivities
  ! abstot  - Non-adjacent layer absorptivites
  ! emstot  - Total emissivity

  ! Those need to be saved in output file
  real(rk8), pointer, contiguous, dimension(:,:,:) :: o3prof => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:,:)  :: gasabsnxt => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:,:)  :: gasabstot => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: gasemstot => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: taucldsp => null( )
  logical, save :: doabsems, dolw, dosw
  integer(ik4), save :: ichso4, ichbc, ichoc
  integer(ik4), save :: kclimh, kth, ktf, ksf, kclimf
  real(rk8) :: chfrovrradfr ! chfrq/rafrq

  real(rk8), dimension(nspi), parameter :: wavmin = &
    [ 0.200_rk8, 0.245_rk8, 0.265_rk8, 0.275_rk8, 0.285_rk8 ,&
      0.295_rk8, 0.305_rk8, 0.350_rk8, 0.640_rk8, 0.700_rk8 ,&
      0.701_rk8, 0.701_rk8, 0.701_rk8, 0.701_rk8, 0.702_rk8 ,&
      0.702_rk8, 2.630_rk8, 4.160_rk8, 4.160_rk8 ]

  real(rk8), dimension(nspi), parameter :: wavmax = &
    [ 0.245_rk8, 0.265_rk8, 0.275_rk8, 0.285_rk8, 0.295_rk8 ,&
      0.305_rk8, 0.350_rk8, 0.640_rk8, 0.700_rk8, 5.000_rk8 ,&
      5.000_rk8, 5.000_rk8, 5.000_rk8, 5.000_rk8, 5.000_rk8 ,&
      5.000_rk8, 2.860_rk8, 4.550_rk8, 4.550_rk8 ]

  real(rk8), dimension(14), parameter :: wavnm1 = &
    [ 2600.0_rk8, 3250.0_rk8, 4000.0_rk8, 4650.0_rk8, 5150.0_rk8, &
      6150.0_rk8, 7700.0_rk8, 8050.0_rk8,12850.0_rk8,16000.0_rk8, &
     22650.0_rk8,29000.0_rk8,38000.0_rk8,  820.0_rk8 ]

  real(rk8), dimension(14), parameter :: wavnm2 = &
    [  3250.0_rk8, 4000.0_rk8, 4650.0_rk8, 5150.0_rk8, 6150.0_rk8, &
       7700.0_rk8, 8050.0_rk8,12850.0_rk8,16000.0_rk8,22650.0_rk8, &
      29000.0_rk8,38000.0_rk8,50000.0_rk8, 2600.0_rk8 ]

end module mod_rad_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
