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
  real(rkx), pointer, contiguous, dimension(:,:,:) :: o3prof => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:,:)  :: gasabsnxt => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:,:)  :: gasabstot => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: gasemstot => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: taucldsp => null( )
  logical, save :: doabsems, dolw, dosw
  integer(ik4), save :: ichso4, ichbc, ichoc
  integer(ik4), save :: kclimh, kth, ktf, ksf, kclimf
  real(rkx) :: chfrovrradfr ! chfrq/rafrq

  real(rkx), dimension(nspi), parameter :: wavmin = &
    [ 0.200_rkx, 0.245_rkx, 0.265_rkx, 0.275_rkx, 0.285_rkx ,&
      0.295_rkx, 0.305_rkx, 0.350_rkx, 0.640_rkx, 0.700_rkx ,&
      0.701_rkx, 0.701_rkx, 0.701_rkx, 0.701_rkx, 0.702_rkx ,&
      0.702_rkx, 2.630_rkx, 4.160_rkx, 4.160_rkx ]

  real(rkx), dimension(nspi), parameter :: wavmax = &
    [ 0.245_rkx, 0.265_rkx, 0.275_rkx, 0.285_rkx, 0.295_rkx ,&
      0.305_rkx, 0.350_rkx, 0.640_rkx, 0.700_rkx, 5.000_rkx ,&
      5.000_rkx, 5.000_rkx, 5.000_rkx, 5.000_rkx, 5.000_rkx ,&
      5.000_rkx, 2.860_rkx, 4.550_rkx, 4.550_rkx ]

  real(rkx), dimension(14), parameter :: wavnm1 = &
    [ 2600.0_rkx, 3250.0_rkx, 4000.0_rkx, 4650.0_rkx, 5150.0_rkx, &
      6150.0_rkx, 7700.0_rkx, 8050.0_rkx,12850.0_rkx,16000.0_rkx, &
     22650.0_rkx,29000.0_rkx,38000.0_rkx,  820.0_rkx ]

  real(rkx), dimension(14), parameter :: wavnm2 = &
    [  3250.0_rkx, 4000.0_rkx, 4650.0_rkx, 5150.0_rkx, 6150.0_rkx, &
       7700.0_rkx, 8050.0_rkx,12850.0_rkx,16000.0_rkx,22650.0_rkx, &
      29000.0_rkx,38000.0_rkx,50000.0_rkx, 2600.0_rkx ]

end module mod_rad_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
