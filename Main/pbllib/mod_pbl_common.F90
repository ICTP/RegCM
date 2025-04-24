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

module mod_pbl_common
!
! Storage parameters and constants related to the boundary layer
!
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_regcm_types

  implicit none

  private

  real(rkx), public, pointer, contiguous, dimension(:,:) :: ricr

  integer(ik4), pointer, contiguous, public, dimension(:,:) :: kmxpbl

  !
  ! Pointers to the TCM state variables
  !
  type(tcm_state), public :: uwstate
  real(rkx), public, pointer, contiguous, dimension(:,:,:,:) :: chiuwten

end module mod_pbl_common

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
