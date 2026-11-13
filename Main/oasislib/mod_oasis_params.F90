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

#ifdef OASIS

module mod_oasis_params

  use mod_intkinds

  implicit none

  private

  ! grid information type
  type infogrd
    character(len=4) :: naNM, naWM ! grid name (No Mask, With Mask)
    integer :: id                   ! partition identificator
    integer :: j1, j2, jll, i1, i2, ill, & ! local indexes
               ja, jb, jgl, ia, ib, igl, & ! global indexes
               nc ! number of corner per cell (1:nc)
  end type infogrd
  public :: infogrd

  ! field information type
  type infofld
    character(len=:), allocatable :: na ! field name
    integer :: id                       ! field identificator
    ! field related definition indexes
    type(infogrd), pointer :: grd => null()
  end type infofld
  public :: infofld

  ! Multi-level field: an OASIS bundle of nz fields on the same 2d grid
  ! (layers, species, size bins), exchanged as one array (j,i,nz).
  type info3dfld
    character(len=:), allocatable :: na ! field name
    integer :: id                       ! field identificator
    integer :: nz                       ! number of fields in the bundle
    type(infogrd), pointer :: grd => null()
  end type info3dfld
  public :: info3dfld

  character(len=6), parameter, public :: comp_name = 'REGCM5' ! component name
  integer(ik4), public :: comp_id ! component identification
  integer(ik4), public :: oasis_lag ! model time lag to other components
                                     ! (variable)

  ! before OASIS-related debug statements
  character(len=8), parameter, public :: oasis_prefix = '[OASIS] '

  ! oasisparam namelist general parameters
  integer(ik4), public :: write_restart_option
  logical, public :: l_write_grids
  integer(ik4), public :: oasis_sync_lag ! model time lag to other components
                                          ! (parameter)
#ifdef ECLM
  ! oasisparam namelist: sizes of the eCLM multi-level fields. They must
  ! match the eCLM configuration.
  integer(ik4), public :: oasis_nlevsoi  ! soil layers (RCM_TSOI, RCM_H2O*)
  integer(ik4), public :: oasis_nddvel   ! dry deposition species (RCM_DDVEL)
  ! MEGAN compounds in RCM_FLXVOC, in the eCLM order
  integer(ik4), parameter, public :: oasis_max_megan = 32
  character(len=8), dimension(oasis_max_megan), public :: oasis_megan_cmp
  ! dust size bins in RCM_FLXDST (CLM dust model)
  integer(ik4), parameter, public :: oasis_ndust = 4
#endif

  ! oasis_get / oasis_put return codes as strings
  character(len=11), dimension(-2:14), parameter, public :: getput_status = &
    [ character(len=11) :: 'NotDef', 'VarUncpl', 'Ok', '', '', 'Recvd', &
      'Sent', 'LocTrans', 'ToRest', 'Output', 'SentOut', 'ToRestOut', &
      'FromRest', 'Input', 'RecvOut', 'FromRestOut', 'WaitGroup' ]

end module mod_oasis_params
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
