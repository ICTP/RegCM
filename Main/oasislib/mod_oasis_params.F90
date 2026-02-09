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

  implicit none (type, external)

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
    type(infogrd), pointer :: grd => null()       ! field related definition indexes
  end type infofld
  public :: infofld

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

  character(len=11), dimension(-2:14), parameter, public :: getput_status = & ! getput kinf string
  ['NotDef', 'VarUncpl', 'Ok', '', '', 'Recvd', 'Sent', 'LocTrans', &
    'ToRest', 'Output', 'SentOut', 'ToRestOut', &
    'FromRest', 'Input', 'RecvOut', 'FromRestOut', 'WaitGroup']

end module mod_oasis_params
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
