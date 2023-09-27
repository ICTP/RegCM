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

#ifdef OASIS

module mod_oasis_params
 
  use mod_intkinds

  implicit none

  private
 
  ! grid information type
  type infogrd
    character(len=4) :: naNM , naWM ! grid name (No Mask , With Mask)
    integer :: id                   ! partition identificator
    integer :: j1 , j2 , jll , i1 , i2 , ill , & ! local indexes
               ja , jb , jgl , ia , ib , igl , & ! global indexes
               nc ! number of corner per cell (1:nc)
  end type infogrd
  public :: infogrd

  ! field information type
  type infofld
    character(len=:), allocatable :: na ! field name
    integer :: id                       ! field identificator
    type(infogrd), pointer :: grd       ! field related definition indexes
  end type infofld
  public :: infofld

  character(len=6) , parameter , public :: comp_name = 'REGCM5' ! component name
  integer(ik4) , public :: comp_id ! component identification
  integer(ik4) , public :: oasis_lag ! model time lag to other components
                                     ! (variable)

  ! before OASIS-related debug statements
  character(len=8) , parameter , public :: oasis_prefix = '[OASIS] '

  ! oasisparam namelist general parameters
  integer(ik4) , public :: write_restart_option
  logical , public :: l_write_grids
  integer(ik4) , public :: oasis_sync_lag ! model time lag to other components
                                          ! (parameter)

  character(len=11), dimension(-2:14), parameter, public :: getput_status = & ! getput kinf string
  (/'NotDef', 'VarUncpl', 'Ok', '', '', 'Recvd', 'Sent', 'LocTrans', &
    'ToRest', 'Output', 'SentOut', 'ToRestOut', &
    'FromRest', 'Input', 'RecvOut', 'FromRestOut', 'WaitGroup'/)

end module mod_oasis_params
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
