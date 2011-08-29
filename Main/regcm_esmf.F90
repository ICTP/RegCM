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

program regcm

!
!**********************************************************************
!
! Used module declarations 
!
!**********************************************************************
!
  use ESMF
!
  use mod_couplerr
  use mod_esmf_atm
  use mod_esmf_ocn
!
!**********************************************************************
!
! Local variable declarations 
!
!**********************************************************************
!
  type(ESMF_VM) :: vm
  type(ESMF_GridComp) :: gcomp
  character(ESMF_MAXSTR), parameter :: conv = 'ESMF'
  character(ESMF_MAXSTR), parameter :: purp = 'General'
!
  type(ESMF_Time) :: time
  type(ESMF_TimeInterval) :: timeInt 
!
  integer :: i, j, rc, localPet, petCount, comm
  character (len=80) :: name, importName, exportName
!
!***********************************************************************
!
! Model Initialization
!
!***********************************************************************
!
!-----------------------------------------------------------------------
! Initialize ESMF
!-----------------------------------------------------------------------
!
  call ESMF_Initialize(vm=vm, rc=rc)
  call CheckError("calling ESMF_Initialize", rc)
!
!-----------------------------------------------------------------------
! Get the VM
!-----------------------------------------------------------------------
!
  call ESMF_VMGet(vm, petCount=petCount, localPet=localPet, rc=rc)
  if (localPet .eq. 0) then
    call ESMF_VMPrint(vm, rc=rc)
  end if
  call CheckError("calling ESMF_VMGet", rc)
!
!-----------------------------------------------------------------------
! Allocate coupler variables
!-----------------------------------------------------------------------
!
!  call allocate_cplvars(petCount, rc)
!
!**********************************************************************
!
! Model Finalize 
!
!**********************************************************************
!
!       call RCM_finalize()
!

end program regcm
