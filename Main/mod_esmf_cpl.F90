!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     This file is part of ICTP RegCM.
!
!     ICTP RegCM is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     ICTP RegCM is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      
!
      module mod_esmf_cpl
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use ESMF
      use mod_couplerr
! 
      implicit none
      private
!
!-----------------------------------------------------------------------
!     Public subroutines 
!-----------------------------------------------------------------------
!
      public  :: CPL_SetServices
      public  :: CPL_SetRun
      public  :: CPL_SetFinalize
!
      contains
!
      subroutine CPL_SetServices(comp, rc)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_CplComp), intent(inout) :: comp
      integer, intent(out) :: rc 
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Register "initialize" routine
!-----------------------------------------------------------------------
!
      call ESMF_CplCompSetEntryPoint(comp,                              &
                                     methodflag=ESMF_METHOD_INITIALIZE, &
                                     userRoutine=CPL_SetInitialize,     &
                                     rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Register "run" routine
!-----------------------------------------------------------------------
!
      call ESMF_CplCompSetEntryPoint(comp,                              &
                                     methodflag=ESMF_METHOD_RUN,        &
                                     userRoutine=CPL_SetRun,            &
                                     rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Register "finalize" routine
!-----------------------------------------------------------------------
!
      call ESMF_CplCompSetEntryPoint(comp,                              &
                                     methodflag=ESMF_METHOD_FINALIZE,   &
                                     userRoutine=CPL_SetFinalize,       &
                                     rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set return flag to success 
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine CPL_SetServices
!
      subroutine CPL_SetInitialize(comp, importState, exportState,     &
                                   clock, rc)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_CplComp), intent(inout) :: comp
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      logical :: flag
      type(ESMF_Config) :: config
      integer :: localPet, petCount, comm, ierr
      integer :: i, j, dir
!
      integer :: itemCount
      character(ESMF_MAXSTR), allocatable :: itemNames(:)
      type(ESMF_VM) :: vm_cpl
      type(ESMF_StateItem_Flag), allocatable :: itemTypes(:)
      type(ESMF_Field) :: dstField, srcField
!
      real(ESMF_KIND_R8), pointer :: ptr(:,:)
!
!-----------------------------------------------------------------------
!     Call ROMS initialization routines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!     
      call ESMF_CplCompGet(comp,                                        &
                           vm=vm_cpl,                                   &
                           rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_VMGet(vm_cpl,                                           &
                      localPet=localPet,                                &
                      petCount=petCount,                                &
                      mpiCommunicator=comm,                             &
                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Reconcile import and export states. Consistent view in all PETs.
!-----------------------------------------------------------------------
!
!     Import state
!
      call ESMF_StateReconcile(importState,                             &
                               vm=vm_cpl,                               &
                               attreconflag=ESMF_ATTRECONCILE_ON,       &
                               rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
! 
!     Export state
!
      call ESMF_StateReconcile(exportState,                             &
                               vm=vm_cpl,                               &
                               attreconflag=ESMF_ATTRECONCILE_ON,       &
                               rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get direction of coupling initialization
!-----------------------------------------------------------------------
!       
      call ESMF_AttributeGet(importState,                               &
                             name=trim(DIRECTION),                      &
                             value=dir,                                 &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!     Print coupling direction info (for debugging)
!  
      if (dir == FORWARD_ON) then
        write(*,fmt="(' PET (', I2, ') Direction = Forward ')") localPet
      else
        write(*,fmt="(' PET (', I2, ') Direction = Backward')") localPet
      end if
!
!-----------------------------------------------------------------------
!     Get and save import state field names 
!-----------------------------------------------------------------------
      call ESMF_StateGet(importState,                                   &
                         itemCount=itemCount,                           &
                         rc=rc) 
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Allocate temporary arrays 
!-----------------------------------------------------------------------
!
      if (.not. allocated(itemNames)) allocate(itemNames(itemCount))
      if (.not. allocated(itemTypes)) allocate(itemTypes(itemCount))
! 
!-----------------------------------------------------------------------
!     Get import state item names and types 
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(importState,                                   &
                         itemNameList=itemNames,                        &
                         itemTypeList=itemTypes,                        &
                         rc=rc)      
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get required item count 
!-----------------------------------------------------------------------
!
      j = 0
      do i = 1, ItemCount
        if ((itemTypes(i) == ESMF_STATEITEM_FIELD) .or.                 &
            (itemTypes(i) == ESMF_STATEITEM_ARRAY)) then
          j = j+1
          write(*,30) localPet, j, '>'//trim(itemNames(j))//'<'
        end if
      end do
 30   format(' PET (', I2, ') - Import Item (',I2,') = ',A)
!
!-----------------------------------------------------------------------
!     Save import state field names
!-----------------------------------------------------------------------
!
      if (dir == FORWARD_ON) then
        if (.not. allocated(itemNamesImportF)) then
          allocate(itemNamesImportF(j))
        end if
      else
        if (.not. allocated(itemNamesImportB)) then
          allocate(itemNamesImportB(j))
        end if
      end if
!
      j = 0
      do i = 1, ItemCount
        if ((itemTypes(i) == ESMF_STATEITEM_FIELD) .or.                 &
            (itemTypes(i) == ESMF_STATEITEM_ARRAY)) then
          j = j+1
          if (dir == FORWARD_ON) then          
            itemNamesImportF(j) = trim(itemNames(i))
          else
            itemNamesImportB(j) = trim(itemNames(i))
          end if
        end if
      end do
!
!-----------------------------------------------------------------------
!     Deallocate temporary arrays 
!-----------------------------------------------------------------------
!
      if (allocated(itemNames)) deallocate(itemNames)
      if (allocated(itemTypes)) deallocate(itemTypes)
!
!-----------------------------------------------------------------------
!     Get and save export state field names 
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(exportState,                                   &
                         itemCount=itemCount,                           &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Allocate temporary arrays 
!-----------------------------------------------------------------------
!
      if (.not. allocated(itemNames)) allocate(itemNames(itemCount))
      if (.not. allocated(itemTypes)) allocate(itemTypes(itemCount))
! 
!-----------------------------------------------------------------------
!     Get export state item names and types 
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(exportState,                                   &
                         itemNameList=itemNames,                        &
                         itemTypeList=itemTypes,                        &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get required item count 
!-----------------------------------------------------------------------
!
      j = 0
      do i = 1, ItemCount
        if ((itemTypes(i) == ESMF_STATEITEM_FIELD) .or.                 &
            (itemTypes(i) == ESMF_STATEITEM_ARRAY)) then
          j = j+1
          write(*,40) localPet, j, '>'//trim(itemNames(j))//'<'
        end if
      end do
 40   format(' PET (', I2, ') - Export Item (',I2,') = ',A)
!
!-----------------------------------------------------------------------
!     Save export state field names
!-----------------------------------------------------------------------
!
      if (dir == FORWARD_ON) then
        if (.not. allocated(itemNamesExportF)) then
          allocate(itemNamesExportF(j))
        end if
      else
        if (.not. allocated(itemNamesExportB)) then
          allocate(itemNamesExportB(j))
        end if
      end if
!
      j = 0
      do i = 1, ItemCount
        if ((itemTypes(i) == ESMF_STATEITEM_FIELD) .or.                 &
            (itemTypes(i) == ESMF_STATEITEM_ARRAY)) then
          j = j+1
          if (dir == FORWARD_ON) then
            itemNamesExportF(j) = trim(itemNames(i))
          else
            itemNamesExportB(j) = trim(itemNames(i))
          end if
        end if
      end do
!
!-----------------------------------------------------------------------
!     Deallocate temporary arrays 
!-----------------------------------------------------------------------
!
      if (allocated(itemNames)) deallocate(itemNames)
      if (allocated(itemTypes)) deallocate(itemTypes)
!
!-----------------------------------------------------------------------
!     Forward coupling initialization
!-----------------------------------------------------------------------
!
      if (dir == FORWARD_ON) then
!
!-----------------------------------------------------------------------
!     Compute weight matrix between gridded component grids. It creates 
!     a sparse matrix operation (stored in routehandle) that contains 
!     the calculations and communications necessary to interpolate from 
!     source field to destination field.
!-----------------------------------------------------------------------
!                
      do i = 1, size(itemNamesImportF, dim=1)
!
!-----------------------------------------------------------------------
!     Get import field
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(importState,                                   &
                         trim(itemNamesImportF(i)),                     &
                         srcField,                                      &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get export field
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(exportState,                                   &
                         trim(itemNamesExportF(i)),                     &
                         dstField,                                      &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create ESMF routhandle 
!-----------------------------------------------------------------------
!
      call ESMF_FieldRegridStore (srcField=srcField,                    &
                                dstField=dstField,                      &
                                routeHandle=routeHandleF,               &
!                                indices=indices,                        &
!                                weights=weights,                        &
                                regridmethod=ESMF_REGRIDMETHOD_BILINEAR,&
                                rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Regrid fields 
!-----------------------------------------------------------------------
!
!      call ESMF_FieldRegrid(srcField, dstField, routeHandleF, rc=rc)
!      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Write field to NetCDF (debug)
!-----------------------------------------------------------------------
!
!      flag = .false.
!      if (any(models(Iatmos)%petList == localpet)) flag = .true.
!      if (flag) then
!      call ESMF_FieldWrite(srcField, 'src_field2.nc', rc=rc)
!      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!      else
!      call ESMF_FieldWrite(dstField, 'dst_field2.nc', rc=rc)
!      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!      end if
      end do
!
!-----------------------------------------------------------------------
!     Backward coupling initialization
!-----------------------------------------------------------------------
!
      else
!
!-----------------------------------------------------------------------
!     Compute weight matrix between gridded component grids. It creates 
!     a sparse matrix operation (stored in routehandle) that contains 
!     the calculations and communications necessary to interpolate from 
!     source field to destination field.
!-----------------------------------------------------------------------
!                
      do i = 1, size(itemNamesImportB, dim=1)
!
!-----------------------------------------------------------------------
!     Get import field
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(importState,                                   &
                         trim(itemNamesImportB(i)),                     &
                         srcField,                                      &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get export field
!-----------------------------------------------------------------------
!
      call ESMF_StateGet(exportState,                                   &
                         trim(itemNamesExportB(i)),                     &
                         dstField,                                      &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create ESMF routhandle 
!
!     The ESMF_UNMAPPEDACTION_IGNORE flag is set because ocean model 
!     grid is smaller than atmosphere model grid.
!-----------------------------------------------------------------------
!
      call ESMF_FieldRegridStore (srcField=srcField,                    &
                              dstField=dstField,                        &
                              unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                              routeHandle=routeHandleB,                 &
!                              indices=indices,                          &
!                              weights=weights,                          &
                              regridmethod=ESMF_REGRIDMETHOD_BILINEAR,  &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Regrid fields 
!-----------------------------------------------------------------------
!
!      call ESMF_FieldRegrid(srcField, dstField, routeHandleB, rc=rc)
!      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Write field to NetCDF (debug)
!-----------------------------------------------------------------------
!
!      flag = .false.
!      if (any(models(Iatmos)%petList == localpet)) flag = .true.
!      if (flag) then
!      call ESMF_FieldWrite(srcField, 'src_field2.nc', rc=rc)
!      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!      else
!      call ESMF_FieldWrite(dstField, 'dst_field2.nc', rc=rc)
!      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!      end if
      end do
      end if
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine CPL_SetInitialize
!
      subroutine CPL_SetRun(comp, importState, exportState, clock, rc)
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_CplComp), intent(inout) :: comp
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Call CPL initialization routines
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine CPL_SetRun
!
      subroutine CPL_SetFinalize(comp, importState, exportState,       &
                                  clock, rc)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_CplComp), intent(inout) :: comp
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Initialize return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!     Terminate CPL execution.  Close all NetCDF files.
!-----------------------------------------------------------------------
!
      end subroutine CPL_SetFinalize
!
      end module mod_esmf_cpl 
