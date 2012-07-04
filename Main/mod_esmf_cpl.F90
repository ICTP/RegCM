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
      integer :: localPet, petCount, comm
      integer :: i, j, n, itype, vid, sid, did
!
      integer :: itemCount
      character(ESMF_MAXSTR), allocatable :: itemNames(:)
      type(ESMF_StateItem_Flag), allocatable :: itemTypes(:)
      type(ESMF_Field) :: dstField, srcField
      type(ESMF_Grid) :: dstGrid, srcGrid
      type(ESMF_ArraySpec) :: arrSpec  
      type(ESMF_RegridMethod_Flag) :: regridMethod
!
!-----------------------------------------------------------------------
!     Call ROMS initialization routines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!     
      call ESMF_VMGet(cplVM,                                            &
                      localPet=localPet,                                &
                      petCount=petCount,                                &
                      mpiCommunicator=comm,                             &
                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Reconcile import and export states. Consistent view in all PETs.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Import state
!-----------------------------------------------------------------------
!
      call ESMF_StateReconcile(importState,                             &
                               vm=cplVM,                                &
                               rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
! 
!-----------------------------------------------------------------------
!     Export state
!-----------------------------------------------------------------------
!
      call ESMF_StateReconcile(exportState,                             &
                               vm=cplVM,                                &
                               rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: print coupling direction
!-----------------------------------------------------------------------
!  
      if ((cpl_dbglevel > 0) .and. (localPet == 0)) then
      if (DIRECTION == FORWARD_ON) then
        write(*,fmt="(' PET (', I2, ') Direction = Forward ')") localPet
      else
        write(*,fmt="(' PET (', I2, ') Direction = Backward')") localPet
      end if
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
          if (cpl_dbglevel > 0 .and. localPet == 0) then
            write(*,30) localPet, j, '>'//trim(itemNames(j))//'<'
          end if
        end if
      end do
!
!-----------------------------------------------------------------------
!     Save import state field names
!-----------------------------------------------------------------------
!
      if (DIRECTION == FORWARD_ON) then
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
          if (DIRECTION == FORWARD_ON) then          
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
          if (cpl_dbglevel > 0 .and. localPet == 0) then
            write(*,40) localPet, j, '>'//trim(itemNames(j))//'<'
          end if
        end if
      end do
!
!-----------------------------------------------------------------------
!     Save export state field names
!-----------------------------------------------------------------------
!
      if (DIRECTION == FORWARD_ON) then
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
          if (DIRECTION == FORWARD_ON) then
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
      if (DIRECTION == FORWARD_ON) then
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
!     Get import field mesh id
!-----------------------------------------------------------------------
!
      vid = getVarID(models(Iatmos)%dataExport(:,1),itemNamesImportF(i))
      if (vid < 0) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      sid = models(Iatmos)%dataExport(vid,1)%gtype
!
!-----------------------------------------------------------------------
!     Get interpolation type 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (srcField,                                 &
                              name='interpolation_type',                &
                              value=itype,                              &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create frac field if the interpolation type is conservative 
!-----------------------------------------------------------------------
!
      if (itype == Iconsv) then   
      call ESMF_FieldGet(srcField,                                      &
                         arrayspec=arrSpec,                             &
                         grid=srcGrid,                                  &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      fracFieldFS = ESMF_FieldCreate(srcGrid,                           &
                                     arrSpec,                           &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, &
                                     rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
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
!     Get export field mesh id
!-----------------------------------------------------------------------
!
      vid = getVarID(models(Iocean)%dataImport(:,1),itemNamesExportF(i))
      if (vid < 0) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      did = models(Iocean)%dataImport(vid,1)%gtype
!
      if ((cpl_dbglevel > 0) .and. (localPet == 0)) then
        print*, 'Init Forward - ',                                      &
                trim(models(Iocean)%dataImport(vid,1)%name), ' / ',     &
                trim(itemNamesExportF(i)), ' / ',                       &
                trim(GRIDDES(sid)), ' --> ',                            &
                trim(GRIDDES(did))
      end if
!
!-----------------------------------------------------------------------
!     Create frac field if the interpolation type is conservative 
!-----------------------------------------------------------------------
!
      if (itype == Iconsv) then
      call ESMF_FieldGet(dstField,                                      &
                         arrayspec=arrSpec,                             &
                         grid=dstGrid,                                  &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      fracFieldFD = ESMF_FieldCreate(dstGrid,                           &
                                     arrSpec,                           &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, &
                                     rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
!-----------------------------------------------------------------------
!     Create ESMF routhandle (atm --> ocn) 
!-----------------------------------------------------------------------
!
      if (itype == Iconsv) then
      regridMethod = ESMF_REGRIDMETHOD_CONSERVE
      call ESMF_FieldRegridStore (srcField=srcField,                    &
                                  dstField=dstField,                    &
                                  srcFracField=fracFieldFS,             &
                                  dstFracField=fracFieldFD,             &
                                  routeHandle=routeHandleFC(sid,did),   &
                                  regridmethod=regridMethod,            &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      else
      regridMethod = ESMF_REGRIDMETHOD_BILINEAR
      call ESMF_FieldRegridStore (srcField=srcField,                    &
                                  dstField=dstField,                    &
                                  routeHandle=routeHandleFB(sid,did),   &
                                  regridmethod=regridMethod,            &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
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
!     Get import field mesh id
!-----------------------------------------------------------------------
!
      vid = getVarID(models(Iocean)%dataExport(:,1),itemNamesImportB(i))
      sid = models(Iocean)%dataExport(vid,1)%gtype
!
!-----------------------------------------------------------------------
!     Get interpolation type 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (srcField,                                 &
                              name='interpolation_type',                &
                              value=itype,                              &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create frac field if the interpolation type is conservative 
!-----------------------------------------------------------------------
!
      if (itype == Iconsv) then
      call ESMF_FieldGet(srcField,                                      &
                         arrayspec=arrSpec,                             &
                         grid=srcGrid,                                  &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      fracFieldBS = ESMF_FieldCreate(srcGrid,                           &
                                     arrSpec,                           &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, &
                                     rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
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
!     Get export field mesh id
!-----------------------------------------------------------------------
!
      vid = getVarID(models(Iatmos)%dataImport(:,1),itemNamesExportB(i))
      did = models(Iatmos)%dataImport(vid,1)%gtype
!
      if ((cpl_dbglevel > 0) .and. (localPet == 0)) then
        print*, 'Init Backward - ',                                     &
        trim(models(Iatmos)%dataImport(vid,1)%name), ' / ',             &
        trim(itemNamesExportB(i)), ' / ',                               &
        trim(GRIDDES(did)), ' --> ',                                    &
        trim(GRIDDES(sid))
      end if
!     
!-----------------------------------------------------------------------
!     Create frac field if the interpolation type is conservative 
!-----------------------------------------------------------------------
!
      if (itype == Iconsv) then
      call ESMF_FieldGet(dstField,                                      &
                         arrayspec=arrSpec,                             &
                         grid=dstGrid,                                  &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      fracFieldBD = ESMF_FieldCreate(dstGrid,                           &
                                     arrSpec,                           &
                                     staggerloc=ESMF_STAGGERLOC_CENTER, &
                                     rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
!-----------------------------------------------------------------------
!     Create ESMF routhandle (ocn --> atm)
!
!     The ESMF_UNMAPPEDACTION_IGNORE flag is set because ocean model 
!     grid is smaller than atmosphere model grid.
!-----------------------------------------------------------------------
!
      if (itype == Iconsv) then
      regridMethod = ESMF_REGRIDMETHOD_CONSERVE
      call ESMF_FieldRegridStore (srcField=srcField,                    &
                                  dstField=dstField,                    &
                                  srcFracField=fracFieldBS,             &
                                  dstFracField=fracFieldBD,             &
                                  routeHandle=routeHandleBC(did,sid),   &
                                  regridmethod=regridMethod,            &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      else
      regridMethod = ESMF_REGRIDMETHOD_BILINEAR 
      call ESMF_FieldRegridStore (srcField=srcField,                    &
                              dstField=dstField,                        &
                              srcMaskValues=(/0/),                      &
                              unmappedaction=ESMF_UNMAPPEDACTION_IGNORE,&
                              routeHandle=routeHandleBB(did,sid),       &
                              regridmethod=regridMethod,                &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
      end do
      end if
!
!-----------------------------------------------------------------------
!     Formats 
!-----------------------------------------------------------------------
!
 30   format(' PET (', I2, ') - Import Item (',I2,') = ',A)
 40   format(' PET (', I2, ') - Export Item (',I2,') = ',A)
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
!
      logical :: flag
      integer :: i, j, itype, vid, sid, did
      integer :: localPet, petCount, comm
!
      type(ESMF_Field) :: dstField, srcField
!
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!     
      call ESMF_VMGet(cplVM,                                            &
                      localPet=localPet,                                &
                      petCount=petCount,                                &
                      mpiCommunicator=comm,                             &
                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Reconcile import and export states. Consistent view in all PETs.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Import state
!-----------------------------------------------------------------------
!
      call ESMF_StateReconcile(importState,                             &
                               vm=cplVM,                                &
                               rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
! 
!-----------------------------------------------------------------------
!     Export state
!-----------------------------------------------------------------------
!
      call ESMF_StateReconcile(exportState,                             &
                               vm=cplVM,                                &
                               rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: print coupling direction
!-----------------------------------------------------------------------
!  
      if (cpl_dbglevel > 0 .and. localPet == 0) then
      if (DIRECTION == FORWARD_ON) then
        write(*,fmt="(' PET (', I2, ') Direction = Forward ')") localPet
      else
        write(*,fmt="(' PET (', I2, ') Direction = Backward')") localPet
      end if
      end if
!
!-----------------------------------------------------------------------
!     Forward coupling run 
!-----------------------------------------------------------------------
!
      if (DIRECTION == FORWARD_ON) then
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
!     Get export field mesh id
!-----------------------------------------------------------------------
!
      vid = getVarID(models(Iatmos)%dataExport(:,1),itemNamesImportF(i))
      sid = models(Iatmos)%dataExport(vid,1)%gtype
!
!-----------------------------------------------------------------------
!     Get interpolation type 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (srcField,                                 &
                              name='interpolation_type',                &
                              value=itype,                              &
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
!     Get export field mesh id
!-----------------------------------------------------------------------
!
      vid = getVarID(models(Iocean)%dataImport(:,1),itemNamesExportF(i))
      did = models(Iocean)%dataImport(vid,1)%gtype
!
      if ((cpl_dbglevel > 0) .and. (localPet == 0)) then
        print*, 'Run Forward  - ',                                      &
                trim(models(Iocean)%dataImport(vid,1)%name), ' / ',     &
                trim(itemNamesImportF(i)), ' / ',                       &
                trim(GRIDDES(sid)), ' --> ',                            &
                trim(GRIDDES(did))
      end if
!
!-----------------------------------------------------------------------
!     Regrid fields 
!-----------------------------------------------------------------------
!
      if (itype == Iconsv) then
      call ESMF_FieldRegrid(srcField,                                   &
                            dstField,                                   &
                            routeHandleFC(sid,did),                     &
                            rc=rc)
      else
      call ESMF_FieldRegrid(srcField,                                   &
                            dstField,                                   &
                            routeHandleFB(sid,did),                     &
                            rc=rc)
      end if
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      end do
!
!-----------------------------------------------------------------------
!     Backward coupling run 
!-----------------------------------------------------------------------
!
      else
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
!     Get import field mesh id
!-----------------------------------------------------------------------
!
      vid = getVarID(models(Iocean)%dataExport(:,1),itemNamesImportB(i))
      sid = models(Iocean)%dataExport(vid,1)%gtype
!
!-----------------------------------------------------------------------
!     Get interpolation type 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (srcField,                                 &
                              name='interpolation_type',                &
                              value=itype,                              &
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
!     Get export field mesh id
!-----------------------------------------------------------------------
!
      vid = getVarID(models(Iatmos)%dataImport(:,1),itemNamesExportB(i))
      did = models(Iatmos)%dataImport(vid,1)%gtype
!
      if ((cpl_dbglevel > 0) .and. (localPet == 0)) then
        print*, 'Run Backward  - ',                                     &
        trim(models(Iatmos)%dataImport(vid,1)%name), ' / ',             &
        trim(itemNamesImportB(i)), ' / ',                               &
        trim(GRIDDES(did)), ' --> ',                                    &
        trim(GRIDDES(sid))
      end if
!
!-----------------------------------------------------------------------
!     Regrid fields 
!-----------------------------------------------------------------------
!
      if (itype == Iconsv) then
      call ESMF_FieldRegrid(srcField,                                   &
                            dstField,                                   &
                            routeHandleBC(did,sid),                     &
                            zeroregion=ESMF_REGION_SELECT,              &
                            rc=rc)
      else
      call ESMF_FieldRegrid(srcField,                                   &
                            dstField,                                   &
                            routeHandleBB(did,sid),                     &
                            zeroregion=ESMF_REGION_SELECT,              &
                            rc=rc)
      end if
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      end do
      end if
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine CPL_SetRun
!
      subroutine CPL_SetFinalize (comp, importState, exportState,       &
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
      integer :: i, j
!
!-----------------------------------------------------------------------
!     Call ESMF finalize routines
!-----------------------------------------------------------------------
!
      call ESMF_ClockDestroy(clock, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompDestroy(comp, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      do i = 1, ubound(models(Iatmos)%mesh, dim=1)
        do j = 1, ubound(models(Iocean)%mesh, dim=1)
          call ESMF_FieldBundleSMMRelease(routehandleFB(i,j), rc=rc)
          call ESMF_FieldBundleSMMRelease(routehandleFC(i,j), rc=rc)
          call ESMF_FieldBundleSMMRelease(routehandleBB(i,j), rc=rc)
          call ESMF_FieldBundleSMMRelease(routehandleBC(i,j), rc=rc)
        end do
      end do
!
      call ESMF_FieldDestroy(fracFieldFS, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_FieldDestroy(fracFieldFD, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_FieldDestroy(fracFieldBS, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_FieldDestroy(fracFieldBD, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine CPL_SetFinalize
!
      end module mod_esmf_cpl 
