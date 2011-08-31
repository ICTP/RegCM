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
      module mod_esmf_atm
!
!***********************************************************************
!
!     Used module declarations 
!
!***********************************************************************
!
      use ESMF
      use mod_regcm_interface
      use mod_couplerr
!
      implicit none
      private
!
!***********************************************************************
!
!     Public subroutines 
!
!***********************************************************************
!
      public  :: RCM_SetServices
      public  :: RCM_SetRun
      public  :: RCM_SetFinalize
!
      contains
!
      subroutine RCM_SetServices(comp, rc)
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_GridComp), intent(inout) :: comp
      integer, intent(out) :: rc
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
!
!***********************************************************************
!
!     Register "initialize" routine
!
!***********************************************************************
!
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_INITIALIZE,&
                                      userRoutine=RCM_SetInitialize,    &
                                      rc=rc)
!
!***********************************************************************
!
!     Register "run" routine
!
!***********************************************************************
!
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_RUN,       &
                                      userRoutine=RCM_SetRun,           &
                                      rc=rc)
!
!***********************************************************************
!
!     Register "finalize" routine
!
!***********************************************************************
!
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_FINALIZE,  &
                                      userRoutine=RCM_SetFinalize,      &
                                      rc=rc)
!
!***********************************************************************
!
!     Set return flag to success 
!
!***********************************************************************
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetServices
!
      subroutine RCM_SetInitialize(comp, importState, exportState,      &
                                   clock, rc)
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_GridComp), intent(inout) :: comp
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out) :: rc
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
      integer :: localPet, petCount, comm, ierr
!
!***********************************************************************
!
!     Call RCM initialization routines
!
!***********************************************************************
!
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
! 
      call ESMF_GridCompGet(comp, vm=models(Iatmos)%vm, rc=rc)
!
      call ESMF_VMGet(models(Iatmos)%vm,                                &
                      localPet=localPet,                                &
                      petCount=petCount,                                &
                      mpiCommunicator=comm,                             &
                      rc=rc)  
!
!-----------------------------------------------------------------------
!     Initialize the gridded component 
!-----------------------------------------------------------------------
!
      call MPI_Comm_dup(comm, models(Iatmos)%comm, ierr)
      call RCM_initialize(mpiCommunicator=models(Iatmos)%comm)
!
!-----------------------------------------------------------------------
!     Set-up ESMF internal clock for gridded component
!-----------------------------------------------------------------------
!
      call RCM_SetClock(clock, rc)
!
!-----------------------------------------------------------------------
!     Set-up excgange grid for gridded component 
!-----------------------------------------------------------------------
!
      call RCM_SetGridArrays(comp, rc)
!
!-----------------------------------------------------------------------
!     Load ROMS exchange grid arrays
!-----------------------------------------------------------------------
!
      call RCM_PutGridData(localPet, rc)
!
!-----------------------------------------------------------------------
!     Set-up import/export states
!-----------------------------------------------------------------------
!
      call RCM_SetStates(ImportState, ExportState, localPet, rc)
!
!-----------------------------------------------------------------------
!     Load export initial conditions data.
!-----------------------------------------------------------------------
!
      call RCM_PutExportData(localPet, rc)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetInitialize
!
      subroutine RCM_SetRun(comp, importState, exportState,             &
                            clock, rc)
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_GridComp), intent(inout) :: comp
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out) :: rc
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
      logical :: first
      integer :: localPet, petCount, comm
!
!***********************************************************************
!
!     Call RCM run routines
!
!***********************************************************************
!
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
! 
      call ESMF_VMGet(models(Iatmos)%vm,                                &
                      localPet=localPet,                                &
                      petCount=petCount,                                &
                      mpiCommunicator=comm,                             &
                      rc=rc)
!
!-----------------------------------------------------------------------
!     Get RCM internal clock current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (models(Iatmos)%clock,                         &
                          currTime=models(Iatmos)%curTime,              &
                          rc=rc)
!
      CALL ESMF_TimeGet (models(Iatmos)%curTime,                        &
                         yy=models(Iatmos)%time%year,                   &
                         mm=models(Iatmos)%time%month,                  &
                         dd=models(Iatmos)%time%day,                    &
                         h=models(Iatmos)%time%hour,                    &
                         m=models(Iatmos)%time%minute,                  &
                         s=models(Iatmos)%time%second,                  &
                         timeZone=models(Iatmos)%time%zone,             &
                         timeStringISOFrac=models(Iatmos)%time%stamp,   &
                         dayOfYear=models(Iatmos)%time%yday,            &
                         rc=rc)
!
!-----------------------------------------------------------------------
!     Write current time (debug)
!-----------------------------------------------------------------------
!
      write(*, 20) localPet, 'Current Time',                            &
                   trim(models(Iatmos)%time%stamp)
!
!-----------------------------------------------------------------------
!     Set start and end time 
!-----------------------------------------------------------------------
!
!      first = .true.
!      timestr = d_zero
!      tdif = idate2 - idate1
!      timeend = tohours(tdif) * secph
!
!      call RCM_run(timestr, timeend, first)
!
!-----------------------------------------------------------------------
!     Formats 
!-----------------------------------------------------------------------
!
 20   format(' PET (', I2, ') - RCM Model ', A, ' = ', A)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetRun
!
      subroutine RCM_SetFinalize(comp, importState, exportState,        &
                                 clock, rc)
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_GridComp), intent(inout) :: comp
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
!     Call RCM finalize routines
!-----------------------------------------------------------------------
!
      call RCM_finalize()
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetFinalize
!    
      subroutine RCM_SetClock(clock, rc)
!
!***********************************************************************
!
!     Imported modules 
!
!***********************************************************************
!
      use mod_runparams
      use mod_date
!
      implicit none
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      TYPE(ESMF_Clock), intent(inout) :: clock
      integer, intent(inout) :: rc 
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
      integer :: ref_year,   str_year,   end_year
      integer :: ref_month,  str_month,  end_month
      integer :: ref_day,    str_day,    end_day
      integer :: ref_hour,   str_hour,   end_hour
      integer :: ref_minute, str_minute, end_minute
      integer :: ref_second, str_second, end_second
      character (len=80) :: str
!
!***********************************************************************
!
!     Create gridded component clock 
!
!***********************************************************************
!
!-----------------------------------------------------------------------
!     Create ESMF calendar
!-----------------------------------------------------------------------
!
      str = 'Mixed Gregorian/Julian calendar'
      models(Iatmos)%cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,  &
                                        name=TRIM(str),                 &
                                        rc=rc)
!
!-----------------------------------------------------------------------
!     Set Reference time.
!-----------------------------------------------------------------------
!
      call split_idate(idate0, ref_year, ref_month, ref_day, ref_hour)
      ref_minute = 0
      ref_second = 0
!
      call ESMF_TimeSet (models(Iatmos)%refTime,                        &
                         yy=ref_year,                                   &
                         mm=ref_month,                                  &
                         dd=ref_day,                                    &
                         h=ref_hour,                                    &
                         m=ref_minute,                                  &
                         s=ref_second,                                  &
                         calendar=models(Iatmos)%cal,                   &
                         rc=rc)
!
!-----------------------------------------------------------------------
!     Set start time
!-----------------------------------------------------------------------
!
      call split_idate(idate1, str_year, str_month, str_day, str_hour)
      str_minute = 0
      str_second = 0
!
      call ESMF_TimeSet (models(Iatmos)%strTime,                        &
                         yy=str_year,                                   &
                         mm=str_month,                                  &
                         dd=str_day,                                    &
                         h=str_hour,                                    &
                         m=str_minute,                                  &
                         s=str_second,                                  &
                         calendar=models(Iatmos)%cal,                   &
                         rc=rc)
!
!-----------------------------------------------------------------------
!     Set stop time
!-----------------------------------------------------------------------
!
      call split_idate(idate2, end_year, end_month, end_day, end_hour)
      end_minute = 0
      end_second = 0
!
      call ESMF_TimeSet (models(Iatmos)%endTime,                        &
                         yy=end_year,                                   &
                         mm=end_month,                                  &
                         dd=end_day,                                    &
                         h=end_hour,                                    &
                         m=end_minute,                                  &
                         s=end_second,                                  &
                         calendar=models(Iatmos)%cal,                   &
                         rc=rc)
!
!-----------------------------------------------------------------------
!     Set time interval
!-----------------------------------------------------------------------
!
      call ESMF_TimeIntervalSet (models(Iatmos)%dtsec,                  &
                                 s_r8=dtsec,                            &
                                 rc=rc)  
!
!-----------------------------------------------------------------------
!     Create time clock.
!-----------------------------------------------------------------------
!
      str = 'Model clock (Atmosphere)'
      models(Iatmos)%clock = ESMF_ClockCreate (name=TRIM(str),          &
                                  refTime=models(Iatmos)%refTime,       &
                                  timeStep=models(Iatmos)%dtsec,        &
                                  startTime=models(Iatmos)%strTime,     &
                                  stopTime=models(Iatmos)%endTime,      &
                                  rc=rc)
!
!-----------------------------------------------------------------------
!     Copy clock
!-----------------------------------------------------------------------
!
      clock = ESMF_ClockCreate (models(Iatmos)%clock, rc=rc)
!
!-----------------------------------------------------------------------
!     Validate time clock
!-----------------------------------------------------------------------
!
      call ESMF_ClockValidate (models(Iatmos)%clock,                    &
                               rc=rc)
!
!-----------------------------------------------------------------------
!     Get gridded component internal clock current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (models(Iatmos)%clock,                         &
                          currTime=models(Iatmos)%curTime,              &
                          rc=rc)
!
!-----------------------------------------------------------------------
!     Put current time into ESM_Time variable
!-----------------------------------------------------------------------
!
      CALL ESMF_TimeGet (models(Iatmos)%curTime,                        &
                         yy=models(Iatmos)%time%year,                   &
                         mm=models(Iatmos)%time%month,                  &
                         dd=models(Iatmos)%time%day,                    &
                         h=models(Iatmos)%time%hour,                    &
                         m=models(Iatmos)%time%minute,                  &
                         s=models(Iatmos)%time%second,                  &
                         timeZone=models(Iatmos)%time%zone,             &
                         timeStringISOFrac=models(Iatmos)%time%stamp,   &
                         dayOfYear=models(Iatmos)%time%yday,            &
                         rc=rc)
!
!-----------------------------------------------------------------------
!     Set return flag to success
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
      end subroutine RCM_SetClock
!
     subroutine RCM_SetGridArrays (comp, rc)
!
!***********************************************************************
!
!     Used module declarations 
!
!***********************************************************************
!
      use mod_dynparam, only : nproc, iy, iym2, jx, jxp, jendx, jendl
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_GridComp), intent(inout) :: comp
      integer, intent(out) :: rc
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
      integer :: ng, tile, i, j, n
      integer :: localPet, petCount, comm
      integer, allocatable :: deBlockList(:,:,:)
      integer, allocatable :: TLWidth(:,:), TUWidth(:,:)
      integer, allocatable :: CLW_c(:,:), CUW_c(:,:)
      integer, dimension(2) :: CLW, CUW, TLW, TUW
      integer, dimension(2) :: deCount, minIndex, maxIndex
      TYPE (ESMF_ARRAY) :: GrdArray
       real(ESMF_KIND_R8), pointer :: farrayPtr(:,:)
!
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!  
      call ESMF_VMGet (models(Iatmos)%vm,                               &
                       localPet=localPet,                               &
                       petCount=petCount,                               &
                       mpiCommunicator=comm,                            &
                       rc=rc)
!
!-----------------------------------------------------------------------
!     Set RCM domain decomposition variables
!-----------------------------------------------------------------------
!
!     Loop over number of nested/composed grids.
      do n = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Set ESMF Layout and distribution objects
!-----------------------------------------------------------------------
!
      deCount = (/ 1, nproc /)       
      models(Iatmos)%deLayout(n) = ESMF_DELayoutCreate (                &
                                        models(Iatmos)%vm,              &
                                        deCountList=deCount,            &
                                        petList=models(Iatmos)%petList, &
                                        rc=rc)
!
!-----------------------------------------------------------------------
!     Validate and print DELayout
!-----------------------------------------------------------------------
!
!      call ESMF_DELayoutValidate(models(Iatmos)%deLayout(n), rc=rc)
!      call ESMF_DELayoutPrint(models(Iatmos)%deLayout(n), rc=rc)
!
!-----------------------------------------------------------------------
!     Create ESMF DistGrid based on model domain decomposition
!-----------------------------------------------------------------------
!     
      if (.not.allocated(deBlockList)) then
        allocate (deBlockList(2, 2, nproc))
      end if
!
      do tile = 1, nproc 
        deBlockList(1,1,tile) = 1
        deBlockList(1,2,tile) = iym2 
        deBlockList(2,1,tile) = (tile-1)*jxp+1
        deBlockList(2,2,tile) = tile*jxp 
      end do
!      if (localPet == 0) then
        print 10, 'Istr = ',(deBlockList(1,1,tile),tile=1,nproc)
        print 10, 'Iend = ',(deBlockList(1,2,tile),tile=1,nproc)
        print 10, 'Jstr = ',(deBlockList(2,1,tile),tile=1,nproc)
        print 10, 'Jend = ',(deBlockList(2,2,tile),tile=1,nproc)
!      end if
10    format(1x,a,64i5)
!
!     Cordinates of the lower and upper corner of the patch
!
      minIndex = (/ 1, 1 /)
      maxIndex = (/ iym2, jx /)
!
      models(Iatmos)%distGrid(n) = ESMF_DistGridCreate (minIndex,       &
                                   maxIndex,                            &
                                   deBlockList=deBlockList,             &
                                   deLayout=models(Iatmos)%deLayout(n), &
                                   vm=models(Iatmos)%vm,                &
                                   rc=rc)
!
!-----------------------------------------------------------------------
!     Validate and print DistGrid
!-----------------------------------------------------------------------
!
!      call ESMF_DistGridValidate(models(Iatmos)%distGrid(n), rc=rc)
!      call ESMF_DistGridPrint(models(Iatmos)%distGrid(n), rc=rc)
!
!-----------------------------------------------------------------------
!     Set array descriptor
!-----------------------------------------------------------------------
!
      call ESMF_ArraySpecSet (models(Iatmos)%arrSpec(n),                &
                              typekind=ESMF_TYPEKIND_R8,                &
                              rank=2,                                   &
                              rc=rc)
!
!-----------------------------------------------------------------------    
!     Set computational region widths
!-----------------------------------------------------------------------    
!
      if (.not. allocated(CLW_c)) then
        allocate(CLW_c(2, 0:nproc-1))
        allocate(CUW_c(2, 0:nproc-1))
        allocate(TLWidth(2, 0:nproc-1))
        allocate(TUWidth(2, 0:nproc-1))
      end if
!
      do tile = 0, nproc-1
        TLWidth(1,tile) = 1 
        TLWidth(2,tile) = 1
        TUWidth(1,tile) = 1
        TUWidth(2,tile) = 1
        CLW_c(1,tile) = 1       
        CLW_c(2,tile) = 1       
        CUW_c(1,tile) = 1       
        CUW_c(2,tile) = 1
      end do
!
      TLW = (/TLWidth(1,localPet), TLWidth(2,localPet)/)
      TUW = (/TUWidth(1,localPet), TUWidth(2,localPet)/)
      CLW = (/CLW_c(1,localPet), CLW_c(2,localPet)/)
      CUW = (/CUW_c(1,localPet), CUW_c(2,localPet)/)
!
!-----------------------------------------------------------------------    
!     Create exchange arrays
!-----------------------------------------------------------------------    
!
      do i = 1, ubound(models(Iatmos)%mesh, dim=1)
!
!       Create array
!
        grdArray = ESMF_ArrayCreate(arrayspec=models(Iatmos)%arrSpec(n),&
                                    distgrid=models(Iatmos)%distGrid(n),&
!                                    computationalLWidth=CLW,            &
!                                    computationalUWidth=CUW,            &
!                                    totalLWidth=TLW,                    &
!                                    totalUWidth=TUW,                    &
                                    indexflag=ESMF_INDEX_GLOBAL,        &
                                    rc=rc)
        models(Iatmos)%mesh(i,n)%lat%array = grdArray
!
!       Get data pointer from array
!
        call ESMF_ArrayGet (models(Iatmos)%mesh(i,n)%lon%array,           &
                        farrayPtr=models(Iatmos)%mesh(i,n)%lon%field,&
                        rc=rc)
        print*, ubound(models(Iatmos)%mesh(i,n)%lon%field, dim=1),      &
                ubound(models(Iatmos)%mesh(i,n)%lon%field, dim=2)
!
!       Set adjustable settings of array object  
!
        call ESMF_ArraySet (array=models(Iatmos)%mesh(i,n)%lat%array,     &
                           name=trim(models(Iatmos)%mesh(i,n)%lat%name),&
                            rc=rc)
!
!       Create array
!
        grdArray = ESMF_ArrayCreate (arrayspec=models(Iatmos)%arrSpec(n),& 
                               distgrid=models(Iatmos)%distGrid(n),     &
                               indexflag=ESMF_INDEX_GLOBAL,             &
                               rc=rc)
        models(Iatmos)%mesh(i,n)%lon%array = grdArray
!
!       Get data pointer from array
!
        call ESMF_ArrayGet (models(Iatmos)%mesh(i,n)%lon%array,           &
                        farrayPtr=models(Iatmos)%mesh(i,n)%lon%field,&
                        rc=rc)          
!
!       Set adjustable settings of array object  
!
        call ESMF_ArraySet (array=models(Iatmos)%mesh(i,n)%lon%array,     &
                        name=trim(models(Iatmos)%mesh(i,n)%lon%name),&
                        rc=rc)
!  
!       Add array to export state
!   
        call ESMF_StateAdd (models(Iatmos)%stateExport,                 &
                           (/ models(Iatmos)%mesh(i,n)%lat%array,         &
                              models(Iatmos)%mesh(i,n)%lon%array /),      &
                           rc=rc)
!
!       Initialize the grid array
!
      print*, "1L ", i,n, localPet, lbound(models(Iatmos)%mesh(i,n)%lat%field, dim=1)
      print*, "1U ", i,n, localPet, ubound(models(Iatmos)%mesh(i,n)%lat%field, dim=1)
      print*, "2L ", i,n, localPet, lbound(models(Iatmos)%mesh(i,n)%lat%field, dim=2)
      print*, "2U ", i,n, localPet, ubound(models(Iatmos)%mesh(i,n)%lat%field, dim=2)
        models(Iatmos)%mesh(i,n)%lat%field = 0.0d0
        models(Iatmos)%mesh(i,n)%lon%field = 0.0d0
      end do 
      end do
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetGridArrays
!
      subroutine RCM_PutGridData (localPet, rc)
!
!**********************************************************************
!
!     Used module declarations 
!
!**********************************************************************
!
      use mod_dynparam, only : iy, iym2, jxp
      use mod_atm_interface, only : mddom 
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      integer, intent(in) :: localPet
      integer, intent(inout) :: rc
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
      integer :: i, j, n, jj
!     
!-----------------------------------------------------------------------
!     Load grid data for cross points 
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iatmos)
      print*, "1L ", n, localPet, lbound(models(Iatmos)%mesh(1,n)%lat%field, dim=1)
      print*, "1U ", n, localPet, ubound(models(Iatmos)%mesh(1,n)%lat%field, dim=1)
      print*, "2L ", n, localPet, lbound(models(Iatmos)%mesh(1,n)%lat%field, dim=2)
      print*, "2U ", n, localPet, ubound(models(Iatmos)%mesh(1,n)%lat%field, dim=2)
      print*, "*1L ", n, localPet, lbound(mddom%xlon)
      print*, "*1U ", n, localPet, ubound(mddom%xlon)
      print*, "*2L ", n, localPet, lbound(mddom%xlon)
      print*, "*2U ", n, localPet, ubound(mddom%xlon)
      do j = 1, jxp
        do i = 1, iym2
          !jj = jxp*localPet+j
          models(Iatmos)%mesh(1,n)%lat%field(i,j) = mddom%xlat(i,j)
          models(Iatmos)%mesh(1,n)%lon%field(i,j) = mddom%xlon(i,j)
        end do
      end do 
      end do
!      call ESMF_ArrayPrint(models(Iatmos)%grid(1)%lat%array, rc=status)
!      call ESMF_ArrayPrint(models(Iatmos)%grid(1)%lon%array, rc=status)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_PutGridData
!
      subroutine RCM_SetStates (importState, exportState, localPet,     &
                                rc)
!
!**********************************************************************
!
!     Used module declarations 
!
!**********************************************************************
!
      use mod_dynparam, only : nproc, iy, jx, jxp, jendl
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      integer, intent(in) :: localPet
      integer, intent(inout) :: rc
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
      integer :: i, j, n
!
!-----------------------------------------------------------------------
!     Initialize the import and export fields 
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Create export state arrays.
!-----------------------------------------------------------------------
!
      do i = 1, ubound(models(Iatmos)%dataExport(:,n), dim=1)
!
!       Create ESMF array
!
        models(Iatmos)%dataExport(i,n)%array = ESMF_ArrayCreate (       &
                                    arrayspec=models(Iatmos)%arrSpec(n),&
                                    distgrid=models(Iatmos)%distGrid(n),&
                                    rc=rc)
!
!       Get data pointer from ESMF array 
!
        call ESMF_ArrayGet (models(Iatmos)%dataExport(i,n)%array,       &
                       farrayPtr=models(Iatmos)%dataExport(i,n)%field,  &
                       rc=rc)
!
!       Set array name
        call ESMF_ArraySet (array=models(Iatmos)%dataExport(i,n)%array, &
                        name=trim(models(Iatmos)%dataExport(i,n)%name), &
                        rc=rc)
!
!       Add array to export state
!
        call ESMF_StateAdd (models(Iatmos)%stateExport,                 &
                           (/ models(Iatmos)%dataExport(i,n)%array /),  &
                           rc=rc)
!
!       Initialize export field to zero to avoid infinities or NaNs.
!
        models(Iatmos)%dataExport(i,n)%field = 0.0d0
      end do
!     
!-----------------------------------------------------------------------
!     Create import state arrays.
!-----------------------------------------------------------------------
!
!      print*, "**", ubound(models(Iatmos)%dataImport(:,n), dim=1), "**"
!      print*, "**", ubound(models(Iatmos)%dataExport(:,n), dim=1), "**"
      do i = 1, ubound(models(Iatmos)%dataImport(:,n), dim=1)
!
!       Create ESMF array
!        
        models(Iatmos)%dataImport(i,n)%array = ESMF_ArrayCreate (       &
                                    arrayspec=models(Iatmos)%arrSpec(n),&
                                    distgrid=models(Iatmos)%distGrid(n),&
                                    rc=rc)
!
!       Get data pointer from ESMF array 
!
        call ESMF_ArrayGet (models(Iatmos)%dataImport(i,n)%array,       &
                       farrayPtr=models(Iatmos)%dataImport(i,n)%field,  &
                       rc=rc)
!
!       Set array name
        call ESMF_ArraySet (array=models(Iatmos)%dataImport(i,n)%array, &
                        name=trim(models(Iatmos)%dataImport(i,n)%name), &
                        rc=rc)
!
!       Add array to export state
!
        call ESMF_StateAdd (models(Iatmos)%stateImport,                 &
                           (/ models(Iatmos)%dataImport(i,n)%array /),  &
                           rc=rc)
!
!       Initialize export field to zero to avoid infinities or NaNs.
!
        models(Iatmos)%dataImport(i,n)%field = 0.0d0
      end do
      end do
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetStates
!
      subroutine RCM_PutExportData (localPet, rc)
!
!**********************************************************************
!
!     Used module declarations 
!
!**********************************************************************
!
      use mod_bats_common !, only : fbat
      use mod_dynparam, only : iy, iym2, jxp
!
!***********************************************************************
!
!     Imported variable declarations 
!     
!***********************************************************************
!
      integer, intent(in) :: localPet
      integer, intent(inout) :: rc
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
      character (len=80) :: name
      integer :: i, j, k, n
!
!-----------------------------------------------------------------------
!     Put data into ESMF arrays 
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Load export fields
!-----------------------------------------------------------------------
!
      do k = 1, ubound(models(Iatmos)%dataExport(:,n), dim=1)
        name = models(Iatmos)%dataExport(i,n)%name
        if (trim(adjustl(name)) == "Pair") then
          do j = 1 , jxp
            do i = 1 , iym2
              models(Iatmos)%dataExport(k,n)%field = (sfps(i,j)+ptop)*  &
                                                      d_1000 
            end do        
          end do
        else if (trim(adjustl(name)) == "Tair") then
          do j = 1 , jxp
            do i = 1 , iym2
              models(Iatmos)%dataExport(k,n)%field = thatm(i,kz,j) 
            end do
          end do
        else if (trim(adjustl(name)) == "Qair") then
          do j = 1 , jxp
            do i = 1 , iym2
              models(Iatmos)%dataExport(k,n)%field = qvatm(i,kz,j)/     &
                                                  (d_one+qvatm(i,kz,j)) 
            end do
          end do
        else if (trim(adjustl(name)) == "Uwind") then
          do j = 1 , jxp
            do i = 1 , iym2
              models(Iatmos)%dataExport(k,n)%field = uatm(i,kz,j) 
            end do
          end do
        else if (trim(adjustl(name)) == "Vwind") then
          do j = 1 , jxp
            do i = 1 , iym2
              models(Iatmos)%dataExport(k,n)%field = vatm(i,kz,j)
            end do
          end do
        else if (trim(adjustl(name)) == "rain") then
          do j = 1 , jxp
            do i = 1 , iym2
              models(Iatmos)%dataExport(k,n)%field = pptnc(i,j) +       &
                                                     pptc(i,j)
            end do
          end do
        else if (trim(adjustl(name)) == "swrad") then
          do j = 1 , jxp
            do i = 1 , iym2
              models(Iatmos)%dataExport(k,n)%field = fsw2d(i,j) 
            end do
          end do
        else if (trim(adjustl(name)) == "lwrad_down") then
          do j = 1 , jxp
            do i = 1 , iym2
              models(Iatmos)%dataExport(k,n)%field = flw2d(i,j) 
            end do
          end do
        end if
      end do
      end do
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_PutExportData
!
      end module mod_esmf_atm
