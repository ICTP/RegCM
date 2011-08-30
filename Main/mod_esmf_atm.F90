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
      integer :: localPet, petCount, comm
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
      call RCM_initialize(mpiCommunicator=comm)
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
!      call RCM_SetStates(ImportState, ExportState, MyRank, status)
!
!-----------------------------------------------------------------------
!     Load export initial conditions data.
!-----------------------------------------------------------------------
!
!      call RCM_PutExportData(Myrank, status)
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
!     
!
!***********************************************************************
!
!     Call RCM initialization routines
!
!***********************************************************************
!
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
      call split_idate(idate1, end_year, end_month, end_day, end_hour)
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
      call ESMF_TimeIntervalSet (models(Iatmos)%dt,                     &
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
                                  timeStep=models(Iatmos)%dt,           &
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
      integer :: tile, i, j
      integer :: localPet, petCount, comm
      integer, allocatable :: deBlockList(:,:,:)
      integer, allocatable :: TLWidth(:,:), TUWidth(:,:)
      integer, allocatable :: CLW_c(:,:), CUW_c(:,:)
      integer, allocatable :: CLW_d(:,:), CUW_d(:,:)
      integer, dimension(2) :: CLW, CUW, TLW, TUW
      integer, dimension(2) :: deCount, minCorner, maxCorner
      TYPE (ESMF_ARRAY) :: GrdArray
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
!     Set ESMF Layout and distribution objects
!-----------------------------------------------------------------------
!
      deCount = (/ 1, nproc /)       
      models(Iatmos)%deLayout = ESMF_DELayoutCreate (models(Iatmos)%vm, &
                                         deCountList=deCount,           &
                                         petList=models(Iatmos)%petList,&
                                         rc=rc)
!
!-----------------------------------------------------------------------
!     Validate and print DELayout
!-----------------------------------------------------------------------
!
!      call ESMF_DELayoutValidate(models(Iatmos)%deLayout, rc=status)
!      call ESMF_DELayoutPrint(models(Iatmos)%deLayout, rc=status)
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
      if (localPet == 0) then
        print 10, 'Istr = ',(deBlockList(1,1,tile),tile=1,nproc)
        print 10, 'Iend = ',(deBlockList(1,2,tile),tile=1,nproc)
        print 10, 'Jstr = ',(deBlockList(2,1,tile),tile=1,nproc)
        print 10, 'Jend = ',(deBlockList(2,2,tile),tile=1,nproc)
      end if
10    format(1x,a,64i5)
!
!     Cordinates of the lower and upper corner of the patch
!
      minCorner = (/ 1, 1 /)
      maxCorner = (/ iym2, jx /)
!
      models(Iatmos)%distGrid = ESMF_DistGridCreate (minCorner,         &
                                     maxCorner,                         &
                                     deBlockList=deBlockList,           &
                                     deLayout=models(Iatmos)%deLayout,  &
                                     vm=models(Iatmos)%vm,              &
                                     rc=rc)
!
!-----------------------------------------------------------------------
!     Validate and print DistGrid
!-----------------------------------------------------------------------
!
!      call ESMF_DistGridValidate(models(Iatmos)%distGrid, rc=rc)
!      call ESMF_DistGridPrint(models(Iatmos)%distGrid, rc=rc)
!
!-----------------------------------------------------------------------
!     Set array descriptor
!-----------------------------------------------------------------------
!
      call ESMF_ArraySpecSet (models(Iatmos)%arrSpec,                   &
                              rank=2,                                   &
                              typekind=ESMF_TYPEKIND_R8,                &
                              rc=rc)
!
!-----------------------------------------------------------------------    
!     Create exchange arrays
!-----------------------------------------------------------------------    
!
      do i = 1, models(Iatmos)%nGrid
        grdArray = ESMF_ArrayCreate (arrayspec=models(Iatmos)%arrSpec,  &
                               distgrid=models(Iatmos)%distGrid,        &
                               indexflag=ESMF_INDEX_GLOBAL,             &
                               rc=rc)
        models(Iatmos)%grid(i)%lat%array = grdArray
!
        grdArray = ESMF_ArrayCreate (arrayspec=models(Iatmos)%arrSpec,  &
                               distgrid=models(Iatmos)%distGrid,        &
                               indexflag=ESMF_INDEX_GLOBAL,             &
                               rc=rc)
        models(Iatmos)%grid(i)%lon%array = grdArray
!
!       Get array pointer
!
        call ESMF_ArrayGet (models(Iatmos)%grid(i)%lat%array,           &
                        farrayPtr=models(Iatmos)%grid(i)%lat%data%field,&
                        rc=rc)          
        call ESMF_ArrayGet (models(Iatmos)%grid(i)%lon%array,           &
                        farrayPtr=models(Iatmos)%grid(i)%lon%data%field,&
                        rc=rc)          
!
!       Set adjustable settings of array object  
!
        call ESMF_ArraySet (array=models(Iatmos)%grid(i)%lat%array,     &
                        name=trim(models(Iatmos)%grid(i)%lat%data%name),&
                        rc=rc)
        call ESMF_ArraySet (array=models(Iatmos)%grid(i)%lon%array,     &
                        name=trim(models(Iatmos)%grid(i)%lon%data%name),&
                        rc=rc)
!  
!       Add array to export state
!   
        call ESMF_StateAdd (models(Iatmos)%stateExport,                 &
                           (/ models(Iatmos)%grid(i)%lat%array,         &
                              models(Iatmos)%grid(i)%lon%array /),      &
                           rc=rc)
!
!       Initialize the grid array
!
        models(Iatmos)%grid(i)%lat%data%field = 0.0d0
        models(Iatmos)%grid(i)%lon%data%field = 0.0d0
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
      integer :: i, j, jj
!     
!-----------------------------------------------------------------------
!     Load grid data for cross points 
!-----------------------------------------------------------------------
!
      do j = 1 , jxp
        do i = 1 , iym2
          jj = jxp*localPet+j
          models(Iatmos)%grid(1)%lat%data%field(i,jj) = mddom%xlon(i,j)
          models(Iatmos)%grid(1)%lon%data%field(i,jj) = mddom%xlon(i,j)
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
      end module mod_esmf_atm
