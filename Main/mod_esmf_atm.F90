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
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use ESMF
!
      use mod_regcm_interface
      use mod_atm_interface
      use mod_couplerr
      use mod_runparams
      use mod_dynparam
!
      implicit none
      private
!
!-----------------------------------------------------------------------
!     Public subroutines 
!-----------------------------------------------------------------------
!
      public  :: RCM_SetServices
      public  :: RCM_SetRun
      public  :: RCM_SetFinalize
!
      contains
!
      subroutine RCM_SetServices(comp, rc)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp), intent(inout) :: comp
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Register "initialize" routine     
!-----------------------------------------------------------------------
! 
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_INITIALIZE,&
                                      userRoutine=RCM_SetInitialize,    &
                                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Register "run" routine    
!-----------------------------------------------------------------------
! 
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_RUN,       &
                                      userRoutine=RCM_SetRun,           &
                                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Register "finalize" routine    
!-----------------------------------------------------------------------
! 
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_FINALIZE,  &
                                      userRoutine=RCM_SetFinalize,      &
                                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set return flag to success 
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetServices
!
      subroutine RCM_SetInitialize(comp, importState, exportState,      &
                                   clock, rc)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_constants, only : d_zero
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp), intent(inout) :: comp
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: localPet, petCount, ierr
!
      type(ESMF_TimeInterval) :: dtrun     
!
!-----------------------------------------------------------------------
!     Call RCM initialization routines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
! 
      call ESMF_GridCompGet(comp, vm=models(Iatmos)%vm, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_VMGet(models(Iatmos)%vm,                                &
                      localPet=localPet,                                &
                      petCount=petCount,                                &
                      mpiCommunicator=models(Iatmos)%comm,              &
                      rc=rc)  
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Initialize the gridded component 
!-----------------------------------------------------------------------
!
      call RCM_initialize(mpiCommunicator=models(Iatmos)%comm)
!
!-----------------------------------------------------------------------
!     Set-up ESMF internal clock for gridded component
!-----------------------------------------------------------------------
!
      call RCM_SetClock(clock)
!
!-----------------------------------------------------------------------
!     Run atmospheric model for only one time step (dt) to get initial
!     output from BATS to feed ocean model
!-----------------------------------------------------------------------
!
      call RCM_run(d_zero, dtsec)
!
!-----------------------------------------------------------------------
!     Update model clock 
!-----------------------------------------------------------------------
!
      call ESMF_TimeIntervalSet(dtrun, s_r8=dtsec, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_ClockAdvance(models(Iatmos)%clock,                      &
                             timeStep=dtrun,                            &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set-up grid and load coordinate data 
!-----------------------------------------------------------------------
!
      call RCM_SetGridArrays()
!
!-----------------------------------------------------------------------
!     Set-up import/export states and load initial data
!-----------------------------------------------------------------------
!
      call RCM_SetStates(localPet)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetInitialize
!
      subroutine RCM_SetRun(comp, importState, exportState, clock, rc)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp), intent(inout) :: comp
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: localPet, petCount, comm
      real*8 :: timestr1, timestr2, timeend, timepass
!
      type(ESMF_Time) :: currTime
      type(ESMF_TimeInterval) :: dt1, dt2, dt3, dt4
!
!-----------------------------------------------------------------------
!     Call RCM run routines
!-----------------------------------------------------------------------
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
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get RCM internal clock current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (models(Iatmos)%clock,                         &
                          currTime=models(Iatmos)%curTime,              &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_TimeGet (models(Iatmos)%curTime,                        &
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
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: write current time
!-----------------------------------------------------------------------
!
      if ((cpldbglevel > 0) .and. (localPet == 0)) then
      write(*,20)localPet,'Current Time',trim(models(Iatmos)%time%stamp)
      end if
!
!-----------------------------------------------------------------------
!     Get coupler current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (cplClock,                                     &
                          currTime=currTime,                            &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set RCM start time (timestr1)
!-----------------------------------------------------------------------
!
      dt1 = models(Iatmos)%curTime-models(Iatmos)%strTime
      call ESMF_TimeIntervalGet (dt1,                                   &
                                 s_r8=timestr1,                         &
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get coupler component start time (timestr2)
!-----------------------------------------------------------------------
!
      dt2 = currTime-cplStartTime
      call ESMF_TimeIntervalGet (dt2,                                   &
                                 s_r8=timestr2,                         &
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set RCM end time (timeend) and simulation lenght (timepass) 
!-----------------------------------------------------------------------
!
      if (timestr1 > timestr2) then
         dt3 = cplTimeStep
      else
         dt3 = dt1+cplTimeStep 
      end if
!
      call ESMF_TimeIntervalGet (dt3,                                   &
                                 s_r8=timeend,                          &
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      timepass = timeend-timestr1
!
!-----------------------------------------------------------------------
!     Get import data
!-----------------------------------------------------------------------
!
      if (dt3 .gt. cplTimeStep) then 
        call RCM_GetImportData()
      end if
!
!-----------------------------------------------------------------------
!     Run RCM
!-----------------------------------------------------------------------
!
      if ((cpl_dbglevel > 0) .and. (localPet == 0)) then
        write(*, fmt="(A28,4F15.2)") '[debug] -- Run ATM component',    &
              timestr1, timeend, timepass, timestr2
      end if
!
      call RCM_run(timestr1, timeend)
!
!-----------------------------------------------------------------------
!     Update model clock 
!-----------------------------------------------------------------------
!
      call ESMF_TimeIntervalSet(dt4, s_r8=timepass, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_ClockAdvance(models(Iatmos)%clock, timeStep=dt4, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Sync the PETs 
!-----------------------------------------------------------------------
!
      call ESMF_VMBarrier(models(Iatmos)%vm, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Put export data
!-----------------------------------------------------------------------
!
      call RCM_PutExportData(localPet)
!
!-----------------------------------------------------------------------
!     Formats 
!-----------------------------------------------------------------------
!
 20   format(' PET (', I3, ') - ATM Model ', A, ' = ', A)
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
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp), intent(inout) :: comp
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      type(ESMF_Clock), intent(inout) :: clock
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
! 
      integer :: i, n
!
!-----------------------------------------------------------------------
!     Call RCM finalize routines
!-----------------------------------------------------------------------
!
      call RCM_finalize()
!
!-----------------------------------------------------------------------
!     Call ESMF finalize routines
!-----------------------------------------------------------------------
!
      call ESMF_ClockDestroy(clock, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_StateDestroy(importState, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_StateDestroy(exportState, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      do n = 1, nNest(Iatmos) 
      do i = 1, ubound(models(Iatmos)%dataImport(:,n), dim=1)
      call ESMF_FieldDestroy(models(Iatmos)%dataImport(i,n)%field, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do
!
      do i = 1, ubound(models(Iatmos)%dataExport(:,n), dim=1)
      call ESMF_FieldDestroy(models(Iatmos)%dataExport(i,n)%field, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do
!
      call ESMF_GridDestroy(models(Iatmos)%grid(n), rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do
!
      call ESMF_StateDestroy(importState, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_StateDestroy(exportState, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompDestroy(comp, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetFinalize
!    
      subroutine RCM_SetClock(clock)
!
!-----------------------------------------------------------------------
!     Imported modules 
!-----------------------------------------------------------------------
!
      use mod_date 
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Clock), intent(inout) :: clock
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: ref_year,   str_year,   end_year
      integer :: ref_month,  str_month,  end_month
      integer :: ref_day,    str_day,    end_day
      integer :: ref_hour,   str_hour,   end_hour
      integer :: ref_minute, str_minute, end_minute
      integer :: ref_second, str_second, end_second
      character (len=80) :: name
      integer :: rc
!
!-----------------------------------------------------------------------
!     Create gridded component clock 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Create ESMF calendar
!-----------------------------------------------------------------------
!
      name = 'Mixed Gregorian/Julian calendar'
      models(Iatmos)%cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,  &
                                        name=trim(name),                &
                                        rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
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
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
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
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
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
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set time interval
!-----------------------------------------------------------------------
!
      call ESMF_TimeIntervalSet (models(Iatmos)%dtsec,                  &
                                 s_r8=dtsec,                            &
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create time clock.
!-----------------------------------------------------------------------
!
      name = 'Model clock (Atmosphere)'
      models(Iatmos)%clock = ESMF_ClockCreate (name=trim(name),         &
                                  refTime=models(Iatmos)%refTime,       &
                                  timeStep=models(Iatmos)%dtsec,        &
                                  startTime=models(Iatmos)%strTime,     &
                                  stopTime=models(Iatmos)%endTime,      &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Copy clock
!-----------------------------------------------------------------------
!
      clock = ESMF_ClockCreate (models(Iatmos)%clock, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Validate time clock
!-----------------------------------------------------------------------
!
      call ESMF_ClockValidate (models(Iatmos)%clock, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get gridded component internal clock current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (models(Iatmos)%clock,                         &
                          currTime=models(Iatmos)%curTime,              &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
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
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Attach time information to export state as attribute
!-----------------------------------------------------------------------
!
      name = 'start time'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             valueList=(/ str_year  , str_month,        &
                                          str_day   , str_hour ,        &
                                          str_minute, str_second /),    &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'stop time'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             valueList=(/ end_year  , end_month,        &
                                          end_day   , end_hour ,        &
                                          end_minute, end_second /),    &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'time step'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=int(dtsec),                          &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set return flag to success
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetClock
!
      subroutine RCM_SetGridArrays ()
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, n, rc
      integer :: localPet, petCount, comm, localDECount
      character (len=40) :: name
!
      type(ESMF_Field) :: grdField
      type(ESMF_StaggerLoc) :: staggerLoc
      real(ESMF_KIND_R8), pointer :: ptrX(:,:), ptrY(:,:)
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
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set RCM domain decomposition variables
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Create ESMF DistGrid based on model domain decomposition
!-----------------------------------------------------------------------
!
      models(Iatmos)%distGrid(n) = ESMF_DistGridCreate (                &
                                        minIndex=(/ 1, 1 /),            &
                                        maxIndex=(/ jxm1, iym1 /),      &
                                        regDecomp=(/nproc,1/),          &
                                        rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: print DistGrid
!-----------------------------------------------------------------------
!
      if (cpldbglevel > 1) then
      call ESMF_DistGridValidate(models(Iatmos)%distGrid(n), rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_DistGridPrint(models(Iatmos)%distGrid(n), rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
!-----------------------------------------------------------------------
!     Set array descriptor
!-----------------------------------------------------------------------
!
      call ESMF_ArraySpecSet (models(Iatmos)%arrSpec(n),                &
                              typekind=ESMF_TYPEKIND_R8,                &
                              rank=2,                                   &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      do i = 1, ubound(models(Iatmos)%mesh, dim=1) 
!
!-----------------------------------------------------------------------
!     Set staggering type 
!-----------------------------------------------------------------------
!
      if (models(Iatmos)%mesh(i,n)%gtype == Icross) then
        staggerLoc = ESMF_STAGGERLOC_CENTER
      else if (models(Iatmos)%mesh(i,n)%gtype == Idot) then
        staggerLoc = ESMF_STAGGERLOC_CORNER
      end if
!
!-----------------------------------------------------------------------
!     Create ESMF Grid
!-----------------------------------------------------------------------
!
      if (i == 1) then
      models(Iatmos)%grid(n) = ESMF_GridCreate (                        &
                                    distgrid=models(Iatmos)%distGrid(n),&
                                    indexflag=ESMF_INDEX_GLOBAL,        &
                                    name="atm_grid",                    &
                                    rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
!-----------------------------------------------------------------------
!     Allocate coordinates 
!-----------------------------------------------------------------------
!
      call ESMF_GridAddCoord (models(Iatmos)%grid(n),                   &
                              staggerLoc=staggerLoc,                    &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iatmos)%grid(n),                        &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointers and set coordinates for the grid 
!-----------------------------------------------------------------------
! 
      do j = 0, localDECount-1
        call ESMF_GridGetCoord (models(Iatmos)%grid(n),                 &
                                localDE=j,                              &
                                staggerLoc=staggerLoc,                  &
                                coordDim=1,                             &
                                farrayPtr=ptrX,                         &
                                rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
!
        call ESMF_GridGetCoord (models(Iatmos)%grid(n),                 &
                                localDE=j,                              &
                                staggerLoc=staggerLoc,                  &
                                coordDim=2,                             &
                                farrayPtr=ptrY,                         &
                                rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)       
        end if
!
!-----------------------------------------------------------------------
!     Debug: write size of pointers    
!-----------------------------------------------------------------------
!
        write(*,30) localPet, j, "PTR-"//                               &
                    trim(GRIDDES(models(Iatmos)%mesh(i,n)%gtype)),      &
                    lbound(ptrX, dim=1), ubound(ptrX, dim=1),           &
                    lbound(ptrX, dim=2), ubound(ptrX, dim=2)
!
!-----------------------------------------------------------------------
!     Fill the pointers    
!-----------------------------------------------------------------------
!
        if (models(Iatmos)%mesh(i,n)%gtype == Idot) then
          write(*,30) localPet, j, "ATM-"//                             &
                   trim(GRIDDES(models(Iatmos)%mesh(i,n)%gtype)),       &
                   lbound(mddom%dlon, dim=1), ubound(mddom%dlon, dim=1),&
                   lbound(mddom%dlon, dim=2), ubound(mddom%dlon, dim=2)
!
          ptrX = mddom%dlon
          ptrY = mddom%dlat
        else if (models(Iatmos)%mesh(i,n)%gtype == Icross) then
          if (localPet .eq. models(Iatmos)%nproc-1) then 
            write(*,30) localPet, j, "ATM-"//                           &
                        trim(GRIDDES(models(Iatmos)%mesh(i,n)%gtype)),  &
                        lbound(mddom%xlon(1:jxp-1,1:iym1), dim=1),      &
                        ubound(mddom%xlon(1:jxp-1,1:iym1), dim=1),      &
                        lbound(mddom%xlon(1:jxp-1,1:iym1), dim=2),      &
                        ubound(mddom%xlon(1:jxp-1,1:iym1), dim=2)
!
              ptrX = mddom%xlon(1:jxp-1,1:iym1)
              ptrY = mddom%xlat(1:jxp-1,1:iym1)
          else
            write(*,30) localPet, j, "ATM-"//                           &
                        trim(GRIDDES(models(Iatmos)%mesh(i,n)%gtype)),  &
                        lbound(mddom%xlon(:,1:iym1), dim=1),            &
                        ubound(mddom%xlon(:,1:iym1), dim=1),            &
                        lbound(mddom%xlon(:,1:iym1), dim=2),            &
                        ubound(mddom%xlon(:,1:iym1), dim=2)
!
            ptrX = mddom%xlon(:,1:iym1)
            ptrY = mddom%xlat(:,1:iym1)
          end if
        end if
!
!-----------------------------------------------------------------------
!     Nullify pointers 
!-----------------------------------------------------------------------
!
        if (associated(ptrY)) then
          nullify(ptrY)
        end if
        if (associated(ptrX)) then
          nullify(ptrX)
        end if
      end do
!
!-----------------------------------------------------------------------
!     Validate Grid 
!-----------------------------------------------------------------------
!
!      call ESMF_GridValidate(models(Iatmos)%grid(n), rc=rc)
!      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Write ESMF Grid in VTK format (debug) 
!-----------------------------------------------------------------------
!
      if (cpldbglevel .gt. 1) then
      print*, '[debug] -- write grid information to file '//            &
      '>atmos_'//trim(GRIDDES(models(Iatmos)%mesh(i,n)%gtype))//'point<'
      call ESMF_GridWriteVTK(models(Iatmos)%grid(n),                    &
                         filename="atmos_"//                            &
                         trim(GRIDDES(models(Iatmos)%mesh(i,n)%gtype))//&
                         "point",                                       &
                         staggerLoc=staggerLoc,                         &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
      end do
      end do
!
!-----------------------------------------------------------------------
!     Format definition 
!-----------------------------------------------------------------------
!
 30   format(" PET(",I3,") - DE(",I2,") - ", A10, " : ", 4I8)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetGridArrays
!
      subroutine RCM_SetStates (localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_bats_common
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      integer, intent(in) :: localPet
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, id, n, rc, localDECount
      logical :: flag
      character (len=40) :: name
      type(ESMF_StaggerLoc) :: staggerLoc
!
!-----------------------------------------------------------------------
!     Initialize the import and export fields 
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Create export state fields 
!-----------------------------------------------------------------------
!
      do i = 1, size(models(Iatmos)%dataExport(:,n), dim=1)
!
!-----------------------------------------------------------------------
!     Set staggering type 
!-----------------------------------------------------------------------
!
      if (models(Iatmos)%dataExport(i,n)%gtype == Icross) then
        staggerLoc = ESMF_STAGGERLOC_CENTER
        id = getMeshID(models(Iatmos)%mesh(:,n), Icross)
      else if (models(Iatmos)%dataExport(i,n)%gtype == Idot) then
        staggerLoc = ESMF_STAGGERLOC_CORNER
        id = getMeshID(models(Iatmos)%mesh(:,n), Idot)
      end if
!
!-----------------------------------------------------------------------
!     Create field 
!-----------------------------------------------------------------------
!
      name = models(Iatmos)%dataExport(i,n)%name
      models(Iatmos)%dataExport(i,n)%field = ESMF_FieldCreate (         &
                                  models(Iatmos)%grid(n),               &
                                  models(Iatmos)%arrSpec(n),            &
                                  indexflag=ESMF_INDEX_GLOBAL,          &
                                  staggerloc=staggerLoc,                &
                                  name=trim(name),                      &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Attach interpolation type attribute of field 
!-----------------------------------------------------------------------
!
      name = 'interpolation_type'
      call ESMF_AttributeSet(models(Iatmos)%dataExport(i,n)%field,      &
                             name=trim(name),                           &
                             value=models(Iatmos)%dataExport(i,n)%itype,&
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iatmos)%grid(n),                        &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Put data into state 
!-----------------------------------------------------------------------
! 
      do j = 0, localDECount-1
!
!-----------------------------------------------------------------------
!     Get pointer from ESMF Field 
!-----------------------------------------------------------------------
!
      call ESMF_FieldGet (models(Iatmos)%dataExport(i,n)%field,         &
                          localDe=j,                                    &
                          farrayPtr=models(Iatmos)%dataExport(i,n)%ptr, &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Initialize pointer 
!-----------------------------------------------------------------------
!
      models(Iatmos)%dataExport(i,n)%ptr = 0.0d0
!
      end do
!
!-----------------------------------------------------------------------
!     Add fields to export state
!-----------------------------------------------------------------------
!
      call ESMF_StateAdd (models(Iatmos)%stateExport,                   &
                         (/ models(Iatmos)%dataExport(i,n)%field /),    &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do
!     
!-----------------------------------------------------------------------
!     Create import state arrays
!-----------------------------------------------------------------------
!
      do i = 1, ubound(models(Iatmos)%dataImport(:,n), dim=1)
!
!-----------------------------------------------------------------------
!     Set staggering type 
!-----------------------------------------------------------------------
!
      if (models(Iatmos)%dataImport(i,n)%gtype == Icross) then
        staggerLoc = ESMF_STAGGERLOC_CENTER
        id = getMeshID(models(Iatmos)%mesh(:,n), Icross)        
      else if (models(Iatmos)%dataImport(i,n)%gtype == Idot) then
        staggerLoc = ESMF_STAGGERLOC_CORNER
        id = getMeshID(models(Iatmos)%mesh(:,n), Idot)        
      end if
!
!-----------------------------------------------------------------------
!     Create field 
!-----------------------------------------------------------------------
!
      name = models(Iatmos)%dataImport(i,n)%name
      models(Iatmos)%dataImport(i,n)%field = ESMF_FieldCreate (         &
                                  models(Iatmos)%grid(n),               &
                                  models(Iatmos)%arrSpec(n),            &
                                  indexflag=ESMF_INDEX_GLOBAL,          &
                                  staggerloc=staggerLoc,                &
                                  name=trim(name),                      &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iatmos)%grid(n),                        &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Put data into state 
!-----------------------------------------------------------------------
! 
      do j = 0, localDECount-1
!
!-----------------------------------------------------------------------
!     Get pointer from ESMF Field 
!-----------------------------------------------------------------------
!
      call ESMF_FieldGet (models(Iatmos)%dataImport(i,n)%field,         &
                          localDe=j,                                    &
                          farrayPtr=models(Iatmos)%dataImport(i,n)%ptr, &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Initialize pointer
!-----------------------------------------------------------------------
!      
      models(Iatmos)%dataImport(i,n)%ptr = MISSING_R8 
!
      end do
!
!-----------------------------------------------------------------------
!     Add fields to import state
!-----------------------------------------------------------------------
!
      call ESMF_StateAdd (models(Iatmos)%stateImport,                   &
                         (/ models(Iatmos)%dataImport(i,n)%field /),    &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do
      end do
!
!-----------------------------------------------------------------------
!     Attach coupled model parameters to export state as attribute 
!-----------------------------------------------------------------------
!
      name = 'coupler time step'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=int(cpldt),                          &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'exchange variable mode'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=cplexvars,                           &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'interpolation type'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=cplinterp,                           &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'boundary smoothing'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=cplbdysmooth,                        &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'debug level'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=cpldbglevel,                         &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine RCM_SetStates
!
      subroutine RCM_PutExportData (localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_bats_common
      use mod_dynparam, only : kz, ici1, ici2, jci1 , jci2
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      integer, intent(in) :: localPet
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, ii, jj, id, n, rc, localDECount
      integer, dimension(2) :: dims
      logical :: flag
      character (len=40) :: name
      character (len=100) :: outfile
      real(sp), allocatable, dimension(:,:) :: urot, vrot
!
      type(ESMF_StaggerLoc) :: staggerLoc
!
!-----------------------------------------------------------------------
!     Initialize the import and export fields 
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Rotate wind components (u and v) to earth coordinates 
!-----------------------------------------------------------------------
!
!      call calc_uvmet (u10m_o, v10m_o, urot, vrot, localPet) 
!
!-----------------------------------------------------------------------
!     Create export state fields 
!-----------------------------------------------------------------------
!
      do i = 1, size(models(Iatmos)%dataExport(:,n), dim=1)
      name = models(Iatmos)%dataExport(i,n)%name
!
!-----------------------------------------------------------------------
!     Debug: write size of pointers    
!-----------------------------------------------------------------------
!
      if (cpl_dbglevel > 1) then
        write(*,40) localPet, "PTR-"//trim(adjustl(name)),              &
                    lbound(models(Iatmos)%dataExport(i,n)%ptr, dim=1),  &
                    ubound(models(Iatmos)%dataExport(i,n)%ptr, dim=1),  &
                    lbound(models(Iatmos)%dataExport(i,n)%ptr, dim=2),  &
                    ubound(models(Iatmos)%dataExport(i,n)%ptr, dim=2)
      end if
!
!-----------------------------------------------------------------------
!     Fill pointers with data 
!-----------------------------------------------------------------------
!
      if (trim(adjustl(name)) == "Pair") then ! mb
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound((sfps+ptop)*d_10, dim=1),                  &
                      ubound((sfps+ptop)*d_10, dim=1),                  &
                      lbound((sfps+ptop)*d_10, dim=2),                  &
                      ubound((sfps+ptop)*d_10, dim=2)
        end if
        models(Iatmos)%dataExport(i,n)%ptr = (sfps+ptop)*d_10
      else if (trim(adjustl(name)) == "Tair") then ! Kelvin 
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound(t2m, dim=2), ubound(t2m, dim=2), &
                      lbound(t2m, dim=3), ubound(t2m, dim=3)
        end if
        if (localPet .eq. 0) then
          models(Iatmos)%dataExport(i,n)%ptr(2:,2:iym2) = t2m(1,:,:) 
        else if (localPet .eq. models(Iatmos)%nproc-1) then
          models(Iatmos)%dataExport(i,n)%ptr(:jxm2,2:iym2) = t2m(1,:,:)
        else
          models(Iatmos)%dataExport(i,n)%ptr(:,2:iym2) = t2m(1,:,:)
        end if
      else if (trim(adjustl(name)) == "Qair") then ! kg/kg
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound(q2m, dim=2), ubound(q2m, dim=2),           &
                      lbound(q2m, dim=3), ubound(q2m, dim=3)
        end if
        if (localPet .eq. 0) then
          models(Iatmos)%dataExport(i,n)%ptr(2:,2:iym2) = q2m(1,:,:)/   &
                                                   (d_one+q2m(1,:,:))
        else if (localPet .eq. models(Iatmos)%nproc-1) then
          models(Iatmos)%dataExport(i,n)%ptr(:jxm2,2:iym2) = q2m(1,:,:)/&
                                                   (d_one+q2m(1,:,:)) 
        else
          models(Iatmos)%dataExport(i,n)%ptr(:,2:iym2) = q2m(1,:,:)/    &
                                                   (d_one+q2m(1,:,:)) 
        end if
      else if (trim(adjustl(name)) == "swrad") then ! W/m2
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound(fsw, dim=1), ubound(fsw, dim=1),           &
                      lbound(fsw, dim=2), ubound(fsw, dim=2)
        end if
        models(Iatmos)%dataExport(i,n)%ptr = fsw
      else if (trim(adjustl(name)) == "lwrad_down") then ! W/m2
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound(flwd, dim=1), ubound(flwd, dim=1),         &
                      lbound(flwd, dim=2), ubound(flwd, dim=2)
        end if
        models(Iatmos)%dataExport(i,n)%ptr = flwd
      else if (trim(adjustl(name)) == "lwrad") then ! W/m2
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound(flw, dim=1), ubound(flw, dim=1),           &
                      lbound(flw, dim=2), ubound(flw, dim=2)
        end if
        models(Iatmos)%dataExport(i,n)%ptr = flw
      else if (trim(adjustl(name)) == "rain") then ! mm/s
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound(totpr, dim=1), ubound(totpr, dim=1),       &
                      lbound(totpr, dim=2), ubound(totpr, dim=2)
        end if
        if (localPet .eq. 0) then
          models(Iatmos)%dataExport(i,n)%ptr(2:,2:iym2) = totpr 
        else if (localPet .eq. models(Iatmos)%nproc-1) then
          models(Iatmos)%dataExport(i,n)%ptr(:jxm1,2:iym2) = totpr 
        else
          models(Iatmos)%dataExport(i,n)%ptr(:,2:iym2) = totpr
        end if
      else if (trim(adjustl(name)) == "Uwind") then ! m/s
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound(u10m, dim=2), ubound(u10m, dim=2),         &
                      lbound(u10m, dim=3), ubound(u10m, dim=3)
        end if
        if (localPet .eq. 0) then
          models(Iatmos)%dataExport(i,n)%ptr(2:,2:iym2) = u10m(1,:,:)
        else if (localPet .eq. models(Iatmos)%nproc-1) then
          models(Iatmos)%dataExport(i,n)%ptr(:jxm2,2:iym2) = u10m(1,:,:)
        else
          models(Iatmos)%dataExport(i,n)%ptr(:,2:iym2) = u10m(1,:,:)
        end if
      else if (trim(adjustl(name)) == "Vwind") then ! m/s
        if (cpl_dbglevel > 1) then
          write(*,40) localPet, "ATM-"//trim(adjustl(name)),            &
                      lbound(v10m, dim=2), ubound(v10m, dim=2),         &
                      lbound(v10m, dim=3), ubound(v10m, dim=3)
        end if
        if (localPet .eq. 0) then
          models(Iatmos)%dataExport(i,n)%ptr(2:,2:iym2) = v10m(1,:,:)
        else if (localPet .eq. models(Iatmos)%nproc-1) then
          models(Iatmos)%dataExport(i,n)%ptr(:jxm2,2:iym2) = v10m(1,:,:)
        else
          models(Iatmos)%dataExport(i,n)%ptr(:,2:iym2) = v10m(1,:,:)
        end if
      end if
!
!-----------------------------------------------------------------------
!     Debug: write field to ASCII file    
!-----------------------------------------------------------------------
!
      if (cpldbglevel > 3) then
      write(outfile,                                                    &
            fmt='(A10,"_",A3,"_",I4,"-",I2.2,"-",                       &
            I2.2,"_",I2.2,"_",I2.2,".txt")')                            &
            'atm_export',                                               &
            trim(adjustl(name)),                                        &
            models(Iatmos)%time%year,                                   &
            models(Iatmos)%time%month,                                  &
            models(Iatmos)%time%day,                                    &
            models(Iatmos)%time%hour,                                   &
            localPet
!
      open (unit=99, file = trim(outfile)) 
      call print_matrix_r8(models(Iatmos)%dataExport(i,n)%ptr,          &
                           1, 1, localPet, 99, "ATM_PTR")
      close(99)
      end if
!
!-----------------------------------------------------------------------
!     Debug: write field to file
!-----------------------------------------------------------------------
!
      if (cpldbglevel > 2) then
      write(outfile,                                                    &
            fmt='(A10,"_",A3,"_",I4,"-",I2.2,"-",I2.2,"_",I2.2,".nc")') &
            'atm_export',                                               &
            trim(adjustl(name)),                                        &
            models(Iatmos)%time%year,                                   &
            models(Iatmos)%time%month,                                  &
            models(Iatmos)%time%day,                                    &
            models(Iatmos)%time%hour
!
      call ESMF_FieldWrite(models(Iatmos)%dataExport(i,n)%field,        &
                           trim(adjustl(outfile)),                      &
                           rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
      end do
      end do
!
!-----------------------------------------------------------------------
!     Formats 
!-----------------------------------------------------------------------
!
 40   format(" PET(",I3,") - ", A14, " : ", 4I8) 
!
      end subroutine RCM_PutExportData
!
      subroutine RCM_GetImportData
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
#ifdef ROMSICE
      use mod_bats_romsocn, only : sst2d, hice2d
#else
      use mod_bats_romsocn, only : sst2d
#endif
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, id, n, rc
      integer :: localPet, petCount, comm, localDECount
      character (len=40) :: name
      character (len=100) :: outfile
      real*8 :: scale_factor, add_offset
!
      real(ESMF_KIND_R8), pointer :: ptr(:,:)
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
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Loop over number of nested/composed meshs 
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iatmos)
!
!-----------------------------------------------------------------------
!     Get import fields 
!-----------------------------------------------------------------------
!
      do i = 1, size(models(Iatmos)%dataImport(:,n), dim=1)
      name = models(Iatmos)%dataImport(i,n)%name
!
!-----------------------------------------------------------------------
!     Set staggering type 
!-----------------------------------------------------------------------
!
      if (models(Iatmos)%dataImport(i,n)%gtype == Icross) then
        id = getMeshID(models(Iatmos)%mesh(:,n), Icross)
      else if (models(Iatmos)%dataImport(i,n)%gtype == Idot) then
        id = getMeshID(models(Iatmos)%mesh(:,n), Idot)
      end if
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iatmos)%grid(n),                        &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointer
!-----------------------------------------------------------------------
! 
      do j = 0, localDECount-1
      call ESMF_FieldGet (models(Iatmos)%dataImport(i,n)%field,         &
                          localDE=j,                                    &
                          farrayPtr=ptr,                                &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Put data to RCM variable
!-----------------------------------------------------------------------
!
      scale_factor = models(Iatmos)%dataImport(i,n)%scale_factor
      add_offset = models(Iatmos)%dataImport(i,n)%add_offset
!
      select case (trim(adjustl(name)))
      case('SST')      
          sst2d = (ptr*scale_factor)+add_offset
#ifdef ROMSICE
      case('Hice')      
          hice2d = (ptr*scale_factor)+add_offset
#endif
      end select
      end do
!
!-----------------------------------------------------------------------
!     Debug: write field to ASCII file    
!-----------------------------------------------------------------------
!
      if (cpldbglevel > 3) then
      write(outfile,                                                    &
            fmt='(A10,"_",A3,"_",I4,"-",I2.2,"-",                       &
            I2.2,"_",I2.2,"_",I2.2,".txt")')                            &
            'atm_import',                                               &
            trim(adjustl(name)),                                        &
            models(Iatmos)%time%year,                                   &
            models(Iatmos)%time%month,                                  &
            models(Iatmos)%time%day,                                    &
            models(Iatmos)%time%hour,                                   &
            localPet
!
      open (unit=99, file = trim(outfile))
      call print_matrix_r8(ptr, 1, 1, localPet, 99, "ATM_PTR")
      close(99)
      end if
!
!-----------------------------------------------------------------------
!     Debug: write field to NetCDF file
!-----------------------------------------------------------------------
!
      if (cpldbglevel > 2) then
      write(outfile,                                                    &
            fmt='(A10,"_",A3,"_",I4,"-",I2.2,"-",I2.2,"_",I2.2,".nc")') &
            'atm_import',                                               &
            trim(adjustl(name)),                                        &
            models(Iatmos)%time%year,                                   &
            models(Iatmos)%time%month,                                  &
            models(Iatmos)%time%day,                                    &
            models(Iatmos)%time%hour
!
      call ESMF_FieldWrite(models(Iatmos)%dataImport(i,n)%field,        &
                           trim(adjustl(outfile)),                      &
                           rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
!-----------------------------------------------------------------------
!     Nullify pointer to make sure that it does not point on a random 
!     part in the memory 
!-----------------------------------------------------------------------
!
      if (associated(ptr)) then
        nullify(ptr)
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
      end subroutine RCM_GetImportData
!
      end module mod_esmf_atm
