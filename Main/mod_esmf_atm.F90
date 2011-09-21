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
      use mod_regcm_interface
      use mod_atm_interface
      use mod_couplerr
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
      use mod_runparams, only : dtsec
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
      logical :: first
      integer :: localPet, petCount, comm, ierr
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
                      mpiCommunicator=comm,                             &
                      rc=rc)  
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
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
      call RCM_SetClock(clock)
!
!-----------------------------------------------------------------------
!     Run atmospheric model for only one time step (dt) to get initial
!     output from BATS to feed ocean model
!-----------------------------------------------------------------------
!
      first = .true.
      call RCM_run(d_zero, dtsec, first)
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
      logical :: first
      integer :: localPet, petCount, comm
      real*8 :: timestr, timeend
!
      type(ESMF_Time) :: currTime
      type(ESMF_Time) :: time1, time2
      type(ESMF_TimeInterval) :: dt1
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
!     Write current time (debug)
!-----------------------------------------------------------------------
!
      write(*, 20) localPet, 'Current Time',                            &
                   trim(models(Iatmos)%time%stamp)
!
!-----------------------------------------------------------------------
!     Get coupler current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (clock,                                        &
                          currTime=currTime,                            &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set start and end time 
!-----------------------------------------------------------------------
!
      dt1 = models(Iatmos)%curTime-models(Iatmos)%strTime
      call ESMF_TimeIntervalGet (dt1,                                   &
                                 s_r8=timestr,                          &
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      if (models(Iatmos)%curTime > currTime) then
         dt1 = cplTimeStep-dt1
      else
         dt1 = dt1+cplTimeStep 
      end if
!
      call ESMF_TimeIntervalGet (dt1,                                   &
                                 s_r8=timeend,                          &
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      print*, "** turuncu ** ", timestr, timeend
!
!      first = .false.
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
      subroutine RCM_SetClock(clock)
!
!-----------------------------------------------------------------------
!     Imported modules 
!-----------------------------------------------------------------------
!
      use mod_runparams
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
!     Attach coupled model parameters to export state as attribute 
!-----------------------------------------------------------------------
!
      name = 'VTK on/off'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=vtk_on,                              &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'debug level'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=debug_level,                         &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'coupler time step'
      call ESMF_AttributeSet(models(Iatmos)%stateExport,                &
                             name=trim(name),                           &
                             value=int(dtcpl),                          &
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
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_runparams, only : vtk_on
      use mod_dynparam, only : nproc, iy, iym2, jx, jxp, jendx,         &
                               jendl, debug_level
      use mod_atm_interface, only : mddom
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
                                        maxIndex=(/ iy, jx /),          &
                                        regDecomp=(/1,nproc/),          &
                                        rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: validate and print DistGrid
!-----------------------------------------------------------------------
!
      if ((localPet == 0) .and. (debug_level > 2)) then
        call ESMF_DistGridValidate(models(Iatmos)%distGrid(n), rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
        call ESMF_DistGridPrint(models(Iatmos)%distGrid(n), rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
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
      models(Iatmos)%mesh(i,n)%grid = ESMF_GridCreate (                 &
                                    distgrid=models(Iatmos)%distGrid(n),&
                                    indexflag=ESMF_INDEX_GLOBAL,        &
                                    name="atm_grid",                    &
                                    rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Allocate coordinates 
!-----------------------------------------------------------------------
!
      call ESMF_GridAddCoord (models(Iatmos)%mesh(i,n)%grid,            &
                              staggerLoc=staggerLoc,                    &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iatmos)%mesh(i,n)%grid,                 &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointers and set coordinates for the grid 
!-----------------------------------------------------------------------
! 
      do j = 0, localDECount-1
        call ESMF_GridGetCoord (models(Iatmos)%mesh(i,n)%grid,          &
                                localDE=j,                              &
                                staggerLoc=staggerLoc,                  &
                                coordDim=1,                             &
                                farrayPtr=ptrX,                         &
                                rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
!
        call ESMF_GridGetCoord (models(Iatmos)%mesh(i,n)%grid,          &
                                localDE=j,                              &
                                staggerLoc=staggerLoc,                  &
                                coordDim=2,                             &
                                farrayPtr=ptrY,                         &
                                rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)       
        end if
!
        if (models(Iatmos)%mesh(i,n)%gtype == Icross) then
          ptrX = mddom%xlon
          ptrY = mddom%xlat
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
!     Write ESMF Grid in VTK format (debug) 
!-----------------------------------------------------------------------
!
      if (vtk_on) then
        call ESMF_GridWriteVTK(models(Iatmos)%mesh(i,n)%grid,           &
                               filename="atmos_src",                    &
                               rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
      end if
!
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
      subroutine RCM_SetStates (localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_runparams, only : debug_level
      use mod_dynparam, only : nproc, iy, iym1, iym2, jx, jxp, ptop
      use mod_bats_common
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
      real(ESMF_KIND_R8), pointer :: ptr(:,:)
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
                                  models(Iatmos)%mesh(id,n)%grid,       &
                                  models(Iatmos)%arrSpec(n),            &
                                  staggerloc=staggerLoc,                &
                                  name=trim(name),                      &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iatmos)%mesh(id,n)%grid,                &
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
                          farrayPtr=ptr,                                &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Initialize pointer 
!-----------------------------------------------------------------------
!
      ptr = 0.0d0
!
!-----------------------------------------------------------------------
!     Put data     
!-----------------------------------------------------------------------
!      
      if (trim(adjustl(name)) == "Pair") then
!       unit: millibar       
        ptr = (sfps(:,1:jxp)+ptop)*d_10
      else if (trim(adjustl(name)) == "Tair") then 
!       unit: celsius
!        ptr = sts2%tg(:,1:jxp)-tzero
        ptr(2:iym1,:) = transpose(t2m_o)
        ptr(1,:) = ptr(2,:)
        ptr(iy,:) = ptr(iym1,:)
      else if (trim(adjustl(name)) == "Qair") then
!       unit: kg/kg
        ptr(2:iym1,:) = transpose(q2m_o)
        ptr(1,:) = ptr(2,:)
        ptr(iy,:) = ptr(iym1,:)
      else if (trim(adjustl(name)) == "swrad") then
        ! one less
        ! ptr: 1-88,1-16,17-32,33-48,49-64 - fsw2d:1-87,1-16
        ptr(1:iym1,:) = fsw2d
        ptr(iym1,:) = fsw2d(iym2,:) 
      else if (trim(adjustl(name)) == "lwrad_down") then
        ! one less
        ! ptr: 1-88,1-16,17-32,33-48,49-64 - flw2d:1-87,1-16
        ptr(1:iym1,:) = flw2d
        ptr(iym1,:) = flw2d(iym2,:)
      else if (trim(adjustl(name)) == "rain") then
        ! one less
        ! ptr: 1-88,1-16,17-32,33-48,49-64 - pptc:1-87,1-16
        ptr(1:iym1,:) = pptnc+pptc 
        ptr(iym1,:) = pptnc(iym2,:)+pptc(iym2,:)
      else if (trim(adjustl(name)) == "Uwind") then
!       unit: meter/second
        ptr(2:iym1,:) = transpose(u10m_o)
        ptr(1,:) = ptr(2,:)
        ptr(iy,:) = ptr(iym1,:)
      else if (trim(adjustl(name)) == "Vwind") then
!       unit: meter/second
        ptr(2:iym1,:) = transpose(v10m_o)
      end if
      end do
!
!-----------------------------------------------------------------------
!     Write field to file (debug)     
!-----------------------------------------------------------------------
!  
      call ESMF_FieldWrite(models(Iatmos)%dataExport(i,n)%field,        &
                           'atm_export_'//trim(adjustl(name))//'.nc',   &
                           rc=rc)
!
!-----------------------------------------------------------------------
!     Add fields to export state
!-----------------------------------------------------------------------
!
      call ESMF_StateAdd (models(Iatmos)%stateExport,                   &
                         (/ models(Iatmos)%dataExport(i,n)%field /),    &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
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
      else if (models(Iatmos)%dataImport(i,n)%gtype == Idot) then
        staggerLoc = ESMF_STAGGERLOC_CORNER
      end if
!
!-----------------------------------------------------------------------
!     Create field 
!-----------------------------------------------------------------------
!
      name = models(Iatmos)%dataImport(i,n)%name
      models(Iatmos)%dataImport(i,n)%field = ESMF_FieldCreate (         &
                                  models(Iatmos)%mesh(i,n)%grid,        &
                                  models(Iatmos)%arrSpec(n),            &
                                  staggerloc=staggerLoc,                &
                                  name=trim(name),                      &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iatmos)%mesh(i,n)%grid,                 &
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
                          farrayPtr=ptr,                                &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Put data     
!-----------------------------------------------------------------------
!      
      if (trim(adjustl(name)) == "SST") then
        ptr = 0.0
      end if
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
      end subroutine RCM_SetStates
!
      end module mod_esmf_atm
