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
      module mod_esmf_ocn
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use ESMF
      use mod_regcm_interface
      use mod_couplerr
!
!-----------------------------------------------------------------------
!     Used module declarations (ROMS Component routines)
!-----------------------------------------------------------------------
!
      use ocean_control_mod, only : ROMS_initialize
      use ocean_control_mod, only : ROMS_run
      use ocean_control_mod, only : ROMS_finalize
! 
      implicit none
      private
!
!-----------------------------------------------------------------------
!     Public subroutines 
!-----------------------------------------------------------------------
!
      public  :: ROMS_SetServices
      public  :: ROMS_SetRun
      public  :: ROMS_SetFinalize
!
      contains
!
      subroutine ROMS_SetServices(comp, rc)
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
                                      userRoutine=ROMS_SetInitialize,   &
                                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Register "run" routine    
!-----------------------------------------------------------------------
!
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_RUN,       &
                                      userRoutine=ROMS_SetRun,          &
                                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Register "finalize" routine    
!-----------------------------------------------------------------------
! 
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_FINALIZE,  &
                                      userRoutine=ROMS_SetFinalize,     &
                                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set return flag to success 
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetServices
!
      subroutine ROMS_SetInitialize(comp, importState, exportState,     &
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
      logical :: flag
      integer :: localPet, petCount, comm, ierr 
!
!-----------------------------------------------------------------------
!     Call ROMS initialization routines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!     
      call ESMF_GridCompGet(comp, vm=models(Iocean)%vm, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_VMGet(models(Iocean)%vm,                                &
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
      flag = .TRUE.
      call MPI_Comm_dup(comm, models(Iocean)%comm, ierr)
      call ROMS_initialize(flag, MyCOMM=models(Iocean)%comm)
!
!-----------------------------------------------------------------------
!     Set-up ESMF internal clock for grided component 
!-----------------------------------------------------------------------
!
      call ROMS_SetClock(clock)
!
!-----------------------------------------------------------------------
!     Set-up excgange mesh for meshded component 
!-----------------------------------------------------------------------
!
      call ROMS_SetGridArrays()
!
!-----------------------------------------------------------------------
!     Set-up import/export states
!-----------------------------------------------------------------------
!
      call ROMS_SetStates()
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetInitialize
!
      subroutine ROMS_SetRun(comp, importState, exportState, clock, rc)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_param, only : Ngrids
      use mod_scalars, only : ntstart, ntend
      use mod_stepping, only : kstp, nstp
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
      logical, save :: first = .true.
      integer, save :: tstr(Ngrids)
      integer, save :: tend(Ngrids)
      integer :: localPet, petCount, comm, nsteps, ng, rc2
!
      type(ESMF_Time) :: currTime
!
!-----------------------------------------------------------------------
!     Call ROMS initialization routines
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
! 
      call ESMF_VMGet(models(Iocean)%vm,                                &
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
      call ESMF_ClockGet (models(Iocean)%clock,                         &
                          currTime=models(Iocean)%curTime,              &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_TimeGet (models(Iocean)%curTime,                        &
                         yy=models(Iocean)%time%year,                   &
                         mm=models(Iocean)%time%month,                  &
                         dd=models(Iocean)%time%day,                    &
                         h=models(Iocean)%time%hour,                    &
                         m=models(Iocean)%time%minute,                  &
                         s=models(Iocean)%time%second,                  &
                         timeZone=models(Iocean)%time%zone,             &
                         timeStringISOFrac=models(Iocean)%time%stamp,   &
                         dayOfYear=models(Iocean)%time%yday,            &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Write current time (debug)
!-----------------------------------------------------------------------
!
      if (localPet .eq. models(Iocean)%petList(1)) then
        write(*, 20) localPet, 'Current Time',                          &
                     trim(models(Iocean)%time%stamp)
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
!     Set RCM start and end time steps
!-----------------------------------------------------------------------
!
!      nsteps=int((currTime-models(Iocean)%curTime)/models(Iocean)%dtsec)
      nsteps = int(cplTimeStep/models(Iocean)%dtsec) 
!
      if (first) then
        first = .false.
        do ng = 1, Ngrids
          tstr(ng) = ntstart(ng)
          tend(ng) = tstr(ng)+nsteps-1
        end do
      else
        do ng = 1, Ngrids
          tstr(ng) = tend(ng)+1 
          tend(ng) = tstr(ng)+nsteps-1
        end do
      end if
!
!-----------------------------------------------------------------------
!     Add extra time step at stop time to finalize ROMS IO 
!-----------------------------------------------------------------------
!
      do ng = 1, Ngrids
        if (tend(ng) == ntend(ng)) then
          tend(ng) = tend(ng)+1
        end if
      end do      
!
!-----------------------------------------------------------------------
!     Get import data 
!-----------------------------------------------------------------------
!
!      call ROMS_GetImportData (localPet, kstp, nstp, rc2)
!
!-----------------------------------------------------------------------
!     Run ROMS
!-----------------------------------------------------------------------
!
      if (localPet .eq. models(Iocean)%petList(1)) then
        write(*, fmt="(A2,3I10)") 'TR', tstr, tend, nsteps
      end if
!
      call ROMS_run (tstr, tend)
!
!-----------------------------------------------------------------------
!     Put export data
!-----------------------------------------------------------------------
!
!      call ROMS_PutExportData (localPet, kstp, nstp, rc2)
!
!-----------------------------------------------------------------------
!     Update model clock 
!-----------------------------------------------------------------------
!
      call ESMF_ClockAdvance (models(Iocean)%clock,                     &
                              timeStep=cplTimeStep,                     &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Sync the PETs 
!-----------------------------------------------------------------------
!
      call ESMF_VMBarrier(models(Iocean)%vm, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Formats 
!-----------------------------------------------------------------------
!
 20   format(' PET (', I2, ') - OCN Model ', A, ' = ', A)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetRun
!
      subroutine ROMS_SetFinalize(comp, importState, exportState,       &
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
!     Initialize return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!     Terminate ROMS execution.  Close all NetCDF files.
!-----------------------------------------------------------------------
!
      call ROMS_finalize()
!
      end subroutine ROMS_SetFinalize
!    
      subroutine ROMS_SetClock(clock)
!
!-----------------------------------------------------------------------
!     Imported modules 
!-----------------------------------------------------------------------
!
      use mod_param
      use mod_scalars
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_Clock), intent(inout) :: clock
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
!
      integer :: ng, myTimeStep
      real(r8) :: MyStartTime, MyStopTime
      real(r8) :: hour, minute, yday
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
      if (int(time_ref) == -2) then
        ref_year=1968
        ref_month=5
        ref_day=23
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='Modified Julian day number, Gregorian Calendar'
        models(Iocean)%cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,&
                                                 name=trim(name),       &
                                                 rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
      else if (int(time_ref) == -1) then
        ref_year=1
        ref_month=1
        ref_day=1
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='360-day, 30 days per month'
        models(Iocean)%cal = ESMF_CalendarCreate(ESMF_CALKIND_360DAY,   &
                                                 name=trim(name),       &
                                                 rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
      else if (int(time_ref) == 0) then
        ref_year=1
        ref_month=1
        ref_day=1
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='Julian Calendar, leap year if divisible by 4'
        models(Iocean)%cal = ESMF_CalendarCreate(ESMF_CALKIND_JULIAN,   &
                                                 name=trim(name),       &
                                                 rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
      else if (time_ref > 0.0_r8) then
        ref_year=int(r_date(2))
        ref_month=int(r_date(4))
        ref_day=int(r_date(5))
        ref_hour=int(r_date(6))
        ref_minute=int(r_date(7))
        ref_second=int(r_date(8))
        name='Gregorian Calendar'
        models(Iocean)%cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,&
                                                 name=trim(name),       &
                                                 rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
      end if
!
!-----------------------------------------------------------------------
!     Set Reference time.
!-----------------------------------------------------------------------
!
      call ESMF_TimeSet (models(Iocean)%refTime,                        &
                         yy=ref_year,                                   &
                         mm=ref_month,                                  &
                         dd=ref_day,                                    &
                         h=ref_hour,                                    &
                         m=ref_minute,                                  &
                         s=ref_second,                                  &
                         calendar=models(Iocean)%cal,                   &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set start time
!-----------------------------------------------------------------------
!
      MyStartTime = minval(tdays)
      call caldate (r_date, MyStartTime, str_year, yday, str_month,     &
                    str_day, hour)
      minute=(hour-Aint(hour))*60.0_r8
      str_hour=int(hour)
      str_minute=int(minute)
      str_second=int((minute-Aint(minute))*60.0_r8)
!
      call ESMF_TimeSet (models(Iocean)%strTime,                        &
                         yy=str_year,                                   &
                         mm=str_month,                                  &
                         dd=str_day,                                    &
                         h=str_hour,                                    &
                         m=str_minute,                                  &
                         s=str_second,                                  &
                         calendar=models(Iocean)%cal,                   &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set stop time
!-----------------------------------------------------------------------
!
      MyStopTime=0.0_r8
      do ng = 1, nNest(Iocean)
        MyStopTime=MAX(MyStopTime,                                      &
                       tdays(ng)+(REAL(ntimes(ng),r8)*dt(ng))*sec2day)
      end do
      call caldate (r_date, MyStopTime, end_year, yday, end_month,      &
                    end_day, hour)
      minute=(hour-Aint(hour))*60.0_r8
      end_hour=int(hour)
      end_minute=int(minute)
      end_second=int((minute-Aint(minute))*60.0_r8)
!
      call ESMF_TimeSet (models(Iocean)%endTime,                        &
                         yy=end_year,                                   &
                         mm=end_month,                                  &
                         dd=end_day,                                    &
                         h=end_hour,                                    &
                         m=end_minute,                                  &
                         s=end_second,                                  &
                         calendar=models(Iocean)%cal,                   &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set time interval
!-----------------------------------------------------------------------
!
      myTimeStep = int(minval(dt))
      call ESMF_TimeIntervalSet (models(Iocean)%dtsec,                  &
                                 s=myTimeStep,                          &
                                 rc=rc)  
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create time clock.
!-----------------------------------------------------------------------
!
      name='Model clock (Ocean)'
      models(Iocean)%clock = ESMF_ClockCreate (name=trim(name),         &
                                  refTime=models(Iocean)%refTime,       &
                                  timeStep=models(Iocean)%dtsec,        &
                                  startTime=models(Iocean)%strTime,     &
                                  stopTime=models(Iocean)%endTime,      &
                                  rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Copy clock
!-----------------------------------------------------------------------
!
      clock = ESMF_ClockCreate (models(Iocean)%clock, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Validate time clock
!-----------------------------------------------------------------------
!
      call ESMF_ClockValidate (models(Iocean)%clock, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get meshded component internal clock current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (models(Iocean)%clock,                         &
                          currTime=models(Iocean)%curTime,              &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Put current time into ESM_Time variable
!-----------------------------------------------------------------------
!
      call ESMF_TimeGet (models(Iocean)%curTime,                        &
                         yy=models(Iocean)%time%year,                   &
                         mm=models(Iocean)%time%month,                  &
                         dd=models(Iocean)%time%day,                    &
                         h=models(Iocean)%time%hour,                    &
                         m=models(Iocean)%time%minute,                  &
                         s=models(Iocean)%time%second,                  &
                         timeZone=models(Iocean)%time%zone,             &
                         timeStringISOFrac=models(Iocean)%time%stamp,   &
                         dayOfYear=models(Iocean)%time%yday,            &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Attach time information to export state as attribute
!-----------------------------------------------------------------------
!
      name = 'start time'
      call ESMF_AttributeSet(models(Iocean)%stateExport,                &
                             name=trim(name),                           &
                             valueList=(/ str_year  , str_month,        &
                                          str_day   , str_hour ,        &
                                          str_minute, str_second /),    &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'stop time'
      call ESMF_AttributeSet(models(Iocean)%stateExport,                &
                             name=trim(name),                           &
                             valueList=(/ end_year  , end_month,        &
                                          end_day   , end_hour ,        &
                                          end_minute, end_second /),    &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      name = 'time step'
      call ESMF_AttributeSet(models(Iocean)%stateExport,                &
                             name=trim(name),                           &
                             value=myTimeStep,                          &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Set return flag to success
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetClock
!
      subroutine ROMS_SetGridArrays ()
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_grid , only : GRID
      use mod_param, only : NtileI, NtileJ, BOUNDS, Lm, Mm
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, ii, jj, n, rc 
      integer :: localPet, petCount, comm, localDECount
      integer :: IstrR, IendR, JstrR, JendR
      integer :: IstrU, IendU, JstrU, JendU     
      integer :: IstrV, IendV, JstrV, JendV     
      integer :: staggerEdgeLWidth(2)
      integer :: staggerEdgeUWidth(2)
!
      type(ESMF_Decomp_Flag) :: deCompFlag(2)
      type(ESMF_StaggerLoc) :: staggerLoc
      real(ESMF_KIND_R8), pointer :: ptrX(:,:), ptrY(:,:)
      integer, pointer :: ptrM(:,:)
!
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!  
      call ESMF_VMGet (models(Iocean)%vm,                               &
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
      do n = 1, nNest(Iocean)
!
!-----------------------------------------------------------------------
!     Get limits of the arrays (based on PET and grid id)
!-----------------------------------------------------------------------
!
      IstrR = BOUNDS(n)%IstrR(localPet)
      IendR = BOUNDS(n)%IendR(localPet)
      JstrR = BOUNDS(n)%JstrR(localPet)
      JendR = BOUNDS(n)%JendR(localPet)
!
      IstrU = BOUNDS(n)%Istr(localPet)
      IendU = BOUNDS(n)%IendR(localPet)
      JstrU = BOUNDS(n)%JstrR(localPet)
      JendU = BOUNDS(n)%JendR(localPet)
!
      IstrV = BOUNDS(n)%IstrR(localPet)
      IendV = BOUNDS(n)%IendR(localPet)
      JstrV = BOUNDS(n)%Jstr(localPet)
      JendV = BOUNDS(n)%JendR(localPet)
!
!-----------------------------------------------------------------------
!     Create ESMF DistGrid based on model domain decomposition
!-----------------------------------------------------------------------
!
      deCompFlag = (/ ESMF_DECOMP_RESTFIRST, ESMF_DECOMP_RESTFIRST /)
!
      models(Iocean)%distGrid(n) = ESMF_DistGridCreate (                &
                                   minIndex=(/ 1, 1 /),                 &
                                   maxIndex=(/ Lm(n), Mm(n) /),         &
                                   regDecomp=(/ NtileI(n), NtileJ(n) /),&
                                   decompflag=deCompFlag,               &
                                   rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: validate and print DistGrid
!-----------------------------------------------------------------------
!
      if ((localPet == 0) .and. (cpl_debug_level > 2)) then
        call ESMF_DistGridValidate(models(Iocean)%distGrid(n), rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
        call ESMF_DistGridPrint(models(Iocean)%distGrid(n), rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
      end if
!
!-----------------------------------------------------------------------
!     Set array descriptor
!-----------------------------------------------------------------------
!
      call ESMF_ArraySpecSet(models(Iocean)%arrSpec(n),                 &
                             typekind=ESMF_TYPEKIND_R8,                 &
                             rank=2,                                    &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      do i = 1, ubound(models(Iocean)%mesh, dim=1)
!
!-----------------------------------------------------------------------
!     Set staggering type 
!-----------------------------------------------------------------------
!
      if (models(Iocean)%mesh(i,n)%gtype == Iupoint) then
        staggerLoc = ESMF_STAGGERLOC_EDGE1
        staggerEdgeLWidth = (/0,1/)
        staggerEdgeUWidth = (/1,1/)
      else if (models(Iocean)%mesh(i,n)%gtype == Ivpoint) then
        staggerLoc = ESMF_STAGGERLOC_EDGE2
        staggerEdgeLWidth = (/1,0/)
        staggerEdgeUWidth = (/1,1/)
      else if (models(Iocean)%mesh(i,n)%gtype == Icross) then
        staggerLoc = ESMF_STAGGERLOC_CENTER
        staggerEdgeLWidth = (/1,1/)
        staggerEdgeUWidth = (/1,1/)
      end if
!
!-----------------------------------------------------------------------
!     Create ESMF Grid
!-----------------------------------------------------------------------
!
      models(Iocean)%mesh(i,n)%grid = ESMF_GridCreate (                 &
                                    distgrid=models(Iocean)%distGrid(n),&
                                    indexflag=ESMF_INDEX_GLOBAL,        &
                                    gridEdgeLWidth=(/1,1/),             &
                                    gridEdgeUWidth=(/1,1/),             &
                                    name="ocn_grid",                    &
                                    rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Allocate coordinates 
!-----------------------------------------------------------------------
!
      call ESMF_GridAddCoord (models(Iocean)%mesh(i,n)%grid,            &
                              staggerLoc=staggerLoc,                    &
                              staggerEdgeLWidth=staggerEdgeLWidth,      &
                              staggerEdgeUWidth=staggerEdgeUWidth,      &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Allocate items for masking
!-----------------------------------------------------------------------
!
      call ESMF_GridAddItem (models(Iocean)%mesh(i,n)%grid,             &
                             staggerLoc=staggerLoc,                     &
                             itemflag=ESMF_GRIDITEM_MASK,               &
                             rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iocean)%mesh(i,n)%grid,                 &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointers and set coordinates for the grid 
!-----------------------------------------------------------------------
! 
      do j = 0, localDECount-1
        call ESMF_GridGetCoord (models(Iocean)%mesh(i,n)%grid,          &
                                localDE=j,                              &
                                staggerLoc=staggerLoc,                  &
                                coordDim=1,                             &
                                farrayPtr=ptrX,                         &
                                rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
        call ESMF_GridGetCoord (models(Iocean)%mesh(i,n)%grid,          &
                                localDE=j,                              &
                                staggerLoc=staggerLoc,                  &
                                coordDim=2,                             &
                                farrayPtr=ptrY,                         &
                                rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
        call ESMF_GridGetItem (models(Iocean)%mesh(i,n)%grid,           &
                               localDE=j,                               &
                               staggerLoc=staggerLoc,                   &
                               itemflag=ESMF_GRIDITEM_MASK,             &
                               farrayPtr=ptrM,                          &
                               rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
        if (models(Iocean)%mesh(i,n)%gtype == Icross) then
          do jj = JstrR, JendR
            do ii = IstrR, IendR
              ptrX(ii,jj) = GRID(n)%lonr(ii,jj)
              ptrY(ii,jj) = GRID(n)%latr(ii,jj)
              ptrM(ii,jj) = int(GRID(n)%rmask(ii,jj))
            end do
          end do       
!
          if (cpl_debug_level > 2) then
            write(*,99) localPet, j, 'R-M', IstrR, IendR,  JstrR, JendR
            write(*,99) localPet, j, 'R-E', lbound(ptrY, dim=1),        &
                        ubound(ptrY, dim=1), lbound(ptrY, dim=2),       &
                        ubound(ptrY, dim=2) 
          end if
!
          if (cpl_vtk_on) then
            call ESMF_GridWriteVTK(models(Iocean)%mesh(i,n)%grid,       &
                                   filename="ocean_RHOpoint")
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
        else if (models(Iocean)%mesh(i,n)%gtype == Iupoint) then
          do jj = JstrU, JendU
            do ii = IstrU, IendU
              ptrX(ii,jj) = GRID(n)%lonu(ii,jj)
              ptrY(ii,jj) = GRID(n)%latu(ii,jj)
              ptrM(ii,jj) = int(GRID(n)%rmask(ii,jj))
            end do
          end do
!
          if (cpl_debug_level > 2) then
            write(*,99) localPet, j, 'U-M', IstrU, IendU,  JstrU, JendU
            write(*,99) localPet, j, 'U-E', lbound(ptrY, dim=1),        &
                        ubound(ptrY, dim=1), lbound(ptrY, dim=2),       &
                        ubound(ptrY, dim=2)
          end if
!
          if (cpl_vtk_on) then
            call ESMF_GridWriteVTK(models(Iocean)%mesh(i,n)%grid,       &
                                   filename="ocean_Upoint")
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
        else if (models(Iocean)%mesh(i,n)%gtype == Ivpoint) then
          do jj = JstrV, JendV
            do ii = IstrV, IendV
              ptrX(ii,jj) = GRID(n)%lonv(ii,jj)
              ptrY(ii,jj) = GRID(n)%latv(ii,jj)
              ptrM(ii,jj) = int(GRID(n)%vmask(ii,jj))
            end do
          end do
!
          if (cpl_debug_level > 2) then
            write(*,99) localPet, j, 'V-M', IstrU, IendU,  JstrU, JendU
            write(*,99) localPet, j, 'V-E', lbound(ptrY, dim=1),        &
                        ubound(ptrY, dim=1), lbound(ptrY, dim=2),       &
                        ubound(ptrY, dim=2)
          end if
!        
          if (cpl_vtk_on) then       
            call ESMF_GridWriteVTK(models(Iocean)%mesh(i,n)%grid,       &
                                   filename="ocean_Vpoint")
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
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
        if (associated(ptrM)) then
          nullify(ptrM)
        end if
!
      end do
      end do
      end do
!
!-----------------------------------------------------------------------
!     Format definition 
!-----------------------------------------------------------------------
!
 99   format(" PET(",I1,") - DE(",I1,") - ", A3, " : ", 4I8)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetGridArrays
!
      subroutine ROMS_SetStates ()
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_ocean, only : OCEAN
      use mod_scalars, only : itemp
      use mod_stepping, only : nstp
      use mod_param, only : NtileI, NtileJ, BOUNDS, N, Lm, Mm
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, ii, jj, id, ng, rc
      integer :: localPet, petCount, comm, localDECount
      integer :: IstrR, IendR, JstrR, JendR
      integer :: IstrU, IendU, JstrU, JendU     
      integer :: IstrV, IendV, JstrV, JendV
      integer :: staggerEdgeLWidth(2)
      integer :: staggerEdgeUWidth(2)
      character (len=40) :: name
!
      type(ESMF_StaggerLoc) :: staggerLoc
      real(ESMF_KIND_R8), pointer :: ptr(:,:)
!
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!  
      call ESMF_VMGet (models(Iocean)%vm,                               &
                       localPet=localPet,                               &
                       petCount=petCount,                               &
                       mpiCommunicator=comm,                            &
                       rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Loop over number of nested/composed meshs.
!-----------------------------------------------------------------------
!
      do ng = 1, nNest(Iocean)
!
!-----------------------------------------------------------------------
!     Get limits of the arrays (based on PET and grid id)
!-----------------------------------------------------------------------
!
      IstrR = BOUNDS(ng)%IstrR(localPet)
      IendR = BOUNDS(ng)%IendR(localPet)
      JstrR = BOUNDS(ng)%JstrR(localPet)
      JendR = BOUNDS(ng)%JendR(localPet)
!
      IstrU = BOUNDS(ng)%Istr(localPet)
      IendU = BOUNDS(ng)%IendR(localPet)
      JstrU = BOUNDS(ng)%JstrR(localPet)
      JendU = BOUNDS(ng)%JendR(localPet)
!
      IstrV = BOUNDS(ng)%IstrR(localPet)
      IendV = BOUNDS(ng)%IendR(localPet)
      JstrV = BOUNDS(ng)%Jstr(localPet)
      JendV = BOUNDS(ng)%JendR(localPet)
!
!-----------------------------------------------------------------------
!     Create export state arrays.
!-----------------------------------------------------------------------
!
      do i = 1, ubound(models(Iocean)%dataExport(:,ng), dim=1)
!
!-----------------------------------------------------------------------
!     Set staggering type 
!-----------------------------------------------------------------------
!
      if (models(Iocean)%dataExport(i,ng)%gtype == Iupoint) then
        staggerLoc = ESMF_STAGGERLOC_EDGE1
        staggerEdgeLWidth = (/0,1/)
        staggerEdgeUWidth = (/1,1/)
      else if (models(Iocean)%dataExport(i,ng)%gtype == Ivpoint) then
        staggerLoc = ESMF_STAGGERLOC_EDGE2
        staggerEdgeLWidth = (/1,0/)
        staggerEdgeUWidth = (/1,1/)
      else if (models(Iocean)%dataExport(i,ng)%gtype == Icross) then
        staggerLoc = ESMF_STAGGERLOC_CENTER
        staggerEdgeLWidth = (/1,1/)
        staggerEdgeUWidth = (/1,1/)
      end if
!
!-----------------------------------------------------------------------
!     Create ESMF Fields 
!-----------------------------------------------------------------------
!
      models(Iocean)%dataExport(i,ng)%field = ESMF_FieldCreate (        &
                        models(Iocean)%mesh(i,ng)%grid,                 &
                        models(Iocean)%arrSpec(ng),                     &
                        staggerLoc=staggerLoc,                          &
                        name=trim(models(Iocean)%dataExport(i,ng)%name),&
                        rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iocean)%mesh(i,ng)%grid,                &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointers and put data 
!-----------------------------------------------------------------------
! 
      do j = 0, localDECount-1
      call ESMF_FieldGet (models(Iocean)%dataExport(i,ng)%field,        &
                          localDE=j,                                    &
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
!     Put data in it 
!-----------------------------------------------------------------------
!
      name = models(Iocean)%dataExport(i,ng)%name
!
      if (trim(adjustl(name)) == "SST") then
        do jj = JstrR, JendR
          do ii = IstrR, IendR
            ptr(ii,jj) = OCEAN(ng)%t(ii,jj,N(ng),nstp(ng),itemp) 
          end do
        end do 
        write(*,99) localPet, j, 'R-I', IstrR, IendR,  JstrR, JendR
        write(*,99) localPet, j, 'R-E', lbound(ptr, dim=1),             &
                    ubound(ptr, dim=1), lbound(ptr, dim=2),             &
                    ubound(ptr, dim=2)         
      end if        
      end do
!
!-----------------------------------------------------------------------
!     Add fields to export state
!-----------------------------------------------------------------------
!
      call ESMF_StateAdd (models(Iocean)%stateExport,                   &
                         (/ models(Iocean)%dataExport(i,ng)%field /),   &
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
!     Create import state arrays.
!-----------------------------------------------------------------------
!
      do i = 1, ubound(models(Iocean)%dataImport(:,n), dim=1)
!
!-----------------------------------------------------------------------
!     Set staggering type 
!-----------------------------------------------------------------------
!
      if (models(Iocean)%dataImport(i,ng)%gtype == Iupoint) then
        id = getMeshID(models(Iocean)%mesh(:,ng), Iupoint)
        staggerLoc = ESMF_STAGGERLOC_EDGE1
        staggerEdgeLWidth = (/0,1/)
        staggerEdgeUWidth = (/1,1/)
      else if (models(Iocean)%dataImport(i,ng)%gtype == Ivpoint) then
        id = getMeshID(models(Iocean)%mesh(:,ng), Ivpoint)
        staggerLoc = ESMF_STAGGERLOC_EDGE2
        staggerEdgeLWidth = (/1,0/)
        staggerEdgeUWidth = (/1,1/)
      else if (models(Iocean)%dataImport(i,ng)%gtype == Icross) then
        id = getMeshID(models(Iocean)%mesh(:,ng), Icross)
        staggerLoc = ESMF_STAGGERLOC_CENTER
        staggerEdgeLWidth = (/1,1/)
        staggerEdgeUWidth = (/1,1/)
      end if
!
!-----------------------------------------------------------------------
!     Create ESMF Fields 
!-----------------------------------------------------------------------
!
      models(Iocean)%dataImport(i,ng)%field = ESMF_FieldCreate (        &
                        models(Iocean)%mesh(id,ng)%grid,                &
                        models(Iocean)%arrSpec(ng),                     &
                        staggerLoc=staggerLoc,                          &
                        name=trim(models(Iocean)%dataImport(i,ng)%name),&
                        rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iocean)%mesh(id,ng)%grid,               &
                         localDECount=localDECount,                     &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get pointers and initialize it 
!-----------------------------------------------------------------------
! 
      do j = 0, localDECount-1
      call ESMF_FieldGet (models(Iocean)%dataImport(i,ng)%field,        &
                          localDE=j,                                    &
                          farrayPtr=ptr,                                &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      ptr = 0.0d0
      end do
!
!      call ESMF_FieldWrite(models(Iocean)%dataImport(i,ng)%field,        &
!                           'dst_field.nc',                              &
!                           rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Add fields to import state
!-----------------------------------------------------------------------
!
      call ESMF_StateAdd (models(Iocean)%stateImport,                   &
                         (/ models(Iocean)%dataImport(i,ng)%field /),   &
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
 99   format(" PET(",I1,") - DE(",I1,") - ", A3, " : ", 4I8)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetStates
!
      subroutine ROMS_PutExportData (localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_param, only : BOUNDS, N
      use mod_ocean, only : OCEAN
      use mod_scalars, only : itemp
      use mod_stepping, only : nstp
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
      character (len=80) :: name
      integer :: i, j, k, ng, rc
      integer :: Istr, Iend, Jstr, Jend
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV
!
!-----------------------------------------------------------------------
!     Loop over number of nested/composed meshs.
!-----------------------------------------------------------------------
!
      do ng = 1, nNest(Iocean)
!
!-----------------------------------------------------------------------
!     Load export fields.
!-----------------------------------------------------------------------
!
      Istr =BOUNDS(ng)%Istr (localPet)
      Iend =BOUNDS(ng)%Iend (localPet)
      Jstr =BOUNDS(ng)%Jstr (localPet)
      Jend =BOUNDS(ng)%Jend (localPet)
      IstrR=BOUNDS(ng)%IstrR(localPet)
      IendR=BOUNDS(ng)%IendR(localPet)
      IstrU=BOUNDS(ng)%IstrU(localPet)
      JstrR=BOUNDS(ng)%JstrR(localPet)
      JendR=BOUNDS(ng)%JendR(localPet)
      JstrV=BOUNDS(ng)%JstrV(localPet) 
!
      do k = 1, ubound(models(Iocean)%dataExport(:,ng), dim=1) 
        name = models(Iocean)%dataExport(k,ng)%name
!
        if (trim(adjustl(name)) == "SST") then 
          do j = JstrR, JendR
            do i = IstrR, IendR
              models(Iocean)%dataExport(k,ng)%ptr =                    &
                        OCEAN(ng)%t(i,j,N(ng),nstp(ng),itemp)
            end do
          end do
        end if
      end do
      end do
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_PutExportData

      end module mod_esmf_ocn
