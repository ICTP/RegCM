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
      use mod_couplerr
      use mod_regcm_interface, only : RCM_initialize,                   &
                                      RCM_run,                          &
                                      RCM_finalize
      use mod_stdio, only : stderr
!
      implicit none
      private
!
!-----------------------------------------------------------------------
!     Public subroutines 
!-----------------------------------------------------------------------
!
      public :: RCM_SetServices
      public :: RCM_SetRun
      public :: RCM_SetFinalize
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
      integer(ik4) :: localPet, petCount, ierr
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
      integer(ik4) :: localPet, petCount, comm
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
      if ((cpl_dbglevel > 0) .and. (localPet == 0)) then
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
!
      if (localPet == 0) then
      call ESMF_TimePrint(models(Iatmos)%strTime,                       &
                          options="string isofrac",                     &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_TimePrint(models(Iatmos)%curTime,                       &
                          options="string isofrac",                     &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_TimeIntervalPrint(dt1,                                          &
                          options="string isofrac",                     &
                          rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if

      call ESMF_TimeIntervalGet (dt1,                                   &
                                 s_r8=timestr1,                         &
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      if (localPet == 0) then
        print*, timestr1
      end if
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
      integer(ik4) :: i, n
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
      use mod_runparams, only : idate0, idate1, idate2, dtsec 
      use mod_runparams, only : split_idate
      use mod_runparams, only : cpldt, cpldbglevel
      use mod_dynparam , only : calendar
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
      integer(ik4) :: ref_year,   str_year,   end_year
      integer(ik4) :: ref_month,  str_month,  end_month
      integer(ik4) :: ref_day,    str_day,    end_day
      integer(ik4) :: ref_hour,   str_hour,   end_hour
      integer(ik4) :: ref_minute, str_minute, end_minute
      integer(ik4) :: ref_second, str_second, end_second
      character (len=80) :: name
      integer(ik4) :: rc
!
!-----------------------------------------------------------------------
!     Create gridded component clock 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Create ESMF calendar
!-----------------------------------------------------------------------
!
      if (calendar == 'gregorian') then
        name = "The Gregorian calendar"
        models(Iatmos)%cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,&
                                                 name=trim(name),       &
                                                 rc=rc)
      else if (calendar == 'noleap' .or. calendar == '365_day') then
        name = "The no-leap calendar"
        models(Iatmos)%cal = ESMF_CalendarCreate(ESMF_CALKIND_NOLEAP,   &
                                                 name=trim(name),       &
                                                 rc=rc)
      else if (calendar == '360_day') then
        name = "The 360-day calendar"
        models(Iatmos)%cal = ESMF_CalendarCreate(ESMF_CALKIND_360DAY,   &
                                                 name=trim(name),       &
                                                 rc=rc)
      else
        name = "The Gregorian calendar"
        models(Iatmos)%cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,&
                                                 name=trim(name),       &
                                                 rc=rc)
      end if
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
      subroutine RCM_SetGridArrays()
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_atm_interface, only : mddom
      use mod_dynparam, only : iy, jx, nproc 
      use mod_runparams, only : cpldt, cpldbglevel
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer(ik4) :: i, j, k, n, rc
      integer(ik4) :: localPet, petCount, comm, localDECount
      integer(ik4) :: unmapped(nproc), mapped(nproc)
      integer, dimension(2) :: cpus_per_dim
      character (len=40) :: name
      character (len=100) :: fmt_123
!
      type(ESMF_StaggerLoc) :: staggerLoc
      type(ESMF_Decomp_Flag) :: decompflag(2)
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
!     Calculate number of CPUs in each direction 
!-----------------------------------------------------------------------
!
      if ( nproc < 4 ) then
        cpus_per_dim(2) = 1
        cpus_per_dim(1) = nproc
      else if ( nproc >= 4 ) then
        cpus_per_dim(2) = (nint(sqrt(dble(nproc)))/2)*2
        if ( iy > int(1.5*dble(jx)) ) then
          cpus_per_dim(2) = cpus_per_dim(2) - 1
          do while ( mod(nproc,cpus_per_dim(2)) /= 0 )
            cpus_per_dim(2) = cpus_per_dim(2) - 1
          end do
        else if ( jx > int(1.5*dble(iy)) ) then
          cpus_per_dim(2) = cpus_per_dim(2) + 1
          do while ( mod(nproc,cpus_per_dim(2)) /= 0 )
            cpus_per_dim(2) = cpus_per_dim(2) + 1
          end do
        else
          do while ( mod(nproc,cpus_per_dim(2)) /= 0 )
            cpus_per_dim(2) = cpus_per_dim(2) + 1
          end do
        end if
        cpus_per_dim(1) = nproc/cpus_per_dim(2)
      end if
!
!-----------------------------------------------------------------------
!     Create ESMF DistGrid based on model domain decomposition
!
!     ESMF is basically using a right handed coordinate system, and 
!     using the Fortran way of using the smallest stride to the first 
!     dimension but RegCM not. The order of dimension is reversed
!     because of this limitation. 
!-----------------------------------------------------------------------
!
      decompflag = (/ ESMF_DECOMP_RESTLAST, ESMF_DECOMP_RESTLAST /)
!
      models(Iatmos)%distGrid(n) = ESMF_DistGridCreate (                &
                                        minIndex=(/ 1, 1 /),            &
                                        maxIndex=(/ iy, jx /),          &
                                        regDecomp=cpus_per_dim,         &
                                        decompflag=decompflag,          &
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
                                    gridEdgeLWidth=(/0,0/),             &
                                    gridEdgeUWidth=(/0,0/),             &
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
      call ESMF_GridGetCoord (models(Iatmos)%grid(n),                   &
                              localDE=j,                                &
                              staggerLoc=staggerLoc,                    &
                              coordDim=1,                               &
                              farrayPtr=ptrX,                           &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) then
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
      call ESMF_GridGetCoord (models(Iatmos)%grid(n),                   &
                              localDE=j,                                &
                              staggerLoc=staggerLoc,                    &
                              coordDim=2,                               &
                              farrayPtr=ptrY,                           &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) then
        call ESMF_Finalize(endflag=ESMF_END_ABORT)       
      end if
!
!-----------------------------------------------------------------------
!     Debug: write size of pointers    
!-----------------------------------------------------------------------
!
      name = GRIDDES(models(Iatmos)%mesh(i,n)%gtype)
!
      write(*,30) localPet, j, adjustl("PTR/ATM/GRD/"//name),           &
                  lbound(ptrX, dim=1), ubound(ptrX, dim=1),             &
                  lbound(ptrX, dim=2), ubound(ptrX, dim=2)
!
!-----------------------------------------------------------------------
!     Fill the pointers    
!-----------------------------------------------------------------------
!
      if (models(Iatmos)%mesh(i,n)%gtype == Idot) then
        write(*,30) localPet, j, adjustl("DAT/ATM/GRD/"//name),         &
                 lbound(mddom%dlon, dim=1), ubound(mddom%dlon, dim=1),  &
                 lbound(mddom%dlon, dim=2), ubound(mddom%dlon, dim=2)
!
        ptrX = transpose(mddom%dlon)
        ptrY = transpose(mddom%dlat)
      else if (models(Iatmos)%mesh(i,n)%gtype == Icross) then
        write(*,30) localPet, j, adjustl("DAT/ATM/GRD/"//name),         &
                 lbound(mddom%xlon, dim=1), ubound(mddom%xlon, dim=1),  &
                 lbound(mddom%xlon, dim=2), ubound(mddom%xlon, dim=2)
!
        ptrX = transpose(mddom%xlon)
        ptrY = transpose(mddom%xlat)
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
      if (cpldbglevel > 1) then
      call ESMF_GridValidate(models(Iatmos)%grid(n), rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
!
!-----------------------------------------------------------------------
!     Write ESMF Grid in VTK format (debug) 
!-----------------------------------------------------------------------
!
      if (cpldbglevel > 1) then
      write(stderr,*) '[debug] -- write grid information to file '//    &
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
 30   format(" PET(",I3,") - DE(",I2,") - ", A20, " : ", 4I8)
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
      use mod_runparams, only : cpldt, cpldbglevel
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
      logical :: flag
      character (len=40) :: name
      integer(ik4) :: i, j, id, n, rc, localDECount
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
      else if (models(Iatmos)%dataExport(i,n)%gtype == Idot) then
        staggerLoc = ESMF_STAGGERLOC_CORNER
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
      models(Iatmos)%dataExport(i,n)%ptr = MISSING_R8
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
                                  models(Iatmos)%grid(n),               &
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
      end module mod_esmf_atm
