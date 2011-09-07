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
!
!-----------------------------------------------------------------------
!     Register "run" routine    
!-----------------------------------------------------------------------
!
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_RUN,       &
                                      userRoutine=ROMS_SetRun,          &
                                      rc=rc)
!
!-----------------------------------------------------------------------
!     Register "finalize" routine    
!-----------------------------------------------------------------------
! 
      call ESMF_GridCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_FINALIZE,  &
                                      userRoutine=ROMS_SetFinalize,     &
                                      rc=rc)
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
!
      call ESMF_VMGet(models(Iocean)%vm,                                &
                      localPet=localPet,                                &
                      petCount=petCount,                                &
                      mpiCommunicator=comm,                             &
                      rc=rc)
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
      call ROMS_SetClock(clock, rc)
!
!-----------------------------------------------------------------------
!     Set-up excgange mesh for meshded component 
!-----------------------------------------------------------------------
!
      call ROMS_SetGridArrays(comp, rc)
!
!-----------------------------------------------------------------------
!     Load ROMS exchange mesh arrays
!-----------------------------------------------------------------------
!
!      call ROMS_PutGridData(MyRank, rc)
!
!-----------------------------------------------------------------------
!     Set-up import/export states
!-----------------------------------------------------------------------
!
!      call ROMS_SetStates(ImportState, ExportState, MyRank, rc)
!
!-----------------------------------------------------------------------
!     Load export initial conditions data.
!-----------------------------------------------------------------------
!
!      call ROMS_PutExportData(Myrank, rc)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetInitialize
!
      subroutine ROMS_SetRun(comp, importState, exportState,            &
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
!     Call ROMS initialization routines
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
      end subroutine ROMS_SetRun
!
      subroutine ROMS_SetFinalize(comp, importState, exportState,       &
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
!     Terminate ROMS execution.  Close all NetCDF files.
!-----------------------------------------------------------------------
!
      call ROMS_finalize()
!
      end subroutine ROMS_SetFinalize
!    
      subroutine ROMS_SetClock(clock, rc)
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
      integer, intent(inout) :: rc
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
!
!-----------------------------------------------------------------------
!     Set time interval
!-----------------------------------------------------------------------
!
      myTimeStep = int(minval(dt))
      call ESMF_TimeIntervalSet (models(Iocean)%dtsec,                  &
                                 s=myTimeStep,                          &
                                 rc=rc)  
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
!
!-----------------------------------------------------------------------
!     Copy clock
!-----------------------------------------------------------------------
!
      clock = ESMF_ClockCreate (models(Iocean)%clock, rc=rc)
!
!-----------------------------------------------------------------------
!     Validate time clock
!-----------------------------------------------------------------------
!
      call ESMF_ClockValidate (models(Iocean)%clock,                    &
                               rc=rc)
!
!-----------------------------------------------------------------------
!     Get meshded component internal clock current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (models(Iocean)%clock,                         &
                          currTime=models(Iocean)%curTime,              &
                          rc=rc)
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
      name = 'stop time'
      call ESMF_AttributeSet(models(Iocean)%stateExport,                &
                             name=trim(name),                           &
                             valueList=(/ end_year  , end_month,        &
                                          end_day   , end_hour ,        &
                                          end_minute, end_second /),    &
                             rc=rc)
      name = 'time step'
      call ESMF_AttributeSet(models(Iocean)%stateExport,                &
                             name=trim(name),                           &
                             value=myTimeStep,                          &
                             rc=rc)
!
!-----------------------------------------------------------------------
!     Set return flag to success
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetClock
!
      subroutine ROMS_SetGridArrays (gcomp, rc)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_param, only : NtileI, NtileJ, BOUNDS, Lm, Mm
      use mod_dynparam, only : debug_level
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESMF_GridComp), intent(inout) :: gcomp
      integer, intent(out) :: rc
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, n 
      integer :: localPet, petCount, comm, localDECount
!
      type(ESMF_Field) :: grdField
      type(ESMF_StaggerLoc) :: staggerLoc
      real(ESMF_KIND_R8), pointer :: ptrX(:,:), ptrY(:,:)
      integer, pointer :: ptrMask(:,:)
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
!
!-----------------------------------------------------------------------
!     Set RCM domain decomposition variables
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iocean)
!
!-----------------------------------------------------------------------
!     Create ESMF DistGrid based on model domain decomposition
!-----------------------------------------------------------------------
!
      models(Iocean)%distGrid(n) = ESMF_DistGridCreate (                &
                                   minIndex=(/ 1, 1 /),                 &
                                   maxIndex=(/ Lm(n), Mm(n) /),         &
                                   regDecomp=(/ NtileI(n), NtileJ(n) /),&
                                   rc=rc)
!
!-----------------------------------------------------------------------
!     Debug: validate and print DistGrid
!-----------------------------------------------------------------------
!
      if ((debug_level > 3) .and. (localPet == 0)) then
        call ESMF_DistGridValidate(models(Iocean)%distGrid(n), rc=rc)
        call ESMF_DistGridPrint(models(Iocean)%distGrid(n), rc=rc)
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
!
      do i = 1, ubound(models(Iocean)%mesh, dim=1)
!
!-----------------------------------------------------------------------
!     Set staggering type 
!-----------------------------------------------------------------------
!
      if (models(Iocean)%mesh(i,n)%gtype == Iupoint) then
        staggerLoc = ESMF_STAGGERLOC_EDGE1
      else if (models(Iocean)%mesh(i,n)%gtype == Ivpoint) then
        staggerLoc = ESMF_STAGGERLOC_EDGE2
      else
        staggerLoc = ESMF_STAGGERLOC_CENTER
      end if
!
!-----------------------------------------------------------------------
!     Create ESMF Grid
!-----------------------------------------------------------------------
!
      models(Iocean)%mesh(i,n)%grid = ESMF_GridCreate (                 &
                                    distgrid=models(Iocean)%distGrid(n),&
                                    name="ocn_grid",                    &
                                    rc=rc)
!
!-----------------------------------------------------------------------
!     Allocate coordinates 
!-----------------------------------------------------------------------
!
      call ESMF_GridAddCoord (models(Iocean)%mesh(i,n)%grid,            &
                              staggerLoc=staggerLoc,                    &
                              rc=rc)
!
!-----------------------------------------------------------------------
!     Allocate items for masking
!-----------------------------------------------------------------------
!
      call ESMF_GridAddItem (models(Iocean)%mesh(i,n)%grid,             &
                             staggerLoc=staggerLoc,                     &
                             itemflag=ESMF_GRIDITEM_MASK,               &
                             rc=rc)
!
!-----------------------------------------------------------------------
!     Get number of local DEs
!-----------------------------------------------------------------------
! 
      call ESMF_GridGet (models(Iocean)%mesh(i,n)%grid,                 &
                         localDECount=localDECount,                     &
                         rc=rc)
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
                                farrayPtr=ptrY,                         &
                                rc=rc)
        ptrY = 0.0
!
        call ESMF_GridGetCoord (models(Iocean)%mesh(i,n)%grid,          &
                                localDE=j,                              &
                                staggerLoc=staggerLoc,                  &
                                coordDim=2,                             &
                                farrayPtr=ptrX,                         &
                                rc=rc)
        ptrX = 0.0
!
        call ESMF_GridGetItem (models(Iocean)%mesh(i,n)%grid,           &
                               localDE=j,                               &
                               staggerLoc=staggerLoc,                   &
                               itemflag=ESMF_GRIDITEM_MASK,             &
                               farrayPtr=ptrMask,                       &
                               rc=rc)
        ptrMask = 0
      end do
      end do
      end do
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine ROMS_SetGridArrays
!
      subroutine ROMS_PutGridData (localPet, rc)
!
!**********************************************************************
!
!     Used module declarations 
!
!**********************************************************************
!
      use mod_param, only : BOUNDS
      use mod_grid , only : GRID
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
      integer :: i, j, k, n
      integer :: Istr, Jstr, IstrR, IendR, JstrR, JendR
!
!-----------------------------------------------------------------------
!     Loop over number of nested/composed meshs.
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iocean)
!
!-----------------------------------------------------------------------
!     Get limits of the arrays (based on PET and grid id)
!-----------------------------------------------------------------------
!
      Istr  = BOUNDS(n)%Istr (localPet)
      Jstr  = BOUNDS(n)%Jstr (localPet)
      IstrR = BOUNDS(n)%IstrR(localPet)
      IendR = BOUNDS(n)%IendR(localPet)
      JstrR = BOUNDS(n)%JstrR(localPet)
      JendR = BOUNDS(n)%JendR(localPet)
!     
!-----------------------------------------------------------------------
!     Load grid data for grid points 
!-----------------------------------------------------------------------
! 
      do k = 1, ubound(models(Iocean)%mesh, dim=1)
!       RHO 
        if (models(Iocean)%mesh(k,n)%gtype == Icross) then
          do j = JstrR, JendR
            do i = IstrR, IendR
              models(Iocean)%mesh(k,n)%lon%ptr(i,j)=GRID(n)%lonr(i,j)
              models(Iocean)%mesh(k,n)%lat%ptr(i,j)=GRID(n)%latr(i,j)
!              models(Iocean)%mesh(k,n)%mask%field(i,j)=GRID(n)%rmask(i,j)
            end do
          end do
!       U
        else if (models(Iocean)%mesh(k,n)%gtype == Iupoint) then
          do j = JstrR, JendR
            do i = Istr, IendR
              models(Iocean)%mesh(k,n)%lon%ptr(i,j)=GRID(n)%lonu(i,j)
              models(Iocean)%mesh(k,n)%lat%ptr(i,j)=GRID(n)%latu(i,j)
!              models(Iocean)%mesh(k,n)%mask%field(i,j)=GRID(n)%umask(i,j)
            end do
          end do
!       V
        else if (models(Iocean)%mesh(k,n)%gtype == Ivpoint) then  
          do j = Jstr, JendR
            do i = IstrR, IendR
              models(Iocean)%mesh(k,n)%lon%ptr(i,j)=GRID(n)%lonv(i,j)
              models(Iocean)%mesh(k,n)%lat%ptr(i,j)=GRID(n)%latv(i,j)
!              models(Iocean)%mesh(k,n)%mask%field(i,j)=GRID(n)%vmask(i,j)
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
      end subroutine ROMS_PutGridData
!
      subroutine ROMS_SetStates (importState, exportState, MyRank,      &
                                 rc)
!
!**********************************************************************
!
!     Used module declarations 
!
!**********************************************************************
!
      use mod_param, only : NtileI, NtileJ, BOUNDS, Lm, Mm
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_State), intent(inout) :: importState
      type(ESMF_State), intent(inout) :: exportState
      integer, intent(in) :: MyRank
      integer, intent(inout) :: rc
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
      integer :: i, n
      integer, dimension(2) :: CLW, CUW, TLW, TUW
      integer, dimension(2) :: CLW_r, CUW_r
      integer, dimension(2) :: CLW_u, CUW_u
      integer, dimension(2) :: CLW_v, CUW_v
!
!-----------------------------------------------------------------------
!     Loop over number of nested/composed meshs.
!-----------------------------------------------------------------------
!
      do n = 1, nNest(Iocean)
!
!-----------------------------------------------------------------------
!     Set ROMS domain decomposition variables
!-----------------------------------------------------------------------
!
      TLW(1)=BOUNDS(n)%Istr(MyRank)-BOUNDS(n)%LBi(MyRank)
      TLW(2)=BOUNDS(n)%Jstr(MyRank)-BOUNDS(n)%LBj(MyRank)
      TUW(1)=BOUNDS(n)%UBi(MyRank)-BOUNDS(n)%Iend(MyRank)
      TUW(2)=BOUNDS(n)%UBj(MyRank)-BOUNDS(n)%Jend(MyRank)
!
      CLW_r(1)=BOUNDS(n)%Istr(MyRank)-BOUNDS(n)%IstrR(MyRank)
      CLW_r(2)=BOUNDS(n)%Jstr(MyRank)-BOUNDS(n)%JstrR(MyRank)
      CUW_r(1)=BOUNDS(n)%IendR(MyRank)-BOUNDS(n)%Iend(MyRank)
      CUW_r(2)=BOUNDS(n)%JendR(MyRank)-BOUNDS(n)%Jend(MyRank)
!
      CLW_u(1)=0
      CLW_u(2)=CLW_r(2)
      CUW_u(1)=CUW_r(1)
      CUW_u(2)=CUW_r(2)
!
      CLW_v(1)=CLW_r(1)
      CLW_v(2)=0
      CUW_v(1)=CUW_r(1)
      CUW_v(2)=CUW_r(2)
!
!-----------------------------------------------------------------------
!     Create export state arrays.
!-----------------------------------------------------------------------
!
      do i = 1, ubound(models(Iocean)%dataExport(:,n), dim=1)
        if (models(Iocean)%dataExport(i,n)%gtype == Iupoint) then
          CLW = (/ CLW_u(1), CLW_u(2) /)
          CUW = (/ CUW_u(1), CUW_u(2) /)
        else if (models(Iocean)%dataExport(i,n)%gtype == Ivpoint) then 
          CLW = (/ CLW_v(1), CLW_v(2) /)
          CUW = (/ CUW_v(1), CUW_v(2) /)
        else
          CLW = (/ CLW_r(1), CLW_r(2)/)
          CUW = (/ CUW_r(1), CUW_r(2)/)
        end if
!
!       Create array
!
        models(Iocean)%dataExport(i,n)%array = ESMF_ArrayCreate (       &
                                  arrayspec=models(Iocean)%arrSpec(n),  &
                                  distgrid=models(Iocean)%distGrid(n),  &
                                     computationalLWidth=CLW,           &
                                     computationalUWidth=CUW,           &
                                     totalLWidth=TLW,                   &
                                     totalUWidth=TUW,                   &
                                     indexflag=ESMF_INDEX_GLOBAL,       &
                                     rc=rc)
!
!       Get data pointer from array
!
        call ESMF_ArrayGet (models(Iocean)%dataExport(i,n)%array,       &
                 farrayPtr=models(Iocean)%dataExport(i,n)%ptr,        &
                 rc=rc) 
!
!       Set array name
!
        call ESMF_ArraySet (array=models(Iocean)%dataExport(i,n)%array, &
                 name=trim(models(Iocean)%dataExport(i,n)%long_name),   &
                 rc=rc)
!
!       Add array to export state
!     
        call ESMF_StateAdd (models(Iocean)%stateExport,                 &
                            (/ models(Iocean)%dataExport(i,n)%array /), &
                            rc=rc)
!
!       Initialize export field to zero to avoid infinities or NaNs.
!
        models(Iocean)%dataExport(i,n)%ptr = 0.0d0
      end do
!
!-----------------------------------------------------------------------
!     Create import state arrays.
!-----------------------------------------------------------------------
!
      do i = 1, ubound(models(Iocean)%dataImport(:,n), dim=1)
        if (models(Iocean)%dataImport(i,n)%gtype == Iupoint) then
          CLW = (/ CLW_u(1), CLW_u(2) /)
          CUW = (/ CUW_u(1), CUW_u(2) /)
        else if (models(Iocean)%dataImport(i,n)%gtype == Ivpoint) then
          CLW = (/ CLW_v(1), CLW_v(2) /)
          CUW = (/ CUW_v(1), CUW_v(2) /)
        else
          CLW = (/ CLW_r(1), CLW_r(2)/)
          CUW = (/ CUW_r(1), CUW_r(2)/)
        end if
!
!       Create array
!
        models(Iocean)%dataImport(i,n)%array = ESMF_ArrayCreate (       &
                                  arrayspec=models(Iocean)%arrSpec(n),  &
                                  distgrid=models(Iocean)%distGrid(n),  &
                                     computationalLWidth=CLW,           &
                                     computationalUWidth=CUW,           &
                                     totalLWidth=TLW,                   &
                                     totalUWidth=TUW,                   &
                                     indexflag=ESMF_INDEX_GLOBAL,       &
                                     rc=rc)
!
!       Get data pointer from array
!
        call ESMF_ArrayGet (models(Iocean)%dataImport(i,n)%array,       &
                 farrayPtr=models(Iocean)%dataImport(i,n)%ptr,        &
                 rc=rc)
!
!       Set array name
!
        call ESMF_ArraySet (array=models(Iocean)%dataImport(i,n)%array, &
                 name=trim(models(Iocean)%dataImport(i,n)%long_name),   &
                 rc=rc)
!
!       Add array to export state
!     
        call ESMF_StateAdd (models(Iocean)%stateImport,                 &
                            (/ models(Iocean)%dataImport(i,n)%array /), &
                            rc=rc)
!
!       Initialize export field to zero to avoid infinities or NaNs.
!
        models(Iocean)%dataImport(i,n)%ptr = 0.0d0
      end do
      end do
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc=ESMF_SUCCESS
!
      end subroutine ROMS_SetStates
!
      subroutine ROMS_PutExportData (MyRank, rc)
!
!**********************************************************************
!
!     Used module declarations 
!
!**********************************************************************
!
      use mod_param, only : BOUNDS, N
      use mod_ocean, only : OCEAN
      use mod_scalars, only : itemp
      use mod_stepping, only : nstp
!
!***********************************************************************
!
!     Imported variable declarations 
!     
!***********************************************************************
!
      integer, intent(in) :: MyRank
      integer, intent(inout) :: rc
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
      character (len=80) :: name
      integer :: i, j, k, ng
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
      Istr =BOUNDS(ng)%Istr (MyRank)
      Iend =BOUNDS(ng)%Iend (MyRank)
      Jstr =BOUNDS(ng)%Jstr (MyRank)
      Jend =BOUNDS(ng)%Jend (MyRank)
      IstrR=BOUNDS(ng)%IstrR(MyRank)
      IendR=BOUNDS(ng)%IendR(MyRank)
      IstrU=BOUNDS(ng)%IstrU(MyRank)
      JstrR=BOUNDS(ng)%JstrR(MyRank)
      JendR=BOUNDS(ng)%JendR(MyRank)
      JstrV=BOUNDS(ng)%JstrV(MyRank)      
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
      rc=ESMF_SUCCESS
!
      end subroutine ROMS_PutExportData

      end module mod_esmf_ocn
