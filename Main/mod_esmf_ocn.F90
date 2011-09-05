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
!     ROMS Component routines.
!
      use ocean_control_mod, only : ROMS_initialize
      use ocean_control_mod, only : ROMS_run
      use ocean_control_mod, only : ROMS_finalize
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
      public  :: ROMS_SetServices
      public  :: ROMS_SetRun
      public  :: ROMS_SetFinalize
!
      contains
!
      subroutine ROMS_SetServices(comp, rc)
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
                                      userRoutine=ROMS_SetInitialize,   &
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
                                      userRoutine=ROMS_SetRun,          &
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
                                      userRoutine=ROMS_SetFinalize,     &
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
      end subroutine ROMS_SetServices
!
      subroutine ROMS_SetInitialize(comp, importState, exportState,     &
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
      logical :: flag
      type(ESMF_Config) :: config
      integer :: MyRank, Nnodes, comm, ierr 
!
!***********************************************************************
!
!     Call ROMS initialization routines
!
!***********************************************************************
!
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!     
      call ESMF_GridCompGet(comp,                                       &
                            vm=models(Iocean)%vm,                       &
                            rc=rc)

      call ESMF_VMGet(models(Iocean)%vm,                                &
                      localPet=MyRank,                                  &
                      petCount=Nnodes,                                  &
                      mpiCommunicator=comm,                             &
                      rc=rc)
!
!-----------------------------------------------------------------------
!     Initialize the meshded component
!-----------------------------------------------------------------------
!
      flag = .TRUE.
      call MPI_Comm_dup(comm, models(Iocean)%comm, ierr)
      call ROMS_initialize(flag, MyCOMM=models(Iocean)%comm)
!
!-----------------------------------------------------------------------
!     Set-up ESMF internal clock for meshded component 
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
      call ROMS_PutGridData(MyRank, rc)
!
!-----------------------------------------------------------------------
!     Set-up import/export states
!-----------------------------------------------------------------------
!
      call ROMS_SetStates(ImportState, ExportState, MyRank, rc)
!
!-----------------------------------------------------------------------
!     Load export initial conditions data.
!-----------------------------------------------------------------------
!
      call ROMS_PutExportData(Myrank, rc)
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
      subroutine ROMS_SetClock(clock, status)
!
!**********************************************************************
!
!     Imported modules 
!
!**********************************************************************
!
      use mod_param
      use mod_scalars
!
!**********************************************************************
!
!     Imported variable declarations 
!
!**********************************************************************
!
      TYPE(ESMF_Clock), intent(inout) :: clock
      integer, intent(inout) :: status
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
      integer :: ref_year,   str_year,   end_year
      integer :: ref_month,  str_month,  end_month
      integer :: ref_day,    str_day,    end_day
      integer :: ref_hour,   str_hour,   end_hour
      integer :: ref_minute, str_minute, end_minute
      integer :: ref_second, str_second, end_second
!
      integer :: ng, MyTimeStep
      real(r8) :: MyStartTime, MyStopTime
      real(r8) :: hour, minute, yday
      character (len=80) :: name
!
!**********************************************************************
!
!     Create meshded component clock 
!
!**********************************************************************
!
!-----------------------------------------------------------------------
!     Create ESMF calendar
!-----------------------------------------------------------------------
!
      IF (INT(time_ref).eq.-2) THEN
        ref_year=1968
        ref_month=5
        ref_day=23
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='Modified Julian day number, Gregorian Calendar'
        models(Iocean)%cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,&
                                                 name=TRIM(name),       &
                                                 rc=status)
      ELSE IF (INT(time_ref).eq.-1) THEN
        ref_year=1
        ref_month=1
        ref_day=1
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='360-day, 30 days per month'
        models(Iocean)%cal = ESMF_CalendarCreate(ESMF_CALKIND_360DAY,   &
                                                 name=TRIM(name),       &
                                                 rc=status)
      ELSE IF (INT(time_ref).eq.0) THEN
        ref_year=1
        ref_month=1
        ref_day=1
        ref_hour=0
        ref_minute=0
        ref_second=0
        name='Julian Calendar, leap year if divisible by 4'
        models(Iocean)%cal = ESMF_CalendarCreate(ESMF_CALKIND_JULIAN,   &
                                                 name=TRIM(name),       &
                                                 rc=status)
      ELSE IF (time_ref.gt.0.0_r8) THEN
        ref_year=INT(r_date(2))
        ref_month=INT(r_date(4))
        ref_day=INT(r_date(5))
        ref_hour=INT(r_date(6))
        ref_minute=INT(r_date(7))
        ref_second=INT(r_date(8))
        name='Gregorian Calendar'
        models(Iocean)%cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN,&
                                                 name=TRIM(name),       &
                                                 rc=status)
      END IF
!
!-----------------------------------------------------------------------
!     Set Reference time.
!-----------------------------------------------------------------------
!
      CALL ESMF_TimeSet (models(Iocean)%refTime,                        &
                         yy=ref_year,                                   &
                         mm=ref_month,                                  &
                         dd=ref_day,                                    &
                         h=ref_hour,                                    &
                         m=ref_minute,                                  &
                         s=ref_second,                                  &
                         calendar=models(Iocean)%cal,                   &
                         rc=status)
!
!-----------------------------------------------------------------------
!     Set start time
!-----------------------------------------------------------------------
!
      MyStartTime = MINVAL(tdays)
      CALL caldate (r_date, MyStartTime, str_year, yday, str_month,     &
                    str_day, hour)
      minute=(hour-AINT(hour))*60.0_r8
      str_hour=INT(hour)
      str_minute=INT(minute)
      str_second=INT((minute-AINT(minute))*60.0_r8)
!
      CALL ESMF_TimeSet (models(Iocean)%strTime,                        &
                         yy=str_year,                                   &
                         mm=str_month,                                  &
                         dd=str_day,                                    &
                         h=str_hour,                                    &
                         m=str_minute,                                  &
                         s=str_second,                                  &
                         calendar=models(Iocean)%cal,                   &
                         rc=status)
!
!-----------------------------------------------------------------------
!     Set stop time
!-----------------------------------------------------------------------
!
      MyStopTime=0.0_r8
      DO ng = 1, nNest(Iocean)
        MyStopTime=MAX(MyStopTime,                                      &
                       tdays(ng)+(REAL(ntimes(ng),r8)*dt(ng))*sec2day)
      END DO
      CALL caldate (r_date, MyStopTime, end_year, yday, end_month,      &
                    end_day, hour)
      minute=(hour-AINT(hour))*60.0_r8
      end_hour=INT(hour)
      end_minute=INT(minute)
      end_second=INT((minute-AINT(minute))*60.0_r8)
!
      CALL ESMF_TimeSet (models(Iocean)%endTime,                        &
                         yy=end_year,                                   &
                         mm=end_month,                                  &
                         dd=end_day,                                    &
                         h=end_hour,                                    &
                         m=end_minute,                                  &
                         s=end_second,                                  &
                         calendar=models(Iocean)%cal,                   &
                         rc=status)
!
!-----------------------------------------------------------------------
!     Set time interval
!-----------------------------------------------------------------------
!
      MyTimeStep = INT(MINVAL(dt))
      call ESMF_TimeIntervalSet (models(Iocean)%dtsec,                  &
                                 s=MyTimeStep,                          &
                                 rc=status)  
!
!-----------------------------------------------------------------------
!     Create time clock.
!-----------------------------------------------------------------------
!
      name='ROMS model time clock'
      models(Iocean)%clock = ESMF_ClockCreate (name=TRIM(name),         &
                                  refTime=models(Iocean)%refTime,       &
                                  timeStep=models(Iocean)%dtsec,        &
                                  startTime=models(Iocean)%strTime,     &
                                  stopTime=models(Iocean)%endTime,      &
                                  rc=status)
!
!-----------------------------------------------------------------------
!     Copy clock
!-----------------------------------------------------------------------
!
      clock = ESMF_ClockCreate (models(Iocean)%clock, rc=status)
!
!-----------------------------------------------------------------------
!     Validate time clock
!-----------------------------------------------------------------------
!
      call ESMF_ClockValidate (models(Iocean)%clock,                    &
                               rc=status)
!
!-----------------------------------------------------------------------
!     Get meshded component internal clock current time
!-----------------------------------------------------------------------
!
      call ESMF_ClockGet (models(Iocean)%clock,                         &
                          currTime=models(Iocean)%curTime,              &
                          rc=status)
!
!-----------------------------------------------------------------------
!     Put current time into ESM_Time variable
!-----------------------------------------------------------------------
!
      CALL ESMF_TimeGet (models(Iocean)%curTime,                        &
                         yy=models(Iocean)%time%year,                   &
                         mm=models(Iocean)%time%month,                  &
                         dd=models(Iocean)%time%day,                    &
                         h=models(Iocean)%time%hour,                    &
                         m=models(Iocean)%time%minute,                  &
                         s=models(Iocean)%time%second,                  &
                         timeZone=models(Iocean)%time%zone,             &
                         timeStringISOFrac=models(Iocean)%time%stamp,   &
                         dayOfYear=models(Iocean)%time%yday,            &
                         rc=status)
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
                             rc=status)
      name = 'stop time'
      call ESMF_AttributeSet(models(Iocean)%stateExport,                &
                             name=trim(name),                           &
                             valueList=(/ end_year  , end_month,        &
                                          end_day   , end_hour ,        &
                                          end_minute, end_second /),    &
                             rc=status)
      name = 'time step'
      call ESMF_AttributeSet(models(Iocean)%stateExport,                &
                             name=trim(name),                           &
                             value=MyTimeStep,                          &
                             rc=status)
!
!-----------------------------------------------------------------------
!  Set return flag to success.
!-----------------------------------------------------------------------
!
      status = ESMF_SUCCESS
!
      end subroutine ROMS_SetClock
!
      subroutine ROMS_SetGridArrays (gcomp, status)
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
      type(ESMF_GridComp), intent(inout) :: gcomp
      integer, intent(out) :: status
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
      integer :: i, n, tile 
      integer :: myRank, nNodes, myComm
      integer, allocatable :: deBlockList(:,:,:) 
      integer, allocatable :: TLWidth(:,:), TUWidth(:,:)
      integer, allocatable :: CLW_r(:,:), CUW_r(:,:)
      integer, allocatable :: CLW_u(:,:), CUW_u(:,:)
      integer, allocatable :: CLW_v(:,:), CUW_v(:,:)
      integer, dimension(2) :: CLW, CUW, TLW, TUW
      integer, dimension(2) :: deCount, minIndex, maxIndex
      TYPE (ESMF_ARRAY) :: grdArray
!
!-----------------------------------------------------------------------
!     Query Virtual Machine (VM) environment for the MPI
!     communicator handle     
!-----------------------------------------------------------------------
!  
      call ESMF_GridCompGet (gcomp,                                     &
                             vm=models(Iocean)%vm,                      &
                             rc=status)

      call ESMF_VMGet (models(Iocean)%vm,                               &
                       localPet=myRank,                                 &
                       petCount=nNodes,                                 &
                       mpiCommunicator=myComm,                          &
                       rc=status)
!
!-----------------------------------------------------------------------
!     Set RCM domain decomposition variables
!-----------------------------------------------------------------------
!
!     Loop over number of nested/composed meshs.
      do n = 1, nNest(Iocean)
!
!-----------------------------------------------------------------------
!     Set ESMF Layout and distribution objects
!-----------------------------------------------------------------------
!
      deCount = (/ NtileI(n), NtileJ(n) /)
!
!      models(Iocean)%deLayout(n) = ESMF_DELayoutCreate(models(Iocean)%vm,  &
!                                         deCountList=deCount,           &
!                                         petList=models(Iocean)%petList,&
!                                         rc=status)
!      call ESMF_DELayoutPrint(models(Iocean)%deLayout(n))
!      print*, "** turuncu ** PET LIST =", models(Iocean)%petList
!
!-----------------------------------------------------------------------
!     Create ESMF DistGrid based on model domain decomposition
!-----------------------------------------------------------------------
!
      if (.not.allocated(deBlockList)) then
        allocate (deBlockList(2, 2, NtileI(n)*NtileJ(n)))
      end if 
!
      do tile = 0, NtileI(n)*NtileJ(n)-1
          deBlockList(1,1,tile+1) = BOUNDS(n)%Istr(tile)
          deBlockList(1,2,tile+1) = BOUNDS(n)%Iend(tile)
          deBlockList(2,1,tile+1) = BOUNDS(n)%Jstr(tile)
          deBlockList(2,2,tile+1) = BOUNDS(n)%Jend(tile)
      end do 
!
!     Cordinates of the lower and upper corner of the patch
!
      minIndex = (/ 1, 1 /)
      maxIndex = (/ Lm(n), Mm(n) /)
!
!      models(Iocean)%distGrid(n) = ESMF_DistGridCreate(minCorner,          &
!                                     maxCorner,                         &
!                                     deBlockList=deBlockList,           &
!                                     deLayout=models(Iocean)%deLayout(n),  &
!                                     vm=models(Iocean)%vm,              &
!                                     rc=status)
      models(Iocean)%distGrid(n) = ESMF_DistGridCreate (minIndex=minIndex,&
                                                        maxIndex=maxIndex,&
                                                        regDecomp=(/NtileI(n), NtileJ(n)/),&
                                                        rc=status)
!
!-----------------------------------------------------------------------
!     Set array descriptor
!-----------------------------------------------------------------------
!
      call ESMF_ArraySpecSet(models(Iocean)%arrSpec(n),                    &
                             rank=2,                                    &
                             typekind=ESMF_TYPEKIND_R8,                 &
                             rc=status)
!
!-----------------------------------------------------------------------
!     Define computational and total (memory) regions
!-----------------------------------------------------------------------
!
      if (.not. allocated(TLWidth)) then
        allocate(TLWidth(2,0:NtileI(n)*NtileJ(n)-1))
        allocate(TUWidth(2,0:NtileI(n)*NtileJ(n)-1))
        allocate(CLW_r(2,0:NtileI(n)*NtileJ(n)-1))
        allocate(CUW_r(2,0:NtileI(n)*NtileJ(n)-1))
        allocate(CLW_u(2,0:NtileI(n)*NtileJ(n)-1))
        allocate(CUW_u(2,0:NtileI(n)*NtileJ(n)-1))
        allocate(CLW_v(2,0:NtileI(n)*NtileJ(n)-1))
        allocate(CUW_v(2,0:NtileI(n)*NtileJ(n)-1))
      end if 
!
      do tile = 0, NtileI(n)*NtileJ(n)-1
        TLWidth(1,tile) = BOUNDS(n)%Istr(tile)-BOUNDS(n)%LBi(tile)
        TLWidth(2,tile) = BOUNDS(n)%Jstr(tile)-BOUNDS(n)%LBj(tile)
        TUWidth(1,tile) = BOUNDS(n)%UBi(tile)-BOUNDS(n)%Iend(tile)
        TUWidth(2,tile) = BOUNDS(n)%UBj(tile)-BOUNDS(n)%Jend(tile)
!
        CLW_r(1,tile) = BOUNDS(n)%Istr(tile)-BOUNDS(n)%IstrR(tile)
        CLW_r(2,tile) = BOUNDS(n)%Jstr(tile)-BOUNDS(n)%JstrR(tile)
        CUW_r(1,tile) = BOUNDS(n)%IendR(tile)-BOUNDS(n)%Iend(tile)
        CUW_r(2,tile) = BOUNDS(n)%JendR(tile)-BOUNDS(n)%Jend(tile)
!
        CLW_u(1,tile) = 0
        CLW_u(2,tile) = CLW_r(2,tile)
        CUW_u(1,tile) = CUW_r(1,tile)
        CUW_u(2,tile) = CUW_r(2,tile)
!
        CLW_v(1,tile) = CLW_r(1,tile)
        CLW_v(2,tile) = 0
        CUW_v(1,tile) = CUW_r(1,tile)
        CUW_v(2,tile) = CUW_r(2,tile)
      end do
!
      TLW = (/ TLWidth(1,MyRank), TLWidth(2,MyRank) /)
      TUW = (/ TUWidth(1,MyRank), TUWidth(2,MyRank) /)
!
!-----------------------------------------------------------------------
!     Print out the information
!-----------------------------------------------------------------------
!
      if (myRank == 0) then
        print *, ' '
        print *, 'Horizontal decomposition indices per tile:'
        print *, ' '
        print 10, 'Istr   = ',(BOUNDS(n)%Istr(tile),tile=0,Nnodes-1)
        print 10, 'IstrU  = ',(BOUNDS(n)%IstrU(tile),tile=0,Nnodes-1)
        print 10, 'Iend   = ',(BOUNDS(n)%Iend(tile),tile=0,Nnodes-1)
        print 10, 'Jstr   = ',(BOUNDS(n)%Jstr(tile),tile=0,Nnodes-1)
        print 10, 'JstrV  = ',(BOUNDS(n)%JstrV(tile),tile=0,Nnodes-1)
        print 10, 'Jend   = ',(BOUNDS(n)%Jend(tile),tile=0,Nnodes-1)
        print *, ' '
        print 10, 'LBi    = ',(BOUNDS(n)%LBi(tile),tile=0,Nnodes-1)
        print 10, 'UBi    = ',(BOUNDS(n)%UBi(tile),tile=0,Nnodes-1)
        print 10, 'LBj    = ',(BOUNDS(n)%LBj(tile),tile=0,Nnodes-1)
        print 10, 'UBj    = ',(BOUNDS(n)%UBj(tile),tile=0,Nnodes-1)
        print *, ' '
        print 10, 'TLWi   = ',(TLWidth(1,tile),tile=0,Nnodes-1)
        print 10, 'TLWj   = ',(TLWidth(2,tile),tile=0,Nnodes-1)
        print 10, 'TUWi   = ',(TUWidth(1,tile),tile=0,Nnodes-1)
        print 10, 'TUWj   = ',(TUWidth(2,tile),tile=0,Nnodes-1)
        print *, ' '
        print 10, 'CLWi_r = ',(CLW_r(1,tile),tile=0,Nnodes-1)
        print 10, 'CLWj_r = ',(CLW_r(2,tile),tile=0,Nnodes-1)
        print 10, 'CUWi_r = ',(CUW_r(1,tile),tile=0,Nnodes-1)
        print 10, 'CUWj_r = ',(CUW_r(2,tile),tile=0,Nnodes-1)
        print *, ' '
        print 10, 'CLWi_u = ',(CLW_u(1,tile),tile=0,Nnodes-1)
        print 10, 'CLWj_u = ',(CLW_u(2,tile),tile=0,Nnodes-1)
        print 10, 'CUWi_u = ',(CUW_u(1,tile),tile=0,Nnodes-1)
        print 10, 'CUWj_u = ',(CUW_u(2,tile),tile=0,Nnodes-1)
        print *, ' '
        print 10, 'CLWi_v = ',(CLW_u(1,tile),tile=0,Nnodes-1)
        print 10, 'CLWj_v = ',(CLW_u(2,tile),tile=0,Nnodes-1)
        print 10, 'CUWi_v = ',(CUW_u(1,tile),tile=0,Nnodes-1)
        print 10, 'CUWj_v = ',(CUW_u(2,tile),tile=0,Nnodes-1)
 10     format(1x,a,64i5)
      end if
!
!
!-----------------------------------------------------------------------    
!     Create exchange arrays
!-----------------------------------------------------------------------    
!
      do i = 1, ubound(models(Iocean)%mesh, dim=1)
!
        if (models(Iocean)%mesh(i,n)%gtype == Iupoint) then
          CLW = (/ CLW_u(1,myRank), CLW_u(2,myRank) /)
          CUW = (/ CUW_u(1,myRank), CUW_u(2,myRank) /)
        else if (models(Iocean)%mesh(i,n)%gtype == Ivpoint) THEN
          CLW = (/ CLW_v(1,myRank), CLW_v(2,myRank) /)
          CUW = (/ CUW_v(1,myRank), CUW_v(2,myRank) /)
        else
          CLW = (/ CLW_r(1,myRank), CLW_r(2,myRank)/)
          CUW = (/ CUW_r(1,myRank), CUW_r(2,myRank)/)
        end if
!
        grdArray = ESMF_ArrayCreate (           &
                                  arrayspec=models(Iocean)%arrSpec(n),     &
                                  distgrid=models(Iocean)%distGrid(n),     &
                                  computationalLWidth=CLW,              &
                                  computationalUWidth=CUW,              &
                                  totalLWidth=TLW,                      &
                                  totalUWidth=TUW,                      &
                                  indexflag=ESMF_INDEX_GLOBAL,          &
                                  rc=status)
        models(Iocean)%mesh(i,n)%lat%array = grdArray
!
!       Get array pointer
!
        call ESMF_ArrayGet(models(Iocean)%mesh(i,n)%lat%array,            &
                 farrayPtr=models(Iocean)%mesh(i,n)%lat%field,       &
                 rc=status)
!
        grdArray = ESMF_ArrayCreate (           &
                                  arrayspec=models(Iocean)%arrSpec(n),     &
                                  distgrid=models(Iocean)%distGrid(n),     &
                                  computationalLWidth=CLW,              &
                                  computationalUWidth=CUW,              &
                                  totalLWidth=TLW,                      &
                                  totalUWidth=TUW,                      &
                                  indexflag=ESMF_INDEX_GLOBAL,          &
                                  rc=status)
        models(Iocean)%mesh(i,n)%lon%array = grdArray
!
!       Get array pointer
!
        call ESMF_ArrayGet(models(Iocean)%mesh(i,n)%lon%array,            &
                 farrayPtr=models(Iocean)%mesh(i,n)%lon%field,       &
                 rc=status) 
!
!       Set adjustable settings of array object  
!
        call ESMF_ArraySet(array=models(Iocean)%mesh(i,n)%lat%array,      &
                 name=trim(models(Iocean)%mesh(i,n)%lat%long_name),  &
                 rc=status)
        call ESMF_ArraySet(array=models(Iocean)%mesh(i,n)%lon%array,      &
                 name=trim(models(Iocean)%mesh(i,n)%lon%long_name),  &
                 rc=status)
!  
!       Add array to export state
!        
        call ESMF_StateAdd(models(Iocean)%stateExport,                  &
                           (/ models(Iocean)%mesh(i,n)%lat%array /),      &
                           rc=status)
        call ESMF_StateAdd(models(Iocean)%stateExport,                  &
                           (/ models(Iocean)%mesh(i,n)%lon%array /),      &
                           rc=status)
!
!       Initialize the mesh array
!
        models(Iocean)%mesh(i,n)%lat%field = 0.0d0
        models(Iocean)%mesh(i,n)%lon%field = 0.0d0
!        models(Iocean)%mesh(i,n)%mask%field = 0.0d0
      end do
      end do
!
!     Deallocate temporary arrays
!
      deallocate (TLWidth)
      deallocate (TUWidth)
      deallocate (CLW_r)
      deallocate (CUW_r)
      deallocate (CLW_u)
      deallocate (CUW_u)
      deallocate (CLW_v)
      deallocate (CUW_v)
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      status = ESMF_SUCCESS
!
      end subroutine ROMS_SetGridArrays
!
      subroutine ROMS_PutGridData (localPet, status)
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
      integer, intent(inout) :: status
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
              models(Iocean)%mesh(k,n)%lon%field(i,j)=GRID(n)%lonr(i,j)
              models(Iocean)%mesh(k,n)%lat%field(i,j)=GRID(n)%latr(i,j)
!              models(Iocean)%mesh(k,n)%mask%field(i,j)=GRID(n)%rmask(i,j)
            end do
          end do
!       U
        else if (models(Iocean)%mesh(k,n)%gtype == Iupoint) then
          do j = JstrR, JendR
            do i = Istr, IendR
              models(Iocean)%mesh(k,n)%lon%field(i,j)=GRID(n)%lonu(i,j)
              models(Iocean)%mesh(k,n)%lat%field(i,j)=GRID(n)%latu(i,j)
!              models(Iocean)%mesh(k,n)%mask%field(i,j)=GRID(n)%umask(i,j)
            end do
          end do
!       V
        else if (models(Iocean)%mesh(k,n)%gtype == Ivpoint) then  
          do j = Jstr, JendR
            do i = IstrR, IendR
              models(Iocean)%mesh(k,n)%lon%field(i,j)=GRID(n)%lonv(i,j)
              models(Iocean)%mesh(k,n)%lat%field(i,j)=GRID(n)%latv(i,j)
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
      status = ESMF_SUCCESS
!
      end subroutine ROMS_PutGridData
!
      subroutine ROMS_SetStates (importState, exportState, MyRank,      &
                                 status)
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
      integer, intent(inout) :: status
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
                                     rc=status)
!
!       Get data pointer from array
!
        call ESMF_ArrayGet (models(Iocean)%dataExport(i,n)%array,       &
                 farrayPtr=models(Iocean)%dataExport(i,n)%field,        &
                 rc=status) 
!
!       Set array name
!
        call ESMF_ArraySet (array=models(Iocean)%dataExport(i,n)%array, &
                 name=trim(models(Iocean)%dataExport(i,n)%long_name),   &
                 rc=status)
!
!       Add array to export state
!     
        call ESMF_StateAdd (models(Iocean)%stateExport,                 &
                            (/ models(Iocean)%dataExport(i,n)%array /), &
                            rc=status)
!
!       Initialize export field to zero to avoid infinities or NaNs.
!
        models(Iocean)%dataExport(i,n)%field = 0.0d0
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
                                     rc=status)
!
!       Get data pointer from array
!
        call ESMF_ArrayGet (models(Iocean)%dataImport(i,n)%array,       &
                 farrayPtr=models(Iocean)%dataImport(i,n)%field,        &
                 rc=status)
!
!       Set array name
!
        call ESMF_ArraySet (array=models(Iocean)%dataImport(i,n)%array, &
                 name=trim(models(Iocean)%dataImport(i,n)%long_name),   &
                 rc=status)
!
!       Add array to export state
!     
        call ESMF_StateAdd (models(Iocean)%stateImport,                 &
                            (/ models(Iocean)%dataImport(i,n)%array /), &
                            rc=status)
!
!       Initialize export field to zero to avoid infinities or NaNs.
!
        models(Iocean)%dataImport(i,n)%field = 0.0d0
      end do
      end do
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      status=ESMF_SUCCESS
!
      end subroutine ROMS_SetStates
!
      subroutine ROMS_PutExportData (MyRank, status)
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
      integer, intent(inout) :: status
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
              models(Iocean)%dataExport(k,ng)%field =                    &
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
      status=ESMF_SUCCESS
!
      end subroutine ROMS_PutExportData

      end module mod_esmf_ocn
