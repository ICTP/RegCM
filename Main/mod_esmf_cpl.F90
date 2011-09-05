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
!***********************************************************************
!
!     Used module declarations 
!
!***********************************************************************
!
      use ESMF
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
      public  :: CPL_SetServices
      public  :: CPL_SetRun
      public  :: CPL_SetFinalize
!
      contains
!
      subroutine CPL_SetServices(comp, rc)
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_CplComp), intent(inout) :: comp
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
      call ESMF_CplCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_INITIALIZE,&
                                      userRoutine=CPL_SetInitialize,    &
                                      rc=rc)
!
!***********************************************************************
!
!     Register "run" routine
!
!***********************************************************************
!
      call ESMF_CplCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_RUN,       &
                                      userRoutine=CPL_SetRun,           &
                                      rc=rc)
!
!***********************************************************************
!
!     Register "finalize" routine
!
!***********************************************************************
!
      call ESMF_CplCompSetEntryPoint(comp,                             &
                                      methodflag=ESMF_METHOD_FINALIZE,  &
                                      userRoutine=CPL_SetFinalize,      &
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
      end subroutine CPL_SetServices
!
      subroutine CPL_SetInitialize(comp, importState, exportState,     &
                                   clock, rc)
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_CplComp), intent(inout) :: comp
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
      integer :: MyRank, Nnodes, comm, ierr, i 
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
      call ESMF_CplCompGet(comp,                                        &
                           vm=vm,                                       &
                           rc=rc)

      call ESMF_VMGet(vm,                                               &
                      localPet=MyRank,                                  &
                      petCount=Nnodes,                                  &
                      mpiCommunicator=comm,                             &
                      rc=rc)
!
!-----------------------------------------------------------------------
!     Reconcile import and export states. Consistent view in all PETs.
!-----------------------------------------------------------------------
!
!     Import state
!
      do i = 1, ubound(models, dim=1)
        call ESMF_StateReconcile(models(i)%stateImport,                 &
                                 vm=vm,                                 &
                                 attreconflag=ESMF_ATTRECONCILE_ON,     &
                                 rc=rc)
      end do
! 
!     Export state
!
      do i = 1, ubound(models, dim=1)
        call ESMF_StateReconcile(models(i)%stateExport,                 &
                                 vm=vm,                                 &
                                 attreconflag=ESMF_ATTRECONCILE_ON,     &
                                 rc=rc)
      end do
!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine CPL_SetInitialize
!
      subroutine CPL_SetRun(comp, importState, exportState,            &
                             clock, rc)
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      type(ESMF_CplComp), intent(inout) :: comp
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
!     Call CPL initialization routines
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
      end subroutine CPL_SetRun
!
      subroutine CPL_SetFinalize(comp, importState, exportState,       &
                                  clock, rc)
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
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
!
      end subroutine CPL_SetFinalize
!
      end module mod_esmf_cpl 
