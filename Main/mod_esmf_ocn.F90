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
!     Initialize the gridded component
!-----------------------------------------------------------------------
!
      flag = .TRUE.
      call MPI_Comm_dup(comm, models(Iocean)%comm, ierr)
      call ROMS_initialize(flag, MyCOMM=models(Iocean)%comm)
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
      end module mod_esmf_ocn
