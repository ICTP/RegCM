!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      program regcm
!
!-----------------------------------------------------------------------
!     Used module declarations  
!-----------------------------------------------------------------------
!
      use ESMF
!
      use mod_couplerr
      use mod_esmf_atm
      use mod_esmf_ocn
      use mod_esmf_cpl
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations  
!-----------------------------------------------------------------------
!
      logical :: first
      integer :: i, j, localPet, petCount, comm, rc
      character (len=80) :: str1, str2
!
!-----------------------------------------------------------------------
!     Coupled Model Initialization
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Initialize ESMF framework and get default VM
!-----------------------------------------------------------------------
!
      call ESMF_Initialize(vm=cplVM, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get information from VM (MPI Communicator, number of PETs etc.)
!-----------------------------------------------------------------------
!
      call ESMF_VMGet(cplVM,                                            &
                      petCount=petCount,                                &
                      localPet=localPet,                                &
                      mpiCommunicator=comm,                             &
                      rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Allocate coupler variables
!-----------------------------------------------------------------------
!
      call allocate_cpl(petCount, rc)
!
!-----------------------------------------------------------------------
!     Create gridded components
!-----------------------------------------------------------------------
!
      do i = 1, nModels 
        if (i == Iatmos) then
          str1 = 'Gridded Component I - Atmosphere' 
        else if (i == Iocean) then
          str1 = 'Gridded Component II - Ocean' 
        end if
        models(i)%comp = ESMF_GridCompCreate(name=trim(str1),           &
                                        petList=models(i)%petList,      &
                                        rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do
!
!-----------------------------------------------------------------------
!     Create coupler component. The coupler component run all PETs
!-----------------------------------------------------------------------
!
      str1 = 'Coupler Component'
      cplComp = ESMF_CplCompCreate(name=trim(str1), rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Register gridded components
!-----------------------------------------------------------------------
!
      call ESMF_GridCompSetServices(models(Iatmos)%comp,                &
                                    RCM_SetServices,                    &
                                    rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_GridCompSetServices(models(Iocean)%comp,                &
                                    ROMS_SetServices,                   &
                                    rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompSetServices(cplComp,                             &
                                   CPL_SetServices,                     &
                                   rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Create gridded components export/import state objects 
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        if (i == Iatmos) then
          str1 = 'Import state (Atmosphere)'
          str2 = 'Export state (Atmosphere)'
        else if (i == Iocean) then
          str1 = 'Import state (Ocean)' 
          str2 = 'Export state (Ocean)'
        end if
!
        models(i)%stateExport = ESMF_StateCreate(name=trim(str2),       &
                                stateintent=ESMF_STATEINTENT_EXPORT,    &
                                rc=rc)
!
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
!
        models(i)%stateImport = ESMF_StateCreate(name=trim(str1),       &
                                stateintent=ESMF_STATEINTENT_IMPORT,    &
                                rc=rc)
        if (rc /= ESMF_SUCCESS) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if
      end do
!
      do i = 1, nModels
!
!-----------------------------------------------------------------------
!     Initialize gridded components 
!-----------------------------------------------------------------------
!
      call ESMF_GridCompInitialize(models(i)%comp,                      &
                                   importState=models(i)%stateImport,   &
                                   exportState=models(i)%stateExport,   &
                                   rc=rc) 
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Reconcile time (set start time, stop time and time step) and
!     coupled model parameters. It is called in the loop because we 
!     need to use coupled model parameters in ocean component
!-----------------------------------------------------------------------
!
      first = .false.
      if (i == Iatmos) first = .true.
      call time_reconcile(first) 
      end do
!
!-----------------------------------------------------------------------
!     Initialize coupler component
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Forward coupling
!-----------------------------------------------------------------------
!
      call ESMF_AttributeSet (models(Iatmos)%stateExport,               &
                              name=trim(DIRECTION),                     &
                              value=FORWARD_ON,                         &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompInitialize(cplComp,                              &
                                 importState=models(Iatmos)%stateExport,&
                                 exportState=models(Iocean)%stateImport,&
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!      
      call ESMF_VMBarrier(cplVM, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Backward coupling
!-----------------------------------------------------------------------
!
      call ESMF_AttributeSet (models(Iocean)%stateExport,               &
                              name=trim(DIRECTION),                     &
                              value=FORWARD_OFF,                        &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_CplCompInitialize(cplComp,                              &
                                 importState=models(Iocean)%stateExport,&
                                 exportState=models(Iatmos)%stateImport,&
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!      
      call ESMF_VMBarrier(cplVM, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Run components
!-----------------------------------------------------------------------
!
      do while (.not. ESMF_ClockIsStopTime(cplClock))
!
!-----------------------------------------------------------------------
!     Run gridded components
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        call ESMF_GridCompRun (models(i)%comp,                          &
                               importState=models(i)%stateImport,       &
                               exportState=models(i)%stateExport,       &
                               clock=cplClock,                          &
                               rc=rc)
!
!-----------------------------------------------------------------------
!     Set coupling direction 
!-----------------------------------------------------------------------
!
      if (i == Iatmos) then
      call ESMF_AttributeSet (models(Iocean)%stateExport,               &
                              name=trim(DIRECTION),                     &
                              value=FORWARD_OFF,                        &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      j = Iocean
      else
      call ESMF_AttributeSet (models(Iatmos)%stateExport,               &
                              name=trim(DIRECTION),                     &
                              value=FORWARD_ON,                         &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      j = Iatmos
      end if
!
!-----------------------------------------------------------------------
!     Run coupler component
!-----------------------------------------------------------------------
!
      call ESMF_CplCompRun(cplComp,                                     &
                           importState=models(i)%stateExport,           &
                           exportState=models(j)%stateImport,           &
                           rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do
!
!-----------------------------------------------------------------------
!     Update coupler component clock 
!-----------------------------------------------------------------------
!
      call ESMF_ClockAdvance(cplClock, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Debug: print clock 
!-----------------------------------------------------------------------
!
      call ESMF_ClockPrint (cplClock, "currTime string", rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      end do
!
      call ESMF_VMBarrier(cplVM, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
! 
!-----------------------------------------------------------------------
!     Finalize gridded components 
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        call ESMF_GridCompFinalize (models(i)%comp,                     &
                                    exportState=models(i)%stateExport,  &
                                    importState=models(i)%stateImport,  &
                                    rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do 
!
!-----------------------------------------------------------------------
!     Finalize coupler components.
!-----------------------------------------------------------------------
!
      call ESMF_CplCompFinalize (cplComp,                               &
                                 exportState=models(Iatmos)%stateExport,&
                                 importState=models(Iocean)%stateImport,&
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
! 
!-----------------------------------------------------------------------
!     Terminates all the ESMF/MPI processing 
!-----------------------------------------------------------------------
!
      call ESMF_Finalize(rc=rc)
!
      end program regcm
