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
      use mod_runparams, only : dtcpl
!
!-----------------------------------------------------------------------
!     Local variable declarations  
!-----------------------------------------------------------------------
!
      real*8 :: dt
      integer :: iarr(6)
      integer :: i, j, localPet, petCount, comm, rc
      character (len=80) :: str1, str2
!
!***********************************************************************
!
!     Coupled Model Initialization
!
!***********************************************************************
!
!-----------------------------------------------------------------------
!     Initialize ESMF framework and get default VM
!-----------------------------------------------------------------------
!
      call ESMF_Initialize(vm=vm,                                       &
                           rc=rc)
!
!-----------------------------------------------------------------------
!     Get information from VM (MPI Communicator, number of PETs etc.)
!-----------------------------------------------------------------------
!
      call ESMF_VMGet(vm,                                               &
                      petCount=petCount,                                &
                      localPet=localPet,                                &
                      mpiCommunicator=comm,                             &
                      rc=rc)
!
!-----------------------------------------------------------------------
!     Allocate coupler variables
!-----------------------------------------------------------------------
!
      call allocate_cpl(petCount, rc)
!
!***********************************************************************
!
!     Create gridded components
!
!***********************************************************************
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
      end do
!
!-----------------------------------------------------------------------
!     Create coupler component. The coupler component run all PETs
!-----------------------------------------------------------------------
!
      str1 = 'Coupler Component'
      comp = ESMF_CplCompCreate(name=trim(str1), rc=rc)
!
!**********************************************************************
!
!     Register gridded components
!
!**********************************************************************
!
      call ESMF_GridCompSetServices(models(Iatmos)%comp,                &
                                    RCM_SetServices,                    &
                                    rc=rc)
      call ESMF_GridCompSetServices(models(Iocean)%comp,                &
                                    ROMS_SetServices,                   &
                                    rc=rc)
      call ESMF_CplCompSetServices(comp,                                &
                                   CPL_SetServices,                     &
                                   rc=rc)
!
!**********************************************************************
!
!     Create gridded components export/import state objects 
!
!**********************************************************************
!
      do i = 1, nModels
        if (i == Iatmos) then
          str1 = 'Import state (Atmosphere)'
          str2 = 'Export state (Atmosphere)'
        else if (i == Iocean) then
          str1 = 'Import state (Ocean)' 
          str2 = 'Export state (Ocean)'
        end if
        models(i)%stateExport = ESMF_StateCreate(name=trim(str2),       &
                                stateintent=ESMF_STATEINTENT_EXPORT,    &
                                rc=rc)
        models(i)%stateImport = ESMF_StateCreate(name=trim(str1),       &
                                stateintent=ESMF_STATEINTENT_IMPORT,    &
                                rc=rc)
      end do
!
!**********************************************************************
!
!     Initialize gridded components 
!
!**********************************************************************
!
      do i = 1, nModels
        call ESMF_GridCompInitialize(models(i)%comp,                    &
                                     importState=models(i)%stateImport, &
                                     exportState=models(i)%stateExport, &
                                     rc=rc)  
      end do
!
!**********************************************************************
!
!     Reconcile time (set start time, stop time and time step) 
!
!**********************************************************************
!
      call time_reconcile() 
!
!-----------------------------------------------------------------------
!     Initialize coupler component
!-----------------------------------------------------------------------
!
!     Forward coupling
!
      call ESMF_AttributeSet (models(Iatmos)%stateExport,               &
                              name=FORWARD_INIT,                        &
                              value=FORWARD_ON,                         &
                              rc=rc)
!
      call ESMF_AttributeSet (models(Iatmos)%stateExport,               &
                              name=BACKWARD_INIT,                       &
                              value=BACKWARD_OFF,                       &
                              rc=rc)
!
      call ESMF_CplCompInitialize(comp,                                 &
                                 importState=models(Iatmos)%stateExport,&
                                 exportState=models(Iocean)%stateImport,&
                                 rc=rc)
!      
      call ESMF_VMBarrier(vm, rc=rc)



! 
!**********************************************************************
!
!     Finalize gridded components 
!
!**********************************************************************
!
      do i = 1, nModels
        call ESMF_GridCompFinalize (models(i)%comp,                     &
                                    exportState=models(i)%stateExport,  &
                                    importState=models(i)%stateImport,  &
                                    rc=rc)
      end do 
! 
!**********************************************************************
!
!     Terminates all the ESMF/MPI processing 
!
!**********************************************************************
!
      call ESMF_Finalize(rc=rc)
!
      end program regcm
