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
!**********************************************************************
!
!     Used module declarations 
!
!**********************************************************************
!
      use ESMF
!
      use mod_couplerr
      use mod_esmf_atm
      use mod_esmf_ocn
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
      type(ESMF_VM) :: vm
      type(ESMF_GridComp) :: gcomp
!
      integer :: i, j, rc, localPet, petCount, comm
      character (len=80) :: name, istr, estr
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
          name = 'Gridded Component I - Atmosphere' 
        else if (i == Iocean) then
          name = 'Gridded Component II - Ocean' 
        end if
        models(i)%comp = ESMF_GridCompCreate(name=TRIM(name),           &
                                        petList=models(i)%petList,      &
                                        rc=rc)
      end do
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
!
!**********************************************************************
!
!     Create gridded components export/import state objects 
!
!**********************************************************************
!
      do i = 1, nModels
        if (i == Iatmos) then
          istr = 'Import state (Atmosphere)'
          estr = 'Export state (Atmosphere)'
        else if (i == Iocean) then
          istr = 'Import state (Ocean)' 
          estr = 'Export state (Ocean)'
        end if
        models(i)%stateExport = ESMF_StateCreate(name=trim(estr),       &
                                stateintent=ESMF_STATEINTENT_EXPORT,    &
                                rc=rc)
        models(i)%stateImport = ESMF_StateCreate(name=trim(istr),       &
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
!     Set start time, stop time and time step 
!
!**********************************************************************
!
      do i = 1, Nmodels
        print*, "--", i, localPet, trim(models(i)%time%stamp)
      end do
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
      end program regcm
