!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

program regcm

!
!**********************************************************************
!
! Used module declarations 
!
!**********************************************************************
!
  use ESMF_Mod

  use mod_coupler
  use mod_esmf
!
!**********************************************************************
!
! Local variable declarations 
!
!**********************************************************************
!
  type(ESMF_VM) :: vm
  type(ESMF_GridComp) :: gcomp
  character(ESMF_MAXSTR), parameter :: conv = 'ESMF'
  character(ESMF_MAXSTR), parameter :: purp = 'General'
!
  type(ESMF_Time) :: time
  type(ESMF_TimeInterval) :: timeInt 
!
  integer :: status, rank, npet, comm
  character (len=80) :: name, importName, exportName
!
!**********************************************************************
!
! Model Initialization
!
!**********************************************************************
!
  call ESMF_Initialize(vm=vm,                                       &
     &                     rc=status)
!
  call ESMF_VMGet(vm,                                               &
     &                localPet=rank,                                    &
     &                petCount=npet,                                    &
     &                mpiCommunicator=comm,                             &
     &                rc=status)
!
  call allocate_cplvars(npet)
!
!**********************************************************************
!
! Create gridded components
!
!**********************************************************************
!
  do i = 1, nModels 
    if (i == 1) then
      name = 'atmos' 
    else if (i == 2) then
      name = 'ocean'
    end if
!
    models(i)%gcomp = ESMF_GridCompCreate(name=TRIM(name),          &
     &                                        petList=models(i)%petList,&
     &                                        rc=status)
  end do
!
!**********************************************************************
!
! Add attribute to gridded components and print out
!
!**********************************************************************
!
  do i = 1, nModels
!
!       Add attribute
!
    call ESMF_AttributeAdd(models(i)%gcomp,                         &
     &                         convention=conv,                         &
     &                         purpose=purp,                            &
     &                         rc=status)
!
!       Set attributes
!
    do j = 1, size(models(i)%attName)
      call ESMF_AttributeSet(models(i)%gcomp,                       &
     &                           name=models(i)%attName(j),             &
     &                           value=models(i)%attValue(j),           &
     &                           convention=conv,                       &
     &                           purpose=purp,                          &
     &                           rc=status)
    end do
!
!       Print attributes
!
    call ESMF_AttributeWrite(models(i)%gcomp,                       &
     &                           convention=conv,                       &
     &                           purpose=purp,                          &
     &                           attwriteflag=ESMF_ATTWRITE_XML,        &
     &                           rc=status)
  end do
!
!**********************************************************************
!
! Register gridded components
!
!**********************************************************************
!
  call ESMF_GridCompSetServices(models(Iatmos)%gcomp,              &
     &                              RCM_SetServices,                   &
     &                              status)
!
  call ESMF_GridCompSetServices(models(Iocean)%gcomp,              &
     &                              ROMS_SetServices,                  &
     &                              status)
!
!**********************************************************************
!
! Create gridded components export/import state objects 
!
!**********************************************************************
!
  do i = 1, nModels
    if (i == 1) then
      importName = 'Atmosphere gridded component import state'
      exportName = 'Atmosphere gridded component export state'
    else if (i == 2) then
      importName = 'Ocean gridded component import state'
      exportName = 'Ocean gridded component export state'
    end if
!
    models(i)%stateExport = ESMF_StateCreate(trim(exportName),     &
     &                                           ESMF_STATE_EXPORT,    &
     &                                           rc=status)
!
    models(i)%stateImport = ESMF_StateCreate(trim(importName),     &
     &                                           ESMF_STATE_IMPORT,    &
     &                                           rc=status)
  end do
!
!**********************************************************************
!
! Initialize gridded components 
!
!**********************************************************************
!
  do i = 1, nModels
    call ESMF_GridCompInitialize(models(i)%gcomp,                  &
     &                               importState=models(i)%stateImport,&
     &                               exportState=models(i)%stateExport,&
     &                               rc=status)  
!
!       Print time information (for debug purposes)
!
    if (rank == models(i)%petList(1)) then
      call ESMF_ClockGet(models(i)%clock,                          &
     &                       startTime=time,                           & 
     &                       timeStep=timeInt)
      call ESMF_TimePrint(time,                                    &
     &                        options="string",                        &
     &                        rc=status)
      call ESMF_TimeIntervalPrint(timeInt,                         &
     &                                options="string",                &
     &                                rc=status)
    end if
  end do
!
!**********************************************************************
!
! Set time interval to exchange data between gridded components 
!
!**********************************************************************
!
!      do i = 1, nModels 
!
!      end do
!
!      call ESMF_TimeIntervalSet(TimeStep(0),                            &
!     &                          s_r8=TimeInterval(Iocean,Iatmos),       &
!     &                          rc=status)
!
!
!      first = .TRUE.
!      timestr = 0.0D0
!      timeend = idatediff(idate2,idate1)*3600.0
!      call RCM_run(timestr, timeend, first)
!
!**********************************************************************
!
! Model Finalize 
!
!**********************************************************************
!
!       call RCM_finalize()
!

end program regcm
