      module mod_couplerr
!
!***********************************************************************
!
!     Imported modules 
!
!***********************************************************************
!
      use ESMF
!
      implicit none
!
!***********************************************************************
!
!     User defined data types for coupling interface 
!
!***********************************************************************
!
!-----------------------------------------------------------------------
!     Earth System Model (ESM) time data type 
!-----------------------------------------------------------------------
!
      type ESM_Time
        integer :: year                 
        integer :: month                
        integer :: day                  
        integer :: hour                 
        integer :: minute               
        integer :: second               
        integer :: yday              
        integer :: zone             
        character (len=30) :: stamp
      end type ESM_Time
!
!-----------------------------------------------------------------------
!     Earth System Model (ESM) field data type 
!-----------------------------------------------------------------------
!
      type ESM_Field
        integer :: fid
        integer :: gtype        
        character (len=40) :: name
        character (len=80) :: long_name
        character (len=80) :: units
        type(ESMF_Array) :: array
        real*8, dimension(:,:), pointer :: field
      end type ESM_Field
!
!-----------------------------------------------------------------------
!     Earth System Model (ESM) mesh data type 
!-----------------------------------------------------------------------
!
      type ESM_Mesh
        integer :: gid
        integer :: gtype
        type(ESM_Field) :: lat
        type(ESM_Field) :: lon
        type(ESM_Field) :: mask
      end type ESM_Mesh
!
!-----------------------------------------------------------------------
!     Earth System Model (ESM) high-level generic data type 
!-----------------------------------------------------------------------
!
      type ESM_Model
        integer :: mid 
        integer :: comm
        type(ESMF_VM) :: vm 
        integer, allocatable :: petList(:) 
        type(ESMF_GridComp) :: comp
        type(ESM_Mesh), allocatable :: mesh(:,:)        
        type(ESMF_DELayout), allocatable :: deLayout(:) 
        type(ESMF_DistGrid), allocatable :: distGrid(:)
        type(ESMF_ArraySpec), allocatable :: arrSpec(:)
        type(ESM_Field), allocatable :: dataExport(:,:)
        type(ESM_Field), allocatable :: dataImport(:,:)
        type(ESMF_State) :: stateExport
        type(ESMF_State) :: stateImport
        type(ESMF_Calendar) :: cal
        type(ESMF_Time) :: refTime
        type(ESMF_Time) :: strTime
        type(ESMF_Time) :: endTime
        type(ESMF_Time) :: curTime
        type(ESMF_TimeInterval) :: dtsec
        type(ESMF_Clock) :: clock
        type(ESM_Time) :: time 
      end type ESM_Model
!
!***********************************************************************
!
!     Global module variables 
!
!***********************************************************************
!
      type(ESM_Model), save, allocatable :: models(:) 
!
!-----------------------------------------------------------------------
!     Number of gridded component or model (Atmosphere/Ocean)
!-----------------------------------------------------------------------
!
      integer, parameter :: nModels = 2
!
!-----------------------------------------------------------------------
!     Number of nested grid in each model 
!-----------------------------------------------------------------------
!
      integer :: nNest(nModels)
!
!-----------------------------------------------------------------------
!     Gridded model indices
!-----------------------------------------------------------------------
!
      integer :: Iatmos  = 1
      integer :: Iocean  = 2
!
!-----------------------------------------------------------------------
!     Staggered grid point indices
!     d --------- d   ----- v -----  
!     |           |   |           |
!     |     c     |   u     c     u
!     |           |   |           |
!     d --------- d   ----- v -----     
!     Arakawa - B     Arakawa - C
!     RegCM           ROMS
!-----------------------------------------------------------------------
!
      integer :: Icross  = 1
      integer :: Idot    = 2
      integer :: Iupoint = 3
      integer :: Ivpoint = 4
!
!-----------------------------------------------------------------------
!     Coupler component variables 
!-----------------------------------------------------------------------
!
      type(ESMF_CplComp) :: comp
      type(ESMF_VM) :: vm
      type(ESMF_Time) :: startTime, stopTime
      type(ESMF_TimeInterval) :: timeStep
      type(ESMF_Clock) :: clock
!
      contains
!
      subroutine allocate_cpl(petCount, rc)
      implicit none
!
!**********************************************************************
!
!     Imported variable declarations 
!
!**********************************************************************
!
      integer, intent(in) :: petCount
      integer, intent(inout) :: rc
!
!**********************************************************************
!
!     Local variable declarations 
!
!**********************************************************************
!
      integer :: i, j, k, nPets
!
!-----------------------------------------------------------------------
!     Initialize the coupler variables
!-----------------------------------------------------------------------
!
      ! number of nested grid for each gridded component
      ! the nested grid is not supported in this version (only mother)
      nNest = 1
!
      if (.not. allocated(models)) then
        ! check for number of pets
        ! number of model must be less or equal than number of pets 
        if (nModels > petCount) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if          

        ! allocate user defined data types for models
        allocate(models(nModels))

        do i = 1, nModels
          ! allocate array to store pet list
          if (i .eq. 1) then
            nPets = (petCount/nModels)+mod(petCount, nModels)
          else
            nPets = petCount/nModels
          end if
          if (.not. allocated(models(i)%petList)) then
            allocate(models(i)%petList(nPets))
          end if

          ! assign pets to model (each component has its own pet) 
          ! For example; two models and seven cpu (or pet)
          ! model a - 0, 2, 4, 6
          ! model b - 1, 3, 5
          do j = 1, nPets 
            models(i)%petList(j) = (i+(j-1)*nModels)-1
          end do
        end do
      end if
!
!-----------------------------------------------------------------------
!     Allocate mesh variable (depend on model staggering type) 
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        if (i == Iatmos) then
          k = 2 ! cross and dot points (dot is not used actually)
        else if (i == Iocean) then
          k = 3 ! cross, u and v points
        end if
!
        if (.not. allocated(models(i)%mesh)) then
          allocate(models(i)%mesh(k,nNest(i)))
        end if
      end do
!
!-----------------------------------------------------------------------
!     Set mesh or grid information for each model 
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        if (i == Iatmos) then
          do j = 1, nNest(i)
!           cross (or cell center) points (t, p etc.)
            models(i)%mesh(1,j)%gid = 1
            models(i)%mesh(1,j)%gtype = Icross
!
            models(i)%mesh(1,j)%lon%gtype = Icross
            models(i)%mesh(1,j)%lon%name = 'xlon'
            models(i)%mesh(1,j)%lon%long_name = 'longitude at cross'
            models(i)%mesh(1,j)%lon%units = 'degrees_east'
!
            models(i)%mesh(1,j)%lat%gtype = Icross
            models(i)%mesh(1,j)%lat%name = 'xlat'
            models(i)%mesh(1,j)%lat%long_name = 'latitude at cross'
            models(i)%mesh(1,j)%lat%units = 'degrees_north'
!
            models(i)%mesh(1,j)%mask%gtype = Icross
            models(i)%mesh(1,j)%mask%name = 'mask'
            models(i)%mesh(1,j)%mask%long_name = 'land sea mask'
            models(i)%mesh(1,j)%mask%units = '1'
!
!           dot (or cell corners) points (u and v)
            models(i)%mesh(2,j)%gid = 2
            models(i)%mesh(2,j)%gtype = Idot
!
            models(i)%mesh(2,j)%lon%gtype = Idot
            models(i)%mesh(2,j)%lon%name = 'dlon'
            models(i)%mesh(2,j)%lon%long_name = 'longitude at dot'
            models(i)%mesh(2,j)%lon%units = 'degrees_east'
!
            models(i)%mesh(2,j)%lat%gtype = Idot
            models(i)%mesh(2,j)%lat%name = 'dlat'
            models(i)%mesh(2,j)%lat%long_name = 'latitude at dot'
            models(i)%mesh(2,j)%lat%units = 'degrees_north'
!
            models(i)%mesh(2,j)%mask%gtype = Idot
            models(i)%mesh(2,j)%mask%name = 'mask'
            models(i)%mesh(2,j)%mask%long_name = 'land sea mask'
            models(i)%mesh(2,j)%mask%units = '1'
          end do
        else if (i == Iocean) then
          do j = 1, nNest(i)
!           rho (or cell center) points
            models(i)%mesh(1,j)%gid = 1
            models(i)%mesh(1,j)%gtype = Icross
!
            models(i)%mesh(1,j)%lon%gtype = Icross
            models(i)%mesh(1,j)%lon%name = 'lonr'
            models(i)%mesh(1,j)%lon%long_name = 'longitude at rho'
            models(i)%mesh(1,j)%lon%units = 'degrees_east'
!
            models(i)%mesh(1,j)%lat%gtype = Icross
            models(i)%mesh(1,j)%lat%name = 'latr'
            models(i)%mesh(1,j)%lat%long_name = 'latitude at rho'
            models(i)%mesh(1,j)%lat%units = 'degrees_north'
!
            models(i)%mesh(1,j)%mask%gtype = Icross
            models(i)%mesh(1,j)%mask%name = 'mask_rho'
            models(i)%mesh(1,j)%mask%long_name = 'mask on rho'
            models(i)%mesh(1,j)%mask%units = '1'
!
!           u points
            models(i)%mesh(2,j)%gid = 2
            models(i)%mesh(2,j)%gtype = Iupoint
!
            models(i)%mesh(2,j)%lon%gtype = Iupoint
            models(i)%mesh(2,j)%lon%name = 'lonu'
            models(i)%mesh(2,j)%lon%long_name = 'longitude at u'
            models(i)%mesh(2,j)%lon%units = 'degrees_east'
!
            models(i)%mesh(2,j)%lat%gtype = Iupoint
            models(i)%mesh(2,j)%lat%name = 'latu'
            models(i)%mesh(2,j)%lat%long_name = 'latitude at u'
            models(i)%mesh(2,j)%lat%units = 'degrees_north'
!
            models(i)%mesh(2,j)%mask%gtype = Iupoint
            models(i)%mesh(2,j)%mask%name = 'mask_u'
            models(i)%mesh(2,j)%mask%long_name = 'mask on u'
            models(i)%mesh(2,j)%mask%units = '1'
!
!           v points
            models(i)%mesh(3,j)%gid = 3
            models(i)%mesh(3,j)%gtype = Ivpoint
!
            models(i)%mesh(3,j)%lon%gtype = Ivpoint
            models(i)%mesh(3,j)%lon%name = 'lonv'
            models(i)%mesh(3,j)%lon%long_name = 'longitude at v'
            models(i)%mesh(3,j)%lon%units = 'degrees_east'
!
            models(i)%mesh(3,j)%lat%gtype = Ivpoint
            models(i)%mesh(3,j)%lat%name = 'latv'
            models(i)%mesh(3,j)%lat%long_name = 'latitude at v'
            models(i)%mesh(3,j)%lat%units = 'degrees_north'
!
            models(i)%mesh(3,j)%mask%gtype = Ivpoint
            models(i)%mesh(3,j)%mask%name = 'mask_v'
            models(i)%mesh(3,j)%mask%long_name = 'mask on v'
            models(i)%mesh(3,j)%mask%units = '1'
          end do 
        end if
!
!-----------------------------------------------------------------------
!       Allocate ESMF Layout and distribution objects.
!-----------------------------------------------------------------------
!
        if (.not. allocated(models(i)%deLayout)) then
          allocate(models(i)%deLayout(nNest(i)))
        end if
        if (.not. allocated(models(i)%distGrid)) then
          allocate(models(i)%distGrid(nNest(i)))
        end if
        if (.not. allocated(models(i)%arrSpec)) then
          allocate(models(i)%arrSpec(nNest(i)))
        end if
!
!-----------------------------------------------------------------------
!       Set import and export fields 
!-----------------------------------------------------------------------
!
        if (i == Iatmos) then
          if (.not. allocated(models(i)%dataExport)) then 
            allocate(models(i)%dataExport(8,nNest(i)))
          end if       
          if (.not. allocated(models(i)%dataImport)) then 
            allocate(models(i)%dataImport(1,nNest(i)))
          end if
!
          do j = 1, nNest(i)      
            models(i)%dataImport(1,j)%fid = 1
            models(i)%dataImport(1,j)%gtype = Icross
            models(i)%dataImport(1,j)%name = 'SST'
            models(i)%dataImport(1,j)%long_name = 'Ground Temperature'
            models(i)%dataImport(1,j)%units = 'Celsius'
!
            models(i)%dataExport(1,j)%fid = 1
            models(i)%dataExport(1,j)%gtype = Icross
            models(i)%dataExport(1,j)%name = 'Pair'
            models(i)%dataExport(1,j)%long_name = 'Surface Pressure'
            models(i)%dataExport(1,j)%units = '?'
!
            models(i)%dataExport(2,j)%fid = 2
            models(i)%dataExport(2,j)%gtype = Icross
            models(i)%dataExport(2,j)%name = 'Tair'
            models(i)%dataExport(2,j)%long_name = &
                                              'Surface Air Temperature'
            models(i)%dataExport(2,j)%units = 'Celsius'             
!
            models(i)%dataExport(3,j)%fid = 3
            models(i)%dataExport(3,j)%gtype = Icross
            models(i)%dataExport(3,j)%name = 'Qair'
            models(i)%dataExport(3,j)%long_name = &
                                         'Surface Air Specific Humidity'
            models(i)%dataExport(3,j)%units = 'g/kg'
!
            models(i)%dataExport(4,j)%fid = 4
            models(i)%dataExport(4,j)%gtype = Icross
            models(i)%dataExport(4,j)%name = 'swrad'
            models(i)%dataExport(4,j)%long_name = &
                                        'solar shortwave radiation flux'
            models(i)%dataExport(4,j)%units = 'watt meter-2'
!
            models(i)%dataExport(5,j)%fid = 5
            models(i)%dataExport(5,j)%gtype = Icross
            models(i)%dataExport(5,j)%name = 'lwrad_down'
            models(i)%dataExport(5,j)%long_name = &
                                   'downwelling longwave radiation flux'
            models(i)%dataExport(5,j)%units = 'watt meter-2'
!
            models(i)%dataExport(6,j)%fid = 6
            models(i)%dataExport(6,j)%gtype = Icross
            models(i)%dataExport(6,j)%name = 'rain'
            models(i)%dataExport(6,j)%long_name = &
                                         'rain fall rate'
            models(i)%dataExport(6,j)%units ='kilogram meter-2 second-1'
!
            models(i)%dataExport(7,j)%fid = 7
            models(i)%dataExport(7,j)%gtype = Icross
            models(i)%dataExport(7,j)%name = 'Uwind'
            models(i)%dataExport(7,j)%long_name = &
                                         'surface u-wind component'
            models(i)%dataExport(7,j)%units = 'meter second-1'
!
            models(i)%dataExport(8,j)%fid = 8
            models(i)%dataExport(8,j)%gtype = Icross
            models(i)%dataExport(8,j)%name = 'Vwind'
            models(i)%dataExport(8,j)%long_name = &
                                         'surface v-wind component'
            models(i)%dataExport(8,j)%units = 'meter second-1'
          end do 
        else if (i == Iocean) then
          if (.not. allocated(models(i)%dataExport)) then
            allocate(models(i)%dataExport(1,nNest(i)))
          end if
          if (.not. allocated(models(i)%dataImport)) then
            allocate(models(i)%dataImport(8,nNest(i)))
          end if
!
          do j = 1, nNest(i)
            models(i)%dataExport(1,j)%fid = 1
            models(i)%dataExport(1,j)%gtype = Icross
            models(i)%dataExport(1,j)%name = 'SST'
            models(i)%dataExport(1,j)%long_name = 'Ground Temperature'
            models(i)%dataExport(1,j)%units = 'Celsius'
!
            models(i)%dataImport(1,j)%fid = 1
            models(i)%dataImport(1,j)%gtype = Icross
            models(i)%dataImport(1,j)%name = 'Pair'
            models(i)%dataImport(1,j)%long_name = 'Surface Pressure'
            models(i)%dataImport(1,j)%units = '?'
!
            models(i)%dataImport(2,j)%fid = 2
            models(i)%dataImport(2,j)%gtype = Icross
            models(i)%dataImport(2,j)%name = 'Tair'
            models(i)%dataImport(2,j)%long_name = &
                                              'Surface Air Temperature'
            models(i)%dataImport(2,j)%units = 'Celsius'
!
            models(i)%dataImport(3,j)%fid = 3
            models(i)%dataImport(3,j)%gtype = Icross
            models(i)%dataImport(3,j)%name = 'Qair'
            models(i)%dataImport(3,j)%long_name = &
                                         'Surface Air Specific Humidity'
            models(i)%dataImport(3,j)%units = 'g/kg'
!
            models(i)%dataImport(4,j)%fid = 4
            models(i)%dataImport(4,j)%gtype = Icross
            models(i)%dataImport(4,j)%name = 'swrad'
            models(i)%dataImport(4,j)%long_name = &
                                      'solar shortwave radiation flux'
            models(i)%dataImport(4,j)%units = 'watt meter-2'
!
            models(i)%dataImport(5,j)%fid = 5
            models(i)%dataImport(5,j)%gtype = Icross
            models(i)%dataImport(5,j)%name = 'lwrad_down'
            models(i)%dataImport(5,j)%long_name = &
                                 'downwelling longwave radiation flux'
            models(i)%dataImport(5,j)%units = 'watt meter-2'
!
            models(i)%dataImport(6,j)%fid = 6
            models(i)%dataImport(6,j)%gtype = Icross
            models(i)%dataImport(6,j)%name = 'rain'
            models(i)%dataImport(6,j)%long_name = &
                                       'rain fall rate'
            models(i)%dataImport(6,j)%units ='kilogram meter-2 second-1'
!
            models(i)%dataImport(7,j)%fid = 7
            models(i)%dataImport(7,j)%gtype = Icross
            models(i)%dataImport(7,j)%name = 'Uwind'
            models(i)%dataImport(7,j)%long_name = &
                                       'surface u-wind component'
            models(i)%dataImport(7,j)%units = 'meter second-1'
!
            models(i)%dataImport(8,j)%fid = 8
            models(i)%dataImport(8,j)%gtype = Icross
            models(i)%dataImport(8,j)%name = 'Vwind'
            models(i)%dataImport(8,j)%long_name = &
                                       'surface v-wind component'
            models(i)%dataImport(8,j)%units = 'meter second-1'
          end do
        end if   
      end do

!
!-----------------------------------------------------------------------
!     Set return flag to success.
!-----------------------------------------------------------------------
!
      rc = ESMF_SUCCESS
!
      end subroutine allocate_cpl
!
      subroutine time_reconcile()
      implicit none

      integer :: iarr(6)
      integer :: i, j, petCount, localPet, comm, mysec, rc
      character(len=80) :: timeString, name

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
      do i = 1, nModels
!
!-----------------------------------------------------------------------
!       Reconcile export state of the gridded components 
!-----------------------------------------------------------------------
!
        call ESMF_StateReconcile(models(i)%stateExport,                 &
                                 vm=vm,                                 &
                                 attreconflag=ESMF_ATTRECONCILE_ON,     &
                                 rc=rc)
!
!-----------------------------------------------------------------------
!       Get coupler time step and set time interval to exchange 
!       data between gridded components
!-----------------------------------------------------------------------
!  
        if (i == Iatmos) then
          call ESMF_AttributeGet(models(i)%stateExport,                 &
                                 name='coupler time step',              &
                                 value=j,                               &
                                 rc=rc)
          call ESMF_TimeIntervalSet(timeStep,                           &
                                    s=j,                                &
                                    rc=rc)
        end if
!
!-----------------------------------------------------------------------
!       Get start time
!-----------------------------------------------------------------------
!  
        j = ubound(iarr,dim=1)
        call ESMF_AttributeGet(models(i)%stateExport,                   &
                               name='start time',                       &
                               valueList=iarr,                          &
                               itemCount=j,                             &
                               rc=rc)
        call ESMF_TimeSet (models(i)%strTime,                           &
                           yy=iarr(1),                                  &
                           mm=iarr(2),                                  &
                           dd=iarr(3),                                  &
                           h=iarr(4),                                   &
                           m=iarr(5),                                   &
                           s=iarr(6),                                   &
                           calkindflag=ESMF_CALKIND_GREGORIAN,          &
                           rc=rc)
!
!-----------------------------------------------------------------------
!       Get stop time
!-----------------------------------------------------------------------
!  
        call ESMF_AttributeGet(models(i)%stateExport,                   &
                               name='stop time',                        &
                               valueList=iarr,                          &
                               itemCount=j,                             &
                               rc=rc)
        call ESMF_TimeSet (models(i)%endTime,                           &
                           yy=iarr(1),                                  &
                           mm=iarr(2),                                  &
                           dd=iarr(3),                                  &
                           h=iarr(4),                                   &
                           m=iarr(5),                                   &
                           s=iarr(6),                                   &
                           calkindflag=ESMF_CALKIND_GREGORIAN,          &
                           rc=rc)
      end do
!
!-----------------------------------------------------------------------
!     Create coupler component clock 
!-----------------------------------------------------------------------
!
      startTime = models(Iatmos)%strTime
      stopTime = models(Iatmos)%endTime
      do i = 1, nModels
        if (models(i)%strTime < startTime) then
          startTime = models(i)%strTime
        end if      
        if (models(i)%endTime < stopTime) then
          stopTime = models(i)%endTime
        end if 
      end do     
!
      name = 'Coupler component clock'
      clock = ESMF_ClockCreate (name=trim(name),                        &
                                timeStep=timeStep,                      &
                                startTime=startTime,                    &
                                stopTime=stopTime,                      &
                                rc=rc)
!
!-----------------------------------------------------------------------
!     Validate external time clock.
!-----------------------------------------------------------------------
!
      call ESMF_ClockValidate (clock, rc=rc)
!
!-----------------------------------------------------------------------
!     Print
!-----------------------------------------------------------------------
!
      call ESMF_TimeGet(startTime, timeString=timeString, rc=rc)
      write(*,30) localPet, 'Start Time   ', trim(timeString)
      call ESMF_TimeGet(stopTime, timeString=timeString, rc=rc)
      write(*,30) localPet, 'Stop Time    ', trim(timeString)
      call ESMF_TimeIntervalGet(timeStep, s = mysec, rc=rc)
      write(*,40) localPet, 'Time Interval', mysec
!
 30   format (' PET (', I2, ') - ', A, ' = ', A)
 40   format (' PET (', I2, ') - ', A, ' = ', I10)
!
      end subroutine time_reconcile
!
      subroutine check_err(rc)
      implicit none
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      integer, intent(in) :: rc 
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
      integer :: rclocal
!
!-----------------------------------------------------------------------
!     Terminate execution due to fatal error 
!-----------------------------------------------------------------------
!      
      if (ESMF_LogFoundError(rcToCheck=rc,                              &
                             msg=ESMF_LOGERR_PASSTHRU,                  &
                             line=__LINE__,                             &
                             file=__FILE__,                             &
                             rcToReturn=rclocal)) then
        call ESMF_Finalize(endflag=ESMF_END_ABORT, rc=rclocal)
      end if
!
      end subroutine check_err
!
      end module mod_couplerr
