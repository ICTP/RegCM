      module mod_couplerr
!
!-----------------------------------------------------------------------
!     Imported modules 
!-----------------------------------------------------------------------
!
      use ESMF
!
      use mod_realkinds, only : sp, dp
!
      implicit none
!
!-----------------------------------------------------------------------
!     User defined data types for coupling interface 
!-----------------------------------------------------------------------
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
        real*8 :: scale_factor
        real*8 :: add_offset
        type(ESMF_Array) :: array
        type(ESMF_Field) :: field
        type(ESMF_RouteHandle) :: rhandle
        real*8, dimension(:,:), pointer :: ptr
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
        type(ESMF_Grid) :: grid
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
        integer :: nproc
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
!-----------------------------------------------------------------------
!     Global module variables 
!-----------------------------------------------------------------------
!
      type(ESM_Model), allocatable :: models(:) 
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
!     d --------- d   d --- v --- d  
!     |           |   |           |
!     |     c     |   u     c     u
!     |           |   |           |
!     d --------- d   d --- v --- d     
!     Arakawa - B     Arakawa - C
!     RegCM           ROMS (c = rho, d = psi)
!-----------------------------------------------------------------------
!
      character(len=6) :: GRIDDES(4) = (/ "CROSS", "DOT", "U", "V" /)
      integer :: Icross  = 1
      integer :: Idot    = 2
      integer :: Iupoint = 3
      integer :: Ivpoint = 4
!
!-----------------------------------------------------------------------
!     Variables for coupling (direction, import and export variables)  
!-----------------------------------------------------------------------
!
      character(len=40), parameter :: DIRECTION = 'DIRECTION'
!
      integer, parameter :: FORWARD_ON   = 1
      integer, parameter :: FORWARD_OFF  = 0
!
      character(ESMF_MAXSTR), allocatable :: itemNamesImportF(:)
      character(ESMF_MAXSTR), allocatable :: itemNamesExportF(:)
      character(ESMF_MAXSTR), allocatable :: itemNamesImportB(:)
      character(ESMF_MAXSTR), allocatable :: itemNamesExportB(:)
!
!-----------------------------------------------------------------------
!     Coupled model parameters
!-----------------------------------------------------------------------
!
      real*8, parameter :: MISSING_R8 = 1.0d20
!
      integer :: cpl_dtsec, cpl_exvars, cpl_interp, cpl_dbglevel
      logical :: cpl_bdysmooth
!
!-----------------------------------------------------------------------
!     Coupler component variables 
!-----------------------------------------------------------------------
!
      type(ESMF_CplComp) :: cplComp
      type(ESMF_VM) :: cplVM
      type(ESMF_Time) :: cplStartTime, cplStopTime
      type(ESMF_TimeInterval) :: cplTimeStep
      type(ESMF_Clock) :: cplClock
!
!-----------------------------------------------------------------------
!     Coupler component variables 
!-----------------------------------------------------------------------
!
      type(ESMF_RouteHandle) :: routeHandleF
      type(ESMF_RouteHandle) :: routeHandleB
      integer(ESMF_KIND_I4), pointer :: indices(:,:)
      real(ESMF_KIND_R8), pointer :: weights(:)
!
!-----------------------------------------------------------------------
!     Constants
!     Cp  : Specific heat for seawater (Joules/Kg/degC)
!     rho0: Mean density (Kg/m3) - Boussinesq approximation
!-----------------------------------------------------------------------
!
      real*8, parameter :: Cp = 3985.0d0
      real*8, parameter :: rho0 = 1025.0d0
      real*8, parameter :: Hscale=1.0d0/(rho0*Cp)
!
!-----------------------------------------------------------------------
!     Temporary variables ( atm->ocn )
!-----------------------------------------------------------------------
!
      real*8, allocatable, dimension(:,:) :: zi, z995, u995, v995
      real*8, allocatable, dimension(:,:) :: psurf, tsurf, t995, q995
      real*8, allocatable, dimension(:,:) :: t2, q2, u10, v10
      real*8, allocatable, dimension(:,:) :: taux, tauy
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
      character(100) :: fmt_123 
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
          ! set number of processor
          models(i)%nproc = nPets

          ! assign pets to model (each component has its own pet) 
          ! For example; two models and seven cpu (or pet)
          ! model a - 0, 2, 4, 6
          ! model b - 1, 3, 5
          !do j = 1, nPets 
          !  models(i)%petList(j) = (i+(j-1)*nModels)-1
          !end do

          ! For example; two models and seven cpu (or pet)
          ! model a - 0, 1, 2, 3
          ! model b - 4, 5, 6
          do j = 1, nPets
            if (i .eq. 1) then 
              models(i)%petList(j) = j-1
            else
              k = ubound(models(i-1)%petList, dim=1)
              models(i)%petList(j) = models(i-1)%petList(k)+j 
            end if 
          end do
!
          k = ubound(models(i)%petList, dim=1)
          write(fmt_123, fmt="('(A3, ', I3, 'I4)')") k
          if (i .eq. 1) then
            write(*, fmt=trim(fmt_123)) "ATM", models(i)%petList 
          else
            write(*, fmt=trim(fmt_123)) "OCN", models(i)%petList
          end if
!
        end do
      end if
!
!-----------------------------------------------------------------------
!     Allocate mesh variable (depend on model staggering type) 
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        if (i == Iatmos) then
          k = 2 ! cross and dot (B Grid)
        else if (i == Iocean) then 
          k = 1 ! only cross is used for variable SST (C Grid)
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
!            models(i)%mesh(2,j)%gid = 2
!            models(i)%mesh(2,j)%gtype = Iupoint
!
!            models(i)%mesh(2,j)%lon%gtype = Iupoint
!            models(i)%mesh(2,j)%lon%name = 'lonu'
!            models(i)%mesh(2,j)%lon%long_name = 'longitude at u'
!            models(i)%mesh(2,j)%lon%units = 'degrees_east'
!
!            models(i)%mesh(2,j)%lat%gtype = Iupoint
!            models(i)%mesh(2,j)%lat%name = 'latu'
!            models(i)%mesh(2,j)%lat%long_name = 'latitude at u'
!            models(i)%mesh(2,j)%lat%units = 'degrees_north'
!
!            models(i)%mesh(2,j)%mask%gtype = Iupoint
!            models(i)%mesh(2,j)%mask%name = 'mask_u'
!            models(i)%mesh(2,j)%mask%long_name = 'mask on u'
!            models(i)%mesh(2,j)%mask%units = '1'
!
!           v points
!            models(i)%mesh(3,j)%gid = 3
!            models(i)%mesh(3,j)%gtype = Ivpoint
!
!            models(i)%mesh(3,j)%lon%gtype = Ivpoint
!            models(i)%mesh(3,j)%lon%name = 'lonv'
!            models(i)%mesh(3,j)%lon%long_name = 'longitude at v'
!            models(i)%mesh(3,j)%lon%units = 'degrees_east'
!
!            models(i)%mesh(3,j)%lat%gtype = Ivpoint
!            models(i)%mesh(3,j)%lat%name = 'latv'
!            models(i)%mesh(3,j)%lat%long_name = 'latitude at v'
!            models(i)%mesh(3,j)%lat%units = 'degrees_north'
!
!            models(i)%mesh(3,j)%mask%gtype = Ivpoint
!            models(i)%mesh(3,j)%mask%name = 'mask_v'
!            models(i)%mesh(3,j)%mask%long_name = 'mask on v'
!            models(i)%mesh(3,j)%mask%units = '1'
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
            models(i)%dataImport(1,j)%long_name = &
            'Sea Surface Temperature'
            models(i)%dataImport(1,j)%units = 'Kelvin'
            models(i)%dataImport(1,j)%scale_factor = 1.0d0
            models(i)%dataImport(1,j)%add_offset = 273.15d0
!
            models(i)%dataExport(1,j)%fid = 1
            models(i)%dataExport(1,j)%gtype = Icross
            models(i)%dataExport(1,j)%name = 'Pair'
            models(i)%dataExport(1,j)%long_name = 'Surface Pressure'
            models(i)%dataExport(1,j)%units = 'Pascal'
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
            models(i)%dataExport(7,j)%gtype = Idot
            models(i)%dataExport(7,j)%name = 'Uwind'
            models(i)%dataExport(7,j)%long_name = &
            'surface u-wind component'
            models(i)%dataExport(7,j)%units = 'meter second-1'
!
            models(i)%dataExport(8,j)%fid = 8
            models(i)%dataExport(8,j)%gtype = Idot
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
            models(i)%dataExport(1,j)%long_name = & 
           'Sea Surface Temperature'
            models(i)%dataExport(1,j)%units = 'Celsius'
!
            models(i)%dataImport(1,j)%fid = 1
            models(i)%dataImport(1,j)%gtype = Icross
            models(i)%dataImport(1,j)%name = 'Pair'
            models(i)%dataImport(1,j)%long_name = 'Surface Pressure'
            models(i)%dataImport(1,j)%units = 'milibar'
            models(i)%dataImport(1,j)%scale_factor = 1.0d0
            models(i)%dataImport(1,j)%add_offset = 0.0d0
!
            models(i)%dataImport(2,j)%fid = 2
            models(i)%dataImport(2,j)%gtype = Icross
            models(i)%dataImport(2,j)%name = 'Tair'
            models(i)%dataImport(2,j)%long_name = &
            'Surface Air Temperature'
            models(i)%dataImport(2,j)%units = 'Celsius'
            models(i)%dataImport(2,j)%scale_factor = 1.0d0
            models(i)%dataImport(2,j)%add_offset = -273.15d0
!
            models(i)%dataImport(3,j)%fid = 3
            models(i)%dataImport(3,j)%gtype = Icross
            models(i)%dataImport(3,j)%name = 'Qair'
            models(i)%dataImport(3,j)%long_name = &
            'Surface Air Specific Humidity'
            models(i)%dataImport(3,j)%units = 'g/kg'
            models(i)%dataImport(3,j)%scale_factor = 1.0d3
            models(i)%dataImport(3,j)%add_offset = 0.0d0
!
            models(i)%dataImport(4,j)%fid = 4
            models(i)%dataImport(4,j)%gtype = Icross
            models(i)%dataImport(4,j)%name = 'swrad'
            models(i)%dataImport(4,j)%long_name = &
            'solar shortwave radiation flux'
            models(i)%dataImport(4,j)%units = 'Celsius m/s'
            models(i)%dataImport(4,j)%scale_factor = Hscale 
            models(i)%dataImport(4,j)%add_offset = 0.0d0
!
            models(i)%dataImport(5,j)%fid = 5
            models(i)%dataImport(5,j)%gtype = Icross
            models(i)%dataImport(5,j)%name = 'lwrad_down'
            models(i)%dataImport(5,j)%long_name = &
            'downwelling longwave radiation flux'
            models(i)%dataImport(5,j)%units = 'Celsius m/s'
            models(i)%dataImport(5,j)%scale_factor = Hscale
            models(i)%dataImport(5,j)%add_offset = 0.0d0
!
            models(i)%dataImport(6,j)%fid = 6
            models(i)%dataImport(6,j)%gtype = Icross
            models(i)%dataImport(6,j)%name = 'rain'
            models(i)%dataImport(6,j)%long_name = &
            'rain fall rate'
            models(i)%dataImport(6,j)%units ='kilogram meter-2 second-1'
            models(i)%dataImport(6,j)%scale_factor = 1.0d0
            models(i)%dataImport(6,j)%add_offset = 0.0d0
!
            models(i)%dataImport(7,j)%fid = 7
            models(i)%dataImport(7,j)%gtype = Icross
            models(i)%dataImport(7,j)%name = 'Uwind'
            models(i)%dataImport(7,j)%long_name = &
            'surface u-wind component'
            models(i)%dataImport(7,j)%units = 'meter second-1'
            models(i)%dataImport(7,j)%scale_factor = 1.0d0
            models(i)%dataImport(7,j)%add_offset = 0.0d0
!
            models(i)%dataImport(8,j)%fid = 8
            models(i)%dataImport(8,j)%gtype = Icross
            models(i)%dataImport(8,j)%name = 'Vwind'
            models(i)%dataImport(8,j)%long_name = &
            'surface v-wind component'
            models(i)%dataImport(8,j)%units = 'meter second-1'
            models(i)%dataImport(8,j)%scale_factor = 1.0d0
            models(i)%dataImport(8,j)%add_offset = 0.0d0
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
      subroutine time_reconcile(first)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      logical, intent(in) :: first
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: iarr(6)
      integer :: i, nitems, petCount, localPet, comm, mysec, rc
      character(len=80) :: timeString, name
!
!-----------------------------------------------------------------------
!     Get information from VM (MPI Communicator, number of PETs etc.)
!-----------------------------------------------------------------------
!
      call ESMF_VMGet (cplVM,                                           &
                       petCount=petCount,                               &
                       localPet=localPet,                               &
                       mpiCommunicator=comm,                            &
                       rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      do i = 1, nModels
!
!-----------------------------------------------------------------------
!     Reconcile export state of the gridded components 
!-----------------------------------------------------------------------
!
      call ESMF_StateReconcile (models(i)%stateExport,                  &
                                vm=cplVM,                               &
                                attreconflag=ESMF_ATTRECONCILE_ON,      &
                                rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      if (i == Iatmos) then
!
!-----------------------------------------------------------------------
!     Get parameter for coupler component time step and set time 
!     interval to exchange data between gridded components
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (models(i)%stateExport,                    &
                              name='coupler time step',                 &
                              value=cpl_dtsec,                          &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_TimeIntervalSet (cplTimeStep,                           &
                                 s=cpl_dtsec,                           &
                                 rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get parameter for list of exchange variables
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (models(i)%stateExport,                    &
                              name='exchange variable mode',            &
                              value=cpl_exvars,                         &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get parameter for interpolation type 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (models(i)%stateExport,                    &
                              name='interpolation type',                &
                              value=cpl_interp,                         &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get parameter for boundary smoothing 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (models(i)%stateExport,                    &
                              name='boundary smoothing',                &
                              value=cpl_bdysmooth,                      &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get parameter for debugging 
!-----------------------------------------------------------------------
!
      call ESMF_AttributeGet (models(i)%stateExport,                    &
                              name='debug level',                       &
                              value=cpl_dbglevel,                       &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end if
      end do
!
!-----------------------------------------------------------------------
!     Write coupled model parameters
!-----------------------------------------------------------------------
!
      if (cpl_dbglevel > 0) then
        write(*, 30) localPet, 'COUPLED MODEL TIME STEP', cpl_dtsec
        write(*, 30) localPet, 'EXCHANGE VARIABLE MODE ', cpl_exvars
        write(*, 30) localPet, 'INTERPOLATION MODE     ', cpl_interp
        write(*, 20) localPet, 'BOUNDARY SMOOTHING     ', cpl_bdysmooth
        write(*, 30) localPet, 'DEBUG LEVEL            ', cpl_dbglevel
      end if
!
      if (.not. first) then
      do i = 1, nModels
!
!-----------------------------------------------------------------------
!     Get start time
!-----------------------------------------------------------------------
!  
      nitems = size(iarr, dim=1)
      call ESMF_AttributeGet (models(i)%stateExport,                    &
                              name='start time',                        &
                              valueList=iarr,                           &
                              itemCount=nitems,                         &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_TimeSet (models(i)%strTime,                             &
                         yy=iarr(1),                                    &
                         mm=iarr(2),                                    &
                         dd=iarr(3),                                    &
                         h=iarr(4),                                     &
                         m=iarr(5),                                     &
                         s=iarr(6),                                     &
                         calkindflag=ESMF_CALKIND_GREGORIAN,            &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Get stop time
!-----------------------------------------------------------------------
!  
      call ESMF_AttributeGet (models(i)%stateExport,                    &
                              name='stop time',                         &
                              valueList=iarr,                           &
                              itemCount=nitems,                         &
                              rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
      call ESMF_TimeSet (models(i)%endTime,                             &
                         yy=iarr(1),                                    &
                         mm=iarr(2),                                    &
                         dd=iarr(3),                                    &
                         h=iarr(4),                                     &
                         m=iarr(5),                                     &
                         s=iarr(6),                                     &
                         calkindflag=ESMF_CALKIND_GREGORIAN,            &
                         rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      end do
!
!-----------------------------------------------------------------------
!     Create coupler component clock 
!-----------------------------------------------------------------------
!
      cplStartTime = models(Iatmos)%strTime
      cplStopTime = models(Iatmos)%endTime
      do i = 1, nModels
        if (models(i)%strTime > cplStartTime) then
          cplStartTime = models(i)%strTime
        end if      
        if (models(i)%endTime < cplStopTime) then
          cplStopTime = models(i)%endTime
        end if 
      end do     
!
      name = 'Coupler component clock'
      cplClock = ESMF_ClockCreate (name=trim(name),                     &
                                   timeStep=cplTimeStep,                &
                                   startTime=cplStartTime,              &
                                   stopTime=cplStopTime,                &
                                   rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Validate external time clock.
!-----------------------------------------------------------------------
!
      call ESMF_ClockValidate (cplClock, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!-----------------------------------------------------------------------
!     Print
!-----------------------------------------------------------------------
!
      if (cpl_dbglevel > 1) then
      call ESMF_TimeGet(cplStartTime, timeString=timeString, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      write(*,40) localPet, 'Start Time   ', trim(timeString)
!
      call ESMF_TimeGet(cplStopTime, timeString=timeString, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      write(*,40) localPet, 'Stop Time    ', trim(timeString)
!
      call ESMF_TimeIntervalGet(cplTimeStep, s = mysec, rc=rc)
      if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)
      write(*,50) localPet, 'Time Interval', mysec
      end if
!
      end if
!
!-----------------------------------------------------------------------
!     Formats 
!-----------------------------------------------------------------------
!
 20   format(' PET (', I2, ') - Parameter ', A, ' = ', L)      
 30   format(' PET (', I2, ') - Parameter ', A, ' = ', I8)      
 40   format(' PET (', I2, ') - ', A, ' = ', A)
 50   format(' PET (', I2, ') - ', A, ' = ', I10)
!
      end subroutine time_reconcile
!
      integer function getIndex(list, name)
      implicit none
!
!=======================================================================
!                                                                      !
!  This subroutine is used to get indices of specified variable        !
!                                                                      !
!=======================================================================
!
!  Imported variable definitions.
!
      character(len=*), dimension(:), intent(in) :: list
      character(len=*), intent(in) :: name
!
!  Local variable definitions.
!
      integer :: i

      do i = 1, ubound(list, dim=1)
        if (trim(list(i)) .eq. trim(name)) then
          GetIndex = i
          return
        end if
      end do

      end function getIndex
!
      integer function getMeshID(list, gtype)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESM_Mesh), intent(in) :: list(:)
      integer, intent(in) :: gtype
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i
!     
!-----------------------------------------------------------------------
!     Find index of specified mesh 
!-----------------------------------------------------------------------
!
      do i = 1, size(list, dim=1) 
        if (list(i)%gtype == gtype) then
          getMeshID = i
          return
        end if
      end do
      end function getMeshId
!
      subroutine print_matrix_r8(inp, iskip, jskip, pet, id, header)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      real*8, intent(in) :: inp(:,:)
      integer, intent(in) ::  iskip, jskip, pet, id
      character(len=*), intent(in) :: header
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, j, imin, imax, jmin, jmax
      character(100) :: fmt_123
!
!-----------------------------------------------------------------------
!     Write data 
!-----------------------------------------------------------------------
!
      imin = lbound(inp, dim=1)
      imax = ubound(inp, dim=1)
      jmin = lbound(inp, dim=2)
      jmax = ubound(inp, dim=2)
!
      write(6, fmt="('PET(',I2,') - ',A)") pet, trim(header)
!
      write(fmt_123, fmt="('(/, 5X, ', I3, 'I10)')") (imax-imin)+1
      write(id, fmt=trim(fmt_123))  (i, i=imin, imax, iskip)
!   
      write(fmt_123, fmt="('(I5, ', I3, 'F10.2)')") imax
      do j=jmin, jmax, jskip
        write(id, fmt=trim(fmt_123)) j, (inp(i,j),i=imin, imax, iskip)
      end do
!
      return
      end subroutine print_matrix_r8      
!
      subroutine calc_uvmet (u, v, urot, vrot, localPet)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_dynparam, only : iproj, truelatl, truelath
      use mod_atm_interface, only : mddom
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      real(sp), dimension(:,:), intent(in) :: u, v
      integer, intent(in) :: localPet
      real(sp), dimension(:,:), intent(inout) :: urot, vrot
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      real(dp) :: cone
      real(dp) :: PI, RAD_PER_DEG
!
!-----------------------------------------------------------------------
!     Calculate parameters 
!-----------------------------------------------------------------------
!
      PI = atan(1.0d0)*4.0d0
      RAD_PER_DEG = PI/180.0d0 

      cone = 1.0d0
!
      if (iproj .eq. 'LAMCON') then !  Lambert Conformal mapping
        if (abs(truelatl-truelath) > 0.1d0) then
          cone = (log(cos(truelatl*RAD_PER_DEG))-                      &
                  log(cos(truelath*RAD_PER_DEG)))/                     &
                 (log(tan((90.0d0-abs(truelatl))*RAD_PER_DEG*0.5d0))-  &
                  log(tan((90.0d0-abs(truelath))*RAD_PER_DEG*0.5d0)))
        else
          cone = dsin(abs(truelatl)*RAD_PER_DEG)
        end if
      end if
!
      write(*, fmt="(A3, 5I8)") "U--", localPet,                        &
                            lbound(u, dim=1),                           &
                            ubound(u, dim=1),                           &
                            lbound(u, dim=2),                           &
                            ubound(u, dim=2)
      write(*, fmt="(A3, 5I8)") "D--", localPet,                        &
                            lbound(mddom%xlon, dim=1),                  &
                            ubound(mddom%xlon, dim=1),                  &
                            lbound(mddom%xlon, dim=2),                  &
                            ubound(mddom%xlon, dim=2)

      return
      end subroutine calc_uvmet
!
      subroutine ocnrough(zo,zot,zoq,ustar,visa)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_constants
      use mod_bats_common, only : iocnrough
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      real(dp), dimension(:,:), intent(in) :: ustar, visa
      real(dp), dimension(:,:), intent(out) :: zoq, zot
      real(dp), dimension(:,:), intent(inout) :: zo
!
!     Graziano: make this selectable on iocnrough flag
!
      if ( iocnrough == 2 ) then
        zo = (0.013d0*ustar*ustar)*regrav+0.11d0*visa/ustar
      else
        zo = (0.0065d0*ustar*ustar)*regrav
      end if
!
      zoq = zo/dexp(2.67d0*(((ustar*zo)/visa)**d_rfour)-2.57d0)
      zot = zoq
!
      return
      end subroutine ocnrough
!
      function psi (k, zeta)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_constants
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      integer,  intent(in) :: k
      real(dp), intent(in) :: zeta
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      real(dp) :: chik, psi
!
      chik = (d_one-16.0D0*zeta)**d_rfour
      if ( k == 1 ) then
        psi = d_two*dlog((d_one+chik)*d_half) +                       &
                      dlog((d_one+chik*chik)*d_half) -                  &
                d_two*datan(chik) + d_two*datan(d_one)
      else
        psi = d_two*dlog((d_one+chik*chik)*d_half)
      end if
      end function psi
!
      function psi2d (k, zeta)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_constants
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      integer, intent(in) :: k
      real(dp), dimension(:,:), intent(in) :: zeta
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: jxmin, jxmax, iymin, iymax
      real(dp), allocatable, dimension(:,:) :: chik, psi2d
!
      jxmin = lbound(zeta, dim=1)
      jxmax = ubound(zeta, dim=1)
      iymin = lbound(zeta, dim=2)
      iymax = ubound(zeta, dim=2)

      if (.not.allocated(chik))allocate(chik(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(psi2d))allocate(psi2d(jxmin:jxmax,iymin:iymax))
!
      chik = (d_one-16.0D0*zeta)**d_rfour
      if ( k == 1 ) then
        psi2d = d_two*dlog((d_one+chik)*d_half) +                       &
                      dlog((d_one+chik*chik)*d_half) -                  &
                d_two*datan(chik) + d_two*datan(d_one)
      else
        psi2d = d_two*dlog((d_one+chik*chik)*d_half)
      end if
      end function psi2d
!
      subroutine adjustvars(z, t, q, u, v, zi, ps, ts, t2, q2,           &
                           u10, v10, taux, tauy)
!
!-----------------------------------------------------------------------
!     Used module declarations 
!-----------------------------------------------------------------------
!
      use mod_constants
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      real(dp), dimension(:,:), intent(in) :: z 
      real(dp), dimension(:,:), intent(in) :: t 
      real(dp), dimension(:,:), intent(in) :: q 
      real(dp), dimension(:,:), intent(in) :: u 
      real(dp), dimension(:,:), intent(in) :: v 
      real(dp), dimension(:,:), intent(in) :: zi 
      real(dp), dimension(:,:), intent(in) :: ps 
      real(dp), dimension(:,:), intent(in) :: ts
!
      real(dp), dimension(:,:), intent(inout) :: t2
      real(dp), dimension(:,:), intent(inout) :: q2
      real(dp), dimension(:,:), intent(inout) :: u10
      real(dp), dimension(:,:), intent(inout) :: v10
      real(dp), dimension(:,:), intent(inout) :: taux 
      real(dp), dimension(:,:), intent(inout) :: tauy
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer :: i, jxmin, jxmax, iymin, iymax
      real(dp), allocatable, dimension(:,:) :: th, dth, qs, dqh
      real(dp), allocatable, dimension(:,:) :: thv, dthv, rho, xlv, visa
      real(dp), allocatable, dimension(:,:) :: wc, ustr, wspd, um, zo
      real(dp), allocatable, dimension(:,:) :: rb, zeta, obu, tstr, qstr
      real(dp), allocatable, dimension(:,:) :: zot, zoq, thvstr, facttq 
!
      real(dp), parameter :: zetat = 0.465D0
      real(dp), parameter :: zetam = 1.574D0
!
      real(dp), parameter :: z2  = d_two   ! m  (reference height)
      real(dp), parameter :: z10 = d_10    ! m  (reference height)
      real(dp), parameter :: zbeta = d_one ! -  (in computing W_*)
!
!     Get dimension bounds and allocate variables
!
      jxmin = lbound(t, dim=1)
      jxmax = ubound(t, dim=1)
      iymin = lbound(t, dim=2)
      iymax = ubound(t, dim=2)

      if (.not.allocated(th)) allocate(th(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(dth)) allocate(dth(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(qs)) allocate(qs(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(dqh)) allocate(dqh(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(thv)) allocate(thv(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(dthv)) allocate(dthv(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(rho)) allocate(rho(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(xlv)) allocate(xlv(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(visa)) allocate(visa(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(wc)) allocate(wc(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(ustr)) allocate(ustr(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(wspd)) allocate(wspd(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(um)) allocate(um(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(zo)) allocate(zo(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(rb)) allocate(rb(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(zeta)) allocate(zeta(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(obu)) allocate(obu(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(tstr)) allocate(tstr(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(qstr)) allocate(qstr(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(zot)) allocate(zot(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(zoq)) allocate(zoq(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(thvstr)) allocate(thvstr(jxmin:jxmax,iymin:iymax))
      if (.not.allocated(facttq)) allocate(facttq(jxmin:jxmax,iymin:iymax))
!
!     Compute potential temperature and difference
!
      th = (t+tzero)*(d_1000/ps)**rovcp
      dth = t+0.0098d0*z-ts
!
!     Compute specific humidity at saturation using Tetens' formula
!
      qs = (1.0007d0+3.46d-6*ps)*6.1121d0*dexp(17.502d0*t/(240.97d0+ts))
      qs = qs*0.98d0
      qs = ep2*qs/(ps-0.378d0*qs) 
      dqh = q-qs
!
!     Compute virtual potential temperature and difference
!
      thv = th*(d_one+0.61d0*q)
      dthv = dth*(d_one+0.61d0*q) + 0.61d0*th*dqh
!
!     Compute density and latent heat of vap. (J/kg)
!
      rho = ps*d_100/(rgas*(ts+tzero)*(d_one+0.61d0*qs))
      xlv = (2.501d0-0.00237d0*ts)*1.0d+6
!
!     Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11
!
      visa = 1.326d-5*(d_one+6.542d-3*t+8.301d-6*t*t-4.84d-9*t*t*t)
!
!     Initial values of u* (fiction velocity) and convective velocity
!
      wc = d_half
      ustr = 0.06d0
      wspd = dsqrt(u**d_two+v**d_two)
      where ( dthv >= d_zero )
        um = dmax1(wspd, 0.1d0)
      elsewhere
        um = dsqrt(wspd*wspd+wc*wc)
      end where
!
!     Loop to obtain initial and good ustar and zo
!
      do i = 1, 5
        zo = 0.013d0*ustr*ustr*regrav+0.11d0*visa/ustr
        ustr = vonkar*um/dlog(z/zo)
      end do
!
      rb = egrav*z*dthv/(thv*um*um)
      where ( rb >= d_zero )       ! neutral or stable
        zeta = rb*dlog(z/zo)/(d_one-d_five*dmin1(rb,0.19d0))
        zeta = dmin1(d_two,dmax1(zeta,1.0d-6))
      elsewhere                    ! unstable
        zeta = rb*dlog(z/zo)
        zeta = dmax1(-d_100,dmin1(zeta,-1.0d-6))
      end where
      obu = z/zeta
!
!     Main iteration (2-10 iterations would be fine)
!
      do i = 1, 10
        call ocnrough (zo, zot, zoq, ustr, visa)
!
!       Wind
!
        zeta = z/obu
        where ( zeta < -zetam )        ! zeta < -1
          ustr = vonkar*um/(dlog(-zetam*obu/zo)-psi(1,-zetam)+          &
                  psi2d(1,zo/obu)+1.14D0*((-zeta)**onet-(zetam)**onet))
        elsewhere ( zeta < d_zero )    ! -1 <= zeta < 0
          ustr = vonkar*um/(dlog(z/zo)-psi2d(1,zeta)+psi2d(1,zo/obu))
        elsewhere ( zeta <= d_one )    !  0 <= zeta <= 1
          ustr = vonkar*um/(dlog(z/zo)+d_five*zeta-d_five*zo/obu)
        elsewhere                      !  1 < zeta, phi=5+zeta
          ustr = vonkar*um/(dlog(obu/zo)+d_five-d_five*zo/obu+          &
                  (d_five*dlog(zeta)+zeta-d_one))
        end where  
!
!       Temperature
!
        zeta = z/obu
        where ( zeta < -zetat )        ! zeta < -1
          tstr = vonkar*dth/                                            &
                 (dlog(-zetat*obu/zot)-psi(2,-zetat)+psi2d(2,zot/obu)+  &
                 0.8d0*((zetat)**(-onet)-(-zeta)**(-onet)))
        elsewhere ( zeta < d_zero )    ! -1 <= zeta < 0
          tstr = vonkar*dth/(dlog(z/zot)-                              &
                 psi2d(2,zeta)+psi2d(2,zot/obu))
        elsewhere ( zeta <= d_one )    !  0 <= zeta <= 1
          tstr = vonkar*dth/(dlog(z/zot)+d_five*zeta-d_five*zot/obu)
        elsewhere                      !  1 < zeta, phi=5+zeta
          tstr = vonkar*dth/(dlog(obu/zot)+d_five-d_five*zot/obu+       &
                 (d_five*dlog(zeta)+zeta-d_one))
        end where
!
!       Humidity
!
        zeta = z/obu
        where ( zeta < -zetat )        ! zeta < -1
          qstr = vonkar*dqh/                                            &
                 (dlog(-zetat*obu/zoq)-psi(2,-zetat)+psi2d(2,zoq/obu)+  &
                 0.8d0*((zetat)**(-onet)-(-zeta)**(-onet)))
        elsewhere ( zeta < d_zero )    ! -1 <= zeta < 0
          qstr = vonkar*dqh/(dlog(z/zoq)-                               &
                 psi2d(2,zeta)+psi2d(2,zoq/obu))
        elsewhere ( zeta <= d_one )    !  0 <= zeta <= 1
          qstr = vonkar*dqh/(dlog(z/zoq)+d_five*zeta-d_five*zoq/obu)
        elsewhere                      !  1 < zeta, phi=5+zeta
          qstr = vonkar*dqh/(dlog(obu/zoq)+d_five-d_five*zoq/obu+       &
                 (d_five*dlog(zeta)+zeta-d_one))
        end where
!
        thvstr = tstr*(d_one+0.61d0*q) + 0.61d0*th*qstr
!
        zeta = vonkar*egrav*thvstr*z/(ustr**d_two*thv)
        where ( zeta >= d_zero )   !neutral or stable
          um = dmax1(wspd, 0.1d0)
          zeta = dmin1(d_two,dmax1(zeta,1.0d-6))
        elsewhere                  !unstable
          wc = zbeta*(-egrav*ustr*thvstr*zi/thv)**onet
          um = dsqrt(wspd*wspd+wc*wc)
          zeta = dmax1(-d_100,dmin1(zeta,-1.0d-6))
        end where
        obu = z/zeta
      end do
!
!     Compute surface variables (adjust height of the lowest sigma lev.) 
!
      taux = rho*ustr*ustr*u995/um
      tauy = rho*ustr*ustr*u995/um
!
      zeta = z10/obu
      where ( zeta < d_zero )
        u10 = u+(ustr/vonkar)*(dlog(z10/z)-                             &
              (psi2d(1,zeta)-psi2d(1,z/obu)))  
        v10 = v+(ustr/vonkar)*(dlog(z10/z)-                             &
              (psi2d(1,zeta)-psi2d(1,z/obu)))
      elsewhere
        u10 = u+(ustr/vonkar)*(dlog(z10/z)+                             &
              d_five*zeta-d_five*z/obu)
        v10 = v+(ustr/vonkar)*(dlog(z10/z)+                             &
              d_five*zeta-d_five*z/obu)
      end where
!
      facttq = dlog(z995*d_half)/dlog(z995/zo)
      t2 = t+tzero-dth*facttq
      q2 = q-dqh*facttq
!      
      return
      end subroutine adjustvars
!
      end module mod_couplerr
