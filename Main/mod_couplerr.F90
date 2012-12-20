      module mod_couplerr
!
!-----------------------------------------------------------------------
!     Imported modules 
!-----------------------------------------------------------------------
!
      use ESMF
!
      use mod_intkinds, only : ik4
      use mod_realkinds, only : rk8
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
        integer(ik4) :: year                 
        integer(ik4) :: month                
        integer(ik4) :: day                  
        integer(ik4) :: hour                 
        integer(ik4) :: minute               
        integer(ik4) :: second               
        integer(ik4) :: yday              
        integer(ik4) :: zone             
        character (len=30) :: stamp
      end type ESM_Time
!
!-----------------------------------------------------------------------
!     Earth System Model (ESM) field data type 
!-----------------------------------------------------------------------
!
      type ESM_Field
        integer(ik4) :: fid
        integer(ik4) :: gtype        
        integer(ik4) :: itype
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
        integer(ik4) :: gid
        integer(ik4) :: gtype
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
        integer(ik4) :: mid 
        integer(ik4) :: comm
        type(ESMF_VM) :: vm 
        integer(ik4) :: nproc
        integer, allocatable :: petList(:) 
        type(ESMF_GridComp) :: comp
        type(ESM_Mesh), allocatable :: mesh(:,:)        
        type(ESMF_Grid), allocatable :: grid(:)
        type(ESMF_DELayout), allocatable :: deLayout(:) 
        type(ESMF_DistGrid), allocatable :: distGrid(:)
        type(ESMF_ArraySpec), allocatable :: arrSpec(:)
        integer(ik4) :: nExport
        integer(ik4) :: nImport
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
      integer(ik4) :: nNest(nModels)
!
!-----------------------------------------------------------------------
!     Gridded model indices
!-----------------------------------------------------------------------
!
      integer(ik4) :: Iatmos  = 1
      integer(ik4) :: Iocean  = 2
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
      integer(ik4) :: Icross  = 1
      integer(ik4) :: Idot    = 2
      integer(ik4) :: Iupoint = 3
      integer(ik4) :: Ivpoint = 4
!
!-----------------------------------------------------------------------
!     Interpolation type        
!-----------------------------------------------------------------------
!
      integer(ik4) :: Ibilin = 1 
      integer(ik4) :: Iconsv = 2
!
!-----------------------------------------------------------------------
!     Variables for coupling (direction, import and export variables)  
!-----------------------------------------------------------------------
!
      integer(ik4) :: DIRECTION
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
      character(ESMF_MAXSTR) :: config_fname="regcm.rc"
!
      integer(ik4) :: cpl_dtsec, cpl_dbglevel
      integer(ik4) :: ibulk, iregrid
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
!
!     FB - Forward-Bilinear
!     FC - Forward-Conservative
!     BB - Backward-Bilinear
!     BC - Backward-Conservative
!
!     FS - Forward-Source
!     FD - Forward-Destination
!     BS - Backward-Source
!     BD - Backward-Destination
!-----------------------------------------------------------------------
!
      type(ESMF_RouteHandle), allocatable, dimension(:,:) :: routeHandleFB
      type(ESMF_RouteHandle), allocatable, dimension(:,:) :: routeHandleFC
      type(ESMF_RouteHandle), allocatable, dimension(:,:) :: routeHandleBB
      type(ESMF_RouteHandle), allocatable, dimension(:,:) :: routeHandleBC
      type(ESMF_Field) :: fracFieldFS       
      type(ESMF_Field) :: fracFieldFD       
      type(ESMF_Field) :: fracFieldBS       
      type(ESMF_Field) :: fracFieldBD       
!
!-----------------------------------------------------------------------
!     Constants
!     Cp  : Specific heat for seawater (Joules/Kg/degC)
!     rho0: Mean density (Kg/m3) - Boussinesq approximation
!-----------------------------------------------------------------------
!
      real*8, parameter :: Cp = 3985.0d0
      real*8, parameter :: rho0 = 1025.0d0
      real*8, parameter :: Hscale = rho0*Cp
      real*8, parameter :: Hscale2 = 1.0d0/(rho0*Cp)
      real*8, parameter :: day2s = 1.0d0/86400.0d0
      real*8, parameter :: mm2m = 1.0d0/1000.0d0
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
      integer(ik4) :: i, j, k, nPets, petNum1, petNum2
      character(100) :: fmt_123 
      logical :: file_exists
!
      type(ESMF_Config) :: cf
!
!-----------------------------------------------------------------------
!     Initialize the coupler variables
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!     Number of nested grid or each gridded component
!     It is not supported yet !!! Set ones to all ... 
!-----------------------------------------------------------------------
!
      nNest = 1
!
!-----------------------------------------------------------------------
!     Allocate ESM_Model to store information about gridded comp. and
!     distribute the PETs across components 
!-----------------------------------------------------------------------
!
      if (.not. allocated(models)) then
        ! check for number of PETs
        ! number of model must be less or equal than number of PETs
        if (nModels > petCount) then
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
        end if          

        ! allocate user defined data types for models
        allocate(models(nModels))

        ! check configuration file exist or not?
        inquire(file=trim(config_fname), exist=file_exists)

        if (file_exists) then
          cf = ESMF_ConfigCreate(rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if

          call ESMF_ConfigLoadFile(cf, trim(config_fname), rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if

          ! read PETs option
          call ESMF_ConfigFindLabel(cf, 'PETs:', rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
          call ESMF_ConfigGetAttribute(cf, petNum1, rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
          call ESMF_ConfigGetAttribute(cf, petNum2, rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if

          ! check PETs option
          if (petCount .ne. (petNum1+petNum2)) then
            write(*,*) "Number of PETs must be consistent with regcm.rc"
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if

          ! read Bulk option
          call ESMF_ConfigFindLabel(cf, 'BulkFlux:', rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
          call ESMF_ConfigGetAttribute(cf, ibulk, rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if

          ! read interpolation type option
          call ESMF_ConfigFindLabel(cf, 'ConservativeRegrid:', rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
          call ESMF_ConfigGetAttribute(cf, iregrid, rc=rc)
          if (rc /= ESMF_SUCCESS) then
            call ESMF_Finalize(endflag=ESMF_END_ABORT)
          end if
        end if 

        ! assign PETs to giridded components
        do i = 1, nModels
          ! allocate array to store PET list
          if (file_exists) then
            if (i .eq. Iatmos) then 
              nPets = petNum1
            else
              nPets = petNum2
            end if
          else
            if (i .eq. Iatmos) then
              nPets = (petCount/nModels)+mod(petCount, nModels)
            else
              nPets = petCount/nModels
            end if
          end if

          if (.not. allocated(models(i)%petList)) then
            allocate(models(i)%petList(nPets))
          end if
          ! set number of processor
          models(i)%nproc = nPets

          ! assign PETs to model (each component has its own PET) 
          ! For example; two models and seven cpu (or PET)
          ! model a - 0, 2, 4, 6
          ! model b - 1, 3, 5
          !do j = 1, nPets 
          !  models(i)%petList(j) = (i+(j-1)*nModels)-1
          !end do

          ! For example; two models and seven cpu (or PET)
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
!     Set number of import and export fields 
!-----------------------------------------------------------------------
!
      if (ibulk == 1) then
        models(Iatmos)%nExport = 11
        models(Iatmos)%nImport = 2
        models(Iocean)%nExport = 2
        models(Iocean)%nImport = 11
      else
        models(Iatmos)%nExport = 5
        models(Iatmos)%nImport = 2
        models(Iocean)%nExport = 2
        models(Iocean)%nImport = 5
      end if
!
!-----------------------------------------------------------------------
!     Allocate mesh variable (depend on model staggering type) 
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        if (i == Iatmos) then
          k = 2 ! cross and dot points are used (B Grid)
        else if (i == Iocean) then 
          k = 4 
        end if
!
        if (.not. allocated(models(i)%mesh)) then
          allocate(models(i)%mesh(k,nNest(i)))
        end if
!
        if (.not. allocated(models(i)%grid)) then
          allocate(models(i)%grid(nNest(i)))
        end if
      end do
!
!-----------------------------------------------------------------------
!     Allocate routehandle 
!-----------------------------------------------------------------------
!
      if (.not. allocated(routeHandleFB)) then
        allocate(routeHandleFB(ubound(models(Iatmos)%mesh, dim=1),      &
                               ubound(models(Iocean)%mesh, dim=1)))
        allocate(routeHandleFC(ubound(models(Iatmos)%mesh, dim=1),      &
                               ubound(models(Iocean)%mesh, dim=1)))
        allocate(routeHandleBB(ubound(models(Iatmos)%mesh, dim=1),      &
                               ubound(models(Iocean)%mesh, dim=1)))
        allocate(routeHandleBC(ubound(models(Iatmos)%mesh, dim=1),      &
                               ubound(models(Iocean)%mesh, dim=1)))
      end if
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
!           psi points
            models(i)%mesh(2,j)%gid = 2
            models(i)%mesh(2,j)%gtype = Idot
!
            models(i)%mesh(2,j)%lon%gtype = Idot
            models(i)%mesh(2,j)%lon%name = 'lonp'
            models(i)%mesh(2,j)%lon%long_name = 'longitude at psi'
            models(i)%mesh(2,j)%lon%units = 'degrees_east'
!
            models(i)%mesh(2,j)%lat%gtype = Idot
            models(i)%mesh(2,j)%lat%name = 'latp'
            models(i)%mesh(2,j)%lat%long_name = 'latitude at psi'
            models(i)%mesh(2,j)%lat%units = 'degrees_north'
!
            models(i)%mesh(2,j)%mask%gtype = Idot
            models(i)%mesh(2,j)%mask%name = 'mask_psi'
            models(i)%mesh(2,j)%mask%long_name = 'mask on psi'
            models(i)%mesh(2,j)%mask%units = '1'
!
!           u points
            models(i)%mesh(3,j)%gid = 3
            models(i)%mesh(3,j)%gtype = Iupoint

            models(i)%mesh(3,j)%lon%gtype = Iupoint
            models(i)%mesh(3,j)%lon%name = 'lonu'
            models(i)%mesh(3,j)%lon%long_name = 'longitude at u'
            models(i)%mesh(3,j)%lon%units = 'degrees_east'

            models(i)%mesh(3,j)%lat%gtype = Iupoint
            models(i)%mesh(3,j)%lat%name = 'latu'
            models(i)%mesh(3,j)%lat%long_name = 'latitude at u'
            models(i)%mesh(3,j)%lat%units = 'degrees_north'

            models(i)%mesh(3,j)%mask%gtype = Iupoint
            models(i)%mesh(3,j)%mask%name = 'mask_u'
            models(i)%mesh(3,j)%mask%long_name = 'mask on u'
            models(i)%mesh(3,j)%mask%units = '1'

!           v points
            models(i)%mesh(4,j)%gid = 4
            models(i)%mesh(4,j)%gtype = Ivpoint
!
            models(i)%mesh(4,j)%lon%gtype = Ivpoint
            models(i)%mesh(4,j)%lon%name = 'lonv'
            models(i)%mesh(4,j)%lon%long_name = 'longitude at v'
            models(i)%mesh(4,j)%lon%units = 'degrees_east'
!
            models(i)%mesh(4,j)%lat%gtype = Ivpoint
            models(i)%mesh(4,j)%lat%name = 'latv'
            models(i)%mesh(4,j)%lat%long_name = 'latitude at v'
            models(i)%mesh(4,j)%lat%units = 'degrees_north'
!
            models(i)%mesh(4,j)%mask%gtype = Ivpoint
            models(i)%mesh(4,j)%mask%name = 'mask_v'
            models(i)%mesh(4,j)%mask%long_name = 'mask on v'
            models(i)%mesh(4,j)%mask%units = '1'
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
!       Note: Only import fields need scale and offset
!-----------------------------------------------------------------------
!
        if (.not. allocated(models(i)%dataExport)) then
          allocate(models(i)%dataExport(models(i)%nExport,nNest(i)))
        end if
        if (.not. allocated(models(i)%dataImport)) then
          allocate(models(i)%dataImport(models(i)%nImport,nNest(i)))
        end if
!
        if (i == Iatmos) then
          do j = 1, nNest(i)
            models(i)%dataImport(1,j)%fid = 1
            models(i)%dataImport(1,j)%gtype = Icross
            models(i)%dataImport(1,j)%itype = Ibilin
            models(i)%dataImport(1,j)%name = 'SST'
            models(i)%dataImport(1,j)%long_name = &
            'Sea Surface Temperature'
            models(i)%dataImport(1,j)%units = 'Kelvin'
            models(i)%dataImport(1,j)%scale_factor = 1.0d0
            models(i)%dataImport(1,j)%add_offset = 273.16d0
!
            models(i)%dataImport(2,j)%fid = 2
            models(i)%dataImport(2,j)%gtype = Icross
            models(i)%dataImport(2,j)%itype = Ibilin
            models(i)%dataImport(2,j)%name = 'Hice'
            models(i)%dataImport(2,j)%long_name = &
            'Average Ice Thickness in Cell'
            models(i)%dataImport(2,j)%units = 'mm'
            models(i)%dataImport(2,j)%scale_factor = 1.0d3
            models(i)%dataImport(2,j)%add_offset = 0.0d0
!
            if (ibulk == 1) then
!
            models(i)%dataExport(1,j)%fid = 1
            models(i)%dataExport(1,j)%gtype = Icross
            models(i)%dataExport(1,j)%itype = Ibilin
            models(i)%dataExport(1,j)%name = 'Pair'
            models(i)%dataExport(1,j)%long_name = 'Surface Pressure'
            models(i)%dataExport(1,j)%units = 'Pascal'
!
            models(i)%dataExport(2,j)%fid = 2
            models(i)%dataExport(2,j)%gtype = Icross
            models(i)%dataExport(2,j)%itype = Ibilin
            models(i)%dataExport(2,j)%name = 'Tair'
            models(i)%dataExport(2,j)%long_name = &
            'Surface Air Temperature'
            models(i)%dataExport(2,j)%units = 'Kelvin'             
!
            models(i)%dataExport(3,j)%fid = 3
            models(i)%dataExport(3,j)%gtype = Icross
            models(i)%dataExport(3,j)%itype = Ibilin
            models(i)%dataExport(3,j)%name = 'Qair'
            models(i)%dataExport(3,j)%long_name = &
            'Surface Air Specific Humidity'
            models(i)%dataExport(3,j)%units = 'kg/kg'
!
            models(i)%dataExport(4,j)%fid = 4
            models(i)%dataExport(4,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(4,j)%itype = Ibilin
            else
            models(i)%dataExport(4,j)%itype = Iconsv
            end if 
            models(i)%dataExport(4,j)%name = 'Swrad'
            models(i)%dataExport(4,j)%long_name = &
            'solar shortwave radiation flux'
            models(i)%dataExport(4,j)%units = 'watt meter-2'
!
            models(i)%dataExport(5,j)%fid = 5
            models(i)%dataExport(5,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(5,j)%itype = Ibilin
            else
            models(i)%dataExport(5,j)%itype = Iconsv
            end if
            models(i)%dataExport(5,j)%name = 'Lwrad'
            models(i)%dataExport(5,j)%long_name = &
            'net longwave radiation flux'
            models(i)%dataExport(5,j)%units = 'watt meter-2'
!
            models(i)%dataExport(6,j)%fid = 6
            models(i)%dataExport(6,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(6,j)%itype = Ibilin
            else
            models(i)%dataExport(6,j)%itype = Iconsv
            end if
            models(i)%dataExport(6,j)%name = 'Lwrad_down'
            models(i)%dataExport(6,j)%long_name = &
            'downwelling longwave radiation flux'
            models(i)%dataExport(6,j)%units = 'watt meter-2'
!
            models(i)%dataExport(7,j)%fid = 7
            models(i)%dataExport(7,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(7,j)%itype = Ibilin
            else
            models(i)%dataExport(7,j)%itype = Iconsv
            end if
            models(i)%dataExport(7,j)%name = 'Lhflx'
            models(i)%dataExport(7,j)%long_name = &
            'latent heat flux'
            models(i)%dataExport(7,j)%units = 'watt meter-2'
!
            models(i)%dataExport(8,j)%fid = 8
            models(i)%dataExport(8,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(8,j)%itype = Ibilin
            else
            models(i)%dataExport(8,j)%itype = Iconsv
            end if
            models(i)%dataExport(8,j)%name = 'Shflx'
            models(i)%dataExport(8,j)%long_name = &
            'sensible heat flux'
            models(i)%dataExport(8,j)%units = 'watt meter-2'
!
            models(i)%dataExport(9,j)%fid = 9
            models(i)%dataExport(9,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(9,j)%itype = Ibilin
            else
            models(i)%dataExport(9,j)%itype = Iconsv
            end if
            models(i)%dataExport(9,j)%name = 'Rain'
            models(i)%dataExport(9,j)%long_name = &
            'rain fall rate'
            models(i)%dataExport(9,j)%units ='kilogram meter-2 second-1'
!
            models(i)%dataExport(10,j)%fid = 10
            models(i)%dataExport(10,j)%gtype = Icross
            models(i)%dataExport(10,j)%itype = Ibilin
            models(i)%dataExport(10,j)%name = 'Uwind'
            models(i)%dataExport(10,j)%long_name = &
            'surface u-wind component'
            models(i)%dataExport(10,j)%units = 'meter second-1'
!
            models(i)%dataExport(11,j)%fid = 11
            models(i)%dataExport(11,j)%gtype = Icross
            models(i)%dataExport(11,j)%itype = Ibilin
            models(i)%dataExport(11,j)%name = 'Vwind'
            models(i)%dataExport(11,j)%long_name = &
            'surface v-wind component'
            models(i)%dataExport(11,j)%units = 'meter second-1'
!
            else
!
            models(i)%dataExport(1,j)%fid = 1
            models(i)%dataExport(1,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(1,j)%itype = Ibilin
            else
            models(i)%dataExport(1,j)%itype = Iconsv
            end if
            models(i)%dataExport(1,j)%name = 'EminP'
            models(i)%dataExport(1,j)%long_name = &
            'surface freshwater (E-P) flux'
            models(i)%dataExport(1,j)%units ='kilogram meter-2 second-1'
!
            models(i)%dataExport(2,j)%fid = 2
            models(i)%dataExport(2,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(2,j)%itype = Ibilin
            else
            models(i)%dataExport(2,j)%itype = Iconsv
            end if
            models(i)%dataExport(2,j)%name = 'NHeat'
            models(i)%dataExport(2,j)%long_name = &
            'surface net heat flux'
            models(i)%dataExport(2,j)%units = 'watt meter-2'
!
            models(i)%dataExport(3,j)%fid = 3
            models(i)%dataExport(3,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataExport(3,j)%itype = Ibilin
            else
            models(i)%dataExport(3,j)%itype = Iconsv
            end if
            models(i)%dataExport(3,j)%name = 'Swrad'
            models(i)%dataExport(3,j)%long_name = &
            'solar shortwave radiation flux'
            models(i)%dataExport(3,j)%units = 'watt meter-2'
!
            models(i)%dataExport(4,j)%fid = 4
            models(i)%dataExport(4,j)%gtype = Icross
            models(i)%dataExport(4,j)%itype = Ibilin
            models(i)%dataExport(4,j)%name = 'Ustr'
            models(i)%dataExport(4,j)%long_name = &
            'u component of surface momentum flux'
            models(i)%dataExport(4,j)%units = 'm2/s2'
!
            models(i)%dataExport(5,j)%fid = 5
            models(i)%dataExport(5,j)%gtype = Icross
            models(i)%dataExport(5,j)%itype = Ibilin
            models(i)%dataExport(5,j)%name = 'Vstr'
            models(i)%dataExport(5,j)%long_name = &
            'v component of surface momentum flux'
            models(i)%dataExport(5,j)%units = 'm2/s2'
            end if
          end do 
        else if (i == Iocean) then
          do j = 1, nNest(i)
            models(i)%dataExport(1,j)%fid = 1
            models(i)%dataExport(1,j)%gtype = Icross
            models(i)%dataExport(1,j)%itype = Ibilin
            models(i)%dataExport(1,j)%name = 'SST'
            models(i)%dataExport(1,j)%long_name = & 
           'Sea Surface Temperature'
            models(i)%dataExport(1,j)%units = 'Celsius'
!
            models(i)%dataExport(2,j)%fid = 2
            models(i)%dataExport(2,j)%gtype = Icross
            models(i)%dataExport(2,j)%itype = Ibilin
            models(i)%dataExport(2,j)%name = 'Hice'
            models(i)%dataExport(2,j)%long_name = &
            'Average Ice Thickness in Cell'
            models(i)%dataExport(2,j)%units = 'm'
!
            if (ibulk == 1) then
!
            models(i)%dataImport(1,j)%fid = 1
            models(i)%dataImport(1,j)%gtype = Icross
            models(i)%dataImport(1,j)%itype = Ibilin
            models(i)%dataImport(1,j)%name = 'Pair'
            models(i)%dataImport(1,j)%long_name = 'Surface Pressure'
            models(i)%dataImport(1,j)%units = 'hPa'
            models(i)%dataImport(1,j)%scale_factor = 1.0d0
            models(i)%dataImport(1,j)%add_offset = 0.0d0
!
            models(i)%dataImport(2,j)%fid = 2
            models(i)%dataImport(2,j)%gtype = Icross
            models(i)%dataImport(2,j)%itype = Ibilin
            models(i)%dataImport(2,j)%name = 'Tair'
            models(i)%dataImport(2,j)%long_name = &
            'Surface Air Temperature'
            models(i)%dataImport(2,j)%units = 'Celsius'
            models(i)%dataImport(2,j)%scale_factor = 1.0d0
            models(i)%dataImport(2,j)%add_offset = -273.16d0
!
            models(i)%dataImport(3,j)%fid = 3
            models(i)%dataImport(3,j)%gtype = Icross
            models(i)%dataImport(3,j)%itype = Ibilin
            models(i)%dataImport(3,j)%name = 'Qair'
            models(i)%dataImport(3,j)%long_name = &
            'Surface Air Specific Humidity'
!           SPECIFIC_HUMIDITY is defined in ROMS configuration 
            models(i)%dataImport(3,j)%units = 'kg/kg'
            models(i)%dataImport(3,j)%scale_factor = 1.0d0
!           SPECIFIC_HUMIDITY is not defined !!!
!            models(i)%dataImport(3,j)%units = 'g/kg'
!            models(i)%dataImport(3,j)%scale_factor = 1.0d3
            models(i)%dataImport(3,j)%add_offset = 0.0d0
!
            models(i)%dataImport(4,j)%fid = 4
            models(i)%dataImport(4,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(4,j)%itype = Ibilin
            else
            models(i)%dataImport(4,j)%itype = Iconsv
            end if
            models(i)%dataImport(4,j)%name = 'Swrad'
            models(i)%dataImport(4,j)%long_name = &
            'solar shortwave radiation flux'
            models(i)%dataImport(4,j)%units = 'Celsius m/s'
            models(i)%dataImport(4,j)%scale_factor = Hscale2 
            models(i)%dataImport(4,j)%add_offset = 0.0d0
!
            models(i)%dataImport(5,j)%fid = 5
            models(i)%dataImport(5,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(5,j)%itype = Ibilin
            else
            models(i)%dataImport(5,j)%itype = Iconsv
            end if
            models(i)%dataImport(5,j)%name = 'Lwrad'
            models(i)%dataImport(5,j)%long_name = &
            'net longwave radiation flux'
            models(i)%dataImport(5,j)%units = 'Celsius m/s'
            models(i)%dataImport(5,j)%scale_factor = -Hscale2
            models(i)%dataImport(5,j)%add_offset = 0.0d0
!
            models(i)%dataImport(6,j)%fid = 6
            models(i)%dataImport(6,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(6,j)%itype = Ibilin
            else
            models(i)%dataImport(6,j)%itype = Iconsv
            end if
            models(i)%dataImport(6,j)%name = 'Lwrad_down'
            models(i)%dataImport(6,j)%long_name = &
            'downwelling longwave radiation flux'
            models(i)%dataImport(6,j)%units = 'Celsius m/s'
            models(i)%dataImport(6,j)%scale_factor = Hscale2
            models(i)%dataImport(6,j)%add_offset = 0.0d0
!
            models(i)%dataImport(7,j)%fid = 7
            models(i)%dataImport(7,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(7,j)%itype = Ibilin
            else
            models(i)%dataImport(7,j)%itype = Iconsv
            end if
            models(i)%dataImport(7,j)%name = 'Lhflx'
            models(i)%dataImport(7,j)%long_name = &
            'latent heat flux'
            models(i)%dataImport(7,j)%units = 'Celsius m/s'
            models(i)%dataImport(7,j)%scale_factor = Hscale2
            models(i)%dataImport(7,j)%add_offset = 0.0d0
!
            models(i)%dataImport(8,j)%fid = 8
            models(i)%dataImport(8,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(8,j)%itype = Ibilin
            else
            models(i)%dataImport(8,j)%itype = Iconsv
            end if
            models(i)%dataImport(8,j)%name = 'Shflx'
            models(i)%dataImport(8,j)%long_name = &
            'sensible heat flux'
            models(i)%dataImport(8,j)%units = 'Celsius m/s'
            models(i)%dataImport(8,j)%scale_factor = Hscale2
            models(i)%dataImport(8,j)%add_offset = 0.0d0
!
            models(i)%dataImport(9,j)%fid = 9
            models(i)%dataImport(9,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(9,j)%itype = Ibilin
            else
            models(i)%dataImport(9,j)%itype = Iconsv
            end if
            models(i)%dataImport(9,j)%name = 'Rain'
            models(i)%dataImport(9,j)%long_name = &
            'rain fall rate'
            models(i)%dataImport(9,j)%units ='kilogram meter-2 second-1'
            models(i)%dataImport(9,j)%scale_factor = 1.0d0
            models(i)%dataImport(9,j)%add_offset = 0.0d0
!
            models(i)%dataImport(10,j)%fid = 10
            models(i)%dataImport(10,j)%gtype = Icross
            models(i)%dataImport(10,j)%itype = Ibilin
            models(i)%dataImport(10,j)%name = 'Uwind'
            models(i)%dataImport(10,j)%long_name = &
            'surface u-wind component'
            models(i)%dataImport(10,j)%units = 'meter second-1'
            models(i)%dataImport(10,j)%scale_factor = 1.0d0
            models(i)%dataImport(10,j)%add_offset = 0.0d0
!
            models(i)%dataImport(11,j)%fid = 11
            models(i)%dataImport(11,j)%gtype = Icross
            models(i)%dataImport(11,j)%itype = Ibilin
            models(i)%dataImport(11,j)%name = 'Vwind'
            models(i)%dataImport(11,j)%long_name = &
            'surface v-wind component'
            models(i)%dataImport(11,j)%units = 'meter second-1'
            models(i)%dataImport(11,j)%scale_factor = 1.0d0
            models(i)%dataImport(11,j)%add_offset = 0.0d0
!
            else
!
            models(i)%dataImport(1,j)%fid = 1
            models(i)%dataImport(1,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(1,j)%itype = Ibilin
            else
            models(i)%dataImport(1,j)%itype = Iconsv
            end if
            models(i)%dataImport(1,j)%name = 'EminP'
            models(i)%dataImport(1,j)%long_name = &
            'surface freshwater (E-P) flux'
            models(i)%dataImport(1,j)%units = 'm/s'
            models(i)%dataImport(1,j)%scale_factor = mm2m 
            models(i)%dataImport(1,j)%add_offset = 0.0d0
!
            models(i)%dataImport(2,j)%fid = 2
            models(i)%dataImport(2,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(2,j)%itype = Ibilin
            else
            models(i)%dataImport(2,j)%itype = Iconsv
            end if
            models(i)%dataImport(2,j)%name = 'NHeat'
            models(i)%dataImport(2,j)%long_name = &
            'surface net heat flux'
            models(i)%dataImport(2,j)%units = 'degC m/s'
            models(i)%dataImport(2,j)%scale_factor = Hscale2
            models(i)%dataImport(2,j)%add_offset = 0.0d0
!
            models(i)%dataImport(3,j)%fid = 3
            models(i)%dataImport(3,j)%gtype = Icross
            if (iregrid == 0) then
            models(i)%dataImport(3,j)%itype = Ibilin
            else
            models(i)%dataImport(3,j)%itype = Iconsv
            end if
            models(i)%dataImport(3,j)%name = 'Swrad'
            models(i)%dataImport(3,j)%long_name = &
            'solar shortwave radiation flux'
            models(i)%dataImport(3,j)%units = 'degC m/s'
            models(i)%dataImport(3,j)%scale_factor = Hscale2
            models(i)%dataImport(3,j)%add_offset = 0.0d0
!
            models(i)%dataImport(4,j)%fid = 4
            models(i)%dataImport(4,j)%gtype = Iupoint !Icross
            models(i)%dataImport(4,j)%itype = Ibilin
            models(i)%dataImport(4,j)%name = 'Ustr'
            models(i)%dataImport(4,j)%long_name = &
            'u component of surface momentum flux'
            models(i)%dataImport(4,j)%units = 'Pa'
            models(i)%dataImport(4,j)%scale_factor = 1.0d0/rho0
            models(i)%dataImport(4,j)%add_offset = 0.0d0
!
            models(i)%dataImport(5,j)%fid = 5
            models(i)%dataImport(5,j)%gtype = Ivpoint !Icross
            models(i)%dataImport(5,j)%itype = Ibilin
            models(i)%dataImport(5,j)%name = 'Vstr'
            models(i)%dataImport(5,j)%long_name = &
            'v component of surface momentum flux'
            models(i)%dataImport(5,j)%units = 'Pa'
            models(i)%dataImport(5,j)%scale_factor = 1.0d0/rho0
            models(i)%dataImport(5,j)%add_offset = 0.0d0
            end if
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
      integer(ik4) :: iarr(6)
      integer(ik4) :: i, nitems, petCount, localPet, comm, mysec, rc
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
      integer function getVarID(list, name)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      type(ESM_Field), intent(in) :: list(:)
      character (len=*) :: name
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      integer(ik4) :: i
!     
!-----------------------------------------------------------------------
!     Find index of specified field
!-----------------------------------------------------------------------
!
      do i = 1, size(list, dim=1)
        if (trim(list(i)%name) == trim(name)) then
          getVarID = i
          return
        end if
      end do
      getVarID = -1 
      return
      end function getVarID
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
      integer(ik4) :: i
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
      getMeshID = -1 
      return
      end function getMeshID
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
      integer(ik4) :: i, j, imin, imax, jmin, jmax
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
      subroutine print_size_r8 (field, localPet, header)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations 
!-----------------------------------------------------------------------
!
      real(rk8), intent(in) :: field(:,:) 
      integer, intent(in) :: localPet
      character(len=*), intent(in) :: header
!
!-----------------------------------------------------------------------
!     Local variable declarations 
!-----------------------------------------------------------------------
!
      write(*,40) localPet, trim(adjustl(header)),                      &
                  lbound(field, dim=1), ubound(field, dim=1),           &
                  lbound(field, dim=2), ubound(field, dim=2)
!
!-----------------------------------------------------------------------
!     Formats 
!-----------------------------------------------------------------------
!
 40   format(" PET(",I3,") - ", A20, " : ", 4I8)
      end subroutine print_size_r8

      end module mod_couplerr
