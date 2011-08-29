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
!     Earth System Model (ESM) field type 
!-----------------------------------------------------------------------
!
      type ESM_Field
        character (len=40) :: name
        character (len=80) :: long_name
        character (len=80) :: units
        real*8, dimension(:,:), pointer :: field
      end type ESM_Field
!
!-----------------------------------------------------------------------
!     Earth System Model (ESM) data type 
!-----------------------------------------------------------------------
!
      type ESM_Data
        integer :: gridType
        type(ESM_Field) :: data 
        type(ESMF_Array) :: array         
      end type ESM_Data
!
!-----------------------------------------------------------------------
!     Earth System Model (ESM) grid type 
!-----------------------------------------------------------------------
!
      type ESM_Grid
        integer :: gridType
        type(ESM_Data) :: lat
        type(ESM_Data) :: lon
      end type ESM_Grid
!
!-----------------------------------------------------------------------
!     Earth System Model (ESM) high-level generic data type 
!-----------------------------------------------------------------------
!
      type ESM_Model 
!
!-----------------------------------------------------------------------
!       PETs (assigned processors to component)
!-----------------------------------------------------------------------
!
        integer, allocatable :: petList(:) 
!
!-----------------------------------------------------------------------
!       Main variables for component
!-----------------------------------------------------------------------
!
        type(ESMF_VM) :: vm 
        type(ESMF_GridComp) :: comp
        type(ESMF_State) :: stateExport
        type(ESMF_State) :: stateImport
        type(ESMF_DELayout) :: deLayout 
        type(ESMF_DistGrid) :: distGrid
        type(ESMF_ArraySpec) :: arrSpec
!
!-----------------------------------------------------------------------
!       Time, calendar and clock objects
!-----------------------------------------------------------------------
!
        type(ESMF_Calendar) :: cal
        type(ESMF_Time) :: refTime
        type(ESMF_Time) :: strTime
        type(ESMF_Time) :: endTime
        type(ESMF_Time) :: curTime
        type(ESMF_TimeInterval) :: dt
        type(ESMF_Clock) :: clock
        type(ESM_Time) :: time 
!
!-----------------------------------------------------------------------
!       Grid related model variables
!-----------------------------------------------------------------------
!
        integer :: nGrid
        type(ESM_Grid), allocatable :: grid(:)        

      end type ESM_Model
!
!***********************************************************************
!
!     Global module variables 
!
!***********************************************************************
!
      type(ESM_Model), save, allocatable, dimension(:) :: models 
!
!-----------------------------------------------------------------------
!     Number of gridded component or model (Atmosphere/Ocean)
!-----------------------------------------------------------------------
!
      integer :: nModels = 2
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
!     Coupling direction   
!-----------------------------------------------------------------------
!
      integer :: Iexport = 1
      integer :: Iimport = 2
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
      integer :: i, j, k
!
!-----------------------------------------------------------------------
!     Initialize the coupler variables
!-----------------------------------------------------------------------
!
      if (.not. allocated(models)) then
        allocate(models(nModels))
        do i = 1, nModels
          if (.not. allocated(models(i)%petList)) then
            allocate(models(i)%petList(petCount))
          end if
          do j = 1, petCount
            models(i)%petList(j) = j-1 
          end do 
        end do
      end if
!
!-----------------------------------------------------------------------
!     Grid variable (depend on model staggering type) 
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        if (i == Iatmos) then
          k = 2 ! cross and dot points (dot is not used actually)
        else if (i == Iocean) then
          k = 3 ! cross, u and v points
        end if
!
        if (.not. allocated(models(i)%grid)) then
          allocate(models(i)%grid(k))
          models(i)%nGrid = k
        end if
      end do
!
!
!-----------------------------------------------------------------------
!     Set mesh or grid information for each model 
!-----------------------------------------------------------------------
!
      do i = 1, nModels
        if (i == Iatmos) then
!
!         CROSS points (t, p etc.)
!
          models(i)%grid(1)%gridType = Icross
!
          models(i)%grid(1)%lon%data%name = 'xlon'
          models(i)%grid(1)%lon%data%long_name = 'longitude at cross'
          models(i)%grid(1)%lon%data%units = 'degrees_east'
!
          models(i)%grid(1)%lat%data%name = 'xlat'
          models(i)%grid(1)%lat%data%long_name = 'latitude at cross'
          models(i)%grid(1)%lat%data%units = 'degrees_north'
!
!         DOT points (u and v)
!
          models(i)%grid(2)%gridType = Idot
!
          models(i)%grid(2)%lon%data%name = 'dlon'
          models(i)%grid(2)%lon%data%long_name = 'longitude at dot'
          models(i)%grid(2)%lon%data%units = 'degrees_east'
!
          models(i)%grid(2)%lat%data%name = 'dlat'
          models(i)%grid(2)%lat%data%long_name = 'latitude at dot'
          models(i)%grid(2)%lat%data%units = 'degrees_north'
        else if (i == Iocean) then
!
!         RHO points
!
          models(i)%grid(1)%gridType = Icross
!
          models(i)%grid(1)%lon%data%name = 'lonr'
          models(i)%grid(1)%lon%data%long_name = 'longitude at rho'
          models(i)%grid(1)%lon%data%units = 'degrees_east'
!
          models(i)%grid(1)%lat%data%name = 'latr'
          models(i)%grid(1)%lat%data%long_name = 'latitude at rho'
          models(i)%grid(1)%lat%data%units = 'degrees_north'
!
!         U points
!
          models(i)%grid(2)%gridType = Iupoint
!
          models(i)%grid(2)%lon%data%name = 'lonu'
          models(i)%grid(2)%lon%data%long_name = 'longitude at u'
          models(i)%grid(2)%lon%data%units = 'degrees_east'
!
          models(i)%grid(2)%lat%data%name = 'latu'
          models(i)%grid(2)%lat%data%long_name = 'latitude at u'
          models(i)%grid(2)%lat%data%units = 'degrees_north'
!
!         V points
! 
          models(i)%grid(3)%gridType = Ivpoint
!         
          models(i)%grid(3)%lon%data%name = 'lonv'
          models(i)%grid(3)%lon%data%long_name = 'longitude at v'
          models(i)%grid(3)%lon%data%units = 'degrees_east'
!         
          models(i)%grid(3)%lat%data%name = 'latv'
          models(i)%grid(3)%lat%data%long_name = 'latitude at v'
          models(i)%grid(3)%lat%data%units = 'degrees_north'
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
      subroutine check_err(msg, fname, mname, localrc)
      implicit none
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      character(*), intent(in) :: msg
      character(*), intent(in) :: fname 
      character(*), intent(in) :: mname
      integer, intent(in) :: localrc
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
      integer :: rc, rc2, rc3
!
      if (ESMF_LogFoundError(localrc,                                   &
                             msg=trim(msg),                             &
                             rcToReturn=rc)) then
        call ESMF_LogWrite(trim(msg),                                   &
                           ESMF_LOGMSG_INFO,                            &
                           line=__LINE__,                               &
                           file=trim(fname),                            &
                           method=trim(mname),                          &
                           rc=rc2)
        call ESMF_Finalize(rc=rc3, endflag=ESMF_END_ABORT)
      end if
      end subroutine check_err 
!
      end module mod_couplerr
