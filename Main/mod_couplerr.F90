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
        character (len=40) :: name
        character (len=80) :: long_name
        character (len=80) :: units
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
        type(ESMF_VM) :: vm 
        integer :: comm
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
      ! no nested grid (not supported in this version)
      nNest = 1

      if (.not. allocated(models)) then
        ! check number of model greater than number of available pets
        if (nModels > petCount) then
          call abort_all(rc)
        end if          

        ! allocate user defined data types for models
        allocate(models(nModels))

        do i = 1, nModels
          ! allocate array to store pet list
          !if (i .eq. nModels) then
          !  nPets = (petCount/nModels)+mod(petCount, nModels)
          !else
          !  nPets = petCount/nModels
          !end if
          nPets = 2
          if (.not. allocated(models(i)%petList)) then
            allocate(models(i)%petList(nPets))
          end if

          ! assign pets to model (each component has its own pet) 
          ! nModels > petCount - not allowed 
          ! nModels < petCount - last component gets more pet
          ! nModels = petCount - each model gets equal number of pets
          !k = 0
          !do j = ((i-1)*petCount/nModels)+1, i*petCount/nModels
          !  k = k+1
          !  models(i)%petList(k) = j-1 
          !end do 
          !print*, i, models(i)%petList
        end do
          models(Iatmos)%petList(:) = (/ 0, 1 /)
          models(Iocean)%petList(:) = (/ 2, 3 /)
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
            models(i)%mesh(1,j)%name = 'xlon'
            models(i)%mesh(1,j)%long_name = 'longitude at cross'
            models(i)%mesh(1,j)%units = 'degrees_east'
            models(i)%mesh(1,j)%name = 'xlat'
            models(i)%mesh(1,j)%long_name = 'latitude at cross'
            models(i)%mesh(1,j)%units = 'degrees_north'
!           dot (or cell corners) points (u and v)
            models(i)%mesh(2,j)%gid = 2
            models(i)%mesh(2,j)%gtype = Idot
            models(i)%mesh(2,j)%name = 'dlon'
            models(i)%mesh(2,j)%long_name = 'longitude at dot'
            models(i)%mesh(2,j)%units = 'degrees_east'
            models(i)%mesh(2,j)%name = 'dlat'
            models(i)%mesh(2,j)%long_name = 'latitude at dot'
            models(i)%mesh(2,j)%units = 'degrees_north'
          end do
        else if (i == Iocean) then
          do j = 1, nNest(i)
!           rho (or cell center) points
            models(i)%mesh(1,j)%gid = 1
            models(i)%mesh(1,j)%gtype = Icross
            models(i)%mesh(1,j)%name = 'lonr'
            models(i)%mesh(1,j)%long_name = 'longitude at rho'
            models(i)%mesh(1,j)%units = 'degrees_east'
            models(i)%mesh(1,j)%name = 'latr'
            models(i)%mesh(1,j)%long_name = 'latitude at rho'
            models(i)%mesh(1,j)%units = 'degrees_north'
!           u points
            models(i)%mesh(2,j)%gid = 2
            models(i)%mesh(2,j)%gtype = Iupoint
            models(i)%mesh(2,j)%name = 'lonu'
            models(i)%mesh(2,j)%long_name = 'longitude at u'
            models(i)%mesh(2,j)%units = 'degrees_east'
            models(i)%mesh(2,j)%name = 'latu'
            models(i)%mesh(2,j)%long_name = 'latitude at u'
            models(i)%mesh(2,j)%units = 'degrees_north'
!           v points
            models(i)%mesh(3,j)%gid = 3
            models(i)%mesh(3,j)%gtype = Ivpoint
            models(i)%mesh(3,j)%name = 'lonv'
            models(i)%mesh(3,j)%long_name = 'longitude at v'
            models(i)%mesh(3,j)%units = 'degrees_east'
            models(i)%mesh(3,j)%name = 'latv'
            models(i)%mesh(3,j)%long_name = 'latitude at v'
            models(i)%mesh(3,j)%units = 'degrees_north'
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
      function check_err(error_flag, source, routine, message)
      implicit none
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      integer, intent(in) :: error_flag 
      character(*), intent(in) :: source
      character(*), intent(in) :: routine
      character(*), intent(in) :: message
!
!***********************************************************************
!
!     Local variable declarations 
!
!***********************************************************************
!
      logical :: check_err
      integer :: rc, rc2, rc3
!
!-----------------------------------------------------------------------
!     Check error flag
!-----------------------------------------------------------------------
!
      check_err = .false.
      if (error_flag .ne. ESMF_SUCCESS) then
        check_err = .true.
      end if
!
!-----------------------------------------------------------------------
!     Report error message
!-----------------------------------------------------------------------
!


!      if (ESMF_LogFoundError(localrc,                                   &
!                             msg=trim(msg),                             &
!                             rcToReturn=rc)) then
!        call ESMF_LogWrite(trim(msg),                                   &
!                           ESMF_LOGMSG_INFO,                            &
!                           line=__LINE__,                               &
!                           file=trim(fname),                            &
!                           method=trim(mname),                          &
!                           rc=rc2)
!        call ESMF_Finalize(rc=rc3, endflag=ESMF_END_ABORT)
!      end if
      end function check_err 
!
      subroutine abort_all(rc)
      implicit none
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      integer, intent(inout) :: rc 
!
!-----------------------------------------------------------------------
!     Terminate execution due to fatal error 
!-----------------------------------------------------------------------
!      
      call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      stop
      end subroutine abort_all
!
      end module mod_couplerr
