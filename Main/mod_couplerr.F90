      module mod_couplerr
#define ESMF_FILENAME "mod_couplerr.F90"
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

#define ESMF_METHOD "subroutine CheckError"
      subroutine CheckError(msg, localrc)
      implicit none
!
!***********************************************************************
!
!     Imported variable declarations 
!
!***********************************************************************
!
      character(*), intent(in) :: msg
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
                           file=ESMF_FILENAME,                          &
                           method=ESMF_METHOD,                          &
                           rc=rc2)
        call ESMF_Finalize(rc=rc3, endflag=ESMF_END_ABORT)
      end if
      end subroutine CheckError
!
      end module mod_couplerr
