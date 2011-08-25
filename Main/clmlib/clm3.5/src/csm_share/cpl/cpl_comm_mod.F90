!===============================================================================
! SVN $Id: cpl_comm_mod.F90 3380 2007-03-06 05:42:19Z robj $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_comm_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_comm_mod -- Define MPI communication groups and model ID's
!
! !DESCRIPTION:
!     Sets up communicator groups and component ID's.
!     
!     A component ID (CID) is an integer identifying each component in the
!     coupled system.  Valid values are 1 to the total number of models
!     (including the Coupler).  Declaring an integer for each model
!     is a requirement of using MCT.
!
!     Use MPH to define the communicator groups and component ID's.
!
!     This module also declares and defines handy data wrt number of pe's, PID's, CID's
!     relative to world and component communicator groups.
!
! !REVISION HISTORY:
!     2001-Aug-20 - B. Kauffman - new naming convention
!     2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_comm_mod

! !USES:

  use mph_module     ! mph library

  use shr_sys_mod    ! system calls
  use shr_mpi_mod    ! mpi layer

  use cpl_kind_mod   ! kinds
  use mct_mod    ! mct library wrapper
  use cpl_fields_mod ! contains valid component name strings

  implicit none

  private ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_comm_init

! !PUBLIC DATA MEMBERS:

   integer(IN),public :: cpl_comm_wrld         ! = MPI_COMM_WORLD, global comm grp
   integer(IN),public :: cpl_comm_wrld_npe     ! number of pe's in MPI_COMM_WORLD
   integer(IN),public :: cpl_comm_wrld_pid     ! this comp pid in MPI_COMM_WORLD

   integer(IN),public :: cpl_comm_comp         ! this comp communicator group
   integer(IN),public :: cpl_comm_comp_npe     ! number of pe's in comp comm group
   integer(IN),public :: cpl_comm_comp_pid     ! this comp's pid in comp comm group

   integer(IN),public :: cpl_comm_mph_cid      ! MPH component ID, this component
   integer(IN),public :: cpl_comm_mph_cid_atm  ! MPH component ID, atm
   integer(IN),public :: cpl_comm_mph_cid_ice  ! MPH component ID, ice
   integer(IN),public :: cpl_comm_mph_cid_lnd  ! MPH component ID, lnd
   integer(IN),public :: cpl_comm_mph_cid_ocn  ! MPH component ID, ocn
   integer(IN),public :: cpl_comm_mph_cid_cpl  ! MPH component ID, cpl

   integer(IN),public :: cpl_comm_wrld_pe0     ! comm world pe0, this component
   integer(IN),public :: cpl_comm_wrld_pe0_atm ! comm world pe0, atm
   integer(IN),public :: cpl_comm_wrld_pe0_ice ! comm world pe0, ice
   integer(IN),public :: cpl_comm_wrld_pe0_lnd ! comm world pe0, lnd
   integer(IN),public :: cpl_comm_wrld_pe0_ocn ! comm world pe0, ocn
   integer(IN),public :: cpl_comm_wrld_pe0_cpl ! comm world pe0, cpl

!EOP

   save

   character(*),parameter :: modName = 'cpl_comm_mod'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_comm_init -- initialize the coupling/mpi environment.
!
! !DESCRIPTION:
!    This routine calls {\it MPI\_init} for the model with name {\tt name}
!    and returns an {\it MPI\_Communicator} {\tt comm} for use in the
!    calling model.
!    This also sets component ids, and processor ranks relative to world and 
!    component communicator groups.
! 
! !REMARKS:
!    Use cpl_interface_init which calls this routine.
!
! !REVISION HISTORY:
!     2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob -- first prototype
!     2001-Dec-10 - R. Jacob -- switch arguments in cpl_mct_world_init to
!                   to match new version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_comm_init(name,comm)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: name ! name of component name
   integer(IN) ,intent(out) :: comm ! communicator group for component

!EOP

   integer(IN)      :: n    ! generic loop index

   !----- formats -----
   character(*),parameter :: subname = "(cpl_comm_init) "
   character(*),parameter :: F00 = "('(cpl_comm_init) ',4a)"
   character(*),parameter :: F02 = "('(cpl_comm_init) ',a,6i4)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   write(6,F00) "setting up communicators, name = ",trim(name)

   !--- initialize my local comm group via MPH ---
#ifdef SINGLE_EXEC
   comm = mph_local_world(name)
#else
   comm = mph_components(name1=name)
#endif


   !--- query MPH for my local/global comm info, mph cid ---
   cpl_comm_wrld     = mph_global_world
   cpl_comm_wrld_pid = mph_global_proc_id()
   cpl_comm_wrld_npe = mph_global_totprocs()
   cpl_comm_wrld_pe0 = mph_global_id(name,0)

   cpl_comm_comp     = comm
   cpl_comm_mph_cid  = mph_comp_id(name)
   cpl_comm_comp_pid = mph_local_proc_id(cpl_comm_mph_cid)
   cpl_comm_comp_npe = mph_local_totprocs(cpl_comm_mph_cid)

   call shr_mpi_commsize(cpl_comm_comp,n,subName//" MPI comm size")
   write(6,F02) "cpl_comm_comp, size:",cpl_comm_comp,n

   !--- determine mph cid's and comm_world pe0's for all components ---
   do n=1,mph_total_components()
      if (mph_comp_name(n) == cpl_fields_atmname ) then
         cpl_comm_mph_cid_atm  = n
         cpl_comm_wrld_pe0_atm = mph_global_id(cpl_fields_atmname,0)
      elseif (mph_comp_name(n) == cpl_fields_icename ) then
         cpl_comm_mph_cid_ice  = n
         cpl_comm_wrld_pe0_ice = mph_global_id(cpl_fields_icename,0)
      elseif (mph_comp_name(n) == cpl_fields_lndname ) then
         cpl_comm_mph_cid_lnd  = n
         cpl_comm_wrld_pe0_lnd = mph_global_id(cpl_fields_lndname,0)
      elseif (mph_comp_name(n) == cpl_fields_ocnname ) then
         cpl_comm_mph_cid_ocn  = n
         cpl_comm_wrld_pe0_ocn = mph_global_id(cpl_fields_ocnname,0)
      elseif (mph_comp_name(n) == cpl_fields_cplname ) then
         cpl_comm_mph_cid_cpl  = n
         cpl_comm_wrld_pe0_cpl = mph_global_id(cpl_fields_cplname,0)
      else
         write(6,*) subName,'mph_component_name error',n,mph_comp_name(n)
         call shr_sys_abort(subName//'mph_component_name error')
      endif
   enddo

   !--- initialize MCT ---
   call mct_world_init(mph_total_components(),cpl_comm_wrld,cpl_comm_comp,cpl_comm_mph_cid)

   !--- document comm groups, pe0's, mph component ids ---
   write(6,F02) 'comm world    : comm,npe,pid   ', &
       &  cpl_comm_wrld,cpl_comm_wrld_npe,cpl_comm_wrld_pid
   write(6,F02) 'comm component: comm,npe,pid   ', &
       &  cpl_comm_comp,cpl_comm_comp_npe,cpl_comm_comp_pid
   write(6,F02) 'comm world pe0: atm,ice,lnd,ocn,cpl,me ', &
       &  cpl_comm_wrld_pe0_atm,cpl_comm_wrld_pe0_ice,cpl_comm_wrld_pe0_lnd, &
       &  cpl_comm_wrld_pe0_ocn,cpl_comm_wrld_pe0_cpl,cpl_comm_wrld_pe0
   write(6,F02) 'mph cid       : atm,ice,lnd,ocn,cpl,me ', &
       &  cpl_comm_mph_cid_atm ,cpl_comm_mph_cid_ice ,cpl_comm_mph_cid_lnd , &
       &  cpl_comm_mph_cid_ocn ,cpl_comm_mph_cid_cpl ,cpl_comm_mph_cid

end subroutine cpl_comm_init

!===============================================================================
!===============================================================================

end module cpl_comm_mod

