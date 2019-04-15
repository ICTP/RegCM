!===============================================================================
! SVN $Id: cpl_domain_mod.F90 3380 2007-03-06 05:42:19Z robj $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_domain_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_domain_mod -- Grid and decomposition information
!
! !DESCRIPTION:
!     The {\it domain} data type is a fundamental coupler data type.
!     A {\it domain} contains both physical information about a {\it grid}
!     such as latitude and longitude values for the points as well as
!     information about the {\it decomposition}.  
!     The decomposition is described by an MCT {\it GlobalSegMap} (gsMap).
!
!     NOTE:  Currently there is no initialization routine in the module.
!     {\it domains} are initialized during the {\it contract} initialization
!     in {\tt cpl\_contract\_init}.
!
! !REVISION HISTORY:
!     2001-aug-15 - B. Kauffman - created module
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_domain_mod

! !USES:

   use shr_sys_mod    ! shared system call wrappers
   use cpl_kind_mod   ! kinds
   use mct_mod    ! MCT API
   use cpl_comm_mod   ! communicator groups, pids, etc.
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug

   implicit none

   private ! except

! !PUBLIC TYPES:

   public :: cpl_domain

   type cpl_domain
      !--- decomposition-independant data ---
      character(80) :: name   ! = "null" ! name of domain       (eg. "ocean")
      character(80) :: suffix ! = "null" ! netCDF domain suffix (eg. "o")
      integer(IN)   :: n      ! n = ni*nj ~ total number of grid pts (global)
      integer(IN)   :: ni     ! number of 2d array i indicies        (global)
      integer(IN)   :: nj     ! number of 2d array j indicies        (global)

      !--- decomposition-dependant data ---
      type(mct_aVect) :: lGrid ! grid data      
      type(mct_gsMap) :: gsMap ! global seg map (defines decomp)
   end type cpl_domain

! !PUBLIC MEMBER FUNCTIONS:

   public cpl_domain_info     ! print some info about a domain
   public cpl_domain_clean    ! clean/dealloc a domain
   public cpl_domain_compare  ! compare two domains for consistency

! !PUBLIC DATA MEMBERS:

   integer,parameter,public :: cpl_domain_compare_algrid = 9901
   integer,parameter,public :: cpl_domain_compare_oigrid = 9902

!EOP

   character(*),parameter :: modName = "cpl_domain_mod"  ! module name
   integer,parameter      :: pid0 = 0                    ! root PID

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_domain_info -- Write info about domain.
!
! !DESCRIPTION:
!    Write basic information about the input domain {\tt cpl\_domain\_x}
!    to stdout.  This information is useful for debugging.
!
! !REVISION HISTORY:
!     2001-Dec-20 - B. Kauffman -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_domain_info(cpl_domain_x)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_domain)   ,target,intent(in) :: cpl_domain_x  ! domain

!EOP

   !--- local ---

   !--- formats ---
   character(*),parameter :: subName = "(cpl_domain_info) "
   character(*),parameter :: F00 = '("(cpl_domain_info) ",60a)'
   character(*),parameter :: F01 = '("(cpl_domain_info) ", a,i9,2i6)'
   character(*),parameter :: F02 = '("(cpl_domain_info) ", a,i3,2i6)'

!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

   write(6,F00) "<--- domain data type info dump START --->"

   write(6,F00) "name = "    ,trim(cpl_domain_x%name  ),&
                ", suffix = ",trim(cpl_domain_x%suffix)
   write(6,F01) "n,ni,nj ="  ,cpl_domain_x%n,cpl_domain_x%ni,cpl_domain_x%nj
   if (dbug >= 2) then
      write(6,F02) "gsMap comp_id, ngSeg, gSize = ", &
      cpl_domain_x%gsMap%comp_id, &
      cpl_domain_x%gsMap%ngseg,   &
      cpl_domain_x%gsMap%gsize
   end if

   !--- write local grid info ---
   if (dbug == 1) call mct_aVect_info(1,cpl_domain_x%lGrid,cpl_comm_comp,cpl_comm_comp_pid)
!  if (dbug == 2) call mct_aVect_info(2,cpl_domain_x%lGrid,cpl_comm_comp,cpl_comm_comp_pid)
   if (dbug >= 2) call mct_aVect_info(4,cpl_domain_x%lGrid,cpl_comm_comp,cpl_comm_comp_pid)

   write(6,F00) "<--- domain data type info dump END   --->"

   call shr_sys_flush(6)

end subroutine cpl_domain_info

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_domain_clean -- Clean a domain type
!
! !DESCRIPTION:
!    Thie routine deallocates the allocated memory associated with 
!    the input/ouput {\tt dom} argument.
!
! !REVISION HISTORY:
!     2002-Jan-20 - T. Craig -- first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine cpl_domain_clean(dom)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_domain)   ,intent(inout) :: dom  ! domain

!EOP

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------

   dom%name   = "<null>"
   dom%suffix = "<null>"
   dom%n  = 0
   dom%ni = 0
   dom%nj = 0
   call mct_aVect_clean(dom%lGrid)
   call mct_gsMap_clean(dom%gsMap)
   dom%gsMap%comp_id = 0
   dom%gsMap%ngseg   = 0
   dom%gsMap%gsize   = 0
   nullify(dom%gsMap%start)
   nullify(dom%gsMap%length)
   nullify(dom%gsMap%pe_loc)

end subroutine cpl_domain_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_domain_compare - Compare two domains
!
! !DESCRIPTION:
!    Compares two domains and summmarizes the differences.  It compares
!    the size, and also checks that the mask, and both model and mapping
!    area are identical to within a {\tt eps} factor defined below.
!
!    The various {\tt enforce\_*} optional arguments will, if present and true,
!    force this routine to abort if the desired test fails.
!
! !REVISION HISTORY:
!    2004-May-21 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_domain_compare(dom1,dom2,enforce_mask,enforce_grid, &
           &                            enforce_area,enforce_aream, enforce_all, &
           &                            gridtype)

! !USES:
   use cpl_control_mod, only: cpl_control_eps_almask, &
                              cpl_control_eps_algrid, &
                              cpl_control_eps_alarea, &
                              cpl_control_eps_oimask, &
                              cpl_control_eps_oigrid, &
                              cpl_control_eps_oiarea

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_domain),intent(in) :: dom1          ! domain #1
   type(cpl_domain),intent(in) :: dom2          ! domain #2
   logical,optional,intent(in) :: enforce_mask  ! abort if masks differ wrt zero/nonzero
   logical,optional,intent(in) :: enforce_grid  ! abort if grids differ by eps_grid
   logical,optional,intent(in) :: enforce_area  ! abort if area  differ by eps_area
   logical,optional,intent(in) :: enforce_aream ! abort if aream differ by eps_area
   logical,optional,intent(in) :: enforce_all   ! abort for all of the above
   integer,optional,intent(in) :: gridtype      ! gridtype, al or oi

!EOP

   !----- local -----
   character(CL)        :: str      ! generic text string
   integer(IN)          :: n,k      ! generic index
   logical              :: enforce  ! flags abort if inconsistent
   integer(IN)          :: ndiff    ! number of points differing
   integer(IN)          :: npts     ! number of points in a field

   integer(IN)          :: npts1    ! number of points in a dom1 field
   integer(IN)          :: npts2    ! number of points in a dom1 field
   integer(IN)          :: rcode    ! return code
   type(mct_aVect)  :: gGrid1   ! global/gathered bundle data
   type(mct_aVect)  :: gGrid2   ! global/gathered bundle data
   real(R8),allocatable :: data1(:) ! temporary real vector
   real(R8),allocatable :: data2(:) ! temporary real vector
   real(R8),allocatable :: mask (:) ! temporary real vector, domain mask

   real(R8)             :: max_diff             ! maximum diff
   real(R8)             :: diff                 ! average diff
!--now set via namelist, in cpl_control_mod--
!  real(R8),parameter   :: eps_mask = 1.0e-6_R8 ! epsilon for masks
!  real(R8),parameter   :: eps_grid = 1.0e-2_R8 ! epsilon for grid coords
!  real(R8),parameter   :: eps_area = 1.0e-1_R8 ! epsilon for areas
   real(R8)             :: eps_mask             ! epsilon for masks
   real(R8)             :: eps_grid             ! epsilon for grid coords
   real(R8)             :: eps_area             ! epsilon for area
   real(R8)             :: x1,x2                ! temp vars wrt wrap-around lat

   !----- formats -----
   character(*),parameter :: subName = '(cpl_domain_compare) '
   character(*),parameter :: F00 = "('(cpl_domain_compare) ',4a)"
   character(*),parameter :: F01 = "('(cpl_domain_compare) ',a,i6,a)"
   character(*),parameter :: F02 = "('(cpl_domain_compare) ',a,es10.3,a,es10.3)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (present(gridtype)) then
      if (gridtype == cpl_domain_compare_algrid) then
         eps_mask = cpl_control_eps_almask
         eps_grid = cpl_control_eps_algrid
         eps_area = cpl_control_eps_alarea
      elseif (gridtype == cpl_domain_compare_oigrid) then
         eps_mask = cpl_control_eps_oimask
         eps_grid = cpl_control_eps_oigrid
         eps_area = cpl_control_eps_oiarea
      else
         write(6,F00) 'gridtype illegal'
         call shr_sys_abort(subName // "ERROR: gridtype illegal")
      endif
   else
      write(6,F00) 'gridtype required'
      call shr_sys_abort(subName // "ERROR: gridtype required")
   endif

   write(6,F00) "domain #1 name: ",trim(dom1%name)
   write(6,F00) "domain #2 name: ",trim(dom2%name)

   !----------------------------------------------------------------------------
   ! do global gather to create non-distributed aVect (*part* of a bundle)
   !----------------------------------------------------------------------------
   call mct_aVect_gather(dom1%lGrid,gGrid1,dom1%gsMap,pid0,cpl_comm_comp,rcode)
   call mct_aVect_gather(dom2%lGrid,gGrid2,dom2%gsMap,pid0,cpl_comm_comp,rcode)
   
   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0) return

   !----------------------------------------------------------------------------
   ! compare size of domain
   !----------------------------------------------------------------------------
   npts1  = mct_aVect_lsize (gGrid1)
   npts2  = mct_aVect_lsize (gGrid2)
   npts   = npts1

   if (npts1 == npts2) then
      write(6,F01) "the domain size is = ", npts
   else
      write(6,F01) "domain size #1 = ", npts1
      write(6,F01) "domain size #2 = ", npts2
      write(6,F00) "ERROR: domain size mis-match"
      call shr_sys_abort(subName // "ERROR: domain size mis-match")
   end if

   allocate(data1(npts))
   allocate(data2(npts))
   allocate(mask (npts))

   !----------------------------------------------------------------------------
   ! compare domain mask
   !----------------------------------------------------------------------------

   call mct_aVect_getRAttr(gGrid1,"mask",data1,rcode)
   call mct_aVect_getRAttr(gGrid2,"mask",data2,rcode)

   ndiff = 0
   mask = 0.0_r8
   do n=1,npts
      if ( (abs(data1(n)) > eps_mask) .and. (abs(data2(n)) < eps_mask) .or. &
      &    (abs(data1(n)) < eps_mask) .and. (abs(data2(n)) > eps_mask) ) then
          ndiff = ndiff + 1
          if (dbug > 2) write(6,*) 'tcx0 ',n,ndiff,data1(n),data2(n)
      end if
      if ( (abs(data1(n)) > eps_mask) .and. (abs(data2(n)) > eps_mask)) then
          mask(n) = eps_mask * 1.e6_R8
      end if
   end do

   !--- enforce consistency ?? ---
   enforce = .false. 
   if (present(enforce_mask)) then
      if (enforce_mask) enforce = .true. 
   end if
   if (present(enforce_all)) then
      if (enforce_all) enforce = .true. 
   end if
   if (enforce .and. ndiff > 0) then
      write(6,F00) "ERROR: incompatible domain masks"
      call shr_sys_abort(subName // "incompatible domain masks")
   end if

   !----------------------------------------------------------------------------
   ! compare grid points (latitude & longitude)
   !----------------------------------------------------------------------------

   ndiff = 0

   call mct_aVect_getRAttr(gGrid1,"lat",data1,rcode)
   call mct_aVect_getRAttr(gGrid2,"lat",data2,rcode)
   max_diff = 0.0_R8
   do n=1,npts
      if ( abs(mask(n)) > eps_mask ) then
         diff = abs(data1(n)-data2(n))
         max_diff = max(max_diff,diff)
         if ( diff > eps_grid ) then
            ndiff = ndiff + 1
         endif
      end if
   end do
   write(6,F02) "maximum latitude  difference = ",max_diff," eps=",eps_grid

   call mct_aVect_getRAttr(gGrid1,"lon",data1,rcode)
   call mct_aVect_getRAttr(gGrid2,"lon",data2,rcode)
   max_diff = 0.0_R8
   do n=1,npts
      if ( abs(mask(n)) > eps_mask ) then
         x1 = data1(n)
         x2 = data2(n)
         if (x1 > x2) then ! make sure x1 < x2
            x1 = data2(n)
            x2 = data1(n)
         end if
         do while ( (x1+360.0_R8) < (x2+180.0_R8) ) ! longitude is periodic
            x1 = x1 + 360.0_R8
         end do
         diff = abs(x2 - x1)
         max_diff = max(max_diff,diff)
         if ( diff > eps_grid ) then
            ndiff = ndiff + 1
         endif
      end if
   end do
   write(6,F02) "maximum longitude difference = ",max_diff," eps=",eps_grid

   !--- enforce consistency ?? ---
   enforce = .false. 
   if (present(enforce_grid)) then
      if (enforce_grid) enforce = .true. 
   end if
   if (present(enforce_all)) then
      if (enforce_all) enforce = .true. 
   end if
   if (enforce .and. ndiff > 0) then
      write(6,F00) "ERROR: incompatible domain grid coordinates"
      call shr_sys_abort(subName // "incompatible domain grid coordinates")
   end if

   !----------------------------------------------------------------------------
   ! compare area
   !----------------------------------------------------------------------------

   call mct_aVect_getRAttr(gGrid1,"area",data1,rcode)
   call mct_aVect_getRAttr(gGrid2,"area",data2,rcode)

   ndiff = 0
   max_diff = 0.0_R8
   do n=1,npts
      if ( abs(mask(n)) > eps_mask ) then
         if (data2(n)/=0.0_R8) diff = abs((data2(n)-data1(n))/data2(n)) 
         max_diff = max(max_diff,diff)
         if ( diff > eps_area ) then
            ndiff = ndiff + 1
            if (dbug > 2) write(6,*) 'tcx3 ',n,ndiff,mask(n),data1(n),data2(n)
         endif
      end if
   end do
   write(6,F02) "maximum relative error of area (model) = ",max_diff," eps=",eps_area

   !--- enforce consistency ?? ---
   enforce = .false. 
   if (present(enforce_area)) then
      if (enforce_area) enforce = .true. 
   end if
   if (present(enforce_all)) then
      if (enforce_all) enforce = .true. 
   end if
   if (enforce .and. ndiff > 0) then
      write(6,F00) "ERROR: incompatible domain area (model)"
      call shr_sys_abort(subName // "incompatible domain area (model)")
   end if

   !----------------------------------------------------------------------------
   ! compare aream
   !----------------------------------------------------------------------------

   call mct_aVect_getRAttr(gGrid1,"aream",data1,rcode)
   call mct_aVect_getRAttr(gGrid2,"aream",data2,rcode)

   ndiff = 0
   max_diff = 0.0_R8
   do n=1,npts
      if ( abs(mask(n)) > eps_mask ) then
         if (data2(n)/=0.0_R8) diff = abs((data2(n)-data1(n))/data2(n)) 
         max_diff = max(max_diff,diff)
         if ( diff > eps_area ) ndiff = ndiff + 1
      end if
   end do
   write(6,F02) "maximum relative error of area (map)   = ",max_diff," eps=",eps_area

   !--- enforce consistency ?? ---
   enforce = .false. 
   if (present(enforce_area)) then
      if (enforce_area) enforce = .true. 
   end if
   if (present(enforce_all)) then
      if (enforce_all) enforce = .true. 
   end if
   if (enforce .and. ndiff > 0) then
      write(6,F00) "ERROR: incompatible domain area (map)"
      call shr_sys_abort(subName // "incompatible domain area (map)")
   end if

   !----------------------------------------------------------------------------
   ! clean-up, deallocate
   !----------------------------------------------------------------------------
   deallocate(data1)
   deallocate(data2)
   deallocate(mask )
   call mct_aVect_clean(gGrid1)
   call mct_aVect_clean(gGrid2)
   call shr_sys_flush(6)

end subroutine cpl_domain_compare

!===============================================================================
!===============================================================================

end module cpl_domain_mod
