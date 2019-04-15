!===============================================================================
! SVN $Id: cpl_bundle_mod.F90 3380 2007-03-06 05:42:19Z robj $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_bundle_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_bundle_mod -- fundamental data type definition
!
! !DESCRIPTION:
!     The {\it bundle} data type is a fundamental coupler data type.
!     A {\it bundle} consists of one or more fields, all of which
!     share the same {\it domain}.  The field data is stored together in an MCT
!     {\it AttributeVector} ({\tt mct\_aVect}) which provides a flexible and
!     indexible storage type.   Following MCT, the individual fields (T, u, SST,
!     etc.) are referred to as ``attributes''.   Each attribute contains the
!     data for a field in a one-dimensional vector.  The size of the vector is
!     equal to the local size of the {\it domain} on a processor.
!
!     This module defines the {\it bundle} datatype and provides several
!     methods for manipulating bundles.
!
! !REVISION HISTORY:
!     2002-Sep-10 - T. Craig - add cpl_bundle_split, cpl_bundle_gather
!     2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_bundle_mod

! !USES:

   use mct_mod
   use cpl_comm_mod
   use cpl_domain_mod
   use cpl_kind_mod
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug
   use cpl_control_mod, only: bfbflag=>cpl_control_bfbflag
   use cpl_const_mod, only :cpl_const_spval
   use shr_sys_mod
   use shr_timer_mod
   use shr_mpi_mod

   implicit none

   private ! except

! !PUBLIC TYPES:

   public :: cpl_bundle

   type cpl_bundle
      character(80)            :: name  ! id string for bundle
      type(mct_aVect)      :: data  ! attribute vector containing data
      type(cpl_domain),pointer :: dom   ! domain associated with data
      integer(IN)              :: cnt   ! counter for accumulating bundles
   end type cpl_bundle

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_bundle_init
   public :: cpl_bundle_initv
   public :: cpl_bundle_clean
   public :: cpl_bundle_info
   public :: cpl_bundle_fill
   public :: cpl_bundle_dump
   public :: cpl_bundle_copy
   public :: cpl_bundle_fcopy
   public :: cpl_bundle_split
   public :: cpl_bundle_gather
   public :: cpl_bundle_hasAttr
   public :: cpl_bundle_zero
   public :: cpl_bundle_accum
   public :: cpl_bundle_avg
   public :: cpl_bundle_add
   public :: cpl_bundle_mult
   public :: cpl_bundle_divide
   public :: cpl_bundle_gsum
   

! !PUBLIC DATA MEMBERS:

  ! no public data

!EOP

   !--- module variables ---
   character(*),parameter :: modName = "cpl_bundle_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_init - Initialize the bundle data type
!
! !DESCRIPTION:
!     Initialize the bundle data type {\tt bun}.  {\tt bun} will
!     be given the name {\tt name} and associated with the {\it domain}
!     {\tt dom}.  Memory will be allocated to hold the attributes listed
!     in the input {\tt rList} in an MCT {\it AttributeVector}.
!
! !REVISION HISTORY:
!     2001-Mar-20  - T. Craig, B. Kauffman, R. Jacob - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_init(bun,name,rList,dom)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(out)       :: bun   ! bundle to initialize
   character(*)    ,intent(in)        :: name  ! name used in netCDF files
   character(*)    ,intent(in)        :: rList ! aVect real data list string 
   type(cpl_domain),intent(in),target :: dom   ! domain assigned to bundle

!EOP

   !--- local ---
   integer(IN)   :: lSize  ! number of data points (local) in an aVect field
   integer(IN)   :: nFlds  ! number of real fields in aVect

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_init) '
   character(*),parameter :: F00 = "('(cpl_bundle_init) ',8a)"
   character(*),parameter :: F01 = "('(cpl_bundle_init) ',3a,i7,3a,i3,2a)"

!-------------------------------------------------------------------------------
! METHOD:
!  1) assign the domain (an input) to the bundle
!  2) set size of bundle aVect to be the same size as that of the domain
!  3) use input str to define real aVect fields (must be a valid aVect list str)
!  4) hard-coded st the aVect has *no* integer data
! NOTE:
!  o memory is allocated for bun%data, but data values are left undefined
!-------------------------------------------------------------------------------

   !--- init bundle ---
   bun%name =  name
   bun%dom  => dom
   bun%cnt  =  0
   lSize = mct_GSMap_lSize(dom%GSMap,cpl_comm_comp)
   call mct_aVect_init(bun%data,' ',rList,lSize)
   call cpl_bundle_zero(bun)

   if (dbug > 0) then
      nFlds = mct_aVect_nRAttr(bun%data)
      write(6,F01)   "bundle: ",trim(bun%name)    , &
                   ", lSize : ",lSize             , &
                   ", domain: ",trim(bun%dom%name), &
                   ", nFlds :" ,nFlds             , &
                   ", rList : ",trim(rList)
   else if (dbug > 2) then
      call cpl_bundle_info(bun)
      call cpl_domain_info(bun%dom)
   end if

end subroutine cpl_bundle_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_initv - Initialize the bundle data type using another bundle
!
! !DESCRIPTION:
!     This routine will initialize {\tt bun} to have the same attributes as
!     the input {\tt bun2} but will be given the new name {\tt name} and
!     will be associated with the {\it domain} {\tt dom}.  Note that this
!     routine only allocates memory and does not copy the contents of
!     {\tt bun2} to {\tt bun}.
!
! !REVISION HISTORY:
!     2002-jan-15 - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_initv(bun,name,bun2,dom)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(out)       :: bun  ! bundle to initialize
   character(*)    ,intent(in)        :: name ! name used in netCDF files
   type(cpl_bundle),intent(in)        :: bun2 ! bundle to "copy" from
   type(cpl_domain),intent(in),target :: dom  ! domain assigned to bundle

!EOP

   !--- local ---
   integer(IN)  :: lSize  ! number of data points ("local size") in an aVect field
   integer(IN)  :: nFlds  ! number of real fields in aVect
   character(CL):: str    ! string for field

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_initv) '
   character(*),parameter :: F00 = "('(cpl_bundle_initv) ',8a)"
   character(*),parameter :: F01 = "('(cpl_bundle_initv) ',a,a15,a25,i7,i3,a)"

!-------------------------------------------------------------------------------
! METHOD:
!  1) assign the domain (an input) to the bundle
!  2) set size of bundle aVect to be the same size as that of the domain
!  3) use input str to define real aVect fields (must be a valid aVect list str)
!  4) hard-coded st the aVect has *no* integer data
! NOTE:
!  o memory is allocated for bun%data, but data values are left undefined
!-------------------------------------------------------------------------------

   !--- init bundle ---
   bun%name =  name
   bun%dom  => dom
   bun%cnt  =  0
   lSize = mct_GSMap_lsize(dom%GSMap,cpl_comm_comp)
   call mct_aVect_init(bun%data,bun2%data,lSize)
   call cpl_bundle_zero(bun)

!  if (dbug > 1) write(6,F00) "initializing bundle = ",trim(bun%name)
   if (dbug > 2) call cpl_bundle_info(bun)
   if (dbug > 2) call cpl_domain_info(bun%dom)

end subroutine cpl_bundle_initv

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_clean - Clean the bundle data type
!
! !DESCRIPTION:
!     Deallocate all allocated memory associated with the
!     input/output argument {\tt bun}.
!
! !REVISION HISTORY:
!     2002-Jan-20 - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_clean(bun)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle)    ,intent(inout)     :: bun ! bundle to initialize

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   bun%name = "unknown"
   call mct_aVect_clean(bun%data)
   nullify(bun%dom)
   bun%cnt = 0

end subroutine cpl_bundle_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_info - Print out bundle info.
!
! !DESCRIPTION:
!     Print out the following information about the bundle {\tt bun}:
!     \begin{itemize}
!       \item the name
!       \item the name of the associated {\it domain}
!       \item the accumulation count
!       \item the names of the attributes
!     \end{itemize}
!
! !REVISION HISTORY:
!     2002-May-09 - B. Kauffman - make's use of mct_aVect_info routine
!     2001-Jun-14 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_info(bun)

! !USES:

   use cpl_fields_mod

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle)    ,intent(in)     :: bun ! bundle to initialize

!EOP

   !--- local ---
   integer(IN)  :: lSize  ! number of data points (local) in an aVect field
   integer(IN)  :: nFlds  ! number of real fields in aVect
   character(CL):: rList  ! aVect rList string

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_info) '
   character(*),parameter :: F00 = "('(cpl_bundle_info) ',8a)"
   character(*),parameter :: F01 = "('(cpl_bundle_info) ',5a,i6)"

!-------------------------------------------------------------------------------
! NOTE: has hard-coded knowledge of MCT internal data-types :(
!-------------------------------------------------------------------------------

   write(6,F01)   "bundle name = ", trim(bun%name)    , &
                ", domain name = ", trim(bun%dom%name), &
                ", accumulation count =",bun%cnt

   call mct_aVect_info(1,bun%data)

   call shr_sys_flush(6)

end subroutine cpl_bundle_info

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_fill - Fill a bundle with test data.
!
! !DESCRIPTION:
!     Fill each field in the input bundle {\tt bun} with a slightly
!     different sine wave.  Used for debugging.
!     
!
! !REVISION HISTORY:
!     2001-Jun-14 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_fill(bun)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout) :: bun   ! bundle to fill

!EOP

   !--- local ---
   integer(IN)  :: i,j    ! generic indicies
   integer(IN)  :: npts   ! number of points (local) in an aVect field
   integer(IN)  :: nflds  ! number of aVect fields (real)

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_fill) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   nflds = mct_aVect_nRAttr(bun%data)
   npts  = mct_aVect_lsize (bun%data)
!DIR$ CONCURRENT
   do j=1,nflds
!DIR$ PREFERVECTOR
   do i=1,npts 
      bun%data%rattr(j,i)=real(j*10,R8) + sin(real(i,R8)/500._R8)
   enddo
   enddo
   bun%cnt = 1

end subroutine cpl_bundle_fill

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_dump - write bundle contents to a file.
!
! !DESCRIPTION:
!     Write out the contents of {\tt bun} on each processor to a
!     seperate file.  Filename will be fort.({\tt iun} +  coupler processor id).
!     
!
! !REVISION HISTORY:
!     2002-Jan-14 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_dump(iun,bun)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(in) :: bun   ! bundle to write
   integer(IN),intent(in)      :: iun   ! base unit number

!EOP

   !--- local ---
   integer(IN) :: aVnum, aVsiz, id1     ! aV properties
   integer(IN) :: i,j                   ! generic indices

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_dump) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   aVnum = mct_aVect_nrAttr (bun%data)
   aVsiz = mct_aVect_lsize  (bun%data)
   id1   = mct_aVect_indexRA(bun%dom%lGrid,'index')
   write(6,*) subName,' info ',aVnum,aVsiz,id1
   do i=1,aVsiz
      write(iun+cpl_comm_comp_pid,*)   &
        bun%dom%lGrid%rAttr(id1,i),    &
       (bun%data%rAttr(j,i),j=1,aVnum)
   enddo
   call shr_sys_flush(6)

end subroutine cpl_bundle_dump

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_copy - copy data from one bundle to another
!
! !DESCRIPTION:
!     This routine copies from input argument {\tt inbun} into the
!     output argument {\tt outbun} the data of all the attributes shared
!     between the two.  If only a subset of shared attributes should be copied,
!     use the optional arguments {\tt bunrList} and {\tt buniList} to specify
!     which attributes should be copied.  If any attributes in {\tt outbun}
!     have different names for the same quantity, provide a corresponding optional
!     translation list, {\tt bunTrList} and {\tt bunTiList}, in addition
!     to {\tt bunrList} and {\tt buniList}.  The translation list is identical
!     to the input List with translated names in place of the names in {\tt inbun}.
!
!     Optional argument {\tt fcopy} directs this routine  to use {\tt cpl\_bundle\_copy}
!     instead of {\tt cpl\_mct\_aVect\_copy}.  This is obsolete since 
!     {\tt cpl\_mct\_aVect\_copy} now uses the same algorithm.
!
! !REVISION HISTORY:
!     2002-Jul-02 - R. Jacob -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_copy(inbun,bunrList,bunTrList,buniList,bunTiList,outbun,fcopy)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(in)          :: inbun    ! bundle to read
   character(*)    ,intent(in),optional :: buniList
   character(*)    ,intent(in),optional :: bunrList
   character(*)    ,intent(in),optional :: bunTiList
   character(*)    ,intent(in),optional :: bunTrList
   type(cpl_bundle),intent(out)         :: outbun    ! bundle to write to
   logical         ,intent(in),optional :: fcopy     ! use fcopy

!EOP
   !--- local vars ---
   logical :: fcopy_loc                  ! call fcopy instead of mct_copy
   logical :: usevector                  ! use vector-friendly mct_copy

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_copy) '
   character(*),parameter :: F01 = '(  "(cpl_bundle_copy) ",3a,i5)'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   fcopy_loc = .false.
   if (present(fcopy)) fcopy_loc = fcopy

#ifdef CPP_VECTOR
!   fcopy_loc = .false.
!   fcopy_loc = .true.
   usevector = .true.
#else
   usevector = .false.
#endif

   !--- If outbun or inbun has no attributes, return ---
   if( .NOT.cpl_bundle_hasAttr(outbun) .or.   &
       .NOT.cpl_bundle_hasAttr(inbun) ) then
     return
   endif

   if (inbun%cnt /= 1) write(6,F01) &
      "WARNING: bundle ",trim(inbun%name)," has accum count =",inbun%cnt

   !--- copy real attributes if specified ---
   if(present(bunrList)) then
     if(present(bunTrList)) then
       call mct_aVect_copy(inbun%data,rList=bunrList,TrList=bunTrList,aVout=outbun%data,vector=usevector)
     else
       call mct_aVect_copy(inbun%data,rList=bunrList,aVout=outbun%data,vector=usevector)
     endif
   endif

   !--- copy integer attributes if specified ---
   if(present(buniList)) then
     if(present(bunTiList)) then
       call mct_aVect_copy(inbun%data,iList=buniList,TiList=bunTiList,aVout=outbun%data,vector=usevector)
     else
       call mct_aVect_copy(inbun%data,iList=buniList,aVout=outbun%data,vector=usevector)
     endif
   endif

   !--- if no lists given, copy all shared attributes ---
   if(.not.present(buniList).and. .not.present(bunrList)) then
      if (fcopy_loc) then
        call cpl_bundle_fcopy(inbun,outbun)
      else
        call mct_aVect_copy(inbun%data,aVout=outbun%data,vector=usevector)
      endif
   endif

   outbun%cnt = inbun%cnt

end subroutine cpl_bundle_copy

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_fcopy - Fast copy version of cpl_bundle_copy. Obsolete.
!
! !DESCRIPTION:
!     This routine copies from input argument {\tt inbun} into the
!     output argument {\tt outbun} the data of all the attributes shared real
!     between the two.
!
!     NOTE:  This routine will be deleted in future versions since it is 
!     no longer necessary.
!
! !REVISION HISTORY:
!     2003-Aug-26 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_fcopy(inbun,outbun)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(in)          :: inbun    ! bundle to read
   type(cpl_bundle),intent(out)         :: outbun    ! bundle to write to

!EOP

   !--- local vars ---
   integer(IN) :: nfld_i,nfld_o                  ! number of fields
   integer(IN) :: kfld_i,kfld_o                  ! field number
   integer(IN) :: nsize_i,nsize_o                ! size of bundles
   type(mct_string) :: item_i,item_o         ! mct string
   character(CL),allocatable :: itemc_i(:),itemc_o(:)     ! item converted to char
   integer(IN),allocatable :: klist_i(:),klist_o(:)  ! list of copy indices
   integer(IN) :: nklist,ni,no,n,noprev

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_fcopy) '
   character(*),parameter :: F01 = '(  "(cpl_bundle_fcopy) ",3a,i5)'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (inbun%cnt /= 1) write(6,F01) &
      "WARNING: bundle ",trim(inbun%name)," has accum count =",inbun%cnt

   nsize_i = mct_aVect_lsize(inbun%data)
   nsize_o = mct_aVect_lsize(outbun%data)
   if (nsize_i .ne. nsize_o) then
     write(6,*) trim(subName),'ERROR on bundle sizes ',    &
       nsize_i,nsize_o
     call shr_sys_abort(trim(subName))
   endif

   nfld_i = mct_aVect_nRAttr(inbun%data)
   nfld_o = mct_aVect_nRAttr(outbun%data)
   allocate(klist_i(nfld_o),klist_o(nfld_o))
   allocate(itemc_i(nfld_i),itemc_o(nfld_o))
   nklist = 0
   noprev = 0

   do no = 1,nfld_o
     call mct_aVect_getRList(item_o,no,outbun%data)
     itemc_o(no) = trim(mct_string_toChar(item_o))
     call mct_string_clean(item_o)
   enddo

   do ni = 1,nfld_i
     call mct_aVect_getRList(item_i,ni,inbun%data)
     itemc_i(ni) = trim(mct_string_toChar(item_i))
     call mct_string_clean(item_i)
   enddo

   do no = 1,nfld_o
     do ni = 1,nfld_i
       if (itemc_i(ni) == itemc_o(no)) then
         if (no.ne.noprev) nklist = nklist + 1
         klist_i(nklist) = ni
         klist_o(nklist) = no
         noprev = no
       endif
    enddo
   enddo

!   write(6,*) trim(subName),'list1 ',nklist,trim(inbun%name),' ',trim(outbun%name)

   if (nklist > 0) then
!     write(6,*) trim(subName),'list2 ',nklist,klist_i(1:nklist),klist_o(1:nklist)
!DIR$ CONCURRENT
     do n = 1,nsize_o
     do no = 1,nklist
        outbun%data%rAttr(klist_o(no),n) = inbun%data%rAttr(klist_i(no),n)
     enddo
     enddo
   endif

   deallocate(klist_i,klist_o)
   deallocate(itemc_i,itemc_o)
   outbun%cnt = inbun%cnt

end subroutine cpl_bundle_fcopy

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_split - Split bundle into multiple bundles.
!
! !DESCRIPTION:
!     Copy data from input bundle {\tt bun\_X} into multiple output bundles.
!     Can have up to 12 output bundles: {\tt bun1, bun2, bun3} etc.
!     On return, each output bundle with an attribute name that matches an attribute in
!     {\tt bun\_X} will have a copy of the data for that attribute.
!
!     The optional argument {\tt fcopy} will use the fast version of cpl\_bundle\_copy.
!
!     This routine is used to ``split'' data in a bundle received in the Coupler
!     from a model into bundles for each data pathway through the Coupler. 
!
!     NOTE: All arguments must already be initialized.
!
! !REVISION HISTORY:
!     2002-Apr-12 - B. Kauffman - first version
!     2002-Jul-02 - R. Jacob - use bundle copy
!     2002-Jul-15 - T. Craig - generalized
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_split(bun_X,bun1,bun2,bun3,bun4,bun5,bun6,bun7,bun8,bun9,bun10,bun11,bun12,fcopy)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),         intent(in ) :: bun_X ! input bundle
   type(cpl_bundle),optional,intent(out) :: bun1  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun2  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun3  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun4  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun5  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun6  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun7  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun8  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun9  ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun10 ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun11 ! split bundle 
   type(cpl_bundle),optional,intent(out) :: bun12 ! split bundle 
   logical         ,optional,intent(in)  :: fcopy ! use fcopy

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_split) '
   character(*),parameter :: F00 = '(  "(cpl_bundle_split) ",4a  )'
   character(CL)          :: str
   logical                :: fcopy_loc

!-------------------------------------------------------------------------------
! METHOD:
! NOTES:
! o all data is on same grid
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! copy all shared attributes to all bundles that were sent
   !----------------------------------------------------------------------------

   fcopy_loc = .false.
   if (present(fcopy)) fcopy_loc = fcopy

   if (.not.present(bun1)) return
   call cpl_bundle_copy(bun_X,outbun=bun1,fcopy=fcopy_loc)
   if (.not.present(bun2)) return
   call cpl_bundle_copy(bun_X,outbun=bun2,fcopy=fcopy_loc)
   if (.not.present(bun3)) return
   call cpl_bundle_copy(bun_X,outbun=bun3,fcopy=fcopy_loc)
   if (.not.present(bun4)) return
   call cpl_bundle_copy(bun_X,outbun=bun4,fcopy=fcopy_loc)
   if (.not.present(bun5)) return
   call cpl_bundle_copy(bun_X,outbun=bun5,fcopy=fcopy_loc)
   if (.not.present(bun6)) return
   call cpl_bundle_copy(bun_X,outbun=bun6,fcopy=fcopy_loc)
   if (.not.present(bun7)) return
   call cpl_bundle_copy(bun_X,outbun=bun7,fcopy=fcopy_loc)
   if (.not.present(bun8)) return
   call cpl_bundle_copy(bun_X,outbun=bun8,fcopy=fcopy_loc)
   if (.not.present(bun9)) return
   call cpl_bundle_copy(bun_X,outbun=bun9,fcopy=fcopy_loc)
   if (.not.present(bun10)) return
   call cpl_bundle_copy(bun_X,outbun=bun10,fcopy=fcopy_loc)
   if (.not.present(bun11)) return
   call cpl_bundle_copy(bun_X,outbun=bun11,fcopy=fcopy_loc)
   if (.not.present(bun12)) return
   call cpl_bundle_copy(bun_X,outbun=bun12,fcopy=fcopy_loc)

end subroutine cpl_bundle_split

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_gather - Copy data into one bundle from many
!
! !DESCRIPTION:
!     Copy into the output argument {\tt bun\_X} all the data from
!     the input bundles {\tt bun1, bun2, bun3} etc. which have the
!     same attribute names.
!
!     The optional argument {\tt fcopy} will use the fast version of cpl\_bundle\_copy.
!
!     If there is a common attribute name between 2 or more of the input bundles
!     which also exists in the output
!     bundle, the values in the output bundle will that of the last listed
!     input bundle.
!
!     This routine is used to ``gather'' data from multiple sources into one bundle
!     for sending from the Coupler to a model.
!
!     NOTE: All arguments must already be initialized.
!
! !REVISION HISTORY:
!     2002-Jun-22 - B. Kauffman - first version
!     2002-Jul-15 - T. Craig - generalized
!     2003-Sep-01 - T. Craig - add fcopy optional argument
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_gather(bun_X, bun1,bun2,bun3,bun4,bun5,bun6,bun7,bun8,bun9,bun10,bun11,bun12,fcopy)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),         intent(out) :: bun_X ! output bundle
   type(cpl_bundle),optional,intent(in ) :: bun1  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun2  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun3  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun4  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun5  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun6  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun7  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun8  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun9  ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun10 ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun11 ! gather bundle
   type(cpl_bundle),optional,intent(in ) :: bun12 ! gather bundle
   logical         ,optional,intent(in)  :: fcopy ! use fcopy

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_gather) '
   character(*),parameter :: F00 = '(  "(cpl_bundle_gather) ",4a  )'
   character(CL)          :: str
   logical                :: fcopy_loc

!-------------------------------------------------------------------------------
! METHOD:
! NOTES:
! o all data is on same grid
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! copy fields
   !----------------------------------------------------------------------------
   ! copy all shared attributes from all bundles that were sent

   fcopy_loc = .false.
   if (present(fcopy)) fcopy_loc = fcopy

   if (.not.present(bun1)) return
   call cpl_bundle_copy(bun1,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun2)) return
   call cpl_bundle_copy(bun2,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun3)) return
   call cpl_bundle_copy(bun3,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun4)) return
   call cpl_bundle_copy(bun4,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun5)) return
   call cpl_bundle_copy(bun5,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun6)) return
   call cpl_bundle_copy(bun6,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun7)) return
   call cpl_bundle_copy(bun7,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun8)) return
   call cpl_bundle_copy(bun8,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun9)) return
   call cpl_bundle_copy(bun9,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun10)) return
   call cpl_bundle_copy(bun10,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun11)) return
   call cpl_bundle_copy(bun11,outbun=bun_X,fcopy=fcopy_loc)
   if (.not.present(bun12)) return
   call cpl_bundle_copy(bun12,outbun=bun_X,fcopy=fcopy_loc)
 
end subroutine cpl_bundle_gather

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_hasAttr
!
! !DESCRIPTION:
!     Return true if input bundle {\tt bun} has any real or integer attributes.
!
! !REVISION HISTORY:
!     2002-Sep-11  - R. Jacob - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function cpl_bundle_hasAttr(bun)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(in ) :: bun

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  if ((mct_list_nitem(bun%data%iList) > 0) .or. &
      (mct_list_nitem(bun%data%rList) > 0)) then
     cpl_bundle_hasAttr = .TRUE.
  else
     cpl_bundle_hasAttr = .FALSE.
  endif

end function cpl_bundle_hasAttr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_zero - Zero values of fields in bundle
!
! !DESCRIPTION:
!     Set all fields (real and integer) in input bundle {\tt bun} to 0.
!     If the optional {\tt fld} argument is present, only set that field to 0.
!
!     The bundle must already be initialized.
!
! !REVISION HISTORY:
!     2002-Sep-15 - T. Craig -- initial version
!     2002-Sep-17 - T. Craig -- add optional fld argument
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_zero(bun,fld)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle)     ,intent(inout) :: bun   ! bundle to zero
   character(*),optional,intent(in)    :: fld   ! field in bundle to zero

!EOP

   !--- local ---
   integer(IN) :: k,n    ! generic indicies
   integer(IN) :: npts   ! number of points (local) in an aVect field

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_zero) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(fld)) then
      npts  = mct_aVect_lsize(bun%data)
      k     = mct_aVect_indexRA(bun%data,fld,perrWith=subName)
      do n=1,npts
         bun%data%rAttr(k,n) = 0.0_R8
      enddo
   else
!     if (cpl_bundle_hasAttr(bun)) call mct_aVect_zero(bun%data)
!CDIR COLLAPSE
      bun%data%rAttr = 0.0_R8
      bun%cnt = 0
   endif

end subroutine cpl_bundle_zero

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_accum - Accumulate fields in a bundle.
!
! !DESCRIPTION:
!     Accumulate fields in input bundle {\tt inbun} in output bundle {\tt outbun}.
!     This mimics cpl\_bundle\_copy but instead of overwriting the data
!     in {\tt outbun} with the data in {\tt inbun}, the data in {\tt inbun}
!     is added to the corresponding fields in {\tt outbun}.
!
!     It is recommended that a cpl\_bundle\_accum only be called once for each set of fields
!     to be accumulated at each accumulation step because of the primitive 
!     nature of the bundle counter (bun%cnt).  It is incremented each
!     time cpl\_bundle\_accum is called for the output bundle.
!
! !REVISION HISTORY:
!     2002-Sep-15 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_accum(inbun,bunrList,bunTrList,buniList,bunTiList,outbun)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(in)          :: inbun   ! bundle to read
   character(*)    ,intent(in),optional :: buniList
   character(*)    ,intent(in),optional :: bunrList
   character(*)    ,intent(in),optional :: bunTiList
   character(*)    ,intent(in),optional :: bunTrList
   type(cpl_bundle),intent(out)         :: outbun   ! bundle to write to

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_accum) '
   character(*),parameter :: F00 = '(  "(cpl_bundle_accum) ",4a  )'
   character(*),parameter :: F01 = '(  "(cpl_bundle_accum) ",3a,i5)'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !--- If outbun or inbun has no attributes, return ---

   if (.NOT.cpl_bundle_hasAttr(outbun) .or.   &
       .NOT.cpl_bundle_hasAttr(inbun) ) then
      return
   endif

   if (inbun%cnt /= 1) write(6,F01) &
      "WARNING: bundle ",trim(inbun%name)," has accum count =",inbun%cnt

   !--- accum real attributes if specified ---
   if (present(bunrList)) then
      if (present(bunTrList)) then
         call mct_aVect_accum(inbun%data,rList=bunrList,TrList=bunTrList,aVout=outbun%data)
      else
         call mct_aVect_accum(inbun%data,rList=bunrList,aVout=outbun%data)
      endif
   endif

   !--- accum integer attributes if specified ---
   if (present(buniList)) then
      if(present(bunTiList)) then
         call mct_aVect_accum(inbun%data,iList=buniList,TiList=bunTiList,aVout=outbun%data)
      else
         call mct_aVect_accum(inbun%data,iList=buniList,aVout=outbun%data)
      endif
   endif

   !--- if no lists given, accum all shared attributes ---
   if (.not.present(buniList).and. .not.present(bunrList)) then
      call mct_aVect_accum(inbun%data,aVout=outbun%data)
   endif

   outbun%cnt = outbun%cnt + 1

end subroutine cpl_bundle_accum

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_avg - averages a bundle
!
! !DESCRIPTION:
!     Average the data in input bundle {\tt bun}.  Divide all fields in the bundle
!     {\tt bun} by the value of the bundle counter, {\tt bun\%cnt}.
!
! !REVISION HISTORY:
!     2002-Sep-15 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_avg(bun)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout) :: bun   ! bundle to read

!EOP

   !--- local ---
   integer(IN) :: i,j    ! generic indicies
   integer(IN) :: npts   ! number of points (local) in an aVect field
   integer(IN) :: nflds  ! number of aVect fields (real)
   real(R8)    :: ravg   ! accumulation count

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_avg) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (bun%cnt == 0) return

   ravg = 1.0_R8/real(bun%cnt,R8)

   nflds = mct_aVect_nRAttr(bun%data)
   npts  = mct_aVect_lsize (bun%data)
!DIR$ CONCURRENT
!DIR$ PREFERVECTOR
   do i=1,npts 
   do j=1,nflds
      bun%data%rattr(j,i) = bun%data%rattr(j,i)*ravg
   enddo
   enddo

   bun%cnt = 1

end subroutine cpl_bundle_avg

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_add - add product of multiple bundles.
!
! !DESCRIPTION:
!     Form the product of individual fields from the input bundles
!     {\tt bun1, bun2, bun3} etc. and {\em add} it to the field {\tt fld} in the
!     input/output argument {\tt bun}.
!     
!     {\tt fld1} is taken from {\tt bun1}, {\tt fld2} is taken from {\tt bun2}, etc.
!     
!
!     If only {\tt fld1} and {\tt bun1} are present, then {\tt fld1} will be added
!     to {\tt bun}
!
!     If the optional argument {\tt scalar} is present, it will also by part of the
!     product.
!
!     If the optional argument {\tt zero} is present, the input/output argument
!     {\tt bun} will be zeroed before adding the product.
!
!     This routine is used to {\it merge} data in the Coupler.
!
! !REVISION HISTORY:
!     2002-Sep-15 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_add(bun,fld,bun1,fld1,bun2,fld2,bun3,fld3,bun4,fld4,bun5,fld5,bun6,fld6,scalar,zero)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle)         ,intent(inout) :: bun    ! bundle output
   character(*)             ,intent(in)    :: fld    ! bun  field name
   type(cpl_bundle)         ,intent(in)    :: bun1   ! bundle input
   character(*)             ,intent(in)    :: fld1   ! bun1 field name
   type(cpl_bundle),optional,intent(in)    :: bun2   ! bundle input
   character(*)    ,optional,intent(in)    :: fld2   ! bun2 field name
   type(cpl_bundle),optional,intent(in)    :: bun3   ! bundle input
   character(*)    ,optional,intent(in)    :: fld3   ! bun3 field name
   type(cpl_bundle),optional,intent(in)    :: bun4   ! bundle input
   character(*)    ,optional,intent(in)    :: fld4   ! bun4 field name
   type(cpl_bundle),optional,intent(in)    :: bun5   ! bundle input
   character(*)    ,optional,intent(in)    :: fld5   ! bun5 field name
   type(cpl_bundle),optional,intent(in)    :: bun6   ! bundle input
   character(*)    ,optional,intent(in)    :: fld6   ! bun6 field name
   real(R8)        ,optional,intent(in)    :: scalar ! scalar multiplier
   logical         ,optional,intent(in)    :: zero   ! zero lhs initially

!EOP

   !--- local ---
   integer(IN),parameter :: mbun=6  ! maximum number of bundles incoming

   integer(IN) :: n,k               ! generic indicies
   integer(IN) :: npts              ! number of points (local) in an aVect field
   integer(IN) :: nptsx             ! number of points (local) in an aVect field
   integer(IN) :: kfld(0:mbun)      ! field number in bundle
   real(R8)    :: rtmp              ! temporary storage
   real(R8)    :: scalar_loc        ! local scalar
   integer(IN) :: kst,knd           ! kstart and kend indices for rhs fields

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_add) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   scalar_loc = 1.0_R8
   if (present(scalar)) scalar_loc = scalar

   if (present(zero)) then
     if (zero) then
       call cpl_bundle_zero(bun,fld)
     endif
   endif

   npts  = mct_aVect_lsize( bun%data)
   if (fld == 'ALL_FIELDS') then
     kst = 1
     knd = mct_aVect_nRAttr(bun%data)
   else
     kst = mct_aVect_indexRA( bun%data,fld ,perrWith=subName)
     knd = kst
   endif

   nptsx = mct_aVect_lsize(bun1%data)
   if (nptsx /= npts) write(6,*) subName,' ERROR: npts error1 ',npts,nptsx
   kfld(1) = mct_aVect_indexRA(bun1%data,fld1,perrWith=subName)

   if (present(bun2)) then
     nptsx = mct_aVect_lsize(bun2%data)
     if (nptsx /= npts) write(6,*) subName,' ERROR: npts error2 ',npts,nptsx
     kfld(2) = mct_aVect_indexRA(bun2%data,fld2,perrWith=subName)
   endif

   if (present(bun3)) then
     nptsx = mct_aVect_lsize(bun3%data)
     if (nptsx /= npts) write(6,*) subName,' ERROR: npts error3 ',npts,nptsx
     kfld(3) = mct_aVect_indexRA(bun3%data,fld3,perrWith=subName)
   endif

   if (present(bun4)) then
     nptsx = mct_aVect_lsize(bun4%data)
     if (nptsx /= npts) write(6,*) subName,' ERROR: npts error4 ',npts,nptsx
     kfld(4) = mct_aVect_indexRA(bun4%data,fld4,perrWith=subName)
   endif

   if (present(bun5)) then
     nptsx = mct_aVect_lsize(bun5%data)
     if (nptsx /= npts) write(6,*) subName,' ERROR: npts error5 ',npts,nptsx
     kfld(5) = mct_aVect_indexRA(bun5%data,fld5,perrWith=subName)
   endif

   if (present(bun6)) then
     nptsx = mct_aVect_lsize(bun6%data)
     if (nptsx /= npts) write(6,*) subName,' ERROR: npts error6 ',npts,nptsx
     kfld(6) = mct_aVect_indexRA(bun6%data,fld6,perrWith=subName)
   endif

   if (present(bun6)) then
#ifdef CPP_VECTOR
     do k=kst,knd
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do k=kst,knd
#endif
       bun%data%rAttr(k,n) = bun%data%rAttr(k,n) &
          + bun1%data%rAttr(kfld(1),n) * bun2%data%rAttr(kfld(2),n) &
          * bun3%data%rAttr(kfld(3),n) * bun4%data%rAttr(kfld(4),n) &
          * bun5%data%rAttr(kfld(5),n) * bun6%data%rAttr(kfld(6),n) &
          * scalar_loc
     enddo
     enddo
   elseif (present(bun5)) then
#ifdef CPP_VECTOR
     do k=kst,knd
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do k=kst,knd
#endif
       bun%data%rAttr(k,n) = bun%data%rAttr(k,n) &
          + bun1%data%rAttr(kfld(1),n) * bun2%data%rAttr(kfld(2),n) &
          * bun3%data%rAttr(kfld(3),n) * bun4%data%rAttr(kfld(4),n) &
          * bun5%data%rAttr(kfld(5),n) &
          * scalar_loc
     enddo
     enddo
   elseif (present(bun4)) then
#ifdef CPP_VECTOR
     do k=kst,knd
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do k=kst,knd
#endif
       bun%data%rAttr(k,n) = bun%data%rAttr(k,n) &
          + bun1%data%rAttr(kfld(1),n) * bun2%data%rAttr(kfld(2),n) &
          * bun3%data%rAttr(kfld(3),n) * bun4%data%rAttr(kfld(4),n) &
          * scalar_loc
     enddo
     enddo
   elseif (present(bun3)) then
#ifdef CPP_VECTOR
     do k=kst,knd
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do k=kst,knd
#endif
       bun%data%rAttr(k,n) = bun%data%rAttr(k,n) &
          + bun1%data%rAttr(kfld(1),n) * bun2%data%rAttr(kfld(2),n) &
          * bun3%data%rAttr(kfld(3),n) &
          * scalar_loc
     enddo
     enddo
   elseif (present(bun2)) then
#ifdef CPP_VECTOR
     do k=kst,knd
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do k=kst,knd
#endif
       bun%data%rAttr(k,n) = bun%data%rAttr(k,n) &
          + bun1%data%rAttr(kfld(1),n) * bun2%data%rAttr(kfld(2),n) &
          * scalar_loc
     enddo
     enddo
   else
#ifdef CPP_VECTOR
     do k=kst,knd
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
     do n=1,npts
#else
     do n=1,npts
     do k=kst,knd
#endif
       bun%data%rAttr(k,n) = bun%data%rAttr(k,n) &
          + bun1%data%rAttr(kfld(1),n) &
          * scalar_loc
     enddo
     enddo
   endif

!   do n=1,npts
!                          rtmp =        bun1%data%rAttr(kfld(1),n)
!     if (present(bun2))   rtmp = rtmp * bun2%data%rAttr(kfld(2),n)
!     if (present(bun3))   rtmp = rtmp * bun3%data%rAttr(kfld(3),n)
!     if (present(bun4))   rtmp = rtmp * bun4%data%rAttr(kfld(4),n)
!     if (present(bun5))   rtmp = rtmp * bun5%data%rAttr(kfld(5),n)
!     if (present(bun6))   rtmp = rtmp * bun6%data%rAttr(kfld(6),n)
!     if (present(scalar)) rtmp = rtmp * scalar

!     bun%data%rAttr(kfld(0),n) = bun%data%rAttr(kfld(0),n) + rtmp
!   enddo

   bun%cnt = 1

end subroutine cpl_bundle_add

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_mult - multiply a bundle by a field.
!
! !DESCRIPTION:
!     Replace each field in {\tt bun} by the product of that field and the
!     field {\tt fld1} from input argument {\tt bun1}.
!
!     If optional argument {\tt bunlist} is present, only those attributes 
!     in {\tt bun} will be replaced.
!
!     If optional argument {\tt initbun} is present, then the data in {\tt bun}
!     is replaced by the product of the data in {\tt initbun} and {\tt fld1}
!     from {\tt bun1}. NOTE:  this assume {\tt initbun} has the exact same
!     attributes in the same order as {\tt bun}.
!
!
! !REVISION HISTORY:
!     2002-Oct-15 - T. Craig -- initial version
!     2003-Jan-02 - T. Craig -- added bundle sub-list option (bunlist)
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_mult(bun,bun1,fld1,bunlist,initbun)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle)     ,intent(inout) :: bun    ! bundle output
   type(cpl_bundle)     ,intent(in)    :: bun1   ! bundle input
   character(*)         ,intent(in)    :: fld1   ! bun1 field name
   character(*),optional,intent(in)    :: bunlist! sublist of field in bun
   type(cpl_bundle),optional,intent(in):: initbun! optional initialization bun

!EOP

   !--- local ---
!  type(cpl_bundle) :: bun_tmp   ! temporary bundle for subroutine
   integer(IN) :: n,m            ! generic indicies
   integer(IN) :: npts           ! number of points (local) in an aVect field
   integer(IN) :: nfld           ! number of fields (local) in an aVect field
   integer(IN) :: nfldi          ! number of fields (local) in an aVect field
   integer(IN) :: nptsx          ! number of points (local) in an aVect field
   integer(IN) :: nptsi          ! number of points (local) in an aVect field
   integer(IN) :: kfld           ! field number of fld1 in bun1
   integer(IN),dimension(:),allocatable :: kfldin   ! field numbers of bunlist in bun
   type(mct_list) :: blist   !  bunlist as a List
   type(mct_string) :: tattr !  an attribute

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_mult) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   nptsx = mct_aVect_lsize(bun1%data)
   npts  = mct_aVect_lsize(bun%data)
   if (nptsx /= npts) write(6,*) subName,' ERROR: npts error1 ',npts,nptsx

   if (present(initbun)) then
     nptsi = mct_aVect_lsize(initbun%data)
     if (nptsi /= npts) write(6,*) subName,' ERROR: nptsi error1 ',npts,nptsi
   endif

   kfld  = mct_aVect_indexRA(bun1%data,fld1,perrWith=subName)

   if (present(bunlist)) then
     call mct_list_init(blist,bunlist)

     nfld=mct_list_nitem(blist)

     allocate(kfldin(nfld))
     do m=1,nfld
       call mct_list_get(tattr,m,blist)
       kfldin(m) = mct_aVect_indexRA(bun%data,mct_string_toChar(tattr))
       call mct_string_clean(tattr)
     enddo
     call mct_list_clean(blist)


     if (present(initbun)) then
#ifdef CPP_VECTOR
       do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
       do n=1,npts
#else
       do n=1,npts
       do m=1,nfld
#endif
         bun%data%rAttr(kfldin(m),n) = initbun%data%rAttr(kfldin(m),n)*bun1%data%rAttr(kfld,n)
       enddo
       enddo
     else
#ifdef CPP_VECTOR
       do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
       do n=1,npts
#else
       do n=1,npts
       do m=1,nfld
#endif
         bun%data%rAttr(kfldin(m),n) = bun%data%rAttr(kfldin(m),n)*bun1%data%rAttr(kfld,n)
       enddo
       enddo
     endif

     deallocate(kfldin)

   else

     nfld  = mct_aVect_nRAttr(bun%data)

     if (present(initbun)) then
       nfldi  = mct_aVect_nRAttr(initbun%data)
       if (nfldi /= nfld) write(6,*) subName,' ERROR: nfldi error1 ',nfld,nfldi
#ifdef CPP_VECTOR
       do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
       do n=1,npts
#else
       do n=1,npts
       do m=1,nfld
#endif
         bun%data%rAttr(m,n) = initbun%data%rAttr(m,n)*bun1%data%rAttr(kfld,n)
       enddo
       enddo
     else
#ifdef CPP_VECTOR
       do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
       do n=1,npts
#else
       do n=1,npts
       do m=1,nfld
#endif
         bun%data%rAttr(m,n) = bun%data%rAttr(m,n)*bun1%data%rAttr(kfld,n)
       enddo
       enddo
     endif

   endif

end subroutine cpl_bundle_mult

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_divide - Divide a bundle by a field
!
! !DESCRIPTION:
!     Divide each attribute in {\tt bun} by the data in {\tt fld1}
!     from the input bundle {\tt bun1}
!
! !REVISION HISTORY:
!     2002-Oct-15 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_divide(bun,bun1,fld1)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle),intent(inout) :: bun    ! bundle output
   type(cpl_bundle),intent(in)    :: bun1   ! bundle input
   character(*)    ,intent(in)    :: fld1   ! bun1 field name

!EOP

   !--- local ---
   integer(IN) :: n,m           ! generic indicies
   integer(IN) :: npts          ! number of points (local) in an aVect field
   integer(IN) :: nfld          ! number of fields (local) in an aVect field
   integer(IN) :: nptsx         ! number of points (local) in an aVect field
   integer(IN) :: kfld          ! field number of fld1 in bun1
!  real(R8)    :: recip         ! reciprical
   real(R8),allocatable  :: recip(:)       ! reciprical

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_divide) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   npts  = mct_aVect_lsize( bun%data)
   nfld  = mct_aVect_nRAttr(bun%data)

   nptsx = mct_aVect_lsize(bun1%data)
   if (nptsx /= npts) write(6,*) subName,' ERROR: npts error1 ',npts,nptsx

   kfld  = mct_aVect_indexRA(bun1%data,fld1,perrWith=subName)


!  do n=1,npts
!    recip = 0.0_R8
!    if (bun1%data%rAttr(kfld,n) /= 0.0) recip = 1.0_R8/bun1%data%rAttr(kfld,n)
!    do m=1,nfld
!      bun%data%rAttr(m,n) = bun%data%rAttr(m,n)*recip
!    enddo
!  enddo

   allocate(recip(npts))

   do n=1,npts
     recip(n) = 0.0_R8
     if (bun1%data%rAttr(kfld,n) /= 0.0_R8) recip(n) = 1.0_R8/bun1%data%rAttr(kfld,n)
   enddo

#ifdef CPP_VECTOR
   do m=1,nfld
!CDIR SELECT(VECTOR)
!DIR$ CONCURRENT
   do n=1,npts
#else
   do n=1,npts
   do m=1,nfld
#endif
       bun%data%rAttr(m,n) = bun%data%rAttr(m,n)*recip(n)
   enddo
   enddo

   deallocate(recip)


end subroutine cpl_bundle_divide

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_bundle_gsum -  Calculate global sum and output
!
! !DESCRIPTION:
!     Calculate the global sum of all the attributes in input
!     bundle {\tt bun} and output the results to stdout.
!
!     The data in {\tt bun} can be weighted by the optional arguments.
!     If present, the data in {\tt bun} will be multiplied by the
!     data from {\tt fld1} of {\tt AV1} and data from {\tt fld2} of
!     {\tt AV2} and the scalar {\tt scalar}.
!
!     If optional argument {\tt istr} is present, it will be included
!     in the output string between the global sum and the name of the
!     attribute.
!
! !REVISION HISTORY:
!     2003-Jan-2 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bundle_gsum(bun,AV1,fld1,AV2,fld2,scalar,istr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_bundle)            ,intent(in)    :: bun    ! bundle input
   type(mct_aVect),optional,intent(in)    :: AV1    ! weight bundle input
   character(*)       ,optional,intent(in)    :: fld1   ! AV1 field name
   type(mct_aVect),optional,intent(in)    :: AV2    ! weight bundle input
   character(*)       ,optional,intent(in)    :: fld2   ! AV2 field name
   real(R8)           ,optional,intent(in)    :: scalar ! scalar for weights
   character(*)       ,optional,intent(in)    :: istr   ! string for print

!EOP

   !--- local ---
   integer(IN) :: i,j               ! generic indicies
   integer(IN) :: AVsiz             ! number of points (local) in an aVect field
   integer(IN) :: AVnum             ! number of fields (local) in an aVect field
   integer(IN) :: AV1siz            ! number of points in AV1
   integer(IN) :: AV2siz            ! number of points in AV2
   integer(IN) :: AV1fld            ! location of fld1 in AV1
   integer(IN) :: AV2fld            ! location of fld2 in AV2
   integer(IN) :: AVsizg            ! number of points (local) in gData
   integer(IN) :: AVnumg            ! number of fields (local) in gData
   type(mct_string) :: item     ! mct string
   character(CL)        :: itemc    ! item converted to char

   real(R8),allocatable :: gsum(:)    ! array for diagnostics
   real(R8),allocatable :: gsumall(:) ! array for diagnostics
   real(R8),allocatable :: weight(:)  ! weights array
   type(cpl_bundle)     :: bun_local  ! local temporary bundle
   type(mct_aVect)  :: gData      ! global/gathered bundle data for bfb sum

   integer(IN),parameter  :: pid0 = 0 ! root process pid = zero
   integer(IN)            :: rcode    ! error code
 
   !--- formats ---
   character(*),parameter :: subName = '(cpl_bundle_gsum) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

     AVnum  = mct_aVect_nRAttr(bun%data)
     AVsiz  = mct_aVect_lsize(bun%data)

     allocate(gsum(AVnum))
     allocate(gsumall(AVnum))
     allocate(weight(AVsiz))
     gsum = 0.0_R8
     gsumall = 0.0_R8
     weight = 1.0_R8
     if (present(scalar)) then
       weight = scalar
     endif

     if (present(AV1)) then
       AV1siz = mct_aVect_lsize(AV1)
       AV1fld = mct_aVect_indexRA(AV1,fld1,perrWith=subName)
       if (AV1siz /= AVsiz) then
         write(6,*) trim(subName),'ERROR AV1siz,AVsiz ',AV1siz,AVsiz
         call shr_sys_flush(6)
         call shr_sys_abort(subName//" ERROR AV1siz")
       endif
       if (fld1 == 'mask') then
         do i=1,AVsiz
           if (abs(AV1%rattr(AV1fld,i)) <= 1.0e-06_R8) weight(i) = 0.0_R8
         enddo
       else
         do i=1,AVsiz
           weight(i) = weight(i)*AV1%rattr(AV1fld,i)
         enddo
       endif
     endif

     if (present(AV2)) then
       AV2siz = mct_aVect_lsize(AV2)
       AV2fld = mct_aVect_indexRA(AV2,fld2,perrWith=subName)
       if (AV2siz /= AVsiz) then
         write(6,*) trim(subName),'ERROR AV2siz,AVsiz ',AV2siz,AVsiz
         call shr_sys_flush(6)
         call shr_sys_abort(subName//" ERROR AV2siz")
       endif
       if (fld2 == 'mask') then
         do i=1,AVsiz
           if (abs(AV2%rattr(AV2fld,i)) <= 1.0e-06_R8) weight(i) = 0.0_R8
         enddo
       else
         do i=1,AVsiz
           weight(i) = weight(i)*AV2%rattr(AV2fld,i)
         enddo
       endif
     endif

     if (bfbflag) then
       call cpl_bundle_initv(bun_local,'bun_local',bun,bun%dom)
!DIR$ CONCURRENT
       do i=1,AVsiz
       do j=1,AVnum
         if (bun%data%rattr(j,i) > 1.01_R8*cpl_const_spval .or. &
             bun%data%rattr(j,i) < 0.99_R8*cpl_const_spval) then
!          bun_local%data%rattr(j,i) = bun%data%rattr(j,i) * weight(i)
           bun_local%data%rattr(j,i) = bun%data%rattr(j,i)
         else
           bun_local%data%rattr(j,i) = 0.0_R8
         endif
       enddo
       enddo
       call mct_aVect_gather(bun_local%data,gData,bun_local%dom%gsMap,pid0,cpl_comm_comp,rcode)
       AVnumg  = mct_aVect_nRAttr(gData)
       AVsizg  = mct_aVect_lsize(gData)
       do i=1,AVsizg
       do j=1,AVnumg
           gsumall(j) = gsumall(j) + gData%rattr(j,i)
       enddo
       enddo
       call cpl_bundle_clean(bun_local)
       call mct_aVect_clean(gData)
     else
!CDIR NOASSOC
       do i=1,AVsiz
       do j=1,AVnum
!        if (bun%data%rattr(j,i) /= cpl_const_spval) then
         if (bun%data%rattr(j,i) > 1.01_R8*cpl_const_spval .or. &
             bun%data%rattr(j,i) < 0.99_R8*cpl_const_spval) then
           gsum(j) = gsum(j) + bun%data%rattr(j,i) * weight(i)
         endif
       enddo
       enddo
       call shr_mpi_sum(gsum,gsumall,cpl_comm_comp,subName)
     endif

     if (cpl_comm_comp_pid == pid0) then
       do j=1,AVnum
         call mct_aVect_getRList(item,j,bun%data)
         itemc = mct_string_toChar(item)
         call mct_string_clean(item)
         if (present(istr)) then
           write(6,100) 'xxx','sorr',j,gsumall(j),trim(istr),trim(itemc)
         else
           write(6,101) 'xxx','sorr',j,gsumall(j),trim(itemc)
         endif
       enddo
     endif
     call shr_sys_flush(6)
     deallocate(gsum)
     deallocate(gsumall)
     deallocate(weight)

100  format('comm_diag ',a3,1x,a4,1x,i3,es26.19,1x,a,1x,a)
101  format('comm_diag ',a3,1x,a4,1x,i3,es26.19,1x,a)


end subroutine cpl_bundle_gsum

!===============================================================================
! !BOP ===========================================================================
!
! !IROUTINE: cpl_bun_bun_op - 
!
! !DESCRIPTION:
!     Experimental routine to provide to provide add/mult/subtract/div
!     capability for bundles.   Not finished.
!
! !REVISION HISTORY:
!     2002-Oct-15 - T. Craig -- initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_bun_bun_op(op,istate,bun0,fld0,bun1,fld1,bun2,fld2,bun3,fld3, &
                         &          bun4,fld4,bun5,fld5,bun6,fld6,scalar)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)             ,intent(in)    :: op     ! operation type
   integer(IN)              ,intent(in)    :: istate ! use initial state 0=false

   type(cpl_bundle)         ,intent(inout) :: bun0   ! bundle output
   character(*)    ,optional,intent(in)    :: fld0   ! bun  field name
   type(cpl_bundle)         ,intent(in)    :: bun1   ! bundle input
   character(*)    ,optional,intent(in)    :: fld1   ! bun1 field name
   type(cpl_bundle),optional,intent(in)    :: bun2   ! bundle input
   character(*)    ,optional,intent(in)    :: fld2   ! bun2 field name
   type(cpl_bundle),optional,intent(in)    :: bun3   ! bundle input
   character(*)    ,optional,intent(in)    :: fld3   ! bun3 field name
   type(cpl_bundle),optional,intent(in)    :: bun4   ! bundle input
   character(*)    ,optional,intent(in)    :: fld4   ! bun4 field name
   type(cpl_bundle),optional,intent(in)    :: bun5   ! bundle input
   character(*)    ,optional,intent(in)    :: fld5   ! bun5 field name
   type(cpl_bundle),optional,intent(in)    :: bun6   ! bundle input
   character(*)    ,optional,intent(in)    :: fld6   ! bun6 field name
   real(R8)        ,optional,intent(in)    :: scalar ! scalar multiplier

! !EOP

   !--- local ---
   integer(IN),parameter :: mbun=6 ! maximum number of bundles incoming

   integer(IN) :: n,m               ! generic indicies
   integer(IN) :: mstart,mend       ! indices for field start,end
   integer(IN) :: nbun              ! number of RHS bundles sent
   integer(IN) :: npts(0:mbun)      ! number of points (local) in an aVect field
   integer(IN) :: nfld(0:mbun)      ! number of fields
   integer(IN) :: kfld(0:mbun)      ! field number in bundle
   logical     :: onefield(0:mbun)  ! is operation on 1 or all fields
   real(R8)    :: rtmp              ! temporary storage

   !--- formats ---
   character(*),parameter :: subName = '(cpl_bun_bun_op) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   npts = -1
   kfld = -1
   onefield = .true.
   nbun = 1
   if (.not.present(fld0)) onefield(0)=.false.
   if (.not.present(fld1)) onefield(1)=.false.

   npts(0) = mct_aVect_lsize(bun0%data)
   nfld(0) = mct_aVect_nRAttr(bun0%data)
   if (onefield(0)) kfld(0) = mct_aVect_indexRA(bun0%data,fld0, &
                                                    perrWith=subName)

   npts(1) = mct_aVect_lsize(bun1%data)
   nfld(1) = mct_aVect_nRAttr(bun1%data)
   if (onefield(1)) kfld(1) = mct_aVect_indexRA(bun1%data,fld1, &
                                                    perrWith=subName)

   if (present(bun2)) then
     nbun = 2
     npts(2) = mct_aVect_lsize(bun2%data)
     if (onefield(2)) kfld(2) = mct_aVect_indexRA(bun2%data,fld2,  &
                                                      perrWith=subName)
     if (npts(2) /= npts(0)) write(6,*) subName,' ERROR: npts2 ',npts(0),npts(2)
   endif

   if (present(bun3)) then
     nbun = 3
     npts(3) = mct_aVect_lsize(bun3%data)
     if (onefield(3)) kfld(3) = mct_aVect_indexRA(bun3%data,fld3,  &
                                                      perrWith=subName)
     if (npts(3) /= npts(0)) write(6,*) subName,' ERROR: npts3 ',npts(0),npts(3)
   endif

   if (present(bun4)) then
     nbun = 4
     npts(4) = mct_aVect_lsize(bun4%data)
     if (onefield(4)) kfld(4) = mct_aVect_indexRA(bun4%data,fld4,  &
                                                      perrWith=subName)
     if (npts(4) /= npts(0)) write(6,*) subName,' ERROR: npts4 ',npts(0),npts(4)
   endif

   if (present(bun5)) then
     nbun = 5
     npts(5) = mct_aVect_lsize(bun5%data)
     if (onefield(5)) kfld(5) = mct_aVect_indexRA(bun5%data,fld5,  &
                                                      perrWith=subName)
     if (npts(5) /= npts(0)) write(6,*) subName,' ERROR: npts5 ',npts(0),npts(5)
   endif

   if (present(bun6)) then
     nbun = 6
     npts(6) = mct_aVect_lsize(bun6%data)
     if (onefield(6)) kfld(6) = mct_aVect_indexRA(bun6%data,fld6,  &
                                                      perrWith=subName)
     if (npts(6) /= npts(0)) write(6,*) subName,' ERROR: npts6 ',npts(0),npts(6)
   endif

   if (onefield(0)) then
     mstart = kfld(0)
     mend   = kfld(0)
   else
     mstart = 1
     mend   = nfld(0)
   endif

   do n=1,npts(0)
     rtmp = 1.0_R8
     if (onefield(1))     rtmp = rtmp * bun1%data%rAttr(kfld(1),n)
     if (present(bun2))   rtmp = rtmp * bun2%data%rAttr(kfld(2),n)
     if (present(bun3))   rtmp = rtmp * bun3%data%rAttr(kfld(3),n)
     if (present(bun4))   rtmp = rtmp * bun4%data%rAttr(kfld(4),n)
     if (present(bun5))   rtmp = rtmp * bun5%data%rAttr(kfld(5),n)
     if (present(bun6))   rtmp = rtmp * bun6%data%rAttr(kfld(6),n)
     if (present(scalar)) rtmp = rtmp * scalar

     do m=mstart,mend
       if (.not.onefield(1)) rtmp = rtmp * bun1%data%rAttr(m,n)
       if (istate.eq.0) then
         if (op.eq.'+') bun0%data%rAttr(m,n) = rtmp
         if (op.eq.'-') bun0%data%rAttr(m,n) = -rtmp
         if (op.eq.'*') bun0%data%rAttr(m,n) = rtmp
         if (op.eq.'/') bun0%data%rAttr(m,n) = 1.0_R8/(rtmp)
       else
         if (op.eq.'+') bun0%data%rAttr(m,n) = bun0%data%rAttr(m,n) + rtmp
         if (op.eq.'-') bun0%data%rAttr(m,n) = bun0%data%rAttr(m,n) - rtmp
         if (op.eq.'*') bun0%data%rAttr(m,n) = bun0%data%rAttr(m,n) * rtmp
         if (op.eq.'/') bun0%data%rAttr(m,n) = bun0%data%rAttr(m,n) / rtmp
       endif
     enddo

   enddo

   bun0%cnt = 1

end subroutine cpl_bun_bun_op

!===============================================================================

end module cpl_bundle_mod

