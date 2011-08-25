!===============================================================================
! SVN $Id: dshr_bundle.F90 3286 2007-02-23 23:59:03Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_bundle.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_bundle -- bundle data type and associated methods
!
! !DESCRIPTION:
!    Bundle data type and associated methods.
!
! !REVISION HISTORY:
!     2006-Oct-30 - B. Kauffman - add support for distributed bcast/gather/scatter
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_bundle

! !USES:

   use shr_sys_mod    ! shared system calls
   use shr_string_mod ! string methods
   use shr_cal_mod    ! calendar methods
   use shr_mpi_mod    ! mpi wrappers

   use dshr_kind      ! F90 kinds
   use dshr_domain    ! domain type and methods
   use dshr_const     ! constants

   implicit none
   private

! !PUBLIC TYPES:

   public :: dshr_bundle_bundleType  ! bundle data type

   type dshr_bundle_bundleType  ! bundle of fields with same domain & date
      private
      character(CS)                   :: initialized ! T/F char string
      character(CS)                   :: dateset     ! T/F char string
      character(CS)                   :: name        ! bundle name
      integer(IN)                     :: ni          ! size of i dimension
      integer(IN)                     :: nj          ! size of j dimension
      integer(IN)                     :: nf          ! number of fields
      integer(IN)                     :: cDate       ! data date, yyyymmdd
      integer(IN)                     :: sec         ! elapsed secs on date
      real(R8),pointer                :: data(:,:,:) ! field data (ni,nj,nf)
      character(CX)                   :: fieldList   ! list of data fields
      type(dshr_domain_domainType),pointer :: domain ! associated domain
   end type dshr_bundle_bundleType

! !PUBLIC MEMBER FUNCTIONS:

   public :: dshr_bundle_create          ! Initialize/allocate bundle
   public :: dshr_bundle_clean           ! Clean/deallocate bundle
   public :: dshr_bundle_fill            ! Fill bundle with internal data
   public :: dshr_bundle_checkInit       ! Check if bundle has been initialized
   public :: dshr_bundle_checkDate       ! Check if bundle date has been set

   public :: dshr_bundle_getName         ! Get bundle name
   public :: dshr_bundle_getNumFields    ! Get number of fields in bundle
   public :: dshr_bundle_getFieldList    ! Copy field list out of bundle
   public :: dshr_bundle_getFieldIndex   ! Get index of field, given name
   public :: dshr_bundle_getFieldName    ! Get name of field, given index

   public :: dshr_bundle_getField        ! Copy data out of bundle
   public :: dshr_bundle_putField        ! Copy data into bundle
   public :: dshr_bundle_getDate         ! Copy date out of bundle
   public :: dshr_bundle_putDate         ! Copy date into bundle
   public :: dshr_bundle_getDims         ! Return dims of data in bundle

   public :: dshr_bundle_assignPtr       ! Pointer into bundle data
   public :: dshr_bundle_domainPtr       ! Pointer into bundle domain
   public :: dshr_bundle_copyFields      ! copy all fields shared by two bundles
   public :: dshr_bundle_print           ! Print bundle info
   public :: dshr_bundle_writeCDF        ! Write bundle to CDF file

   public :: dshr_bundle_bcastData       ! COLLECTIVE: b-cast a bundle data
   public :: dshr_bundle_gather          ! COLLECTIVE: gather  local  to global
   public :: dshr_bundle_scatter         ! COLLECTIVE: scatter global to local
   public :: dshr_bundle_extractLocal    ! extract local bundle from global

   public :: dshr_bundle_setAbort        ! set local abort flag
   public :: dshr_bundle_setDebug        ! set local debug flag
   public :: dshr_bundle_getDebug        ! get local debug flag

! !PUBLIC DATA MEMBERS:

   ! no public data members

!EOP

   interface dshr_bundle_assignPtr ; module procedure &
      dshr_bundle_assignPtrField, &
      dshr_bundle_assignPtrIndex, &
      dshr_bundle_assignPtrAll3D
   end interface 

   character(CS),parameter :: dshr_bundle_setTrue  = 'TRUE'
   character(CS),parameter :: dshr_bundle_setFalse = 'FALSE'
   real(R8)     ,parameter :: pi                   = dshr_const_pi

   logical      ,save      :: doabort              = .true.
   integer(IN)  ,save      :: debug                = 0

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_create -- Initialize bundle, allocate memory
!
! !DESCRIPTION:
!     Initialize bundle, allocate memory, init data to special-value
!     \newline
!     call dshr\_bundle\_init(bundle,datm\_data\_domain,"bun\_ocn","u:v:s:t")
!
! !REVISION HISTORY:
!     2005-Mar-3 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_create(bundle,domain,name,fieldList,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType)       ,intent(out) :: bundle    ! bundle to allocate
   type(dshr_domain_domainType),target,intent(in)  :: domain    ! bundle's domain
   character(*)                       ,intent(in)  :: name      ! bundle's name
   character(*)                       ,intent(in)  :: fieldList ! list of data fields
   integer(IN),optional               ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   integer(IN)            :: ni,nj    ! size of bundle dimensions
   integer(IN)            :: nf       ! number of fields
   integer(IN)            :: rcode    ! return code
   integer(IN)            :: lrc      ! local rc

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_create)"
   character(*),parameter :: F00     = "('(dshr_bundle_create) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_create) ',2a,3i6,2x,a)"
   character(*),parameter :: F02     = "('(dshr_bundle_create) ',a,a15,i6,a,2i5,2a)"

!-------------------------------------------------------------------------------
! Notes: allocates memory for a bundle
!-------------------------------------------------------------------------------

   lrc = 0

   !--- verify field list is valid ---
   if ( .not. shr_string_listIsValid(fieldList)) then
      write(6,F00) "bundle ",trim(name)," has invalid fieldList: ",trim(fieldList)
      call dshr_bundle_abort(subName//" invalid field list")
   end if

   !--- get spatial dims and number of fields ---
   call dshr_domain_getDims(domain,ni,nj,rc)
   nf = shr_string_listGetNum(fieldList)
   write(6,F02) "name =",trim(name),ni," x",nj,nf," fields ",trim(fieldList)
   
   !--- allocate memory, clean (dealloc) if bundle is previously initialized ---
   if (dshr_bundle_checkInit(bundle)) call dshr_bundle_clean(bundle)
   allocate(bundle%data(ni,nj,nf))

   !--- set initial values ---
   bundle%initialized =  dshr_bundle_setTrue
   bundle%dateset     =  dshr_bundle_setFalse
   bundle%name        =  name
   bundle%ni          =  ni
   bundle%nj          =  nj
   bundle%nf          =  nf
   bundle%cDate       =   0
   bundle%sec         =   0
   bundle%data        =  dshr_const_spval
   bundle%fieldList   =  fieldList
   bundle%domain      => domain

   if (present(rc)) rc = lrc

end subroutine dshr_bundle_create

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_clean -- clean, deallocate bundle
!
! !DESCRIPTION:
!     Clean, deallocate bundle
!     \newline
!     call dshr\_bundle\_clean(bun,rc)
!
! !REVISION HISTORY:
!     2005-Mar-28 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_clean(bundle,rc)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(inout) :: bundle ! bundle
   integer(IN),optional        ,intent(out)   :: rc     ! return code

!EOP

   !--- local ---
   integer :: rCode  ! return code

   !--- formats ---
   character(*),parameter :: subName =   "(dshr_bundle_clean)"
   character(*),parameter :: F00     = "('(dshr_bundle_clean) ',a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rCode = 0

   deallocate(bundle%data,stat=rCode)
   if (debug > 0.and. rCode > 0) write(6,F00) 'ERROR deallocating bundle data'

   !--- set initial values ---
   bundle%initialized =  dshr_bundle_setFalse 
   bundle%dateset     =  dshr_bundle_setFalse 
   bundle%name        =  "uninitiallized"
   bundle%ni          =  -1
   bundle%nj          =  -1
   bundle%nf          =  -1
   bundle%cDate       =  -1
   bundle%sec         =  -1
!  bundle%data        =  dshr_const_spval ! this data is deallocated
   bundle%fieldList   =  " "
   bundle%domain      => null()

   if (present(rc)) rc = rCode

end subroutine dshr_bundle_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_fill -- Fill bundle data with test data
!
! !DESCRIPTION:
!     Fill bundle data with test data
!     \newline
!     call dshr\_bundle\_fill(bundle,type,rc)
!
! !REVISION HISTORY:
!     2005-Mar-23 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_fill(bundle,type,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(inout) :: bundle ! bundle to allocate/fill
   character(*),optional       ,intent(in)    :: type   ! type of fill 'index','sincos'
   integer(IN),optional        ,intent(out)   :: rc     ! return code

!EOP

   !----- local -----
   integer(IN)            :: ni,nj    ! size of bundle's i and j dimensions
   integer(IN)            :: nf       ! number of fields
   integer(IN)            :: n,i,j    ! loop indexes
   integer(IN)            :: rcode    ! return code
   character(CS)          :: ltype    ! local copy of type
   logical                :: firsterr ! error string

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_fill)"
   character(*),parameter :: F00     = "('(dshr_bundle_fill) ',a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   firsterr = .true.
   rcode = 0
   ltype = 'index'
   if (present(type)) ltype = type

   if (.not.dshr_bundle_checkInit(bundle)) then
     call dshr_bundle_abort(subName//' bundle not intialized '//trim(bundle%name))
     rcode = 1
     if (present(rc)) rc = rcode
     return
   endif

   ni = bundle%ni
   nj = bundle%nj
   nf = bundle%nf

   do n=1,nf
   do j=1,nj
   do i=1,ni
     if (ltype(1:5) == 'index') then
       bundle%data(i,j,n) = real(i*1000+j,R8) + real(n,R8)/100._R8 
     elseif (ltype(1:6) == 'sincos') then
       bundle%data(i,j,n) = cos(real(i,R8)/real(ni,R8)*20)*sin(real(j,R8)/real(nj,R8)*20._R8)*real(n,R8) 
     elseif (ltype(1:9) == 'vecttest1') then
       bundle%data(i,j,n) = 0._R8
       bundle%data(i,j,2) = 1._R8
     elseif (ltype(1:9) == 'vecttest2') then
       bundle%data(i,j,n) = cos(real(i,R8)/real(ni,R8)*2*pi)*sin((real(j-1,R8)/real(nj,R8))*20)*real(n,R8) 
       bundle%data(i,j,n) = 0._R8
       bundle%data(i,j,1) = cos(real(i,R8)/real(ni,R8)*2*pi-pi/2._R8)*real(n,R8) 
       bundle%data(i,j,2) = cos(real(i,R8)/real(ni,R8)*2*pi)*real(n,R8) 
     elseif (ltype(1:8) == 'vecttest') then
       bundle%data(i,j,n) = 1.0_R8
       if (n == 1) bundle%data(i,j,n) =                   real(j-1,R8)/real(nj-1,R8)
       if (n == 2) bundle%data(i,j,n) = 1.0_R8 - real(j-1,R8)/real(nj-1,R8)
     elseif (ltype(1:6) == 'nfield') then
       bundle%data(i,j,n) = n
     else
       if (firsterr) write(6,F00) ' WARNING type not allowed, using index '//trim(ltype)
       firsterr = .false.
       bundle%data(i,j,n) = real(i*1000+j,R8) + real(n,R8)/100._R8        
     endif
   enddo
   enddo
   enddo

   if (present(rc)) rc = rcode

end subroutine dshr_bundle_fill

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_checkInit -- check if bundle has been initialized
!
! !DESCRIPTION:
!  Check if bundle has been initialized
!
! !REVISION HISTORY:
!     2005-Mar-28 - T. Craig - First version
!
! !INTERFACE: ------------------------------------------------------------------

logical function dshr_bundle_checkInit(bun)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType) :: bun

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_checkInit)"
   character(*),parameter :: F00     = "('(dshr_bundle_checkInit) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   dshr_bundle_checkInit = .false.
   if (trim(bun%initialized) == trim(dshr_bundle_setTrue)) then
      dshr_bundle_checkInit = .true.
   endif

end function dshr_bundle_checkInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_checkDate -- check if bundle date has been set
!
! !DESCRIPTION:
!  Check if bundle date has been set
!
! !REVISION HISTORY:
!     2005-May-28 - T. Craig - First version
!
! !INTERFACE: ------------------------------------------------------------------

logical function dshr_bundle_checkDate(bun)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType) :: bun

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_checkDate)"
   character(*),parameter :: F00     = "('(dshr_bundle_checkDate) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   dshr_bundle_checkDate = .false.
   if (trim(bun%dateset) == trim(dshr_bundle_setTrue)) then
      dshr_bundle_checkDate = .true.
   endif

end function dshr_bundle_checkDate

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getName -- Get bundle name
!
! !DESCRIPTION:
!     Get bundle name
!     \newline
!     name = dshr\_bundle\_getName(bun)
!
! !REVISION HISTORY:
!     2005-Aug-01 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

character(CS) function dshr_bundle_getName(bundle)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to put data into

!EOP

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_getName)"
   character(*),parameter :: F00     = "('(dshr_bundle_getName) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   dshr_bundle_getName = bundle%name

end function dshr_bundle_getName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getNumFields -- Get number of fields in bundle
!
! !DESCRIPTION:
!     Get number of fields in bundle
!     \newline
!     call dshr\_bundle\_getNumFields(bun,k,rc)
!
! !REVISION HISTORY:
!     2005-May-05 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_getNumFields(bundle,k,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to put data into
   integer(IN)                 ,intent(out) :: k         ! index of field
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   integer(IN)            :: rCode            ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_getNumFields)"
   character(*),parameter :: F00     = "('(dshr_bundle_getNumFields) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   k = shr_string_listGetNum(bundle%fieldList)

   if (present(rc)) rc = 0

end subroutine dshr_bundle_getNumFields

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getFieldList -- Get field list of bundle datatype
!
! !DESCRIPTION:
!     Get field name list out of bundle
!     \newline
!     call dshr\_bundle\_getfieldList(bun,fldlist,rc)
!
! !REVISION HISTORY:
!     2005-May-15 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_getFieldList (bundle,fldlist,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to get date from
   character(*)                ,intent(out) :: fldlist   ! fields
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   integer(IN) :: lrc

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_getFieldList)"
   character(*),parameter :: F00     = "('(dshr_bundle_getFieldList) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_getFieldList) ',a,i6)"

!-------------------------------------------------------------------------------
!  Notes:
!-------------------------------------------------------------------------------

   lrc = 0

   if (.not. dshr_bundle_checkInit(bundle)) then
      write(6,F00) "ERROR: bundle not initialized"
      call dshr_bundle_abort(subName//"ERROR: bundle not initialized")
   else if (len(fldlist) < len(trim(bundle%fieldList))) then
      write(6,F01) "ERROR: fldList var is too small to recieve fldList: ",len(fldlist)
      call dshr_bundle_abort(subName//" fldlist variable is too small")
   end if

   fldlist = trim(bundle%fieldList)

   if (present(rc)) rc = lrc

end subroutine dshr_bundle_getFieldList

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getFieldIndex -- Get index of field in bundle
!
! !DESCRIPTION:
!     Get index of field in bundle
!     \newline
!     call dshr\_bundle\_getFieldIndex(bun,"taux",k,rc)
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_getFieldIndex(bundle,fldStr,k,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to put data into
   character(*)                ,intent(in)  :: fldStr    ! name of field
   integer(IN)                 ,intent(out) :: k         ! index of field
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   integer(IN)            :: lrc              ! local rc
   logical                :: print            ! print fields not found?

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_getFieldIndex)"
   character(*),parameter :: F00     = "('(dshr_bundle_getFieldIndex) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   lrc = 0

   print = .false.
   if (debug>0) print = .true.

   if (.not. dshr_bundle_checkInit(bundle)) then
      write(6,F00) "ERROR: bundle not initialized"
      call dshr_bundle_abort(subName//" ERROR: bundle not initialized")
   end if

   call shr_string_listGetIndex(bundle%fieldList,fldstr,k,print,rc=lrc)

   if (present(rc)) rc = lrc

end subroutine dshr_bundle_getFieldIndex

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getFieldName -- Get name of field in bundle
!
! !DESCRIPTION:
!     Get name of field in bundle
!     \newline
!     call dshr\_bundle\_getFieldName(bun,k,name,rc)
!
! !REVISION HISTORY:
!     2005-May-05 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_getFieldName(bundle,k,name,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to put data into
   integer(IN)                 ,intent(in)  :: k         ! index of field
   character(*)                ,intent(out) :: name      ! name  of field
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   integer(IN)            :: rCode            ! return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_getFieldName)"
   character(*),parameter :: F00     = "('(dshr_bundle_getFieldName) ',4a)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   call shr_string_listGetName(bundle%fieldList,k,name,rCode)

   if (present(rc)) rc = rCode

end subroutine dshr_bundle_getFieldName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getField -- Copy data out of bundle
!
! !DESCRIPTION:
!     Copy data out of bundle
!     \newline
!     call dshr\_bundle\_getField (bun,taux(:,:),"taux") 
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_getField (bundle,data,fldStr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to get data from
   character(*)                ,intent(in)  :: fldStr    ! name of field to get
   real(R8)                    ,intent(out) :: data(:,:) ! data copied from bundle
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   integer(IN)            :: ni,nj,k  ! size of bundle's i, j & k dimensions
   character(CL)          :: idStr    ! bundle's name

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_getField)"
   character(*),parameter :: F00     = "('(dshr_bundle_getField) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_getField) ',a,i2,a)"
   character(*),parameter :: F02     = "('(dshr_bundle_getField) ',a,2i7)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0
   if (.not. dshr_bundle_checkInit(bundle)) &
      call dshr_bundle_abort(subName//" bundle not initialized")

   !--- confirm proper size of input data ---
   ni    = bundle%ni
   nj    = bundle%nj
   if (size(data,1) /= ni .or. size(data,2) /= nj) then
      write(6,F00) "ERROR: field=",trim(fldStr),", bundle=",trim(bundle%name)
      write(6,F00) "ERROR: input field is not the same size as the bundle"
      write(6,F02) "ERROR: bundle ni,nj= ",ni,nj
      write(6,F02) "ERROR: field  ni,nj= ",size(data,1),size(data,2)
      call dshr_bundle_abort(subName//" data size mismatch")
   end if

   !--- get data from bundle ---
   call dshr_bundle_getFieldIndex(bundle,fldStr,k,rc)
   if (k<1) then
      write(6,F00) "ERROR: field ",trim(fldStr)," is not in bundle ",trim(bundle%name)
      call dshr_bundle_abort(subName//" field not in bundle")
   endif
   data = bundle%data(:,:,k)

end subroutine dshr_bundle_getField

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_putField -- Copy data into bundle
!
! !DESCRIPTION:
!     Copy data into bundle
!     \newline
!     call dshr\_bundle\_putField(bun,taux(:,:),"taux")   ! copy data into bundle
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_putField (bundle,data,fldStr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(inout) :: bundle    ! bundle to put data into
   real(R8)                    ,intent(in)    :: data(:,:) ! data to get from bundle
   character(*)                ,intent(in)    :: fldStr    ! name of field
   integer(IN),optional        ,intent(out)   :: rc        ! return code

!EOP

   !----- integer -----
   integer(IN)            :: ni,nj,k  ! size of bundle's i, j & k dimensions

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_putField)"
   character(*),parameter :: F00     = "('(dshr_bundle_putField) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_putField) ',a,i2,a)"
   character(*),parameter :: F02     = "('(dshr_bundle_putField) ',a,2i7)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0
   if (.not. dshr_bundle_checkInit(bundle)) &
      call dshr_bundle_abort(subName//" bundle not initialized")

   !--- confirm proper size of input data ---
   ni = bundle%ni
   nj = bundle%nj
   if (size(data,1) /= ni .or. size(data,2) /= nj) then
      write(6,F00) "ERROR: field=",trim(fldStr),", bundle=",trim(bundle%name)
      write(6,F00) "ERROR: input field is not the same size as the bundle"
      write(6,F02) "ERROR: bundle ni,nj= ",ni,nj
      write(6,F02) "ERROR: field  ni,nj= ",size(data,1),size(data,2)
      call dshr_bundle_abort(subName//" data size mismatch")
   end if

   !--- copy data into bundle ---
   call dshr_bundle_getFieldIndex(bundle,fldStr,k,rc)
   if (k<1) then
      write(6,F00) "ERROR: field ",trim(fldStr)," is not in bundle ",trim(bundle%name)
      call dshr_bundle_abort(subName//" field not in bundle")
   endif
   bundle%data(:,:,k) = data

end subroutine dshr_bundle_putField

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getDate -- Copy date out of bundle
!
! !DESCRIPTION:
!     Copy date out of bundle
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_getDate (bundle,cDate,sec,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to get date from
   integer(IN)                 ,intent(out) :: cDate     ! date copied from bundle
   integer(IN)                 ,intent(out) :: sec       ! sec copied from bundle
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_getDate)"
   character(*),parameter :: F00     = "('(dshr_bundle_getDate) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_getDate) ',a,i8)"

!-------------------------------------------------------------------------------
!  Notes:
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0
   if (.not. dshr_bundle_checkInit(bundle)) &
      call dshr_bundle_abort(subName//" bundle not initialized")

   cDate = bundle%cDate
   sec   = bundle%sec 

   if (.not.shr_cal_validDate(cDate)) then
     write(6,F01) ' ERROR invalid date = ',cDate
     write(6,F01) ' ERROR bundle name  = ',trim(bundle%name)
     call dshr_bundle_abort(trim(subName)//' illegal date')
   endif

end subroutine dshr_bundle_getDate

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_putDate -- Copy date into bundle
!
! !DESCRIPTION:
!     Copy date into bundle
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_putDate (bundle,cDate,sec,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(inout) :: bundle    ! bundle to get date from
   integer(IN)                 ,intent(in)    :: cDate     ! date copied from bundle
   integer(IN)                 ,intent(in)    :: sec       ! sec copied from bundle
   integer(IN),optional        ,intent(out)   :: rc        ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_putDate)"
   character(*),parameter :: F00     = "('(dshr_bundle_putDate) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_putDate) ',a,i8)"

!-------------------------------------------------------------------------------
! Notes
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0
   if (.not. dshr_bundle_checkInit(bundle)) &
      call dshr_bundle_abort(subName//" bundle not initialized")

   if (.not.shr_cal_validDate(cDate)) then
     write(6,F01) ' ERROR invalid date ',cDate
     call dshr_bundle_abort(trim(subName)//' illegal date')
   endif

   bundle%dateSet = dshr_bundle_setTrue
   bundle%cDate  = cDate
   bundle%sec    = sec  

end subroutine dshr_bundle_putDate

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getDims -- Get dimensions of bundle data
!
! !DESCRIPTION:
!     Get dimensions of bundle data
!     \newline
!     call dshr\_bundle\_getDims(bun,ni,nj,nf,rc)
!
! !REVISION HISTORY:
!     2005-Apr-15 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_getDims (bundle,ni,nj,nf,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to get date from
   integer(IN)                 ,intent(out) :: ni        ! i size
   integer(IN)                 ,intent(out) :: nj        ! j size
   integer(IN)                 ,intent(out) :: nf        ! field size
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   integer(IN) :: lrc

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_getDims)"
   character(*),parameter :: F00     = "('(dshr_bundle_getDims) ',4a)"

!-------------------------------------------------------------------------------
!  Notes:
!-------------------------------------------------------------------------------

   lrc = 0
   if (.not. dshr_bundle_checkInit(bundle)) &
      call dshr_bundle_abort(subName//" bundle not initialized")

   ni = bundle%ni
   nj = bundle%nj
   nf = bundle%nf

   if (present(rc)) rc = lrc

end subroutine dshr_bundle_getDims

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_assignPtrField -- Assign a pointer to a field in a bundle
!
! !DESCRIPTION:
!     Assign a pointer to a field in a bundle
!     \newline
!     call dshr\_bundle\_assignPtrF(bun,"taux",taux\_ptr,rc)
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_assignPtrField(bundle,fldStr,fldPtr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)    :: bundle      ! bundle to point into
   character(*)                ,intent(in)    :: fldStr      ! name of field
   real(R8)                    ,pointer       :: fldPtr(:,:) ! pointer to field
   integer(IN),optional        ,intent(out)   :: rc          ! return code

!EOP

   !----- local -----
   integer(IN)            :: k    ! bundle's field dimension
   integer(IN)            :: lrc  ! local return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_assignPtrField)"
   character(*),parameter :: F00     = "('(dshr_bundle_assignPtrField) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lrc = 0
   if (.not. dshr_bundle_checkInit(bundle)) then
      call dshr_bundle_abort(subName//" bundle not initialized")
      lrc = 1
      fldPtr => null()
      if (present(rc)) rc = lrc
      return
   endif

   !--- assign pointer ---
   call dshr_bundle_getFieldIndex(bundle,fldStr,k,lrc)
   if (k<1) then
      write(6,F00) "ERROR: field ",trim(fldStr)," is not in bundle ",trim(bundle%name)
      fldPtr => null()
      lrc = 1
      call dshr_bundle_abort(subName//" field not in bundle")
   else
      fldPtr => bundle%data(:,:,k)
   endif

   if (present(rc)) rc = lrc

end subroutine dshr_bundle_assignPtrField

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_assignPtrAll3D -- Assign a pointer to the 3d data in a bundle
!
! !DESCRIPTION:
!     Assign a pointer to the 3d data in the bundle
!     \newline
!     call dshr\_bundle\_assignPtrAll3D(bun,data\_ptr,rc)
!
! !REVISION HISTORY:
!     2005-Apr-15 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_assignPtrAll3D(bundle,A3dPtr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle        ! bundle to point into
   real(R8)                    ,pointer     :: A3dPtr(:,:,:) ! pointer to array
   integer(IN),optional        ,intent(out) :: rc            ! return code

!EOP

   !----- local -----
   integer(IN)            :: lrc  ! local return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_assignPtrAll3D)"
   character(*),parameter :: F00     = "('(dshr_bundle_assignPtrAll3D) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lrc = 0
   if (.not. dshr_bundle_checkInit(bundle)) then
      call dshr_bundle_abort(subName//" bundle not initialized")
      A3dPtr => null()
      lrc = 1
   else
      A3dPtr => bundle%data(:,:,:)
   endif
   if (present(rc)) rc = lrc

end subroutine dshr_bundle_assignPtrAll3D

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_assignPtrIndex -- Assign a pointer to an index in a bundle
!
! !DESCRIPTION:
!     Assign a pointer to a field in a bundle via index
!     \newline
!     call dshr\_bundle\_assignPtrI(bun,2,data\_ptr,rc)
!
! !REVISION HISTORY:
!     2005-Mar-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_assignPtrIndex(bundle,index,fldPtr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)    :: bundle      ! bundle to point into
   integer(IN)                 ,intent(in)    :: index       ! index in bundle
   real(R8)                    ,pointer       :: fldPtr(:,:) ! pointer to field
   integer(IN),optional        ,intent(out)   :: rc          ! return code

!EOP

   !----- local -----
   integer(IN)            :: lrc  ! local return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_assignPtrIndex)"
   character(*),parameter :: F00     = "('(dshr_bundle_assignPtrIndex) ',a,2i8)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lrc = 0
   if (.not. dshr_bundle_checkInit(bundle)) then
      call dshr_bundle_abort(subName//" bundle not initialized")
      lrc = 1
      fldPtr => null()
      if (present(rc)) rc = lrc
      return
   endif

   !--- assign pointer ---
   if (index > bundle%nf) then
      write(6,F00) "ERROR: index is too large for bundle "//trim(bundle%name),index,bundle%nf
      lrc = 2
      fldPtr => null()
      call dshr_bundle_abort(subName//" ERROR index not in bundle")
   else if (index < 1) then
      write(6,F00) "ERROR: index must be > 0"
      lrc = 3
      fldPtr => null()
      call dshr_bundle_abort(subName//" ERROR: index must be > 0")
   else
      fldPtr => bundle%data(:,:,index)
   endif

   if (present(rc)) rc = lrc

end subroutine dshr_bundle_assignPtrIndex

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_domainPtr -- Assign a pointer to the domain in a bundle
!
! !DESCRIPTION:
!     Assign a pointer to the domain in a bundle
!     \newline
!     call dshr\_bundle\_domainPtr(bun,domain\_ptr,rc)
!
! !REVISION HISTORY:
!     2005-Feb-28 - B. Kauffman and J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_domainPtr(bundle,domPtr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle  ! bundle to point into
   type(dshr_domain_domainType),pointer     :: domPtr  ! pointer to domain
   integer(IN),optional        ,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(IN)            :: lrc  ! local return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_domainPtr)"
   character(*),parameter :: F00     = "('(dshr_bundle_domainPtr) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lrc = 0
   if (.not. dshr_bundle_checkInit(bundle)) then
      call dshr_bundle_abort(subName//" bundle not initialized")
      lrc = 1
      if (present(rc)) rc = lrc
      return
   endif

   !--- assign pointer ---
   domPtr => bundle%domain

   if (present(rc)) rc = lrc

end subroutine dshr_bundle_domainPtr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_copyFields -- copy all common fields from one bundle to another
!
! !DESCRIPTION:
!    copy all matching fields from input to output bundle
!
! !REVISION HISTORY:
!    2005-Nov-01  - B. Kauffman - relocated from datm_phys to dshr_bundle
!    2005-Jun-07  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_copyFields(bunSrc,bunDst,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bunSrc ! ptr to input  bundle
   type(dshr_bundle_bundleType),intent(out) :: bunDst ! ptr to output bundle
   integer(IN),optional        ,intent(out) :: rc     ! return code

!EOP

   !--- local ---
   real(R8),pointer :: dataSrc(:,:)    ! data from source      bundle
   real(R8),pointer :: dataDst(:,:)    ! data from destination bundle
   integer(IN)      :: kFlds           ! number of fields in a bundle
   integer(IN)      :: k_src,k_dst     ! field index in src & dest bundles
   integer(IN)      :: ni_src,nj_src   ! dimensions of source bundle
   integer(IN)      :: ni_dst,nj_dst   ! dimensions of dest   bundle
   character(CS)    :: fldName         ! name of a field in a bundle
   character(CS)    :: str             ! generic char string

   !--- formats ---
   character(*),parameter :: subName =   "(dshr_bundle_copyFields) "
   character(*),parameter :: F00     = "('(dshr_bundle_copyFields) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0
   if (debug>0) write(6,F00) "copy all matching fields"

   call dshr_bundle_getDims(bunSrc,ni_src,nj_src,kFlds,rc)
   call dshr_bundle_getDims(bunDst,ni_dst,nj_dst,kFlds,rc)
   if (ni_src /= ni_dst .or. nj_src /= nj_dst ) then
      str = "ERROR: source and dest bundles are have different dimensions"
      write(6,F00) trim(str)
      call shr_sys_abort(subName//trim(str))
   end if

   call dshr_bundle_getNumFields(bunDst,kFlds)
   do k_dst=1,kFlds
      call dshr_bundle_getFieldName (bunDst,k_dst,fldName)
      call dshr_bundle_getFieldIndex(bunSrc,fldName,k_src)
      if (k_src > 0) then
         call dshr_bundle_assignPtr(bunDst,fldName,dataDst)
         call dshr_bundle_assignPtr(bunSrc,fldName,dataSrc)
!dir$ concurrent
!cdir nodep
         dataDst = dataSrc ! copy data from one bundle to another
         if (debug>1) write(6,F00) "copy: bunSrc -> bunDst, field = ",trim(fldName)
      else
         if (debug>1) write(6,F00) "copy: bunSrc -> bunDst, field = ",trim(fldName), &
         & " FYI: field not found, not copied"
      end if
   end do

end subroutine dshr_bundle_copyFields

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_print -- Print info about bundle
!
! !DESCRIPTION:
!     Print info about bundle
!     \newline
!     call dshr\_bundle\_print(bun,rc)
!
! !REVISION HISTORY:
!     2005-Mar-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_print(bun,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bun       ! input bundle
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !--- local ---
   character(CL)   :: str     ! generic string
   integer(IN)     :: rCode   ! return code

   !--- formats ---
   character(*),parameter :: subName =   "(dshr_bundle_print)"
   character(*),parameter :: F00     = "('(dshr_bundle_print) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_print) ',a,3i8)"
   character(*),parameter :: F02     = "('(dshr_bundle_print) ',a,2i8)"
   character(*),parameter :: F03     = "('(dshr_bundle_print) ',a,g20.13)"
   character(*),parameter :: F04     = "('(dshr_bundle_print) ',a,2g20.13)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rCode = 0

   if (dshr_bundle_checkInit(bun)) then
      write(6,F00) ' '
      write(6,F00) '        name  : ',   trim(bun%name)
      write(6,F00) '        init  : ',   trim(bun%initialized)
      write(6,F00) '      dateset : ',   trim(bun%dateset)
      write(6,F01) '        size  : ',        bun%ni, bun%nj, bun%nf
      write(6,F02) '        date  : ',        bun%cDate, bun%sec
      write(6,F00) '        field : ',   trim(bun%fieldList)
      write(6,F03) '  data(1,1,1) : ',        bun%data(1,1,1)
      write(6,F04) ' data min/max : ', minval(bun%data), maxval(bun%data)
      call dshr_domain_getName(bun%domain,str)
      write(6,F00) ' domain name  : ',   trim(str)
   else
      write(6,*) 'bundle not initialized, no data to print'
   end if

   call shr_sys_flush(6)

   if (present(rc)) rc = rCode

end subroutine dshr_bundle_print

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_writeCDF -- Write netCDF file for bundle
!
! !DESCRIPTION:
!     Write bundle to netCDF file
!     \newline
!     call dshr\_bundle\_writeCDF(bun,rc)
!
! !REVISION HISTORY:
!     2005-Mar-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_writeCDF(bundle,rc)

   use netcdf
   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! input bundle
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !--- local ---
   logical,save            :: first_call = .true. ! check if first time in routine
   integer(IN),save        :: fnum                ! file counter number
   integer(IN)             :: lrc                 ! local rc, -> rc
   integer(IN)             :: irc                 ! internal rc
   integer(IN)             :: n
   integer(IN)             :: fid,dimid(2)
   integer(IN),allocatable :: varid(:)
   character(CL)           :: fname
   character(CL)           :: vname
   real(R8),pointer        :: darray(:,:)
   integer(IN)             :: rcode

   !--- formats ---
   character(*),parameter :: subName =   "(dshr_bundle_writeCDF)"
   character(*),parameter :: F00     = "('(dshr_bundle_writeCDF) ',a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lrc = 0

   if (.not.first_call) then
      first_call = .false.
      fnum = 0
   endif
 
   fnum = fnum + 1
   if (fnum > 999) call dshr_bundle_abort(subName//' ERROR fnum gt 999')

   write(fname,'(a2,i3.3,a3)') 'bf',fnum,'.nc'

   rcode = nf90_create(trim(fname),NF90_CLOBBER,fid)
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))
   rcode = nf90_def_dim(fid,'x',bundle%ni,dimid(1))
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))
   rcode = nf90_def_dim(fid,'y',bundle%nj,dimid(2))
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))

   if (bundle%nf > 999) call dshr_bundle_abort(subName//' ERROR nf gt 999')
   allocate(varid(bundle%nf))
   do n = 1,bundle%nf
      write(vname,'(a1,i3.3)') 'f',n
      rcode = nf90_def_var(fid,vname,NF90_DOUBLE,dimid,varid(n))
      if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))
   enddo

   rcode = nf90_put_att(fid,NF90_GLOBAL,'bundle_name',trim(bundle%name))
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))
   rcode = nf90_put_att(fid,NF90_GLOBAL,'field_names',trim(bundle%fieldList))
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))
   rcode = nf90_put_att(fid,NF90_GLOBAL,'bundle_date',bundle%cDate)
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))
   rcode = nf90_put_att(fid,NF90_GLOBAL,'bundle_sec',bundle%sec)
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))
   rcode = nf90_enddef(fid)
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))

   do n = 1,bundle%nf
     call dshr_bundle_assignPtr(bundle,n,darray)
     rcode = nf90_put_var(fid,varid(n),darray)
     if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))
   enddo

   rcode = nf90_close(fid)
   if (rcode /= NF90_NOERR) write(6,F00) trim(nf90_strerror(rcode))

   deallocate(varid)

   write(6,F00) ' wrote '//trim(fname)

   if (present(rc)) rc = lrc

end subroutine dshr_bundle_writeCDF

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_bcastData -- broadcast bundle data from root to all
!
! !DESCRIPTION:
!     Broadcast bundle data from root process to all other processes.
!     \newline
!     Assumes all bundles are initialized on all processes.
!     \newline
!     This is a collective opertation.
!     \newline
!     call dshr\_bundle\_bcastData(bundle,comm_group,debug_string)
!
! !REVISION HISTORY:
!     2006-Oct-24 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_bcastData(bun,comm,str,rc)  

! !INPUT/OUTPUT PARAMETERS:

   !----- arguments -----
   type(dshr_bundle_bundleType),intent(inout) :: bun   ! bundle to b-cast
   integer(IN)                 ,intent(in)    :: comm  ! MPI comm group
   character(*),optional       ,intent(in)    :: str   ! debug string (routine name?)
   integer(IN) ,optional       ,intent(out)   :: rc    ! return code

!EOP

   !----- local -----
   integer(IN) :: commSize    ! comm size
   integer(IN) :: commPID     ! comm rank
   integer(IN) :: nijf(3)     ! size of my    bundle (ni,nj,nf)
   integer(IN) :: ni,nj,nf    ! size of bcast bundle
   integer(IN) :: dateSec(2)  ! date & sec packed into one vector
   integer(IN) :: lrc         ! local return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_bundle_bcastData)"
   character(*),parameter :: F00     = "('(dshr_bundle_bcastData) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_bcastData) ',a,3i9)"

!-------------------------------------------------------------------------------
! Note: 
! - broadcast of bundle date & data, assumes field list is the same
! - bundle must initialized and have the same dimensions
!-------------------------------------------------------------------------------

   if (debug > 0) write(6,F00) 'enter'

   lrc = 0

   !--- verify bundle is initialized (allocated) ------------
   if (.not. dshr_bundle_checkInit(bun)) then
      write(6,F00) 'ERROR: bundle not initialized '
      if (present(str)) write(6,F00) 'string arg: ',trim(str)
      call shr_mpi_commSize(comm,commSize,subName)
      call shr_mpi_commRank(comm,commPID ,subName)
      write(6,F01) ' comm size & rank = ',commSize,commPID 
      call dshr_bundle_abort(trim(subName)//'ERROR: bundle not initialized')
   end if

   !--- verify bundle is appropriately sized ----------------
   call dshr_bundle_getDims(bun,nijf(1),nijf(2),nijf(3))
   ni = nijf(1)
   nj = nijf(2)
   nf = nijf(3)
   call shr_mpi_bcast(nijf,comm,subName)
   if (nijf(1) /= ni .or. nijf(2) /= nj .or. nijf(3) /= nf) then
      write(6,F00) 'ERROR bcast dim mismatch '
      write(6,F01) 'my    ni,nj,nf dims: ',ni     ,nj     ,nf
      write(6,F01) 'bcast ni,nj,nf dims: ',nijf(1),nijf(2),nijf(3) 
      if (present(str)) write(6,F00) 'string arg: ',trim(str)
      call shr_mpi_commSize(comm,commSize,subName)
      call shr_mpi_commRank(comm,commPID ,subName)
      write(6,F01) 'comm size & rank = ',commSize,commPID 
      call dshr_bundle_abort(trim(subName)//' ERROR: bcast dim mismatch')
   end if

   !--- bcast date & sec ------------------------------------
   call shr_mpi_commRank(comm,commPID ,subName)
   if (commPID == 0  .and. .not. dshr_bundle_checkDate(bun)) then
      write(6,F00) 'ERROR: date not set on master PID'
      call dshr_bundle_abort(trim(subName)//' ERROR: date not set on master PID')
   end if
   if (commPID == 0) call dshr_bundle_getDate(bun,dateSec(1),dateSec(2))
   call shr_mpi_bcast(dateSec  ,comm,subName)
   if (commPID /= 0) call dshr_bundle_putDate(bun,dateSec(1),dateSec(2))
   if (debug > 0) write(6,F01) 'date & sec =',bun%cDate,bun%sec

   !--- bcast real data -------------------------------------
   call shr_mpi_bcast(bun%data,comm,subName)

   if (present(rc)) rc = lrc

   if (debug > 0) write(6,F00) 'exit'

end subroutine dshr_bundle_bcastData

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_gather -- root pid gathers tiles from all processors
!
! !DESCRIPTION:
!     Root pid gathers tiles from all processors
!     \newline
!     This is a collective opertation.
!     \newline
!     call dshr\_bundle\_gather(bun|_global,bun\_local,comm_group,debug_string)
!
! !REVISION HISTORY:
!     2006-Oct-30 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_gather(bunG,bunL,comm,str,rc)  

! !INPUT/OUTPUT PARAMETERS:

   !----- arguments -----
   type(dshr_bundle_bundleType),intent(inout) :: bunG  ! bundle, global
   type(dshr_bundle_bundleType),intent(in)    :: bunL  ! bundle, local
   integer(IN)                 ,intent(in)    :: comm  ! MPI comm group
   character(*),optional       ,intent(in)    :: str   ! debug string (routine name?)
   integer(IN) ,optional       ,intent(out)   :: rc    ! return code

!EOP

   !----- local -----
   integer(IN)      :: commSize      ! comm size
   integer(IN)      :: commPID       ! comm rank
   integer(IN)      :: n             ! loop index (thru all PID's in this comm)
   integer(IN)      :: tag1,tag2     ! mpi msg tags
   integer(IN)      :: pid           ! loop index, loops over each pid
   integer(IN)      :: ni  ,nj  ,nf  ! size of local  tile
   integer(IN)      :: ngi ,ngj ,ngf ! size of global tile

   type(dshr_domain_domainType),pointer :: domL ! pointer to local  domain

   integer(IN)      :: i0  ,j0       ! start of local tile wrt global domain
   integer(IN)      :: i1  ,j1       ! end   of local tile wrt global domain
   integer(IN)      :: ndi ,ndj      ! number of decomp tiles in i & j
   integer(IN)      :: tileSpec(7)   ! bun/dom info packed into one msg vector
   real(R8),pointer :: dataL(:,:,:)  ! local  bundle data array
   real(R8),pointer :: dataG(:,:,:)  ! global bundle data array
   real(R8),allocatable :: data(:,:,:)  ! local bundle data array for pid=0

   integer(IN)      :: rCode         ! local return code

   !----- formats -----
   character(*),parameter :: subName =   '(dshr_bundle_gather) '
   character(*),parameter :: F00     = "('(dshr_bundle_gather) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_gather) ',a,8i7)"
   character(*),parameter :: F02     = "('(dshr_bundle_gather) ',a,2es13.3)"

!-------------------------------------------------------------------------------
! Notes: 
! - gather global bundle data on root PID from all local PID's
! - local bundles must be initialized on all PID's
! - global bundle must be initialized on root PID
!-------------------------------------------------------------------------------

   if (debug > 0) write(6,F00) 'enter'
   call shr_sys_flush(6)

   rCode = 0

   call shr_mpi_commSize(comm,commSize,subName)
   call shr_mpi_commRank(comm,commPID ,subName)

   !----------------------------------------------------------------------------
   if (debug>1) write(6,F00) "verify bundles are initialized (allocated)"
   call shr_sys_flush(6)
   !----------------------------------------------------------------------------
   if (.not. dshr_bundle_checkInit(bunL)) then
      write(6,F00) 'ERROR: local bundle not initialized '
      if (present(str)) write(6,F00) 'string arg: ',trim(str)
      write(6,F01) ' comm size & rank = ',commSize,commPID 
      call dshr_bundle_abort(subName//'ERROR: local bundle not initialized')
   end if

   if (commPID == 0) then ! only pid = 0 needs and initialized global bundle
      if (.not. dshr_bundle_checkInit(bunG)) then
         write(6,F00) 'ERROR: global bundle not initialized '
         if (present(str)) write(6,F00) 'string arg: ',trim(str)
         write(6,F01) ' comm size & rank = ',commSize,commPID 
         call dshr_bundle_abort(subName//'ERROR: global bundle not initialized')
      end if
   end if

   !----------------------------------------------------------------------------
   if (debug>1) write(6,F00) " gather the data (msg passing is involved)"
   call shr_sys_flush(6)
   !----------------------------------------------------------------------------

   !--- get local bundle dims ---
   call dshr_bundle_assignPtr(bunL,dataL,rCode)
   call dshr_bundle_getDims  (bunL,ni,nj,nf)

   !--- get local decomp info ---
   call dshr_bundle_domainPtr    (bunL,domL)    ! points to local domain
   call dshr_domain_getDecompInfo(domL,i0,j0,ngi,ngj,ndi,ndj)
   tileSpec(1) = i0    ! start of local tile
   tileSpec(2) = j0
   tileSpec(3) = ni    ! size of local tile
   tileSpec(4) = nj
   tileSpec(5) = ngi   ! size of global domain
   tileSpec(6) = ngj
   tileSpec(7) = nf    ! number of fields

   !--- sanity check ---
   if (ni>ngi .or. nj>ngj ) then
      write(6,F01) "local  ni,nj =",ni ,nj
      write(6,F01) "global ni,nj =",ngi,ngj
      call dshr_bundle_abort(subName//'ERROR: local & global bundles not consistent?')
   end if

   !--- get global bundle dims ---
   if (commPID == 0) then
      call dshr_bundle_assignPtr(bunG,dataG,rCode) ! points to global bundle data
      call dshr_bundle_getDims  (bunG,ngi,ngj,ngf) ! get global bundle dims
   end if

   !--- set default value on gathered data ---
   if (commPID == 0) dataG = dshr_const_spval                  

   !--- send msg from all pid's to root pid ---
   do pid=0,commSize-1
      call  shr_mpi_barrier(comm,subName)
      if (pid == 0) then ! src PID == dest PID == root PID => local copy
         if (commPID == 0) then
            if (ngi /= tileSpec(5) .or. ngj /= tileSpec(6) .or. ngf /= tileSpec(7) ) then
               write(6,F01) "ERROR: local & global bundles not consistent"
               write(6,F01) "local  task,ni,nj,nf = ", pid,tileSpec(5),tileSpec(6),tileSpec(7)
               write(6,F01) "global task,ni,nj,nf = ",   0,ngi        ,ngj        ,ngf
               call dshr_bundle_abort(subName//'ERROR: local & global bundles not consistent')
            end if
            i1 = i0 + ni - 1
            j1 = j0 + nj - 1
            dataG(i0:i1,j0:j1,:) = dataL
         end if
      else               ! src PID /= dest PID
         tag1 = 10
         tag2 = 11
         if (pid == commPID) then
            call shr_mpi_send(tileSpec,  0,tag1,comm,subName) ! specify tile
            if (ni*nj*nf > 0) then ! non-zero data size
               call shr_mpi_send(dataL   ,  0,tag2,comm,subName) ! send    tile
               if (debug>1) then
                  write(6,F01) "send: tileSpec = ",tileSpec
                  write(6,F02) "send: dataL min/max =",minVal(dataL),maxVal(dataL)
               end if
            end if
         else if (commPID == 0) then
            call shr_mpi_recv(tileSpec,pid,tag1,comm,subName)
            i0 = tileSpec(1)
            j0 = tileSpec(2)
            ni = tileSpec(3)
            nj = tileSpec(4)
            i1 = i0 + ni - 1
            j1 = j0 + nj - 1
            if (ngi /= tileSpec(5) .or. ngj /= tileSpec(6) .or. ngf /= tileSpec(7) ) then
               write(6,F01) "ERROR: local & global bundles not consistent"
               write(6,F01) "local  task,ni,nj,nf = ", pid,tileSpec(5),tileSpec(6),tileSpec(7)
               write(6,F01) "global task,ni,nj,nf = ",   0,ngi        ,ngj        ,ngf
               call dshr_bundle_abort(subName//'ERROR: local & global bundles not consistent')
            end if
            if (ni*nj*nf > 0) then ! non-zero data size
               allocate(data(ni,nj,ngf))
               call shr_mpi_recv(data,pid,tag2,comm,subName)
               dataG(i0:i1,j0:j1,:) = data(:,:,:) ! insert tile
               if (debug>1) then
                  write(6,F01) "recv: tileSpec = ",tileSpec
                  write(6,F02) "recv: data  min/max =",minVal(data ),maxVal(data )
               end if
               deallocate(data)
            end if
         end if
      end if
   end do

   if (present(rc)) rc = rCode

   if (debug > 0) write(6,F00) 'exit'
   call shr_sys_flush(6)

end subroutine dshr_bundle_gather

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_scatter -- root pid scatters tiles to all processors
!
! !DESCRIPTION:
!     Root pid scatters tiles to all processors
!     \newline
!     This is a collective opertation.
!     \newline
!     call dshr\_bundle\_scatter(bun|_global,bun\_local,comm_group,debug_string)
!
! !REVISION HISTORY:
!     2006-Oct-30 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_scatter(bunG,bunL,comm,str,rc)  

! !INPUT/OUTPUT PARAMETERS:

   !----- arguments -----
   type(dshr_bundle_bundleType),intent(in)    :: bunG  ! bundle, global
   type(dshr_bundle_bundleType),intent(inout) :: bunL  ! bundle, local
   integer(IN)                 ,intent(in)    :: comm  ! MPI comm group
   character(*),optional       ,intent(in)    :: str   ! debug string (routine name?)
   integer(IN) ,optional       ,intent(out)   :: rc    ! return code

!EOP

   !----- local -----
   integer(IN)      :: commSize      ! comm size
   integer(IN)      :: commPID       ! comm rank
   integer(IN)      :: n             ! loop index (thru all PID's in this comm)
   integer(IN)      :: tag1,tag2     ! mpi msg tags
   integer(IN)      :: pid           ! loop index, loops over each pid
   integer(IN)      :: ni  ,nj  ,nf  ! size of local  tile
   integer(IN)      :: ngi ,ngj ,ngf ! size of global tile

   type(dshr_domain_domainType),pointer :: domL ! pointer to local  domain

   integer(IN)          :: i0  ,j0       ! start of local tile wrt global domain
   integer(IN)          :: i1  ,j1       ! end   of local tile wrt global domain
   integer(IN)          :: ndi ,ndj      ! number of decomp tiles in i & j
   integer(IN)          :: tileSpec(7)   ! bun/dom info packed into one msg vector
   real(R8),pointer     :: dataL(:,:,:)  ! local  bundle data array
   real(R8),pointer     :: dataG(:,:,:)  ! global bundle data array
   real(R8),allocatable :: data (:,:,:)  ! temp local bundle data for pid=0

   integer(IN)      :: rCode         ! local return code

   !----- formats -----
   character(*),parameter :: subName =   '(dshr_bundle_scatter) '
   character(*),parameter :: F00     = "('(dshr_bundle_scatter) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_scatter) ',a,8i7)"
   character(*),parameter :: F02     = "('(dshr_bundle_scatter) ',a,2es13.3)"

!-------------------------------------------------------------------------------
! Notes: 
! - scatter global bundle data on root PID to all local PID's
! - local bundles must be initialized on all PID's
! - global bundle must be initialized on root PID
!-------------------------------------------------------------------------------

   if (debug > 0) write(6,F00) 'enter'
   call shr_sys_flush(6)

   rCode = 0

   call shr_mpi_commSize(comm,commSize,subName)
   call shr_mpi_commRank(comm,commPID ,subName)

   !----------------------------------------------------------------------------
   if (debug>1) write(6,F00) "verify bundles are initialized (allocated)"
   call shr_sys_flush(6)
   !----------------------------------------------------------------------------
   if (.not. dshr_bundle_checkInit(bunL)) then
      write(6,F00) 'ERROR: local bundle not initialized '
      if (present(str)) write(6,F00) 'string arg: ',trim(str)
      write(6,F01) ' comm size & rank = ',commSize,commPID 
      call dshr_bundle_abort(subName//'ERROR: local bundle not initialized')
   end if

   if (commPID == 0) then ! only pid = 0 needs and initialized global bundle
      if (.not. dshr_bundle_checkInit(bunG)) then
         write(6,F00) 'ERROR: global bundle not initialized '
         if (present(str)) write(6,F00) 'string arg: ',trim(str)
         write(6,F01) ' comm size & rank = ',commSize,commPID 
         call dshr_bundle_abort(subName//'ERROR: global bundle not initialized')
      end if
   end if

   !----------------------------------------------------------------------------
   if (debug>1) write(6,F00) "scatter the data (msg passing is involved)"
   call shr_sys_flush(6)
   !----------------------------------------------------------------------------

   !--- get local bundle dims ---
   call dshr_bundle_assignPtr(bunL,dataL,rCode)
   call dshr_bundle_getDims  (bunL,ni,nj,nf)

   !--- get local decomp info ---
   call dshr_bundle_domainPtr    (bunL,domL)    ! points to local domain
   call dshr_domain_getDecompInfo(domL,i0,j0,ngi,ngj,ndi,ndj)
   tileSpec(1) = i0    ! start of local tile
   tileSpec(2) = j0
   tileSpec(3) = ni    ! size of local tile
   tileSpec(4) = nj
   tileSpec(5) = ngi   ! size of global domain
   tileSpec(6) = ngj
   tileSpec(7) = nf    ! number of fields

   !--- sanity check ---
   if (ni>ngi .or. nj>ngj ) then
      write(6,F01) "local  ni,nj =",ni ,nj
      write(6,F01) "global ni,nj =",ngi,ngj
      call dshr_bundle_abort(subName//'ERROR: local & global bundles not consistent?')
   end if

   !--- get global bundle dims ---
   if (commPID == 0) then
      call dshr_bundle_assignPtr(bunG,dataG,rCode) ! points to global bundle data
      call dshr_bundle_getDims  (bunG,ngi,ngj,ngf) ! get global bundle dims
   end if

   !--- set default value on scattered data ---
   dataL = dshr_const_spval                  

   !--- each task requests a tile from root, root extracts & sends the data ---
   do pid=0,commSize-1
      call  shr_mpi_barrier(comm,subName)
      if (pid == 0) then ! src PID == dest PID == root PID => local copy
         if (commPID == 0) then
            if (ngi /= tileSpec(5) .or. ngj /= tileSpec(6) .or. ngf /= tileSpec(7) ) then
               write(6,F01) "ERROR: local & global bundles not consistent"
               write(6,F01) "local  task,ni,nj,nf = ", pid,tileSpec(5),tileSpec(6),tileSpec(7)
               write(6,F01) "global task,ni,nj,nf = ",   0,ngi        ,ngj        ,ngf
               call dshr_bundle_abort(subName//'ERROR: local & global bundles not consistent')
            end if
            i1 = i0 + ni - 1
            j1 = j0 + nj - 1
            dataL = dataG(i0:i1,j0:j1,:)
         end if
      else               ! src PID /= dest PID
         tag1 = 10
         tag2 = 11
         if (pid == commPID) then
            call shr_mpi_send(tileSpec,  0,tag1,comm,subName) ! specify tile
            if (ni*nj*nf > 0) then ! non-zero data size
               call shr_mpi_recv(dataL   ,  0,tag2,comm,subName) ! recv    tile
               if (debug>1) then
                  write(6,F01) "send: tileSpec = ",tileSpec
                  write(6,F02) "recv: dataL min/max =",minVal(dataL),maxVal(dataL)
               end if
            end if
         else if (commPID == 0) then 
            call shr_mpi_recv(tileSpec,pid,tag1,comm,subName)
            i0 = tileSpec(1)
            j0 = tileSpec(2)
            ni = tileSpec(3)
            nj = tileSpec(4)
            i1 = i0 + ni - 1
            j1 = j0 + nj - 1
            if (ngi /= tileSpec(5) .or. ngj /= tileSpec(6) .or. ngf /= tileSpec(7) ) then
               write(6,F01) "ERROR: local & global bundles not consistent"
               write(6,F01) "local  task,ni,nj,nf = ", pid,tileSpec(5),tileSpec(6),tileSpec(7)
               write(6,F01) "global task,ni,nj,nf = ",   0,ngi        ,ngj        ,ngf
               call dshr_bundle_abort(subName//'ERROR: local & global bundles not consistent')
            end if
            if (ni*nj*nf > 0) then ! non-zero data size
               allocate(data(ni,nj,ngf))
               data(:,:,:) = dataG(i0:i1,j0:j1,:) ! extract tile
               call shr_mpi_send(data,pid,tag2,comm,subName)
               if (debug>1) then
                  write(6,F01) "recv: tileSpec = ",tileSpec
                  write(6,F02) "send: data  min/max =",minVal(data ),maxVal(data )
               end if
               deallocate(data)
            end if
         end if
      end if
   end do

   if (present(rc)) rc = rCode

   if (debug > 0) write(6,F00) 'exit'
   call shr_sys_flush(6)

end subroutine dshr_bundle_scatter

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_extractLocal -- extract local tile from global bundle
!
! !DESCRIPTION:
!     Extract local tile from global bundle.
!     \newline
!     call dshr\_bundle\_extractLocal(bun|_global,bun\_local)
!
! !REVISION HISTORY:
!     2006-Nov-06 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_extractLocal(bunG,bunL,rc)  

! !INPUT/OUTPUT PARAMETERS:

   !----- arguments -----
   type(dshr_bundle_bundleType),intent(in)    :: bunG  ! bundle, global
   type(dshr_bundle_bundleType),intent(inout) :: bunL  ! bundle, local
   integer(IN) ,optional       ,intent(out)   :: rc    ! return code

!EOP

   !----- local -----
   integer(IN)      :: ni  ,nj  ,nf   ! size of local  tile
   integer(IN)      :: ngi ,ngj ,ngf  ! size of bunG (according to bunL)
   integer(IN)      :: ngi2,ngj2,ngf2 ! size of bunG (according to bunG)

   integer(IN)      :: i0  ,j0        ! start of local tile wrt global domain
   integer(IN)      :: i1  ,j1        ! end   of local tile wrt global domain
   integer(IN)      :: ndi ,ndj       ! number of decomp tiles in i & j
   real(R8),pointer :: dataL(:,:,:)   ! local  bundle data array
   real(R8),pointer :: dataG(:,:,:)   ! global bundle data array

   type(dshr_domain_domainType),pointer :: domL ! pointer to local domain

   integer(IN)      :: rCode         ! local return code

   !----- formats -----
   character(*),parameter :: subName =   '(dshr_bundle_extractLocal) '
   character(*),parameter :: F00     = "('(dshr_bundle_extractLocal) ',4a)"
   character(*),parameter :: F01     = "('(dshr_bundle_extractLocal) ',a,8i7)"
   character(*),parameter :: F02     = "('(dshr_bundle_extractLocal) ',a,2es13.3)"

!-------------------------------------------------------------------------------
! Notes: 
!-------------------------------------------------------------------------------

   debug = 1
   if (debug > 0) write(6,F00) 'enter'
   call shr_sys_flush(6)

   rCode = 0

   !----------------------------------------------------------------------------
   if (debug>1) write(6,F00) "verify bundles are initialized (allocated)"
   call shr_sys_flush(6)
   !----------------------------------------------------------------------------
   if (.not. dshr_bundle_checkInit(bunL)) then
      write(6,F00) 'ERROR: local bundle not initialized '
      call dshr_bundle_abort(subName//'ERROR: local bundle not initialized')
   end if
   if (.not. dshr_bundle_checkInit(bunG)) then
      write(6,F00) 'ERROR: global bundle not initialized '
      call dshr_bundle_abort(subName//'ERROR: global bundle not initialized')
   end if

   !----------------------------------------------------------------------------
   if (debug>1) write(6,F00) "determine local & global domain info"
   call shr_sys_flush(6)
   !----------------------------------------------------------------------------

   !--- get local bundle dims & decomp info---
   call dshr_bundle_assignPtr    (bunL,dataL,rCode) ! points to local data
   call dshr_bundle_domainPtr    (bunL,domL)        ! points to local domain
   call dshr_bundle_getDims      (bunL,ni,nj,nf)    ! get local bundle dims
   call dshr_domain_getDecompInfo(domL,i0,j0,ngi,ngj,ndi,ndj)

   !--- sanity check ---
   if (ni>ngi .or. nj>ngj ) then
      write(6,F01) "local  ni,nj =",ni ,nj
      write(6,F01) "global ni,nj =",ngi,ngj
      call dshr_bundle_abort(subName//'ERROR: local bundle not self-consistent?')
   end if

   !--- get global bundle dims ---
   call dshr_bundle_assignPtr(bunG,dataG,rCode)    ! points to global bundle data
   call dshr_bundle_getDims  (bunG,ngi2,ngj2,ngf2) ! get global bundle dims

   !--- sanity check ---
   if (ngi/=ngi2 .or. ngj/=ngj2 .or. ngf2/=nf ) then
      write(6,F01) "bunL global ni,nj,nf =",ngi ,ngj ,nf
      write(6,F01) "bunG global ni,nj,nf =",ngi2,ngj2,ngf2
      call dshr_bundle_abort(subName//'ERROR: local & global bundles not consistent?')
   end if

   !----------------------------------------------------------------------------
   if (debug>1) write(6,F00) "extract local data from global bundle"
   call shr_sys_flush(6)
   !----------------------------------------------------------------------------
   dataL = dshr_const_spval                  
   i1 = i0 + ni - 1
   j1 = j0 + nj - 1
   dataL(:,:,:) = dataG(i0:i1,j0:j1,:)

   bunL%dateset   = bunG%dateset 
   bunL%cDate     = bunG%cDate
   bunL%sec       = bunG%sec
   bunL%fieldList = bunG%fieldList

   if (present(rc)) rc = rCode

   if (debug > 0) write(6,F00) 'exit'
   call shr_sys_flush(6)
   debug = 0

end subroutine dshr_bundle_extractLocal

!===============================================================================
!===============================================================================

subroutine dshr_bundle_abort(string)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),optional,intent(IN) :: string

!EOP

   !--- local ---
   character(CL)          :: lstring

   !--- formats ---
   character(*),parameter :: subName =   "(dshr_bundle_abort)"
   character(*),parameter :: F00     = "('(dshr_bundle_abort) ',a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lstring = ''
   if (present(string)) lstring = string

   if (doabort) then
      call shr_sys_abort(lstring)
   else
      write(6,F00) ' no abort:'//trim(lstring)
   endif

end subroutine dshr_bundle_abort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_setAbort -- Set local dshr_bundle abort flag
!
! !DESCRIPTION:
!     Set local abort flag, true = abort, false = print and continue
!     \newline
!     call dshr\_bundle\_setAbort(.false.)
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_setAbort(flag)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   logical,intent(in) :: flag

!EOP

   !--- local ---

   !--- formats ---
   character(*),parameter :: subName = "('dshr_bundle_setAbort') "
   character(*),parameter :: F00     = "('(dshr_bundle_setAbort) ',a) "

!-------------------------------------------------------------------------------

   doabort = flag

end subroutine dshr_bundle_setAbort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_setDebug -- Set local debug level
!
! !DESCRIPTION:
!     Set local debug level, 0 = production
!     \newline
!     call dshr\_bundle\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_setDebug(level)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer,intent(in) :: level

!EOP

   !--- local ---

   !--- formats ---
   character(*),parameter :: subName =  "('dshr_bundle_setDebug') "
   character(*),parameter :: F00     = "('(dshr_bundle_setDebug) ',a) "
   character(*),parameter :: F01     = "('(dshr_bundle_setDebug) ',a,i4) "

!-------------------------------------------------------------------------------

   debug = level
   write(6,F01) "debug level reset to ",level

end subroutine dshr_bundle_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_bundle_getDebug -- Return local debug level
!
! !DESCRIPTION:
!     Return local debug level, 0 = production
!     \newline
!     call dshr\_bundle\_getDebug(level)
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_bundle_getDebug(level)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer,intent(out) :: level

!EOP

   !--- local ---

   !--- formats ---
   character(*),parameter :: subName =  "('dshr_bundle_getDebug') "
   character(*),parameter :: F00     = "('(dshr_bundle_getDebug) ',a) "

!-------------------------------------------------------------------------------

   level = debug

end subroutine dshr_bundle_getDebug

!===============================================================================
!===============================================================================

end module dshr_bundle
