!===============================================================================
! SVN $Id: dshr_domain.F90 3456 2007-03-09 23:14:30Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_domain.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_domain -- domain data type and associated methods
!
! !DESCRIPTION:
!    Domain data type and associated methods.
!     
! !REMARKS:
!    1) should domain include z coordinates?
!       integer(IN)      :: nk           ! size of k dimension
!       real(R8),pointer :: z    (:,:,:) ! height coordinate
!    2) should domain store data the same way that bundles store data?
!       integer(IN)      :: nf            ! number of fields
!       real             :: data(:,:,:)   ! coord data (ni,nj,nf)
!       character(CL)    :: fieldNames    ! = "lat:lon:area:mask:index"
!     
! !REVISION HISTORY:
!    2006-Nov-xx - B. Kauffman - add support for parallel data models
!    2006-May-22 - B. Kauffman - add cell fraction data
!    2005-Jan-27 - J. Schramm and B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_domain 

! !USES:

   use shr_sys_mod    ! system calls
   use shr_string_mod ! string/list operators
   use shr_ncread_mod ! netCDF reader
   use shr_mpi_mod    ! mpi wrappers

   use dshr_kind      ! kinds for strong typing
   use dshr_const     ! constants

   implicit none
   private            ! everything is default private

! !PUBLIC TYPES:

   public :: dshr_domain_domainType  ! domain data type

   type dshr_domain_domainType ! domain definition, grid + mask + area
      private                         ! no public access to internal components
      character(CS)    :: name        ! domain name / ID string
      character(CL)    :: Fn          ! source filename, if applicable
      character(CS)    :: initialized ! use char string instead of logical
      character(CS)    :: filled      ! use char string instead of logical
      integer(IN)      :: ni          ! size of i dimension
      integer(IN)      :: nj          ! size of j dimension
      
      integer(IN)      ::  i0 , j0    ! 2d decomp: local start index wrt global domain
      integer(IN)      :: ngi ,ngj    ! 2d decomp: size of global domain
      integer(IN)      :: ndi ,ndj    ! 2d decomp: number of tiles wrt i,j dims

      real(R8),pointer :: lat  (:,:)  ! latitude
      real(R8),pointer :: lon  (:,:)  ! longitude
      real(R8),pointer :: area (:,:)  ! cell area
      real(R8),pointer :: frac (:,:)  ! cell fraction
      real(R8),pointer :: mask (:,:)  ! domain mask
      real(R8),pointer :: index(:,:)  ! global cell index (wrt decomp)
   end type dshr_domain_domainType

! !PUBLIC MEMBER FUNCTIONS:

   public :: dshr_domain_decomp2d      ! compute a 2d domain decomposition
   public :: dshr_domain_bcastInfo     ! COLLECTIVE: bcast info needed by domain_create
   public :: dshr_domain_bcastData     ! COLLECTIVE: bcast data into existing domains
   public :: dshr_domain_getDecompInfo ! get decomposition info
   public :: dshr_domain_putDecompInfo ! put decomposition info
   public :: dshr_domain_extractLocal  ! extract local domain from global domain

   public :: dshr_domain_create     ! allocate memory data arrays
   public :: dshr_domain_clean      ! clean domain, delete memory
   public :: dshr_domain_readData   ! read domain data from a file
   public :: dshr_domain_fill       ! fill domain with bogus test data
   public :: dshr_domain_compare    ! returns true if 2 domains are the same
   public :: dshr_domain_assignPtr  ! pointer into domain data

   public :: dshr_domain_getDims    ! get domain dimensions from data type
   public :: dshr_domain_getData    ! copy data out of domain
   public :: dshr_domain_putData    ! copy data into   domain
   public :: dshr_domain_getName    ! get domain name from data type
   public :: dshr_domain_putName    ! put domain name into data type
   public :: dshr_domain_getFn      ! get domain name from data type
   public :: dshr_domain_putFn      ! put domain name into data type

   public :: dshr_domain_print      ! print domain information
   public :: dshr_domain_setFilled  ! indicates domain is filled with data
   public :: dshr_domain_setAbort   ! set local abort flag 
   public :: dshr_domain_setDebug   ! set internal shr_stream debug level
   public :: dshr_domain_getDebug   ! get internal shr_stream debug level

! !PUBLIC DATA MEMBERS:

   integer(IN),parameter,public :: dshr_domain_compareXYabs      = 1 ! X,Y  relative error
   integer(IN),parameter,public :: dshr_domain_compareXYrel      = 2 ! X,Y  absolute error
   integer(IN),parameter,public :: dshr_domain_compareAreaAbs    = 3 ! area relative error
   integer(IN),parameter,public :: dshr_domain_compareAreaRel    = 4 ! area absolute error
   integer(IN),parameter,public :: dshr_domain_compareMaskIdent  = 5 ! masks are identical
   integer(IN),parameter,public :: dshr_domain_compareMaskZeros  = 6 ! masks have same zeros
   integer(IN),parameter,public :: dshr_domain_compareMaskSubset = 7 ! mask is subset of other

!EOP

   character(CS),parameter :: dshr_domain_setTru = 'TRUE dom'
   character(CS),parameter :: dshr_domain_setFal = 'FALSE d '
   character(CS),parameter :: dshr_domain_Fnnone = 'FNnone  '
   real(R8)     ,parameter :: pi = dshr_const_pi ! locally available pi
   logical      ,save      :: doabort = .true.   ! local abort flag
   integer(IN)  ,save      :: debug = 0          ! local debug level

   !--- manages locally generated domain names, matches them with filenames ---
   integer(IN)  ,parameter :: dshr_domain_nmax = 100
   integer(IN)  ,save      :: dshr_domain_nloc = 0
   character(CL),save      :: dshr_domain_dname(dshr_domain_nmax) ! locally generated domain names
   character(CL),save      :: dshr_domain_fname(dshr_domain_nmax) ! associated filenames
   character(CS),parameter :: dshr_domain_dnameh = 'dom'  ! domain name header

!===============================================================================
contains
!===============================================================================

subroutine dshr_domain_decomp2d(comm,pid, ngi,ngj, ndi,ndj, i0,j0, ni,nj)
 

   !----- arguments -----
   integer(IN),intent(in)  :: comm    ! communicator group in question
   integer(IN),intent(in)  :: pid     ! pid within the comm group
   integer(IN),intent(in)  :: ngi,ngj ! size of global domain
   integer(IN),intent(out) :: ndi,ndj ! number of tiles in i & j dimensions
   integer(IN),intent(out) :: i0,j0   ! start indicies of tile wrt global domain
   integer(IN),intent(out) :: ni,nj   ! size of local tile for this PID

   !----- local -----
   integer(IN) :: commSize  ! MPI communicator group size
   integer(IN) :: i1,j1     ! ending indicies of tile wrt global domain
   integer(IN) :: id,jd     ! indicies of 2d decomp tile for this PID
   integer(IN) :: tSize     ! tile size
   integer(IN) :: i,j,n     ! generic indicies
   logical(IN) :: found     ! T <=> was able to compute a 2d decomposition

   !----- formats -----
   character(*),parameter :: subName = '(dshr_domain_decomp2d) '
   character(*),parameter :: F00   = "('(dshr_domain_decomp2d) ',4a)"
   character(*),parameter :: F01   = "('(dshr_domain_decomp2d) ',a,2i6)"
   character(*),parameter :: F02   = "('(dshr_domain_decomp2d) ',a,i6,4(3x,2i5))"

!-------------------------------------------------------------------------------
! ASSUMPTIONS:
! - integers ndi & ndj exist st ndi*ndj=commSize, ngi >= dni  and  ngj >= dnj
! - pid is in {0,1,2,...,commSize-1}
!-------------------------------------------------------------------------------

   call shr_mpi_commSize(comm,commSize,subName)

   !----------------------------------------------------------------------------
   ! find ndi & ndj st ndi*ndj=commSize (determine number of tiles in i & j)
   !----------------------------------------------------------------------------
   found = .false.
   do ndi=1,ngi
      if (mod(commSize,ndi) == 0) then
         ndj = commSize/ndi
         if (ngj >= ndj) found = .true.
      end if
      if (found) exit
   end do
   if (.not. found) call dshr_domain_abort(trim(subName)//'no valid 2d decomp')

   !----------------------------------------------------------------------------
   ! confirm 2d decomposition is valid:
   ! * commSize = ndi*ndj
   ! * pid is in 0,1,2,...,(commSize-1)
   ! * ngi >= dni  and  ngj >= dnj
   !----------------------------------------------------------------------------
   if (ndi*ndj /= commSize) then
      write(6,F01) "ERROR:           ndi,ndj = ",          ndi,ndj
      write(6,F01) "ERROR: commSize, ndi*ndj = ",commSize, ndi*ndj
      write(6,F00) "ERROR: commSize /= ndi*ndj"
      call dshr_domain_abort(subName//"commSize /= ndi*ndj")
   else if (pid < 0  .or.  commSize <= pid) then
      write(6,F01) "ERROR: bad pid  = ",pid
      write(6,F01) "ERROR: commSize = ",commSize
      call dshr_domain_abort(subName//"bad pid")
   else if (ndi > ngi) then
      write(6,F01) "ERROR: ndi > ngi,  ndi,ngi = ",ndi,ngi
      call dshr_domain_abort(subName//"ndi > ngi")
   else if (ndj > ngj) then
      write(6,F01) "ERROR: ndj > ngj,  ndj,ngj = ",ndj,ngj
      call dshr_domain_abort(subName//"ndj > ngj")
   end if

   !----------------------------------------------------------------------------
   ! determine which 2d tile belongs to this PID and its size & location
   !
   ! get (id,jd) ~ tile's indecies wrt 2d decomposition
   ! note: pid is in [0,commSize-1]
   ! note: id is in [1,ndi], jd is in [1,ndj]
   !----------------------------------------------------------------------------

   !--- get i & j index of tile in 2d decomposition ---
   id = mod(pid,ndi) + 1
   jd = (pid + 1 - id)/ndi + 1

   if ( ((jd-1)*ndi + id - 1) /= PID ) then
      write(6,F01) "ERROR: computing tile index, PID = ",PID
      call dshr_domain_abort(subName//"ERROR: computing tile index")
   end if

   !--- for tile (id,jd), get 1st & last global i index ---
   tSize = ngi/ndi        ! nominal tile size...
   n     = mod(ngi,ndi)   ! ... except 1st n tiles are one cell larger
   if (  id-1 < n ) then  ! all preceeding tiles are size (tSize+1)
      i0 = (id-1)*(tSize+1) + 1
      i1 = i0 + (tSize+1) - 1
   else                   ! n preceeding tiles are size (tSize+1), rest are tSize
      i0 = (id-1)*(tSize) + n + 1
      i1 = i0 +   (tSize) - 1
   end if
   ni = i1 - i0 + 1

   !--- for tile (id,jd), get 1st & last global j index ---
   tSize = ngj/ndj        ! nominal tile size...
   n     = mod(ngj,ndj)   ! ... except 1st n tiles are one cell larger
   if (  jd-1 < n ) then  ! all preceeding tiles are size (tSize+1)
      j0 = (jd-1)*(tSize+1) + 1
      j1 = j0 +   (tSize+1) - 1
   else                   ! n preceeding tiles are size (tSize+1), rest are tSize
      j0 = (jd-1)*(tSize) + n + 1
      j1 = j0 +   (tSize) - 1
   end if
   nj = j1 - j0 + 1

   !--------------------------------------------------------------
   ! document what happened
   !--------------------------------------------------------------
   write(6,F02) "commSize, ngi,ngj, ndi,ndj = ",commSize,ngi,ngj,ndi,ndj
   write(6,F02) "PID, id,jd, ni,nj, i0,i1, j0,j1 =",PID ,id,jd, ni,nj, i0,i1, j0,j1

end subroutine dshr_domain_decomp2d

!===============================================================================
!===============================================================================

subroutine dshr_domain_bcastInfo(domain,name,ni,nj,comm,str,rc)  

   !----- arguments -----
   type(dshr_domain_domainType),intent(in)    :: domain ! domain with info
   character(*)                ,intent(out)   :: name   ! domain's name
   integer(IN)                 ,intent(out)   :: ni,nj  ! domain's dimension
   integer(IN)                 ,intent(in)    :: comm   ! MPI comm group
   character(*),optional       ,intent(in)    :: str    ! debug string (routine name?)
   integer(IN) ,optional       ,intent(out)   :: rc     ! return code

   !----- local -----
   integer(IN) :: PID         ! mpi process ID
   integer(IN) :: rCode       ! local return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_domain_bcastInfo) "
   character(*),parameter :: F00     = "('(dshr_domain_bcastInfo) ',4a)"
   character(*),parameter :: F01     = "('(dshr_domain_bcastInfo) ',a,3i9)"

!-------------------------------------------------------------------------------
! Note: 
!-------------------------------------------------------------------------------

!  if (debug > 0) write(6,F00) 'enter'
   if (.true.   ) write(6,F00) 'enter'
   call shr_sys_flush(6)

   rCode = 0

   !--- get info from source domain ---
   call shr_mpi_commRank( comm,PID ,subName)
   if (PID == 0) then
      call dshr_domain_getDims(domain,ni,nj,rCode)
      call dshr_domain_getName(domain,name ,rCode)
   end if

   !--- bcast info from source domain ---
   call     shr_mpi_bcast(   ni,comm,subName)
   call     shr_mpi_bcast(   nj,comm,subName)
   call     shr_mpi_bcast(name ,comm,subName)

   if (present(rc)) rc = rCode

!  if (debug > 0) write(6,F00) 'exit'
   if (.true.   ) write(6,F00) 'exit'

end subroutine dshr_domain_bcastInfo

!===============================================================================
!===============================================================================

subroutine dshr_domain_bcastData(domain,comm,str,rc)  

   !----- arguments -----
   type(dshr_domain_domainType),intent(inout) :: domain ! domain to b-case
   integer(IN)                 ,intent(in)    :: comm   ! MPI comm group
   character(*),optional       ,intent(in)    :: str    ! debug string (routine name?)
   integer(IN) ,optional       ,intent(out)   :: rc     ! return code

   !----- local -----
   integer(IN) :: commSize    ! comm size
   integer(IN) :: commPID     ! comm rank
   integer(IN) :: nij(2)      ! size of my    domain (ni,nj)
   integer(IN) :: ni,nj       ! size of bcast domain
   integer(IN) :: rCode       ! local return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_domain_bcastData) "
   character(*),parameter :: F00     = "('(dshr_domain_bcastData) ',4a)"
   character(*),parameter :: F01     = "('(dshr_domain_bcastData) ',a,3i9)"

!-------------------------------------------------------------------------------
! Note: 
! - broadcast of domain data
! - domains must initialized and have the same dimensions
!-------------------------------------------------------------------------------

!  if (debug > 0) write(6,F00) 'enter'
   if (.true.   ) write(6,F00) 'enter'

   rCode = 0

   !--- verify domain is initialized (allocated) ---
   if (.not. dshr_domain_checkInit(domain)) then
      write(6,F00) 'ERROR: domain not initialized '
      if (present(str)) write(6,F00) 'string arg: ',trim(str)
      call shr_mpi_commSize(comm,commSize,subName)
      call shr_mpi_commRank(comm,commPID ,subName)
      write(6,F01) ' comm size & rank = ',commSize,commPID 
      call dshr_domain_abort(trim(subName)//'domain not initialized')
   end if

   !--- verify master domain is filled ---
   call shr_mpi_commRank(comm,commPID ,subName)
   if (commPID == 0  .and. .not. dshr_domain_checkFilled(domain)) then
      write(6,F00) 'ERROR: domain not filled '
      if (present(str)) write(6,F00) 'string arg: ',trim(str)
      call shr_mpi_commSize(comm,commSize,subName)
      call shr_mpi_commRank(comm,commPID ,subName)
      write(6,F01) ' comm size & rank = ',commSize,commPID 
      call dshr_domain_abort(trim(subName)//'domain not filled')
   end if

   !--- verify domain is appropriately sized ----------------
   call dshr_domain_getDims(domain,nij(1),nij(2))
   ni = nij(1)
   nj = nij(2)
   call shr_mpi_bcast(nij,comm,subName)
   if (nij(1) /= ni .or. nij(2) /= nj ) then
      write(6,F00) 'ERROR bcast dim mismatch '
      write(6,F01) 'my    ni,nj dims: ',ni    ,nj
      write(6,F01) 'bcast ni,nj dims: ',nij(1),nij(2)
      if (present(str)) write(6,F00) 'string arg: ',trim(str)
      call shr_mpi_commSize(comm,commSize,subName)
      call shr_mpi_commRank(comm,commPID ,subName)
      write(6,F01) 'comm size & rank = ',commSize,commPID 
      call dshr_domain_abort(trim(subName)//'bcast dim mismatch')
   end if

   !--- bcast the data --------------------------------------
   call shr_mpi_bcast(domain%name ,comm,subName)
   call shr_mpi_bcast(domain%lat  ,comm,subName)
   call shr_mpi_bcast(domain%lon  ,comm,subName)
   call shr_mpi_bcast(domain%area ,comm,subName)
   call shr_mpi_bcast(domain%frac ,comm,subName)
   call shr_mpi_bcast(domain%mask ,comm,subName)
   call shr_mpi_bcast(domain%index,comm,subName)

   call dshr_domain_setInit  (domain,.true.)
   call dshr_domain_setFilled(domain,.true.)

   if (present(rc)) rc = rCode

!  if (debug > 0) write(6,F00) 'exit'
   if (.true.   ) write(6,F00) 'exit'

end subroutine dshr_domain_bcastData

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_getDecompInfo -- returns domain dimension info
!
! !DESCRIPTION:
!    Returns domain docomposition info
!
! !REVISION HISTORY:
!    2006-Oct-19 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_getDecompInfo(domain,i0,j0,ngi,ngj,ndi,ndj,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)  :: domain   ! domain in question
   integer(IN)                 ,intent(out) :: i0 ,j0   ! start index of tile
   integer(IN)                 ,intent(out) :: ngi,ngj  ! size of global domain
   integer(IN)                 ,intent(out) :: ndi,ndj  ! number of tiles in i,j
   integer(IN),optional        ,intent(out) :: rc       ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName =   '(dshr_domain_getDecompInfo) '
   character(*),parameter :: F00     = "('(dshr_domain_getDecompInfo) ',4a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

    i0 = domain%i0
    j0 = domain%j0
   ngi = domain%ngi
   ngj = domain%ngj
   ndi = domain%ndi
   ndj = domain%ndj

end subroutine dshr_domain_getDecompInfo

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_putDecompInfo -- returns domain dimension info
!
! !DESCRIPTION:
!    Returns domain docomposition info
!
! !REVISION HISTORY:
!    2006-Oct-19 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_putDecompInfo(domain,i0,j0,ngi,ngj,ndi,ndj,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(out) :: domain   ! domain in question
   integer(IN)                 ,intent(in)  :: i0 ,j0   ! start index of tile
   integer(IN)                 ,intent(in)  :: ngi,ngj  ! size of global domain
   integer(IN)                 ,intent(in)  :: ndi,ndj  ! number of tiles in i,j
   integer(IN),optional        ,intent(out) :: rc       ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName =   '(dshr_domain_putDecompInfo) '
   character(*),parameter :: F00     = "('(dshr_domain_putDecompInfo) ',4a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   domain%i0  = i0
   domain%j0  = j0
   domain%ngi = ngi
   domain%ngj = ngj
   domain%ndi = ndi
   domain%ndj = ndj

end subroutine dshr_domain_putDecompInfo

!===============================================================================
!===============================================================================

subroutine dshr_domain_extractLocal(domainG, domainL, rc)

   !----- arguments -----
   type(dshr_domain_domainType),intent(in)    :: domainG ! global domain
   type(dshr_domain_domainType),intent(inout) :: domainL ! local  domain
   integer(IN)        ,optional,intent(out)   :: rc      ! return code

   !----- local -----
   integer(IN)      :: ngi,ngj    ! size of global domain
   integer(IN)      :: ni ,nj     ! size of local  domain
   integer(IN)      :: i0,j0      ! starting index of local tile
   integer(IN)      :: i1,j1      ! ending   index of local tile
   character(CL)    :: str        ! generic char string
   real(R8),pointer :: dataG(:,:) ! temp pointer into domain data
   integer(IN)      :: rCode      ! local return code

   !----- formats -----
   character(*),parameter :: subName = "(dshr_domain_extractLocal) "
   character(*),parameter :: F00   = "('(dshr_domain_extractLocal) ',4a     )"
   character(*),parameter :: F01   = "('(dshr_domain_extractLocal) ',a,3(2i6,2x))"
   character(*),parameter :: F03   = "('(dshr_domain_extractLocal) ',4(a,'[',i4,',',i4,']'))"

!-------------------------------------------------------------------------------
! ASSUMPTIONS:
!-------------------------------------------------------------------------------

   !--- make sure domains are created (allocated)
   if      (.not. dshr_domain_checkInit(domainL)) then
       write(6,F01) "ERROR: local domain is not initialized (allocated)"
       call dshr_domain_abort(subName//"ERROR: local domain is not initialized (allocated)")
   else if (.not. dshr_domain_checkInit(domainG)) then
       write(6,F01) "ERROR: global domain is not initialized (allocated)"
       call dshr_domain_abort(subName//"ERROR: global domain is not initialized (allocated)")
   else if (.not. dshr_domain_checkFilled(domainG)) then
       write(6,F01) "ERROR: global domain is not filled (with data)"
       call dshr_domain_abort(subName//"ERROR: global domain is not filled (with data)")
   end if

   !--- gather info about local & global domains ---
   ni = domainL%ni
   nj = domainL%nj
   i0 = domainL%i0  ! starting i-index of local tile wrt global domain
   j0 = domainL%j0  ! starting j-index
   i1 = i0 + ni - 1 ! ending   i-index of local tile wrt global domain
   j1 = j0 + nj - 1 ! ending   j-index
   call dshr_domain_getdims(domainG,ngi,ngj,rCode)

   !--- make sure local domain is a subset of global domain ---
   if (  i0 < 1  .or.  ngi < i1  .or.  j0 < 1   .or.  ngj < j1  ) then 
       write(6,F00) "ERROR: local domain not a subset of global domain"
       write(6,F01) "ERROR: size of global domain:  ngi,ngj =",ngi,ngj
       write(6,F01) "ERROR: size of local  domain:   ni, nj =", ni, nj
       write(6,F01) "ERROR: i-range of local tile    i0, i1 =", i0, i1
       write(6,F01) "ERROR: j-range of local tile    j0, j1 =", j0, j1
       call dshr_domain_abort(subName//"ERROR: local domain not subset of global domain")
   end if

   write(6,F03) "extracting tile ",i0,i1,"x",j0,j1," from global domain ",1,ngi,"x",1,ngj

   !--- extract local tile from global domain ---
   call dshr_domain_getFn(domainG,str)
   call dshr_domain_putFn(domainL,str)
   call dshr_domain_assignPtr(domainG,"lat"  ,dataG,rCode)
   call dshr_domain_putData  (domainL, dataG(i0:i1,j0:j1), "lat"  , rCode)
   call dshr_domain_assignPtr(domainG,"lon"  ,dataG,rCode )
   call dshr_domain_putData  (domainL, dataG(i0:i1,j0:j1), "lon"  , rCode)
   call dshr_domain_assignPtr(domainG,"frac" ,dataG,rCode)
   call dshr_domain_putData  (domainL, dataG(i0:i1,j0:j1), "frac" , rCode)
   call dshr_domain_assignPtr(domainG,"area" ,dataG,rCode)
   call dshr_domain_putData  (domainL, dataG(i0:i1,j0:j1), "area" , rCode)
   call dshr_domain_assignPtr(domainG,"mask" ,dataG,rCode)
   call dshr_domain_putData  (domainL, dataG(i0:i1,j0:j1), "mask" , rCode)
   call dshr_domain_assignPtr(domainG,"index",dataG,rCode)
   call dshr_domain_putData  (domainL, dataG(i0:i1,j0:j1), "index", rCode)
   call dshr_domain_setFilled(domainL,.true.)
   
end subroutine dshr_domain_extractLocal

!===============================================================================
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_create -- Initialize/allocate memory for domain info
!
! !DESCRIPTION:
!    Allocate memory for domain information
!
! !REVISION HISTORY:
!    2005-Jan-27 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_create(domain,name,ni,nj,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(inout) :: domain  ! domain to initialize
   character(*)                ,intent(in)    :: name    ! domain name
   integer(IN)                 ,intent(in)    :: ni,nj   ! domain dimensions
   integer(IN),optional        ,intent(out)   :: rc      ! return code

!EOP

   !----- local -----
   integer(IN)   :: rCode   ! return code

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_create')"
   character(*),parameter :: F00     = "('(dshr_domain_create) ',4a)" 
   character(*),parameter :: F01     = "('(dshr_domain_create) ',a,2i6)" 

!-------------------------------------------------------------------------------
! NOTES
! o after this routine, a domain is "initialized" but not "filled"
!
! o we used to dealloc with a return code prior to allocation, but for some
!   unexplained reason, this caused the IBM power4 (bluesky) do dump core,
!   so we switched to an 
!     if (initialized) then allocate
!   approach
!
! o this is no longer done, see above
!   the return code allows the deallocate to fail without stopping program
!   exectution.  This would happen, eg, if the pointer has never been
!   allocated to begin with (dealloc doesn't like undefied pointers).
!   F90 has no way (?) of testing an undefined pointer (eg. never allocated).
!   See Fortran 90/95 Explained, Melcaf, 1998, p 124, 175
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   write(6,F00) "initializing domain, name = ",trim(name)
   write(6,F01) "allocating 2d memory: ni,nj = ",ni,nj

   if (ni <= 0 .or. nj <= 0) then
      write(6,F01) "ERROR: ni,nj = ",ni,nj
      if (present(rc)) rc = 1
      return
   end if

   if (dshr_domain_checkInit(domain)) call dshr_domain_clean(domain,rCode)

   domain%Fn   = dshr_domain_Fnnone
   domain%name = name
   domain%ni   = ni
   domain%nj   = nj
   domain%i0   = -1
   domain%j0   = -1
   domain%ngi  = -1
   domain%ngj  = -1
   domain%ndi  = -1
   domain%ndj  = -1

   !--- allocate ---
   allocate(domain%lon  (ni,nj))
   allocate(domain%lat  (ni,nj))
   allocate(domain%area (ni,nj))
   allocate(domain%frac (ni,nj))
   allocate(domain%mask (ni,nj))
   allocate(domain%index(ni,nj))

   domain%lon   = dshr_const_spval
   domain%lat   = dshr_const_spval
   domain%area  = dshr_const_spval
   domain%frac  = dshr_const_spval
   domain%mask  = dshr_const_spval
   domain%index = dshr_const_spval

   call dshr_domain_setInit(domain,.true.) ! now initialized, but not filled

end subroutine dshr_domain_create

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_clean -- clean domain
!
! !DESCRIPTION:
!    Clean domain including reset values, deallocate memory
!
! !REVISION HISTORY:
!    2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_clean(domain,rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(inout) :: domain  ! domain 
   integer(IN),optional        ,intent(out)   :: rc      ! return code

!EOP

   !----- local -----
   integer(IN)   :: rCode             ! return code

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_clean') "
   character(*),parameter :: F00     = "('(dshr_domain_clean) ',4a)" 

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   if (debug>0) write(6,F00) 'dealloc all domain memory'

   if (present(rc)) rc = 0

   call dshr_domain_setInit  (domain,.false.)
   call dshr_domain_setFilled(domain,.false.)

   domain%name = ''
   domain%Fn   = ''
   domain%ni   = -1
   domain%nj   = -1

   !--- dealloc if allocated ---
   deallocate(domain%lon  ,STAT=rCode)
   if (debug>0 .and. rCode>0) write(6,F00) 'WARNING unable to dealloc domain%lon'
   deallocate(domain%lat  ,STAT=rCode)
   deallocate(domain%area ,STAT=rCode)
   deallocate(domain%frac ,STAT=rCode)
   deallocate(domain%mask ,STAT=rCode)
   deallocate(domain%index,STAT=rCode)


end subroutine dshr_domain_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_readData -- read in domain data from a file
!
! !DESCRIPTION:
!    Read in domain data from a file
!
! !REVISION HISTORY:
!    2005-Jan-27 - J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_readData(domain    ,fileName  ,lonNameIn,latNameIn,maskNameIn, &
                                areaNameIn,fracNameIn,readFrac  ,rc)

! !USES:

   use shr_file_mod  ! shared file functions

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(inout) :: domain     ! domain to put data into
   character(*)                ,intent(in)    :: fileName   ! string 
   character(*),optional       ,intent(in)    ::  lonNameIn ! name of  lon variable
   character(*),optional       ,intent(in)    ::  latNameIn ! name of  lat variable
   character(*),optional       ,intent(in)    :: maskNameIn ! name of mask variable
   character(*),optional       ,intent(in)    :: areaNameIn ! name of area variable
   character(*),optional       ,intent(in)    :: fracNameIn ! name of frac variable
   logical     ,optional       ,intent(in)    :: readFrac   ! T <=> also read frac
   integer (IN),optional       ,intent(out)   :: rc         ! return code

!EOP

   !----- local -----
   integer(IN)   :: i,j,n               ! generic indices
   integer(IN)   :: ni,nj               ! size of i and j dimensions
   integer(IN)   :: rCode               ! return code
   character(CL) :: remoteFn,localFn    ! file name with and without dir path
   character(CL) :: str                 ! data file ID string
   character(CL) :: name                ! name given to domain
   logical       :: readFrac2           ! T <=> read optional frac data

   real   (R8),allocatable ::  lon(:,:) ! temp array for domain lon  info
   real   (R8),allocatable ::  lat(:,:) ! temp array for domain lat  info
   integer(IN),allocatable :: mask(:,:) ! temp array for domain mask info
   real   (R8),allocatable :: area(:,:) ! temp array for domain area info
   real   (R8),allocatable :: frac(:,:) ! temp array for domain frac info

   character(CS) ::  latName            ! name of  lat variable
   character(CS) ::  lonName            ! name of  lon variable
   character(CS) :: maskName            ! name of mask variable
   character(CS) :: areaName            ! name of area variable
   character(CS) :: fracName            ! name of area variable

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_readData')"
   character(*),parameter :: F00    = "('(dshr_domain_readData) ',4a)" 
   character(*),parameter :: F01    = "('(dshr_domain_readData) ',a,2f14.6)"
   character(*),parameter :: F02    = "('(dshr_domain_readData) ',a,2i5)"
   character(*),parameter :: F04    = "('(dshr_domain_readData) ',i6,4a)" 

!----------------------------------------------------------------------------
! Notes:
! o as per shr_file_get(), the file name format is expected to be
!   remoteFn = [location:][directory path]localFn
!   eg. "foobar.nc"  "/home/user/foobar.nc"  "mss:/USER/fobar.nc"
! o assumes a very specific netCDF domain file format wrt var names, etc.
!----------------------------------------------------------------------------

   if (present(rc)) rc = 0
   readFrac2 = .false.  ! true <=> read cell frac data (optional)
   if (present(readFrac)) readFrac2 = readFrac 

   !-------------------------------------------------------------------------
   if (debug>0) write(6,F00) 'get the domain data file into the cwd'
   !-------------------------------------------------------------------------
   remoteFn = fileName
   n = shr_string_lastIndex(remoteFn,"/")
   if (n==0) n = index(remoteFn,":")
   if (n==0) then
      localFn = remoteFn
   else
      localFn = remoteFn(n+1: len_trim(remoteFn) )
   end if
   call shr_file_get(rCode,localFn,remoteFn)

   !-------------------------------------------------------------------------
   if (debug>0) write(6,F00) 'set names of domain data variables'
   !-------------------------------------------------------------------------
    lonName = "xc"  ! default values / standard data model domain file format
    latName = "yc"
   maskName = "mask"
   areaName = "area"
   fracName = "frac"
   if (present( lonNameIn))  lonName =  lonNameIn
   if (present( latNameIn))  latName =  latNameIn
   if (present(maskNameIn)) maskName = maskNameIn
   if (present(areaNameIn)) areaName = areaNameIn
   if (present(fracNameIn)) fracName = fracNameIn

   if (debug>0) then
      write(6,F00) 'longitude name: ',trim(lonName)
      write(6,F00) 'latitude  name: ',trim(latName)
      write(6,F00) 'mask      name: ',trim(maskName)
      write(6,F00) 'area      name: ',trim(areaName)
      write(6,F00) 'frac      name: ',trim(fracName)
   end if

   !-------------------------------------------------------------------------
   if (debug>0) write(6,F00) 'get domain dimension from domain data file'
   !-------------------------------------------------------------------------
!  name = localFn
   name = 'unset'
   do n=1,dshr_domain_nloc
      if (trim(dshr_domain_fname(n)) == trim(localFn)) then
         name = dshr_domain_dname(n)
         exit
      endif
   enddo
   if (trim(name) == 'unset') then
      dshr_domain_nloc = dshr_domain_nloc + 1
      if (dshr_domain_nloc > dshr_domain_nmax) then
         call dshr_domain_abort(subName//"ERROR: dshr_domain_nmax hit")
      endif
      dshr_domain_fname(dshr_domain_nloc) = trim(localFn)
      write(dshr_domain_dname(dshr_domain_nloc),'(a,i3.3)')  &
         trim(dshr_domain_dnameh),dshr_domain_nloc
      name = dshr_domain_dname(dshr_domain_nloc)
   endif
   call shr_ncread_varDimSizes(localFn,maskName,ni,nj)
   if (debug>0) then
      write(6,F00) 'dname = ',trim(name)
      write(6,F00) 'fname = ',trim(localFn)
      write(6,F02) 'ni,nj = ',ni,nj
   end if
   call dshr_domain_create(domain,name,ni,nj,rCode)

   !-------------------------------------------------------------------------
   if (debug>0) write(6,F00) 'open data file, get dims, alloc domain data'
   !-------------------------------------------------------------------------
   allocate(lon(ni,nj),lat(ni,nj),area(ni,nj),mask(ni,nj),frac(ni,nj))

   !--- read data ---
   if (readFrac2) then
      call shr_ncread_domain(localFn,lonName,lon,latName,lat,maskName,mask,areaName,area,fracName,frac)
   else
      call shr_ncread_domain(localFn,lonName,lon,latName,lat,maskName,mask,areaName,area)
      where (mask == 0)
         frac = 0.0_R8
      elsewhere
         frac = 1.0_R8
      end where
   endif

   domain%Fn   = localFn
   domain%lon  = lon
   domain%lat  = lat
   domain%mask = real(mask,R8)
   domain%area = area
   domain%frac = frac

   deallocate(lon,lat,mask,area,frac)

   do i = 1,ni
   do j = 1,nj
      domain%index(i,j) = (j-1)*ni + i
   end do
   end do
   
   !-------------------------------------------------------------------------
   if (debug>0) write(6,F00) 'done.'
   !-------------------------------------------------------------------------
   call dshr_domain_setFilled(domain,.true.)
   call dshr_domain_print    (domain)

end subroutine dshr_domain_readData

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_fill -- create domain and fill with bogus test data
!
! !DESCRIPTION:
!    Create domain and fill with test data (a somewhat reasonable uniform grid)
!
! !REVISION HISTORY:
!    2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_fill(domain,name,ni,nj,rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(inout) :: domain  ! domain 
   character(*),optional       ,intent(in)    :: name    ! domain name
   integer(IN) ,optional       ,intent(in)    :: ni,nj   ! domain dimensions
   integer(IN) ,optional       ,intent(out)   :: rc      ! return code

!EOP

   !----- local -----
   integer(IN)          :: i,j        ! generic indices
   integer(IN)          :: lni,lnj    ! local copies of ni,nj
   real(R8)             :: dx,dy      ! regular lat/lon dx/dy
   real(R8),allocatable :: temp(:,:)  ! temp array for domain info
   integer(IN)          :: rCode      ! return code

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_fill') "
   character(*),parameter :: F01     = "('(dshr_domain_fill) ',a,2f12.6)"

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   if (present(rc)) rc = 0

   !--- assign domain name and allocate domain data ---
   if (present(name) .and. present(ni) .and. present(nj)) then
      call dshr_domain_create(domain,name,ni,nj,rCode)
   endif

   if (dshr_domain_checkInit(domain)) then
      call dshr_domain_getdims(domain,lni,lnj,rCode)
   else
      call dshr_domain_abort(subName//"ERROR: domain not initialized")
      if (present(rc)) rc = 1
      return
   endif

   !-------------------------------------------------------------------------
   ! set domain data
   !-------------------------------------------------------------------------

   allocate(temp(lni,lnj),stat=rCode)
   if (rCode /= 0) then
      call dshr_domain_abort(subName//'problem with temp allocation')
   endif

   !--- longitude data ---
   dx = 360._R8/real(lni,R8)
   do i = 1,lni
      temp(i,:) = (real(i,R8)-0.5_R8)*dx
   enddo
   write(6,F01) 'min/max longitude  value: ',minval(temp),maxval(temp)
   call dshr_domain_putData(domain, temp, "lon", rc)

   !--- latitude data ---
   dy = 180._R8/real(lnj,R8)
   do j = 1,lnj
      temp(:,j) = (real(j,R8)-0.5_R8)*dy - 90._R8
   enddo
   write(6,F01) 'min/max latitude  value: ',minval(temp),maxval(temp)
   call dshr_domain_putData(domain, temp, "lat", rc)

   !--- mask data ---
   temp(:,:) = 1._R8
   if (lni > 10 .and. lnj > 10) then
      temp(2:max(4,lni/3-5),2:max(4,lnj/3-5)) = 0._R8
      temp(2*lni/3,3) = 0._R8
      temp(2*lni/3,5) = 0._R8
      temp(2*lni/3+2,3) = 0._R8
      do j = lnj/2+1,4*lnj/5
      do i = lni/2+1,3*lni/4
         temp(i,j) = 0._R8
      enddo
      enddo
   else
      temp(lni/3,lnj/3) = 0._R8
      temp(lni/2:3*lni/4,3*lnj/4:max(3*lnj/4,lnj-1)) = 0._R8
   endif
   write(6,F01) 'min/max mask      value: ',minval(temp),maxval(temp)
   call dshr_domain_putData(domain, temp, "mask", rc)
   call dshr_domain_putData(domain, temp, "area", rc)

   !--- area data ---
   call dshr_domain_getData(domain, temp, "lat", rc)
   do j = 1,lnj
   do i = 1,lni
      temp(i,j) = dx*pi/180._R8*dy*pi/180._R8*cos(temp(i,j)*pi/180._R8)
   enddo
   enddo
   write(6,F01) 'min/max area      value: ',minval(temp),maxval(temp)
   call dshr_domain_putData(domain, temp, "area", rc)

   deallocate(temp)
   call dshr_domain_setFilled(domain,.true.)

end subroutine dshr_domain_fill

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_compare -- returns TRUE if two domains are the same.
!
! !DESCRIPTION:
!    Returns TRUE if two domains are the the same (within tolerance).
!
! !REVISION HISTORY:
!    2005-May-03 - B. Kauffman - added mulitiple methods
!    2005-May-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function dshr_domain_compare(dom1,dom2,method,eps,rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)  :: dom1   ! 1st domain
   type(dshr_domain_domainType),intent(in)  :: dom2   ! 2nd domain
   integer(IN)                 ,intent(in)  :: method ! selects what to compare
   real(R8)   ,optional        ,intent(in)  :: eps    ! epsilon compare value
   integer(IN),optional        ,intent(out) :: rc     ! return code

!EOP

   !--- local ---
   real(R8)    :: leps         ! local epsilon
   integer(IN) :: i,j          ! counters
   integer(IN) :: rCode        ! error code
 
   !--- formats ---
   character(*),parameter :: subName = '(dshr_domain_compare) '
   character(*),parameter :: F01     = "('(dshr_domain_compare) ',4a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rCode = 0

   !--- set default epsilon, over-ride with user supplied epsilon ---
   leps = 1.0e-6_R8
   if (present(eps)) leps = eps

   !--- assume domains are the same, then look for differences ---
   dshr_domain_compare = .true. 

   !--- verify that domains are initialized ---
   if  (.not. dshr_domain_checkInit(dom1)) then
      call dshr_domain_abort(subName//"ERROR: dom1 not initialized")
      dshr_domain_compare = .false.
      rCode = 1
   end if
   if  (.not. dshr_domain_checkInit(dom2)) then
      call dshr_domain_abort(subName//"ERROR: dom2 not initialized")
      dshr_domain_compare = .false.
      rCode = 1
   end if

   !--- verify that domains have been filled with data ---
   if  (.not.dshr_domain_checkFilled(dom1))  then
      call dshr_domain_abort(subName//"ERROR: dom1 not filled")
      dshr_domain_compare = .false.
      rCode = 1
   end if
   if  (.not.dshr_domain_checkFilled(dom2))  then
      call dshr_domain_abort(subName//"ERROR: dom2 not filled")
      dshr_domain_compare = .false.
      rCode = 1
   end if

   !--- verify that domains have the same array size ---
   if (dom1%ni /= dom2%ni) dshr_domain_compare = .false.
   if (dom1%nj /= dom2%nj) dshr_domain_compare = .false.

   if (.not. dshr_domain_compare ) then
      !--- already failed the comparison test, check no futher ---
   else if (method == dshr_domain_compareXYabs      ) then
      do j=1,dom1%nj
      do i=1,dom1%ni
        if (abs(dom1%lat(i,j)-dom2%lat(i,j)) > leps) dshr_domain_compare = .false.
        if (abs(dom1%lon(i,j)-dom2%lon(i,j)) > leps) dshr_domain_compare = .false.
      enddo
      enddo
   else if (method == dshr_domain_compareXYrel      ) then
      do j=1,dom1%nj
      do i=1,dom1%ni
         if (rdiff(dom1%lat(i,j),dom2%lat(i,j)) > leps) dshr_domain_compare = .false.
         if (rdiff(dom1%lon(i,j),dom2%lon(i,j)) > leps) dshr_domain_compare = .false.
      enddo
      enddo
   else if (method == dshr_domain_compareMaskIdent  ) then
      do j=1,dom1%nj
      do i=1,dom1%ni
         if (dom1%mask(i,j) /= dom2%mask(i,j) ) dshr_domain_compare = .false.
      enddo
      enddo
   else if (method == dshr_domain_compareMaskZeros  ) then
      do j=1,dom1%nj
      do i=1,dom1%ni
         if (dom1%mask(i,j) == 0 .and. dom2%mask(i,j) /= 0) dshr_domain_compare = .false.
         if (dom1%mask(i,j) /= 0 .and. dom2%mask(i,j) == 0) dshr_domain_compare = .false.
      enddo
      enddo
   else if (method == dshr_domain_compareMaskSubset ) then
      do j=1,dom1%nj
      do i=1,dom1%ni
         if (dom1%mask(i,j) /= 0 .and. dom2%mask(i,j) == 0) dshr_domain_compare = .false.
      enddo
      enddo
   else
      write(6,F01) "ERROR: compare method not recognized, method = ",method
      call dshr_domain_abort(subName//"ERROR: compare method not recognized")
      if (present(rc)) rc = 1
   end if

   if (.not.dshr_domain_compare) rCode = 1
   if (present(rc)) rc = rCode

   return

!-------------------------------------------------------------------------------
contains   ! internal subprogram
!-------------------------------------------------------------------------------

   real(R8) function rdiff(v1,v2) ! internal function
      !------------------------------------------
      real(R8),intent(in) :: v1,v2                 ! two values to compare
      real(R8),parameter  :: c0           = 0.0_R8 ! zero
      real(R8),parameter  :: large_number = 1.0e20_R8 ! infinity
      !------------------------------------------
      if (v1 == v2) then
         rdiff = c0
      elseif (v1 == c0 .and. v2 /= c0) then
         rdiff = large_number
      elseif (v2 == c0 .and. v1 /= c0) then
         rdiff = large_number
      else
         rdiff = abs((v2-v1)/v1)
      endif
      !------------------------------------------
   end function rdiff

end function dshr_domain_compare

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_assignPtr -- Assign a pointer to a field in a domain
!
! !DESCRIPTION:
!     Assign a pointer to a field in a domain
!     \newline
!     call dshr\_domain\_assignPtrF(domain,"taux",taux\_ptr,rc)
!
! !REVISION HISTORY:
!     2006-Jun-02 - first version, based on dshr_bundle_assignPtr
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_assignPtr(domain,fldStr,fldPtr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)    :: domain      ! domain to point into
   character(*)                ,intent(in)    :: fldStr      ! name of field
   real(R8)                    ,pointer       :: fldPtr(:,:) ! pointer to field
   integer(IN),optional        ,intent(out)   :: rc          ! return code

!EOP

   !----- local -----
   integer(IN)            :: lrc  ! local return code

   !----- formats -----
   character(*),parameter :: subName =   "(dshr_domain_assignPtr)"
   character(*),parameter :: F00     = "('(dshr_domain_assignPtr) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lrc = 0
   if (.not. dshr_domain_checkInit(domain)) then
      call dshr_domain_abort(subName//"domain not initialized")
      lrc = 1
      fldPtr => null()
      return
   endif

   if      (trim(fldStr) == "lat"   ) then
      fldPtr => domain%lat
   else if (trim(fldStr) == "lon"   ) then
      fldPtr => domain%lon
   else if (trim(fldStr) == "area"  ) then
      fldPtr => domain%area
   else if (trim(fldStr) == "frac"  ) then
      fldPtr => domain%frac
   else if (trim(fldStr) == "mask"  ) then
      fldPtr => domain%mask
   else if (trim(fldStr) == "index" ) then
      fldPtr => domain%index
   else
      write(6,F00) "ERROR: requested field not a vaild part of a domain: ",trim(fldStr)
      fldPtr => null()
      lrc = 1
      call dshr_domain_abort(subName//" field not in domain")
   endif

   if (present(rc)) rc = lrc

end subroutine dshr_domain_assignPtr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_getDims -- returns domain dimension
!
! !DESCRIPTION:
!    Returns domain dimensions
!
! !REVISION HISTORY:
!    2005-Jan-27 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_getDims(domain,ni,nj,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)  :: domain ! domain in question
   integer(IN)                 ,intent(out) :: ni,nj  ! domain dimensions
   integer(IN),optional        ,intent(out) :: rc     ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_getDims')"
   character(*),parameter :: F00     = "('(dshr_domain_getDims) ',4a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   ni = domain%ni
   nj = domain%nj

end subroutine dshr_domain_getDims

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_getData -- copy field data out of domain data type
!
! !DESCRIPTION:
!    Copy field data out of domain data type
!
! !REVISION HISTORY:
!    2005-Jan-27 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_getData(domain,data,fldStr,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)   :: domain    ! domain to get data from
   real(R8)                    ,intent(out)  :: data(:,:) ! data to get from domain
   character(*)                ,intent(in)   :: fldStr    ! string 
   integer(IN),optional        ,intent(out)  :: rc        ! return code

!EOP

   !----- integer -----
   integer(IN)            :: ni,nj ! size of domain's i & j dimensions
   character(CL)          :: idStr ! domain's name

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_getData')"
   character(*),parameter :: F00     = "('(dshr_domain_getData) ',4a)" 
   character(*),parameter :: F01     = "('(dshr_domain_getData) ',a,i2,a)" 
   character(*),parameter :: F02     = "('(dshr_domain_getData) ',a,2i7)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   call dshr_domain_getDims(domain,ni,nj) 
   idStr = domain%name

   !--- confirm proper dimensions of input data ---
   if (size(data,1) /= ni .or. size(data,2) /= nj) then
      write(6,F00) "ERROR: field=",trim(fldStr),", domain=",trim(idStr)
      write(6,F00) "ERROR: input data dims don't match domain dims"
      write(6,F02) "ERROR: domain ni,nj= ",ni,nj
      write(6,F02) "ERROR: field  ni,nj= ",size(data,1),size(data,2)
      if (present(rc)) rc = 1
      return
   end if
 
   !--- copy the data ---
   if     (trim(fldStr) == "lat"  ) then
      data = domain%lat   
   elseif (trim(fldStr) == "lon"  ) then
      data = domain%lon   
   elseif (trim(fldStr) == "mask" ) then
      data = domain%mask  
   elseif (trim(fldStr) == "frac" ) then
      data = domain%frac  
   elseif (trim(fldStr) == "area" ) then
      data = domain%area  
   elseif (trim(fldStr) == "index") then
      data = domain%index 
   else 
      write(6,F00) "WARNING: unknown field=",trim(fldStr),", domain=",trim(idStr)
      if (present(rc)) rc = 1
   endif

end subroutine dshr_domain_getData

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_putData -- copy field data into domain data type
!
! !DESCRIPTION:
!    Copy field data into domain data type
!
! !REVISION HISTORY:
!    2005-Jan-27 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_putData(domain,data,fldStr,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(inout) :: domain    ! domain to put data into
   real(R8)                    ,intent(in)    :: data(:,:) ! data to put into domain
   character(*)                ,intent(in)    :: fldStr    ! field ID string
   integer(IN),optional        ,intent(out)   :: rc        ! return code

!EOP

   !----- integer -----
   integer(IN)            :: ni,nj ! size of domain's i & j dimensions
   character(CL)          :: idStr ! domain's name

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_putData')"
   character(*),parameter :: F00     = "('(dshr_domain_putData) ',4a)" 
   character(*),parameter :: F01     = "('(dshr_domain_putData) ',a,i2,a)" 
   character(*),parameter :: F02     = "('(dshr_domain_putData) ',a,2i7)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   call dshr_domain_getDims(domain,ni,nj)
   idStr = domain%name

   !--- confirm proper size of input data ---
   if (size(data,1) /= ni .or. size(data,2) /= nj) then
      write(6,F00) "ERROR: field=",trim(fldStr),", domain=",trim(idStr)
      write(6,F00) "ERROR: input field is not the same size as the domain"
      write(6,F02) "ERROR: domain ni,nj= ",ni,nj
      write(6,F02) "ERROR: field  ni,nj= ",size(data,1),size(data,2)
      if (present(rc)) rc = 1
      return
   end if

   !--- copy in the data ---
   if     (trim(fldStr) == "lat"  ) then
      domain%lat   = data
   elseif (trim(fldStr) == "lon"  ) then
      domain%lon   = data
   elseif (trim(fldStr) == "mask" ) then
      domain%mask  = data
   elseif (trim(fldStr) == "frac" ) then
      domain%frac  = data
   elseif (trim(fldStr) == "area" ) then
      domain%area  = data
   elseif (trim(fldStr) == "index") then
      domain%index = data
   else 
      write(6,F00) "WARNING: unknown field=",trim(fldStr),", domain=",trim(idStr)
      if (present(rc)) rc = 1
   endif

end subroutine dshr_domain_putData

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_getName -- returns domain name
!
! !DESCRIPTION:
!    Returns domain name
!
! !REVISION HISTORY:
!    2005-Jan-27 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_getName(domain,name,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)   :: domain ! domain in question
   character(*)                ,intent(out)  :: name   ! domain name
   integer(IN),optional        ,intent(out)  :: rc     ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_getName')"
   character(*),parameter :: F00     = "('(dshr_domain_getName) ',4a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   name = domain%name

end subroutine dshr_domain_getName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_putName -- change domain name
!
! !DESCRIPTION:
!    Change domain name
!
! !REVISION HISTORY:
!    2005-Jan-27 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_putName(domain,name,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(out) :: domain ! domain in question
   character(*)                ,intent(in)  :: name   ! domain name
   integer(IN),optional        ,intent(out) :: rc     ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_putName')"
   character(*),parameter :: F00     = "('(dshr_domain_putName) ',4a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   domain%name = name

end subroutine dshr_domain_putName

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_getFn -- returns domain Fn
!
! !DESCRIPTION:
!    Returns domain Fn
!
! !REVISION HISTORY:
!    2005-Jul-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_getFn(domain,Fn,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)   :: domain ! domain in question
   character(*)                ,intent(out)  :: Fn     ! domain Fn
   integer(IN),optional        ,intent(out)  :: rc     ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_getFn')"
   character(*),parameter :: F00     = "('(dshr_domain_getFn) ',4a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   Fn = domain%Fn

end subroutine dshr_domain_getFn

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_putFn -- change domain Fn
!
! !DESCRIPTION:
!    Change domain Fn
!
! !REVISION HISTORY:
!    2005-Jul-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_putFn(domain,Fn,rc)

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(out) :: domain ! domain in question
   character(*)                ,intent(in)  :: Fn     ! domain Fn
   integer(IN),optional        ,intent(out) :: rc     ! return code

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_putFn')"
   character(*),parameter :: F00     = "('(dshr_domain_putFn) ',4a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   domain%Fn = Fn

end subroutine dshr_domain_putFn

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_print -- print domain information
!
! !DESCRIPTION:
!    Print domain information
!
! !REVISION HISTORY:
!    2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_print(domain,rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(inout) :: domain  ! domain 
   integer(IN),optional        ,intent(out)   :: rc      ! return code

!EOP

   !----- local -----
   integer(IN)   :: i,j               ! index
   integer(IN)   :: iCnt              ! counter
   integer(IN)   :: rCode             ! return code

   !----- formats -----
   character(*),parameter :: subName = "('dshr_domain_print') "
   character(*),parameter :: F00     = "('(dshr_domain_print) ',a)" 
   character(*),parameter :: F01     = "('(dshr_domain_print) ',a,2f14.6)"
   character(*),parameter :: F02     = "('(dshr_domain_print) ',a,2i8)"
   character(*),parameter :: F04     = "('(dshr_domain_print) ',i6,4a)" 
   character(*),parameter :: F05     = "('(dshr_domain_print) ',a,i8,a,i8,a)" 

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   if (present(rc)) rc = 0

   write(6,F00) 'name        : '//trim(domain%name)
   write(6,F00) 'filename    : '//trim(domain%Fn)
   write(6,F00) 'init        : '//trim(domain%initialized)
   write(6,F00) 'filled      : '//trim(domain%filled)
   write(6,F02) 'size (ni,nj): ',domain%ni,domain%nj
   write(6,F01) 'min/max lat : ',minval(domain%lat)  ,maxval(domain%lat)
   write(6,F01) 'min/max lon : ',minval(domain%lon)  ,maxval(domain%lon)
   write(6,F01) 'min/max area: ',minval(domain%area) ,maxval(domain%area)
   write(6,F01) 'min/max frac: ',minval(domain%frac) ,maxval(domain%frac)
   write(6,F01) 'min/max mask: ',minval(domain%mask) ,maxval(domain%mask)
   write(6,F01) 'min/max indx: ',minval(domain%index),maxval(domain%index)

   iCnt = 0
   do j=1,domain%nj
   do i=1,domain%ni
      if (nint(domain%mask(i,j)) /= 0) iCnt = iCnt + 1
   enddo
   enddo
   write(6,F05) 'mask: active cells: ',iCnt,' out of ',domain%nj*domain%ni,' pts'

end subroutine dshr_domain_print

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_setInit -- set initialized flag in domain
!
! !DESCRIPTION:
!    Set initialized flag in domain
!
! !REVISION HISTORY:
!    2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_setInit(domain,iFlag,rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(inout):: domain
   logical                     ,intent(in)   :: iFlag  ! setting
   integer(IN),optional        ,intent(out)  :: rc     ! return code

!EOP

   !--- local ---
   character(*),parameter :: subName = '(dshr_domain_setInit) '

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   if (iFlag) then
      domain%initialized = dshr_domain_setTru
   else
      domain%initialized = dshr_domain_setFal
   endif

   if (present(rc)) rc = 0

end subroutine dshr_domain_setInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_checkInit -- check if domain has been initialized
!
! !DESCRIPTION:
!    Check if domain has been initialized
!
! !REVISION HISTORY:
!    2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function dshr_domain_checkInit(domain,rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)  :: domain
   integer(IN),optional        ,intent(out) :: rc     ! return code

!EOP

   !--- local ---
   integer(IN) :: rCode
   character(*),parameter :: subName = '(dshr_domain_checkInit) '

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   dshr_domain_checkInit = .false.
   rCode = 1

   if (trim(domain%initialized) == trim(dshr_domain_setTru)) then
      dshr_domain_checkInit = .true.
      rCode = 0
   endif

   if (present(rc)) rc = rCode

end function dshr_domain_checkInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_setFilled -- set filled flag in domain
!
! !DESCRIPTION:
!    Set filled flag in domain
!
! !REVISION HISTORY:
!    2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_setFilled(domain,iFlag,rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(inout):: domain
   logical                     ,intent(in)   :: iFlag  ! setting
   integer(IN),optional        ,intent(out)  :: rc     ! return code

!EOP

   !--- local ---
   character(*),parameter :: subName = '(dshr_domain_setFilled) '

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   if (iFlag) then
      domain%filled = dshr_domain_setTru
   else
      domain%filled = dshr_domain_setFal
   endif

   if (present(rc)) rc = 0

end subroutine dshr_domain_setFilled

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_checkFilled -- check if domain has been filled with data
!
! !DESCRIPTION:
!    Check if domain has been filled with data
!
! !REVISION HISTORY:
!    2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function dshr_domain_checkFilled(domain,rc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_domain_domainType),intent(in)  :: domain
   integer(IN),optional        ,intent(out) :: rc     ! return code

!EOP

   !--- local ---
   integer(IN) :: rCode
   character(*),parameter :: subName = '(dshr_domain_checkFilled) '

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   dshr_domain_checkFilled = .false.
   rCode = 1

   if (trim(domain%filled) == trim(dshr_domain_setTru)) then
      dshr_domain_checkFilled = .true.
      rCode = 0
   endif

   if (present(rc)) rc = rCode

end function dshr_domain_checkFilled

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_setAbort -- Set local dshr_domain abort flag
!
! !DESCRIPTION:
!    Set local flag, true = abort, false = print and continue
!    \newline
!    call dshr\_domain\_setAbort(.false.)
!
! !REVISION HISTORY:
!    2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_setAbort(flag)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   logical,intent(in) :: flag

!EOP

   !--- local ---

   !--- formats ---
   character(*),parameter :: subName = "('dshr_domain_setAbort') "
   character(*),parameter :: F00     = "('(dshr_domain_setAbort) ',a) "

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   doabort = flag

end subroutine dshr_domain_setAbort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_abort -- abort method
!
! !DESCRIPTION:
!    Abort method for dshr\_domain
!
! !REVISION HISTORY:
!    2005-Mar-27 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_abort(string)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),optional,intent(in) :: string

!EOP

   !--- local ---
   character(CL)          :: lstring
   character(*),parameter :: subName = '(dshr_domain_abort) '

!----------------------------------------------------------------------------
! Notes:
!----------------------------------------------------------------------------

   lstring = ''
   if (present(string)) lstring = string

   write(6,*) trim(subName),' ' ,trim(lstring)
   if (doabort) call shr_sys_abort(lstring)

end subroutine dshr_domain_abort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_setDebug -- Set local shr_domain debug level
!
! !DESCRIPTION:
!    Set local debug level, 0 = production
!    \newline
!    General Usage: call shr\_domain\_setDebug(2)
!
! !REVISION HISTORY:
!    2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_setDebug(level)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer,intent(in) :: level

!EOP

   !--- formats ---
   character(*),parameter :: subName = "('dshr_domain_setDebug') "
   character(*),parameter :: F01     = "('(dshr_domain_setDebug) ',a,i4) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   debug = level
   write(6,F01) "debug level reset to ",level

end subroutine dshr_domain_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_domain_getDebug -- get shr_domain internal debug level
!
! !DESCRIPTION:
!    Get internal debug level, 0 = production
!    \newline
!    General Usage: call shr\_domain\_getDebug(level)
!
! !REVISION HISTORY:
!    2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_domain_getDebug(level)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer,intent(out) :: level

!EOP

   !--- formats ---
   character(*),parameter :: subName = "('dshr_domain_getDebug') "
   character(*),parameter :: F00     = "('(dshr_domain_getDebug) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   level = debug

end subroutine dshr_domain_getDebug

!===============================================================================
!===============================================================================

end module dshr_domain
