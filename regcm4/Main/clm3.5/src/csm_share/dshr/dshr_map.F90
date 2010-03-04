!===============================================================================
! SVN $Id: dshr_map.F90 3488 2007-03-13 00:05:59Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_map.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_map -- map data type and associated methods
!
! !DESCRIPTION:
!     Map data type and associated methods.
!     Built heavily on shr\_map
!     \newline
!     General Usage:
!        type(dshr\_map\_mapType) :: mymap
!        type(dshr\_domain\_domainType) :: Dsrc,Ddst
!        type(dshr\_bundle\_bundleType) :: Bsrc,Bdst
!        call dshr\_map\_mapSet(mymap,Dsrc,Ddst,'map1','remap',shr\_map\_fs\_nn','dstmask')
!        call dshr\_map\_mapData(Bsrc,Bdst,mymap)
!     \newline
!        call dshr\_map\_mapData(Bsrc,Bdst,'nosave','remap','bilinear','nomask','vector')
!     
! !REMARKS:
!     See shr\_map\_mod for more detail information on options
!     
! !REVISION HISTORY:
!     2005-Mar-10 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_map

! !USES:

  use shr_timer_mod   ! timer routines
  use shr_sys_mod     ! system routines
  use shr_string_mod  ! string & list handling
  use shr_map_mod     ! mapping routines

  use dshr_kind
  use dshr_domain
  use dshr_bundle

  implicit none
  private

! !PUBLIC TYPES:

  public :: dshr_map_mapType

  type dshr_map_mapType                          ! like cpl_map datatype
    private
    character(CS)                        :: name ! map name
    type(shr_map_mapType)                :: sMat ! weights and more
    type(dshr_domain_domainType),pointer :: src  ! pointer to src domain
    type(dshr_domain_domainType),pointer :: dst  ! pointer to dst domain
  end type dshr_map_mapType

! !PUBLIC MEMBER FUNCTIONS:

  public :: dshr_map_checkInit         ! check whether map type is set
  public :: dshr_map_checkFilled       ! check whether map wts are set 
  public :: dshr_map_listValidOpts     ! list valid map setting options
  public :: dshr_map_get               ! get stuff out of the datatype
  public :: dshr_map_put               ! put stuff into the datatype
  public :: dshr_map_sMatPtr           ! set a pointer to sMat in map
  public :: dshr_map_mapSet            ! compute weights in map
  public :: dshr_map_mapData           ! map data
  public :: dshr_map_clean             ! clean map datatype
  public :: dshr_map_print             ! print map datatype info
  public :: dshr_map_setDebug          ! set local debug flag
  public :: dshr_map_setTiming         ! set local timing flag
  public :: dshr_map_setAbort          ! set local abort flag

! !PUBLIC DATA MEMBERS:

  character(CS),public,parameter :: dshr_map_fs_name   = shr_map_fs_name

!EOP

  !--- private module parameters ---
  integer(IN),parameter :: dshr_map_ispval  = shr_map_ispval
  real(R8)   ,parameter :: dshr_map_spval   = shr_map_spval
  integer(IN),save      :: debug  = 0           ! debug  level
  integer(IN),save      :: timing = 0           ! timing level
  logical    ,save      :: doabort = .true.

  interface dshr_map_put ; module procedure &
    dshr_map_putCS, &
    dshr_map_putR8, &
    dshr_map_putIN
  end interface

  interface dshr_map_get ; module procedure &
    dshr_map_getCS, &
    dshr_map_getR8, &
    dshr_map_getIN
  end interface

  interface dshr_map_mapSet ; module procedure &
    dshr_map_mapSetD, &
    dshr_map_mapSetB
  end interface

  interface dshr_map_mapData ; module procedure &
    dshr_map_mapDataB, &
    dshr_map_mapDataA, &
    dshr_map_mapDataBnm
  end interface

!===============================================================================
contains
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_checkInit -- Check status of map
!
! !DESCRIPTION:
!     Check init status of map
!     Layer to shr\_map\_checkInit
!     \newline
!     test = dshr\_map\_checkInit(map)
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function dshr_map_checkInit(map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType),intent(in) :: map

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_checkInit) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  dshr_map_checkInit = shr_map_checkInit(map%sMat)

end function dshr_map_checkInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_checkFilled -- Check fill status of map
!
! !DESCRIPTION:
!     Check fill status of map
!     Layer to shr\_map\_checkFilled
!     \newline
!     test = dshr\_map\_checkFilled(map)
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

logical function dshr_map_checkFilled(map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType),intent(in) :: map

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_checkFilled) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  dshr_map_checkFilled = shr_map_checkFilled(map%sMat)

end function dshr_map_checkFilled

!==============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_listValidOpts -- list map options
!
! !DESCRIPTION:
!     List map options
!     Layer to shr\_map\_listValidOpts
!     \newline
!     call dshr\_map\_listValidOpts()
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_listValidOpts()

  implicit none

! !INPUT/OUTPUT PARAMETERS:

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_listValidOpts) '
  character(*),parameter :: F00     = "('(dshr_map_listValidOpts) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  write(6,F00) ':'
  write(6,F00) '  '//trim(dshr_map_fs_name)//'  : any character string'
  call shr_map_listValidOpts()
  call shr_sys_flush(6)

  end subroutine dshr_map_listValidOpts

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_getCS -- get char string from map
!
! !DESCRIPTION:
!     One method of dshr_map_get, get cval in map
!     Layer to shr\_map\_get
!     \newline
!     call dshr\_map\_get(map,'algo',cval)
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_getCS(map,fldStr,cval)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType),intent(in) :: map
  character(*)          ,intent(in) :: fldStr
  character(*)          ,intent(out):: cval

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_getCS) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  if (fldStr == dshr_map_fs_name) then
    cval = map%name
  else
    call shr_map_get(map%sMat,fldStr,cval)
  endif

end subroutine dshr_map_getCS

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_getIN -- get integer from map
!
! !DESCRIPTION:
!     One method of dshr_map_get, get ival in map
!     Layer to shr\_map\_get
!     \newline
!     call dshr\_map\_get(map,shr_map_fs_nwts,ival)
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_getIN(map,fldStr,ival)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType),intent(in) :: map
  character(*)          ,intent(in) :: fldStr
  integer(IN)  ,intent(out):: ival

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_getIN) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  ival = dshr_map_ispval
  call shr_map_get(map%sMat,fldStr,ival)

end subroutine dshr_map_getIN

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_getR8 -- get real from map
!
! !DESCRIPTION:
!     One method of dshr_map_get, get rval in map
!     Layer to shr\_map\_get
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_getR8(map,fldStr,rval)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType),intent(in) :: map
  character(*)          ,intent(in) :: fldStr
  real(R8)  ,intent(out):: rval

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_getR8) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  rval = dshr_map_spval
  call shr_map_get(map%sMat,fldStr,rval)

end subroutine dshr_map_getR8

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_putCS -- put char string to map
!
! !DESCRIPTION:
!     One method of dshr_map_put, put cval in map
!     Layer to shr\_map\_put
!     \newline
!     call dshr\_map\_put(map,'algo','spval')
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_putCS(map,fldStr,cval,verify)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType) ,intent(inout):: map
  character(*)          ,intent(in) :: fldStr
  character(*)          ,intent(in) :: cval
  logical,optional      ,intent(in) :: verify

!EOP

  !--- local ---
  logical :: lverify

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_putCS) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  lverify = .true.
  if (present(verify)) lverify=verify
  if (fldStr == dshr_map_fs_name) then
    map%name = cval
  endif
  call shr_map_put(map%sMat,fldStr,cval,verify=lverify)

end subroutine dshr_map_putCS

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_putIN -- put integer to map
!
! !DESCRIPTION:
!     One method of dshr_map_put, put ival in map
!     Layer to shr\_map\_put
!     \newline
!     call dshr\_map\_put(map,'nwts',-1)
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_putIN(map,fldStr,ival,verify)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType) ,intent(inout):: map
  character(*)          ,intent(in) :: fldStr
  integer(IN)           ,intent(in) :: ival
  logical,optional      ,intent(in) :: verify

!EOP

  !--- local ---
  logical :: lverify

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_putIN) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  lverify = .true.
  if (present(verify)) lverify=verify
  call shr_map_put(map%sMat,fldStr,ival,verify=lverify)

end subroutine dshr_map_putIN

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_putR8 -- put real to map
!
! !DESCRIPTION:
!     One method of dshr_map_put, put rval in map
!     Layer to shr\_map\_put
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_putR8(map,fldStr,rval,verify)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType) ,intent(inout):: map
  character(*)          ,intent(in) :: fldStr
  real(R8)              ,intent(in) :: rval
  logical,optional      ,intent(in) :: verify

!EOP

  !--- local ---
  logical :: lverify

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_putR8) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  lverify = .true.
  if (present(verify)) lverify=verify
  call shr_map_put(map%sMat,fldStr,rval,verify=lverify)

end subroutine dshr_map_putR8

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_sMatPtr -- Assign a pointer to an sMat in a map datatype
!
! !DESCRIPTION:
!     Assign a pointer to an sMat in a map datatype
!     \newline
!     call dshr\_map\_sMatPtr(map,sMat\_ptr)
!
! !REVISION HISTORY:
!     2005-Mar-28 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_sMatPtr(map,sMatPtr,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_map_mapType),target ,intent(in)  :: map     ! map to point into
   type( shr_map_mapType)        ,pointer     :: sMatPtr ! pointer to sMat
   integer(IN),optional          ,intent(out) :: rc      ! return code

!EOP

   !----- local -----
   integer(IN)            :: lrc  ! local return code

   !----- formats -----
   character(*),parameter :: subName =   '(dshr_map_sMatPtr) '
   character(*),parameter :: F00     = "('(dshr_map_sMatPtr) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lrc = 0

   !--- assign pointer ---
   sMatPtr => map%sMat

   if (present(rc)) rc = lrc

end subroutine dshr_map_sMatPtr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_mapSetD - Initialize map method, weights, etc
!
! !DESCRIPTION:
!     One of the methods associated with dshr_map_mapSet.
!     Generate weights based on input domains and some flags
!     (name, type, algo, mask).  Output is a filled map
!     including generated weights.
!     Calls shr\_map\_put and shr\_map\_mapSet
!     \newline
!     call dshr\_map\_mapSet(map,DS,DD,'mymap','remap',shr_map_fs_nn,'dstmask')
!
! !REMARKS:
!     See shr\_map\_mod for more infomation about valid options for switches
!
! !REVISION HISTORY:
!     2005-Mar-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_mapSetD(map,domsrc,domdst,name,type,algo,mask,vect,rc) 

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType)                      ,intent(inout):: map
  type(dshr_domain_domainType)         ,target,intent(in) :: domsrc    
  type(dshr_domain_domainType),optional,target,intent(in) :: domdst   
  character(*)                ,optional       ,intent(in) :: name
  character(*)                ,optional       ,intent(in) :: type
  character(*)                ,optional       ,intent(in) :: algo
  character(*)                ,optional       ,intent(in) :: mask
  character(*)                ,optional       ,intent(in) :: vect
  integer(IN)                 ,optional       ,intent(out):: rc

!EOP

  !--- local ---
  integer(IN)          :: ni,nj          ! dims
  integer(IN)          :: lrc            ! local rc
  real(R8),allocatable :: Xin(:,:)       ! input longitude
  real(R8),allocatable :: Yin(:,:)       ! input latitude
  integer(IN),allocatable:: Min(:,:)     ! input mask
  real(R8),allocatable :: Xout(:,:)      ! output longitude
  real(R8),allocatable :: Yout(:,:)      ! output latitude
  integer(IN),allocatable:: Mout(:,:)    ! output mask

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_mapSetD) '

!-------------------------------------------------------------------------------
! do not clean the map first, map might have attributes already set 
!-------------------------------------------------------------------------------

  map%src => domsrc
  if (present(domdst)) then
    map%dst => domdst
  else
    map%dst => domsrc
  endif
  
  if (present(name)) then
    call dshr_map_put(map,shr_map_fs_name,name)
  endif
  if (present(type)) then
    call shr_map_put(map%sMat,shr_map_fs_type,type)
  endif
  if (present(algo)) then
    call shr_map_put(map%sMat,shr_map_fs_algo,algo)
  endif
  if (present(mask)) then
    call shr_map_put(map%sMat,shr_map_fs_mask,mask)
  endif
  if (present(vect)) then
    call shr_map_put(map%sMat,shr_map_fs_vect,vect)
  endif

  if (.not.shr_map_checkInit(map%sMat)) then
    call dshr_map_abort(subName//'map is not initialized properly')
  endif

  call dshr_domain_getDims(map%src,ni,nj,rc)
  allocate(Xin(ni,nj),Yin(ni,nj),Min(ni,nj))
  call dshr_domain_getData(map%src,Xin,'mask',rc)
  Min = 1
  where (abs(Xin(:,:)) < 0.5_R8) Min(:,:) = 0
  call dshr_domain_getData(map%src,Xin,'lon',rc)
  call dshr_domain_getData(map%src,Yin,'lat',rc)

  call dshr_domain_getDims(map%dst,ni,nj,rc)
  allocate(Xout(ni,nj),Yout(ni,nj),Mout(ni,nj))
  call dshr_domain_getData(map%dst,Xout,'mask',rc)
  Mout = 1
  where (abs(Xout(:,:)) < 0.5_R8) Mout(:,:) = 0
  call dshr_domain_getData(map%dst,Xout,'lon',rc)
  call dshr_domain_getData(map%dst,Yout,'lat',rc)

  call shr_map_mapSet(map%sMat,Xin,Yin,Min,Xout,Yout,Mout,rc=lrc)
  
  if (.not.shr_map_checkFilled(map%sMat)) then
    call dshr_map_abort(subName//'map is not set properly')
  endif

  deallocate(Xin,Yin)
  deallocate(Xout,Yout)
  deallocate(Min,Mout)

  if (present(rc)) rc = lrc

end subroutine dshr_map_mapSetD

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_mapSetB - Initialize map method, weights, etc
!
! !DESCRIPTION:
!     Same as dshr\_map\_mapSetD except with bundle arguments instead of domains
!     \newline
!     call dshr\_map\_mapSet(map,BS,BD,'mymap','remap',shr_map_fs_nn,'dstmask')
!
! !REMARKS:
!     See shr\_map\_mod for more infomation about valid options for switches
!
! !REVISION HISTORY:
!     2005-Mar-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_mapSetB(map,bunsrc,bundst,name,type,algo,mask,vect,rc) 
  ! generate weights based on two input domains contained in bundles

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType)               ,intent(inout):: map
  type(dshr_bundle_bundleType)         ,intent(in)   :: bunsrc    ! use domain
  type(dshr_bundle_bundleType),optional,intent(in)   :: bundst    ! use domain
  character(*)                ,optional,intent(in)   :: name
  character(*)                ,optional,intent(in)   :: type
  character(*)                ,optional,intent(in)   :: algo
  character(*)                ,optional,intent(in)   :: mask
  character(*)                ,optional,intent(in)   :: vect
  integer(IN)                 ,optional,intent(out)  :: rc

!EOP

  !--- local ---
  integer(IN)    :: lrc
  character(CS)  :: lname
  character(CS)  :: ltype
  character(CS)  :: lalgo
  character(CS)  :: lmask
  character(CS)  :: lvect
  type(dshr_domain_domainType),pointer :: dsrc,ddst

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_mapSetB) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  if (present(name)) then
    lname = name
  else
    call dshr_map_get(map,shr_map_fs_name,lname)
  endif

  if (present(type)) then
    ltype = type
  else
    call dshr_map_get(map,shr_map_fs_type,ltype)
  endif

  if (present(algo)) then
    lalgo = algo
  else
    call dshr_map_get(map,shr_map_fs_algo,lalgo)
  endif

  if (present(mask)) then
    lmask = mask
  else
    call dshr_map_get(map,shr_map_fs_mask,lmask)
  endif

  if (present(vect)) then
    lvect = vect
  else
    call dshr_map_get(map,shr_map_fs_vect,lvect)
    if (trim(lvect) /= trim(shr_map_fs_vector)) then
      lvect = shr_map_fs_scalar
    endif
  endif

  call dshr_bundle_domainPtr(bunsrc,dsrc)
  if (present(bundst)) then
    call dshr_bundle_domainPtr(bundst,ddst)
  else
    call dshr_bundle_domainPtr(bunsrc,ddst)
  endif

  call dshr_map_mapSetD(map,dsrc,ddst,lname,ltype,lalgo,lmask,lvect,lrc)

  if (present(rc)) rc = lrc

end subroutine dshr_map_mapSetB

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_mapDataBnm - Map data without having map
!
! !DESCRIPTION:
!     Map data without saving map, similar to mapDataB with no map (nm).
!     Calls mapSet, uses map, then cleans map
!     \newline
!     call dshr\_map\_mapData(BS,BD,'mymap','remap',shr_map_fs_nn,'dstmask')
!
! !REMARKS:
!     See shr\_map\_mod for more infomation about valid options for switches
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_mapDataBnm(bunsrc,bundst,name,type,algo,mask,vect,fsrc,fdst,rc) 
  ! generate weights based on two input domains contained in bundles

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_bundle_bundleType)         ,intent(inout):: bunsrc    ! use domain
  type(dshr_bundle_bundleType),optional,intent(inout):: bundst    ! use domain
  character(*)                ,optional,intent(in)   :: name
  character(*)                         ,intent(in)   :: type
  character(*)                         ,intent(in)   :: algo
  character(*)                         ,intent(in)   :: mask
  character(*)                ,optional,intent(in)   :: vect
  character(*)                ,optional,intent(in)   :: fsrc
  character(*)                ,optional,intent(in)   :: fdst
  integer(IN)                 ,optional,intent(out)  :: rc

!EOP

  !--- local ---
  type(dshr_map_mapType):: map
  integer(IN)    :: lrc
  character(CS)  :: lvect
  character(CS)  :: lname
  logical        :: flists

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_mapDataBnm) '
  character(*),parameter :: F00   = "('(dshr_map_mapdataBnm) ',4a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  lrc = 0

  if (present(fsrc).and.present(fdst)) then
    flists = .true.
    if (.not.shr_string_listIsValid(fsrc).or. &
        .not.shr_string_listIsValid(fdst).or. &
        shr_string_listGetNum(fsrc) /= shr_string_listGetNum(fdst)) then
      write(6,F00) ' ERROR fsrc or fdst invalid'
      write(6,F00) ' fsrc  ',trim(fsrc)
      write(6,F00) ' fdst ',trim(fdst)
      lrc = 1
      call dshr_map_abort(subName//' fsrc and fdst invalid')
    endif
  elseif (.not.present(fsrc).and..not.present(fdst)) then
    flists = .false.
  else
    write(6,F00) ' ERROR fsrc and fdst must both or neither exist'
    lrc = 1
    call dshr_map_abort(subName//' ERROR fsrc, fdst ')
  endif

  if (present(name)) then
    lname = name
  else
    lname = 'none'
  endif

  if (present(vect)) then
    lvect = vect
  else
    lvect = shr_map_fs_scalar
  endif

  if (present(bundst)) then
    call dshr_map_mapSet(map,bunsrc,bundst,name=lname,type=type,algo=algo,mask=mask,vect=lvect,rc=lrc) 
    if (flists) then
      call dshr_map_mapData(bunsrc,bundst,map,fsrc=fsrc,fdst=fdst)
    else
      call dshr_map_mapData(bunsrc,bundst,map)
    endif
  else
    call dshr_map_mapSet(map,bunsrc,name=lname,type=type,algo=algo,mask=mask,vect=lvect,rc=lrc) 
    if (flists) then
      call dshr_map_mapData(bunsrc,bunsrc,map,fsrc=fsrc,fdst=fdst)
    else
      call dshr_map_mapData(bunsrc,bunsrc,map)
    endif
  endif

  call dshr_map_clean(map)

  if (present(rc)) rc = lrc

end subroutine dshr_map_mapDataBnm

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_mapDataB -- map data from bundle to bundle 
!
! !DESCRIPTION:
!     One of the methods associated with dshr\_map\_mapData.
!     Map arrays from source bundle to destination bundle using preset
!     map datatype.  Layer to shr\_map\_mapData for bundles
!     \newline
!     call dshr\_map\_mapData(BS,BD,map)
!
! !REVISION HISTORY:
!     2008-Mar-09  - B. Kauffman - optimize wrt matching bundle lists
!     2005-Mar-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_mapDataB(bunSrc,bunDst,map,fSrc,fDst)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_bundle_bundleType),intent(in)    :: bunSrc  
  type(dshr_bundle_bundleType),intent(inout) :: bunDst 
  type(dshr_map_mapType)      ,intent(in)    :: map
  character(*),optional       ,intent(in)    :: fSrc
  character(*),optional       ,intent(in)    :: fDst

!EOP

  !--- local ---
  integer(IN)          :: n,i,j,k           ! loop indicies
  integer(IN)          :: nis,njs,nfs       ! bundle i,j,k dims, src
  integer(IN)          :: nid,njd,nfd       ! bundle i,j,k dims, dst
  integer(IN)          :: cDate,sec         ! bundle date & sec
  integer(IN)          :: kSrc,kDst         ! bundle field indices
  character(CS)        :: namSrc,namDst     ! field name
  character(CX)        :: bunSrcf,bunDstf   ! field name lists
  character(CX)        :: lfSrc,lfDst       ! field name lists for src/dst
  real(R8),allocatable :: Asrc(:,:)         ! temp array of src fields (k,n)
  real(R8),allocatable :: Adst(:,:)         ! temp array of dst fields (k,n)
  real(R8),pointer     :: dataSrc(:,:)      ! pointer to src field (i,j)
  real(R8),pointer     :: dataDst(:,:)      ! pointer to dst field (i,j)
  integer(IN)          :: lrc               ! local rc
  logical              :: sameList          ! T <=> bunSrcf = bunDstf = lfSrc = lfDst
  integer(IN),save     :: t0=-1,t1,t2,t3,t4 ! timers

  !--- formats ---
  character(*),parameter :: subName =   '(dshr_map_mapDataB) '
  character(*),parameter :: F00     = "('(dshr_map_mapDataB) ',4a) "

!-------------------------------------------------------------------------------
! interpolate data from bunSrc to bunDst
! Note: most common use is two bundles with same lists, fSrc & fDst are absent
!-------------------------------------------------------------------------------

  lrc = 0

  if (timing>0 .and. t0 == -1) then
     call shr_timer_get(t0,subName//"everything")
     call shr_timer_get(t1,subName//"determine field lists")
     call shr_timer_get(t2,subName//"pack temp array")
     call shr_timer_get(t3,subName//"map")
     call shr_timer_get(t4,subName//"unpack temp array")
  end if
  if (timing>0) call shr_timer_start(t0)

  call dshr_bundle_getFieldList(bunSrc,bunSrcf)
  call dshr_bundle_getFieldList(bunDst,bunDstf)
  call dshr_bundle_getDims(bunSrc,nis,njs,nfs)
  call dshr_bundle_getDims(bunDst,nid,njd,nfd)

  !--- determine list of fields to map -----------------------------------------
  if (timing>0) call shr_timer_start(t1)
  sameList = .false.
  if (.not.present(fSrc).and..not.present(fDst)) then
    if ( trim(bunSrcf) == trim(bunDstf) ) then
       sameList = .true.
    else
       call shr_string_listIntersect(bunSrcf,bunDstf,lfSrc)
       lfDst = lfSrc
    end if
  else if (present(fSrc).and.present(fDst)) then
     lfSrc = trim(fSrc)
     lfDst = trim(fDst)
     nfs = shr_string_listGetNum(lfSrc)
     nfd = shr_string_listGetNum(lfDst)
     if (.not.shr_string_listIsValid(lfSrc).or. &
         .not.shr_string_listIsValid(lfDst).or. &
         nfs /= nfd) then
       write(6,F00) ' ERROR lfSrc or lfDst invalid'
       write(6,F00) ' lfSrc ',trim(lfSrc)
       write(6,F00) ' lfDst ',trim(lfDst)
       lrc = 1
       call dshr_map_abort(subName//'lfSrc and lfDst invalid')
     endif
  else
    write(6,F00) ' ERROR fSrc and fDst must both or neither exist'
    lrc = 1
    call dshr_map_abort(subName//'ERROR fSrc, fDst ')
  endif
  if (timing>0) call shr_timer_stop (t1)

  if (lrc /= 0) call dshr_map_abort(subName//'ERROR')

  !--- put data into temp/reshaped array ---------------------------------------
  if (timing>0) call shr_timer_start(t2)
  allocate(Asrc(nfs,nis*njs), Adst(nfd,nid*njd))
  do k = 1,nfs
     !--- get field index ---
     if (sameList) then
        kSrc = k
     else
        call shr_string_listGetName(lfSrc,k,namSrc)
        call shr_string_listGetIndex(bunSrcf,namSrc,kSrc)
        if (kSrc == 0) then
           write(6,F00) 'ERROR in src field for ',trim(namSrc)
           lrc = 1
           call dshr_map_abort(subName//'ERROR in name index')
        endif
     endif
     !--- get field data via pointer ---
     call dshr_bundle_assignPtr(bunSrc,kSrc,dataSrc)
     !--- copy dataSrc into Asrc  ---
     n = 0
     do j=1,njs
     do i=1,nis
        n = n + 1
        Asrc(k,n) = dataSrc(i,j)
     enddo
     enddo
  enddo
  if (timing>0) call shr_timer_stop (t2)

  !--- map the data ------------------------------------------------------------
  if (timing>0) call shr_timer_start(t3)
  call dshr_map_mapDataA(Asrc,Adst,map)
  if (timing>0) call shr_timer_stop (t3)

  !--- get data from temp/reshaped array ---------------------------------------
  if (timing>0) call shr_timer_start(t4)
  do k = 1,nfd
     !--- get field index ---
     if (sameList) then
        kDst = k
     else
        call shr_string_listGetName(lfDst,k,namDst)
        call shr_string_listGetIndex(bunDstf,namDst,kDst)
        if (kDst == 0) then
           write(6,F00) 'ERROR in dst field for ',trim(namDst)
           lrc = 1
           call dshr_map_abort(subName//'ERROR in name index')
        endif
     endif
     !--- get field data via pointer ---
     call dshr_bundle_assignPtr(bunDst,kDst,dataDst)
     !--- copy Adst into dataDst ---
     n = 0
     do j=1,njd
     do i=1,nid
        n = n + 1
        dataDst(i,j) = Adst(k,n)
     enddo
     enddo
  enddo
  deallocate(Asrc,Adst)
  if (timing>0) call shr_timer_stop (t4)

  !--- copy date ---
  call dshr_bundle_getDate(bunSrc,cDate,sec)
  call dshr_bundle_putDate(bunDst,cDate,sec)

  if (timing>0) call shr_timer_stop (t0)

end subroutine dshr_map_mapDataB

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_mapDataA -- map data from array to array
!
! !DESCRIPTION:
!     One of the methods associated with dshr\_map\_mapData.
!     Map arrays from source array to destination array using preset
!     map datatype.  Layer to shr\_map\_mapData.
!     \newline
!     call dshr\_map\_mapData(Ain,Aout,map)
!
! !REVISION HISTORY:
!     2005-Mar-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_mapDataA(arrsrc,arrdst,map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  real(R8)                ,intent(in) :: arrsrc(:,:)   ! (field,index)
  real(R8)                ,intent(out):: arrdst(:,:)   ! (field,index)
  type(dshr_map_mapType)  ,intent(in) :: map

!EOP

  !--- local ---
  integer(IN) :: nssrc,nsdst   ! number of grid cells, src and dst
  integer(IN) :: nfsrc,nfdst   ! number of fields, src and dst

  !--- formats ---
  character(*),parameter :: subName =   '(dshr_map_mapDataA) '
  character(*),parameter :: F01     = "('(dshr_map_mapDataA) ',a,2i8,a,2i8) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  nfsrc = size(arrsrc,1)
  nssrc = size(arrsrc,2)
  nfdst = size(arrdst,1)
  nsdst = size(arrdst,2)

! if (debug > 1) write(6,F01) ' sizes src:',nfsrc,nssrc,' dst:',nfdst,nsdst

  call shr_map_mapData(arrsrc,arrdst,map%sMat)
  
end subroutine dshr_map_mapDataA

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_clean -- cleans map
!
! !DESCRIPTION:
!     Clean map type, reset and deallocate
!     \newline
!     call dshr\_map\_clean(map)
!
! !REVISION HISTORY:
!     2005-Mar-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_clean(map)
  ! clean, zero, deallocate map datatype

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType)      ,intent(inout):: map

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_clean )'

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  nullify(map%src)
  nullify(map%dst)
  call dshr_map_put(map,shr_map_fs_name,'')
  call shr_map_clean(map%sMat)

end subroutine dshr_map_clean

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_print -- print info about map
!
! !DESCRIPTION:
!     print info about map
!     Layer to shr\_map\_listValidOpts
!     \newline
!     call dshr\_map\_print(map)
!
! !REVISION HISTORY:
!     2005-Apr-3 - T. Craig -- first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_print(map)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(dshr_map_mapType),intent(in) :: map

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = '(dshr_map_print) '
  character(*),parameter :: F00   = "('(dshr_map_print) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  write(6,*) ' '
  write(6,F00) ':'//trim(map%name)
  call shr_map_print(map%sMat)
  call shr_sys_flush(6)

end subroutine dshr_map_print

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_setDebug -- Set local dshr_map debug level
!
! !DESCRIPTION:
!     Set local dshr_map debug level, 0 = production
!     call shr\_map\_setDebug with same argument
!     \newline
!     call dshr\_map\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_setDebug(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: iflag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName =   '(dshr_map_setDebug) '
  character(*),parameter :: F00     = "('(dshr_map_setDebug) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  debug = iflag
  call shr_map_setDebug(iflag)

end subroutine dshr_map_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_setTiming -- Set local dshr_map timing level
!
! !DESCRIPTION:
!     Set local dshr_map timing level, 0 = production
!     \newline
!     call dshr\_map\_setTiming(1)
!
! !REVISION HISTORY:
!     2007-Mar-12  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_setTiming(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: iflag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName =   '(dshr_map_setTiming) '
  character(*),parameter :: F00     = "('(dshr_map_setTiming) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  timing = iflag

end subroutine dshr_map_setTiming

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_map_setAbort -- Set local dshr_map abort flag
!
! !DESCRIPTION:
!     Set local dshr_map abort flag, true = abort, false = print and continue
!     call shr\_map\_setAbort with same argument
!     \newline
!     call dshr\_map\_setAbort(.false.)
!
! !REVISION HISTORY:
!     2005-Mar-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_setAbort(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical,intent(in) :: flag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName =   '(dshr_map_setAbort) '
  character(*),parameter :: F00     = "('(dshr_map_setAbort) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  doabort = flag
  call shr_map_setAbort(flag)

end subroutine dshr_map_setAbort

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: dshr_map_abort -- local interface for abort
!
! !DESCRIPTION:
!     Local interface for dshr\_map abort calls
!     \newline
!     call dshr\_map\_sbort(subName//' ERROR illegal option')
!
! !REVISION HISTORY:
!     2005-Mar-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_map_abort(string)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),optional,intent(in) :: string

!XXEOP

  !--- local ---

  !--- formats ---
  character(CL) :: lstring
  character(*),parameter :: subName =   '(dshr_map_abort) '
  character(*),parameter :: F00     = "('(dshr_map_abort) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  lstring = ''
  if (present(string)) lstring = string

  if (doabort) then
    call shr_sys_abort(lstring)
  else
    write(6,F00) trim(lstring)
  endif

end subroutine dshr_map_abort

!===============================================================================
!===============================================================================

end module dshr_map
