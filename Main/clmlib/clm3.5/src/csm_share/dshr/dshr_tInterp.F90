!===============================================================================
! SVN $Id: dshr_tInterp.F90 3322 2007-02-28 00:39:45Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_tInterp.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_tInterp -- time interpolation routines
!
! !DESCRIPTION:
!     dshr shared time interpolation routines
!
! !REVISION HISTORY:
!     2004-Dec-10 - J. Schramm - first version
!     2005-Apr-10 - T. Craig - updated for dshr bundles
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_tInterp
 
! !USES:

   use shr_sys_mod   ! shared system calls
   use shr_cal_mod   ! shared calendar type and methods
   use shr_tInterp_mod ! shared tInterp stuff
   use shr_string_mod  ! strings
   use dshr_kind     ! kinds for strong typing
   use dshr_const    ! constants
   use dshr_bundle   ! bundles

   implicit none

   private ! except

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: dshr_tInterp_data
   public :: dshr_tInterp_setAbort
   public :: dshr_tInterp_setDebug

! !PUBLIC DATA MEMBERS:

   ! no public data

!EOP

   real(R8),parameter :: c0 = dshr_const_c0
   real(R8),parameter :: c1 = dshr_const_c1
   real(R8),parameter :: eps = 1.0E-12_R8
   logical ,save      :: doabort = .true.
   integer ,save      :: debug = 0

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_tInterp_data -- time interpolation of bundles
!
! !DESCRIPTION:
!     Returns time interpolated bundle from two bundle inputs
!     Legal algorithms are (algo):
!       lower   - sets bun to lb data
!       upper   - sets bun to ub data
!       nearest - sets bun to nearest data in time
!       linear  - sets bun to linear interpolation between lb and ub
!     \newline
!     call dshr\_tInterp\_data(bunlb,bunub,bun,algo='linear')
!     \newline
!     call dshr\_tInterp\_data(bunlb,bunub,bun,fsrc='a:b:c',fdst='d:e:f',algo='linear')
!     \newline
!     field dims of all bundles must be same (x,y)
!     field lists in bunlb and bunub are assumed to be be same
!     matches bun field names with bunlb/bunub field names
!     optional field list allows different field names to be related
!     time of ub >= time of lb  for all algos
!     time of ub >= time of model data >= time of lb  for linear
!     if fsrc is present then fdst must be present, must be same number of fields
!
! !REVISION HISTORY:
!     2005-Apr-10 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_tInterp_data(bunlb,bunub,bun,fsrc,fdst,algo,rc)

   implicit none

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(dshr_bundle_bundleType),intent(in)    :: bunlb    ! lower bound
   type(dshr_bundle_bundleType),intent(in)    :: bunub    ! upper bound
   type(dshr_bundle_bundleType),intent(inout) :: bun      ! output bundle
   character(*),optional       ,intent(in)    :: fsrc     ! input field list
   character(*),optional       ,intent(in)    :: fdst     ! output field list
   character(*),optional       ,intent(in)    :: algo     ! algorithm
   integer(IN) ,optional       ,intent(out)   :: rc       ! return code

!EOP

   !----- local  ------
   integer(IN) :: n                      ! counter
   integer(IN) :: k,ks,kd                ! field index: src & dest
   character(CS) :: fldName              ! field name
   integer(IN) :: nfs,nfd                ! number of fields
   integer(IN) :: nilb,njlb,nflb         ! ni,nj,nf for lower bound
   integer(IN) :: niub,njub,nfub         ! ni,nj,nf for upper bound
   integer(IN) :: ni,nj,nf               ! ni,nj,nf for model date
   integer(IN) :: cDatelb,seclb          ! date timestamp in yyyymmdd, seconds
   integer(IN) :: cDateub,secub          ! date timestamp in yyyymmdd, seconds
   integer(IN) :: cDate  ,sec            ! date timestamp in yyyymmdd, seconds
   character(CX) :: bunlbf,bunubf,bunf   ! field name lists
   character(CX) :: lfsrc,lfdst          ! field name lists
   real(R8),pointer :: Alb(:,:)          ! field data for lower bound
   real(R8),pointer :: Aub(:,:)          ! field data for upper bound
   real(R8),pointer :: A(:,:)            ! field data for model date
   real(R8)    :: flb,fub                ! factor for lb and ub
   character(CS) :: lalgo                ! local algo variable
   integer(IN) :: lrc                    ! local rc

   !----- formats -----
   character(*),parameter :: subName = "('dshr_tInterp_data')"
   character(*),parameter :: F00 = "('(dshr_tInterp_data) ',8a)" 
   character(*),parameter :: F01 = "('(dshr_tInterp_data) ',a,2f17.8)" 
   character(*),parameter :: F02 = "('(dshr_tInterp_data) ',a,3f17.8)" 
   character(*),parameter :: F03 = "('(dshr_tInterp_data) ',a,3i8)" 

!-------------------------------------------------------------------------------
! Computes time interpolation factors
!-------------------------------------------------------------------------------

   lrc = 0

   if (present(algo)) then
     lalgo = trim(algo)
   else
     lalgo = 'linear'
   endif

   !--- extract dates from bundles ---
   if (debug>0)  write(6,F00) 'extract dates from bundles'
   call dshr_bundle_getDate(bunlb,cDatelb,seclb)
   call dshr_bundle_getDate(bunub,cDateub,secub)
   call dshr_bundle_getDate(bun,cDate,sec)

   if (debug>0)  write(6,F00) 'call shr_tInterp_getFactors'
   call shr_tInterp_getFactors(cDatelb,seclb,cDateub,secub,cDate,sec, &
     flb,fub,lalgo,lrc)

   !--- check that bundle sizes are mostly consistent ---
   if (debug>0)  write(6,F00) 'check that bundle sizes'
   call dshr_bundle_getDims(bunlb,nilb,njlb,nflb)
   call dshr_bundle_getDims(bunub,niub,njub,nfub)
   call dshr_bundle_getDims(bun  ,ni  ,nj  ,nf  )
   if (nilb /= ni .or. niub /= ni) then
     write(6,F03) ' ERROR i index size ',nilb,ni,niub
     lrc = 1
     call dshr_tInterp_abort(subName//' i index size')
   endif
   if (njlb /= nj .or. njub /= nj) then
     write(6,F03) ' ERROR j index size ',njlb,nj,njub
     lrc = 1
     call dshr_tInterp_abort(subName//' j index size')
   endif
   call dshr_bundle_getFieldList(bunLB,bunLBf)
   call dshr_bundle_getFieldList(bunUB,bunUBf)
   if (trim(bunLBf) /= trim(bunUBf) ) then
      write(6,F03) "ERROR: ub,lb field lists don't match "
      write(6,F03) "ERROR: lb fields: ",trim(bunLBf)
      write(6,F03) "ERROR: ub fields: ",trim(bunUBf)
      lrc = 1
      call dshr_tInterp_abort(subName//"ub,lb field lists don't match")
   endif
   if   (present(fdst) .and. .not. present(fsrc)) then
      write(6,F03) "ERROR: fdst is present and fsrc is not"
      call dshr_tInterp_abort(subName//"fdst is present and fsrc is not")
   endif
   if   (.not. present(fdst) .and. present(fsrc)) then
      write(6,F03) "ERROR: fdst is not present and fsrc is"
      call dshr_tInterp_abort(subName//"fdst is nnot present and fsrc is")
   endif

   if (lrc == 0 .and. .not. present(fdst)) then
      !-------------------------------------------------------------------------
      ! interpolate all fields with matching names
      !-------------------------------------------------------------------------
      do n = 1,nf
         call dshr_bundle_getFieldName (bun,n,fldName)
         call dshr_bundle_getFieldIndex(bunUB,fldName,ks)
         if (ks > 0) then
            !--- get field data via pointer ---
            call dshr_bundle_assignPtr(bunlb,ks,Alb)
            call dshr_bundle_assignPtr(bunub,ks,Aub)
            call dshr_bundle_assignPtr(bun  ,n ,A  )
!dir$ concurrent
!cdir nodep
            A = Alb * flb + Aub * fub           
         end if
      enddo
   else if (lrc == 0 ) then
      !-------------------------------------------------------------------------
      ! interpolate fields as specified in optional input field lists
      !-------------------------------------------------------------------------
      nf = shr_string_listGetNum(fdst)
      do n = 1,nf
         !--- find src field ---
         call shr_string_listGetName   (fsrc,n,fldName)
         call dshr_bundle_getFieldIndex(bunLB ,fldName,ks)
         if (ks<1) then
            write(6,F00) 'ERROR: missing src field = ',trim(fldName)
            call dshr_tInterp_abort(subName//"missing src field = "//trim(fldName))
         end if
         !--- find dst field ---
         call shr_string_listGetName   (fdst,n,fldName)
         call dshr_bundle_getFieldIndex(bun   ,fldName,kd)
         if (kd<1) then
            write(6,F00) 'ERROR: missing src field = ',trim(fldName)
            call dshr_tInterp_abort(subName//"missing src field = "//trim(fldName))
         end if
         !--- interp ---
         call dshr_bundle_assignPtr(bunLB,ks,Alb)
         call dshr_bundle_assignPtr(bunUB,ks,Aub)
         call dshr_bundle_assignPtr(bun  ,kd,A  )
!dir$ concurrent
!cdir nodep
            A = Alb * flb + Aub * fub           
      enddo
   endif

   !--- clean up pointers ---
   nullify(Alb)
   nullify(Aub)
   nullify(A)

   if (debug > 0) then
     write(6,F01) ' algo,flb,fub === '//trim(lalgo),flb,fub
   endif

   if (present(rc)) rc = lrc

end subroutine dshr_tInterp_data

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_tInterp_setAbort -- Set local dshr_tInterp abort flag
!
! !DESCRIPTION:
!     Set local dshr_tInterp abort flag, true = abort, false = print and continue
!     \newline
!     call dshr\_tInterp\_setAbort(.false.)
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_tInterp_setAbort(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical,intent(in) :: flag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('dshr_tInterp_setAbort') "
  character(*),parameter :: F00     = "('(dshr_tInterp_setAbort) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  doabort = flag
  call shr_tInterp_setAbort(doabort)

end subroutine dshr_tInterp_setAbort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_tInterp_setDebug -- Set local dshr_tInterp debug level
!
! !DESCRIPTION:
!     Set local dshr_tInterp debug level, 0 = production
!     \newline
!     call dshr\_tInterp\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_tInterp_setDebug(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: iflag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('dshr_tInterp_setDebug') "
  character(*),parameter :: F00     = "('(dshr_tInterp_setDebug) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  debug = iflag
  call shr_tInterp_setDebug(debug)

end subroutine dshr_tInterp_setDebug

!===============================================================================
!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: dshr_tInterp_abort -- local interface for abort
!
! !DESCRIPTION:
!     Local interface for dshr\_tInterp abort calls
!     \newline
!     call dshr\_tInterp\_abort(subName//' ERROR illegal option')
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_tInterp_abort(string)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),optional,intent(in) :: string

!XXEOP

  !--- local ---

  !--- formats ---
  character(CL) :: lstring
  character(*),parameter :: subName = "('dshr_tInterp_abort') "
  character(*),parameter :: F00     = "('(dshr_tInterp_abort) ',a) "

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

end subroutine dshr_tInterp_abort

!===============================================================================
!===============================================================================

end module dshr_tInterp

