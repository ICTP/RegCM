!===============================================================================
! SVN $Id: dshr_rest.F90 2413 2006-11-08 22:54:08Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_rest.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_rest -- read and write restart files.
!
! !DESCRIPTION:
!    Read and write restart files.
!     
! !REMARKS:
!     
! !REVISION HISTORY:
!    2006-Feb-13 - B. Kauffman - refactored to support branch functionality
!    2005-Mar-02 - B. Kauffman - first place-holder module
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_rest

! !USES:

   use shr_sys_mod     ! system routines
   use shr_file_mod    ! file get & put routines
   use shr_timer_mod   ! timer routines
   use shr_date_mod    ! date data-type and methods
   use shr_string_mod  ! string/list operators
   use shr_stream_mod  ! stream data-type and methods

   use dshr_kind       ! kinds for strong typing
   use dshr_bundle     ! bundle data-type and methods
   use dshr_iocdf      ! cdf routines

   implicit none

   private             ! everything is default private

! !PUBLIC TYPES:

   ! none
 
! !PUBLIC MEMBER FUNCTIONS:

   public :: dshr_rest_write             ! create restart files
   public :: dshr_rest_readPointer       ! read pointer file
   public :: dshr_rest_readBundle        ! read bundle  data file
   public :: dshr_rest_readBundleFldList ! read bundle  field list
   public :: dshr_rest_readStream        ! read stream  data file
   public :: dshr_rest_setDebug          ! set internal debug level
   public :: dshr_rest_getDebug          ! get internal debug level

! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   integer(IN)  ,save      :: debug = 0          ! local debug level

!===============================================================================
contains
!===============================================================================

!BOP ===========================================================================
!
! !IROUTINE: dshr_rest_write -- create a restart files
!
! !DESCRIPTION:
!    Create bundle and stream restart files and an associated restart pointer file.
!
! !REVISION HISTORY:
!    2005-Jul-07 - B. Kauffman - creates file names and pointer file
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_rest_write(bundle,stream,date,modelName,caseName,caseDesc,ptrFn,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType) ,pointer     :: stream(:) ! streams to save
   type(dshr_bundle_bundleType),intent(in)  :: bundle    ! bundle to save
   type(shr_date)              ,intent(in)  :: date      ! model date
   character(*)                ,intent(in)  :: modelName ! model name 
   character(*)                ,intent(in)  :: caseName  ! case name
   character(*)                ,intent(in)  :: caseDesc  ! case description
   character(*)                ,intent(in)  :: ptrFn     ! pointer file name
   integer(IN),optional        ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   character(CL)    :: bundleFn           ! bundle data file name
   character(CL)    :: streamFn           ! stream data file name    
   character(CL)    :: locPtrFn           ! ptrFn with device & path removed
   integer(IN)      :: n                  ! index into char string
   integer(IN)      :: year,month,day,sec ! model date
   character(CS)    :: dateStr            ! date string part of file name
   character(CX)    :: fldList            ! list of fields in bundle
   character(CL)    :: attName            ! cdf attribute containing field list
   character( 8)    :: dStr               ! F90 wall clock date str yyyymmdd                          
   character(10)    :: tStr               ! F90 wall clock time str hhmmss.sss 
   integer(IN)      :: nUnit              ! a file unit number
   integer(IN)      :: fid                ! cdf file id
   integer(IN)      :: rCode              ! return code

   integer(IN),save :: timer = -1         ! timer

   !----- formats -----
   character(*),parameter :: subName = '(dshr_rest_write) '
   character(*),parameter :: F00   = "('(dshr_rest_write) ',4a)" 

!-------------------------------------------------------------------------------
! NOTES
! - local vs. remote restart pointer file names:
!   The full/remote ptrFn is specified by [device:][path/]fileName
!   A local ptr file (in the cwd) with the same fileName is first created 
!   and it is copied to the specified device or location.  
!   See routine shr_file_put()
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   if (timer == -1) call shr_timer_get(timer,"dshr_rest_write")
   call shr_timer_start(timer)

   !--- get date associated with the data ---
   call shr_date_getYMD  (date,year,month,day,sec)

   !----------------------------------------------------------------------------
   ! construct the restart file names
   !----------------------------------------------------------------------------
   dateStr = "123456789+123456"
   dateStr = "yyyy-mm-dd-sssss"
   write(dateStr( 1: 4),'(i4.4)') year
   write(dateStr( 6: 7),'(i2.2)') month
   write(dateStr( 9:10),'(i2.2)') day
   write(dateStr(12:16),'(i5.5)') sec
   streamFn = trim(caseName)//"."//trim(modelName)//".rs."//trim(dateStr)//".bin"
   bundleFn = trim(caseName)//"."//trim(modelName)//".rb."//trim(dateStr)//".nc"

   !---------------------------------------------------------------------------- 
   ! create the restart files
   !---------------------------------------------------------------------------- 
   write(6,F00) 'creating restart bundle  file: ',trim(bundleFn)
   fldList = "null - no bundle fields"
   if (dshr_bundle_checkInit(bundle)) call dshr_bundle_getFieldList(bundle,fldList)
   attName = trim(modelName)//"_restart_fieldList"
   call dshr_iocdf_create(bundleFn,caseDesc)
   call dshr_iocdf_appendAtt(fid,attName,fldList,bundleFn)
   if (fldList(1:5) /= "null ") call dshr_iocdf_append(fid,date,bundle,bundleFn)

   write(6,F00) 'creating restart stream  file: ',trim(streamFn)
   call shr_stream_restWrite(stream,streamFn,caseName,caseDesc)

   !---------------------------------------------------------------------------- 
   ! update the restart pointer file                                      
   !----------------------------------------------------------------------------
   write(6,F00) 'updating restart pointer file: ',trim(ptrFn)                           
   call shr_sys_flush(6)                                                        

   !--- create a local pointer file ---                                         
   n = max(index(ptrFn,":"),index(ptrFn,"/",.true.)) ! rm device & path prefix
   locPtrFn = trim(ptrFn(n+1:len_trim(ptrFn)))       ! from ptrFn
   call date_and_time(dStr,tStr)
   nUnit = shr_sys_ioUnit()
   open (nUnit,file=trim(locPtrFn),form="FORMATTED",status="REPLACE")               
   write(nUnit,'( a)') trim(bundleFn)                                            
   write(nUnit,'( a)') trim(streamFn)                                            
   write(nUnit,'(2a)') "case name           : ",trim(caseName)                                            
   write(nUnit,'(2a)') "case description    : ",trim(caseDesc)                                            
   write(nUnit,'(2a)') 'Pointer file created: ',              &
      &       dStr(1:4)//'-'//dStr(5:6)//'-'//dStr(7:8)//' '  &
      &     //tStr(1:2)//':'//tStr(3:4)//':'//tStr(5:6)                         
   close(nUnit)

   !--- copy pointer file to specfied name & location ---
   call shr_file_put(rCode,trim(locPtrFn),trim(ptrFn))

   !---------------------------------------------------------------------------- 
   ! done
   !---------------------------------------------------------------------------- 
   call shr_timer_stop(timer)

end subroutine dshr_rest_write

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_rest_readPointer -- read restart pointer file
!
! !DESCRIPTION:
!    Read a restart pointer file
!
! !REVISION HISTORY:
!    2005-Nov-21 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_rest_readPointer(ptrFn,bundleFn,streamFn,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)           :: ptrFn    ! pointer file to read
   character(*),intent(out)          :: bundleFn ! bundle restart file name
   character(*),intent(out)          :: streamFn ! stream restart file name
   integer(IN) ,intent(out),optional :: rc       ! return code

!EOP

   !----- local -----
   character(CL)    :: caseName           ! case name
   character(CL)    :: caseDesc           ! case description
   character(CL)    :: timeStamp          ! pointer file time stamp
   integer(IN)      :: nUnit              ! a file unit number
   logical          :: open               ! true if file unit is open 
   character(CL)    :: locPtrFn           ! local pointer file
   integer(IN)      :: n                  ! index into char string
   integer(IN)      :: rCode              ! return code

   !----- formats -----
   character(*),parameter :: subName = '(dshr_rest_readPointer) '
   character(*),parameter :: F00   = "('(dshr_rest_readPointer) ',4a)" 

!-------------------------------------------------------------------------------
! NOTES
! - local vs. remote restart pointer file names:
!   The full/remote ptrFn is specified by [device:][path/]fileName
!   First a local copy of the ptr file (with the same file name) is acquired
!   and then this local copy is read.
!   See routine shr_file_get()
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   !------------------------------------------------------------
   ! aquire and read the restart pointer file
   !------------------------------------------------------------

   !--- aquire file ---
   n = max(index(ptrFn,":"),index(ptrFn,"/",.true.)) ! rm device & path prefix
   locPtrFn = trim(ptrFn(n+1:len_trim(ptrFn)))       ! from ptrFn
   call shr_file_get(rCode,trim(locPtrFn),trim(ptrFn),clobber=.true.)
                                                                                
   !--- read the pointer file ---
   nUnit = shr_sys_ioUnit() ! get an unused unit number
   open (nUnit,file=trim(locPtrFn),form="FORMATTED",status="OLD",action="READ")
   read (nUnit,'(a)') bundleFn                                            
   read (nUnit,'(a)') streamFn
   read (nUnit,'(a)') caseName
   read (nUnit,'(a)') caseDesc
   read (nUnit,'(a)') timeStamp
   close(nUnit)

   write(6,F00)' rpointer file name    : ',trim(ptrFn)
   write(6,F00)' * bundle restart file : ',trim(bundleFn)
   write(6,F00)' * stream restart file : ',trim(streamFn)
   write(6,F00)' * ',trim(caseName)
   write(6,F00)' * ',trim(caseDesc)
   write(6,F00)' * ',trim(timeStamp)
   call shr_sys_flush(6)                                                        

end subroutine dshr_rest_readPointer

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_rest_readBundle -- read restart bundle file
!
! !DESCRIPTION:
!    Read bundle restart file
!
! !REVISION HISTORY:
!    2005-Jul-07 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_rest_readBundle(bundleFn,bundle,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)                ,intent(in)    :: bundleFn ! restart pointer file name
   type(dshr_bundle_bundleType),intent(inout) :: bundle   ! bundle to read
   integer(IN),optional        ,intent(out)   :: rc       ! return code

!EOP

   !----- local -----
   integer(IN)   :: rCode          ! return code

   !----- formats -----
   character(*),parameter :: subName = '(dshr_rest_readBundle) '
   character(*),parameter :: F00   = "('(dshr_rest_readBundle) ',4a)" 

!-------------------------------------------------------------------------------
! NOTES
! o bundle file must have already been initialized/allocated
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   call dshr_iocdf_read(bundleFn,bundle)

end subroutine dshr_rest_readBundle

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_rest_readBundleFldList -- read field list from a bundle file
!
! !DESCRIPTION:
!    Read a restart field list from a bundle restart file
!
! !REVISION HISTORY:
!    2005-Jul-07 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_rest_readBundleFldList(bundleFn,modelName,fldList,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)  :: bundleFn  ! bundle data file name
   character(*)        ,intent(in)  :: modelName ! model name
   character(*)        ,intent(out) :: fldList   ! restart field List
   integer(IN),optional,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   character(CS) :: attName  ! cdf attribute name
   integer(IN)   :: rCode    ! return code

   !----- formats -----
   character(*),parameter :: subName = '(dshr_rest_readBundleFldList) '
   character(*),parameter :: F00   = "('(dshr_rest_readBundleFldList) ',4a)" 

!-------------------------------------------------------------------------------
! NOTES
! o bundle file must have already been initialized/allocated
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   attName = trim(modelName)//"_restart_fieldList"
   call dshr_iocdf_readAtt(bundleFn,attName,fldList)

end subroutine dshr_rest_readBundleFldList

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_rest_readStream -- read stream restart file
!
! !DESCRIPTION:
!    Read stream restart file
!
! !REVISION HISTORY:
!    2005-Jul-07 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_rest_readStream(streamFn,stream,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)               ,intent(in)  :: streamFn  ! stream file name
   type(shr_stream_streamType),pointer     :: stream(:) ! vector of streams to read
   integer(IN),optional       ,intent(out) :: rc        ! return code

!EOP

   !----- local -----
   integer(IN)      :: rCode              ! return code

   !----- formats -----
   character(*),parameter :: subName = '(dshr_rest_readStream) '
   character(*),parameter :: F00   = "('(dshr_rest_readStream) ',4a)" 

!-------------------------------------------------------------------------------
! NOTES
! o stream vector must *not* have been initialized/allocated
!   there's no reasonable way (?) to do this outside of shr_stream_restRead
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   call shr_stream_restRead(stream,streamFn)

end subroutine dshr_rest_readStream

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_rest_setDebug -- Set internal debug level
!
! !DESCRIPTION:
!    Set internal debug level, 0 = production
!    \newline
!    General Usage: call dshr\_rest\_setDebug(2)
!
! !REVISION HISTORY:
!    2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_rest_setDebug(level)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer,intent(in) :: level

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(dshr_rest_setDebug) '
   character(*),parameter :: F01   = "('(dshr_rest_setDebug) ',a,i4) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   debug = level
   write(6,F01) "debug level reset to ",level

end subroutine dshr_rest_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_rest_getDebug -- return internal debug level
!
! !DESCRIPTION:
!    Return internal debug level, 0 = production
!    \newline
!    General Usage: call dshr\_rest\_getDebug(level)
!
! !REVISION HISTORY:
!    2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_rest_getDebug(level)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer,intent(out) :: level

!EOP

   !--- formats ---
   character(*),parameter :: subName = '(dshr_rest_getDebug) '
   character(*),parameter :: F00   = "('(dshr_rest_getDebug) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   level = debug

end subroutine dshr_rest_getDebug

!===============================================================================
!===============================================================================

end module dshr_rest
