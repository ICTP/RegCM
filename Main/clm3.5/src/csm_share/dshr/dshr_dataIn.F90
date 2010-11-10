!===============================================================================
! SVN $Id: dshr_dataIn.F90 2098 2006-10-11 19:11:41Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_dataIn.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_dataIn - reads upper and lower-bound data for a given stream
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!     2005-May-06 - B. Kauffman - first functional version
!     2004-Dec-15 - B. Kauffman - first proto-type version
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_dataIn 

   use shr_sys_mod    ! shared system calls
!  use shr_date_mod   ! shared date type and methods
   use shr_cal_mod    ! shared calendar type and methods
   use shr_stream_mod ! stream type and methods
   use shr_ncread_mod ! shared netCDF file reading module
   use shr_file_mod   ! file get/put routines

   use dshr_kind      ! kinds for strong typing
   use dshr_const     ! constants
   use dshr_bundle    ! bundle type and methods

   implicit none

   private ! default private

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: dshr_dataIn_readLBUB  ! if necessary, read in new LB & UB data
   public :: dshr_dataIn_setDebug  ! set internal dshr_dataIn debug level
   public :: dshr_dataIn_getDebug  ! get internal dshr_dataIn debug level

! !PUBLIC DATA MEMBERS:

  ! none

!EOP

   real(R8)   ,parameter :: spd = DSHR_CONST_CDAY ! seconds per day
   integer(IN),save      :: debug = 0             ! debug level

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_dataIn_readLBUB - update upper and lower bounding data
!
! !DESCRIPTION:
!    Update upper and lower bounding data, if necessary.
!
! !REVISION HISTORY:
!     2006-Oct-11 - B. Kauffman - add copy old UB to new LB functionality
!     2005-Nov-21 - B. Kauffman - add pre-fetch & remove old functionality
!     2005-May-06 - B. Kauffman - first functional version
!     2004-Dec-15 - B. Kauffman - first proto-typ version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_dataIn_readLBUB(stream,mDate,mSec,bunLB,bunUB,newData,rmOldFile,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_stream_streamType) ,intent(inout) :: stream      ! defines data file stream
   integer(IN)                 ,intent(in)    :: mDate,mSec  ! yymmdd & sec of model
   type(dshr_bundle_bundleType),intent(inout) :: bunLB,bunUB ! LB & UB data bundles
   logical                     ,intent(in)    :: rmOldFile   ! rm old/unneeded file
   logical                     ,intent(out)   :: newData     ! T => new data in bundles
   integer(IN),optional        ,intent(out)   :: rc          ! return code

!EOP

   integer(IN)       :: rCode                   ! local rc
   integer(IN)       :: mDateLB,dDateLB,mSecLB  ! yymmdd & sec of model,LB,UB data
   integer(IN)       :: mDateUB,dDateUB,mSecUB  ! yymmdd & sec of model,LB,UB data
   real(R8)          :: rDateM,rDateLB,rDateUB  ! model,LB,UB dates with fractional days
   character(CL)     :: fileName                ! file name (with file path prepended?)
   character(CL)     :: fileLB,fileUB           ! UB & LB files to open
   character(CL)     :: fn_lb,fn_ub             ! UB & LB data files to read
   integer(IN)       ::  n_lb, n_ub             ! UB & LB time index in data files

   logical           :: newLBeqOldUB            ! true iff (new UB == old LB)
   integer(IN)       :: oDateUB,oSecUB          ! UB date *before* reading data ("old" UB)

   integer(IN)       :: k,kFlds                 ! index to fields, number of fields
   character(CL)     :: sFldName                ! field name in data file
   character(CL)     :: bFldName                ! field name in data bundle
   real(R8),pointer  :: dataLB(:,:),dataUB(:,:) ! pointers to data in LB,UB bundles
   character(CL)     :: str                     ! temporary string var
   character(CL)     :: fn_prev,fn_next,path    ! file name: prev,next,path
   logical           :: fileExists              ! a file exists in directory?
   logical           :: newFile                 ! just opened a new file in list?

   character(*),parameter :: subName =  "(dshr_dataIn_readLBUB) "
   character(*),parameter :: F00 =    "('(dshr_dataIn_readLBUB) ',8a)"
   character(*),parameter :: F01 =    "('(dshr_dataIn_readLBUB) ',a,3(i8.8,1x,i5.5,'s  '))"
   character(*),parameter :: F02 =    "('(dshr_dataIn_readLBUB) ',a,3f12.2)"

!-------------------------------------------------------------------------------
! PURPOSE: 
!    1) determine LB & UB input data that brackets the current model time
!    2) update the LB & UB data bundles, if necessary
! Notes:
! o  The newFile logic assumes that when the a stream advances into the next
!    file in the time series, the UB will be in that new file while the LB 
!    remains in the older file.  If this logic fails, the only consequence is 
!    that a new-file message will not be printed and the next file in the 
!    sequence will not be pre-fetched.
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! determine LB & UB of input data that bracket the current model time
   !----------------------------------------------------------------------------

   newData = .false.  ! T => new data samples
   newFile = .false.  ! T => new data files
   rCode = 0

   newLBeqOldUB = .false.
   if (dshr_bundle_checkDate(bunLB) .and. dshr_bundle_checkDate(bunUB)) then
      call dshr_bundle_getDate(bunLB,mDateLB,mSecLB,rCode)
      call dshr_bundle_getDate(bunUB,mDateUB,mSecUB,rCode)
      rDateM  = real(mDate  ,R8) + real(mSec  ,R8)/spd
      rDateLB = real(mDateLB,R8) + real(mSecLB,R8)/spd
      rDateUB = real(mDateUB,R8) + real(mSecUB,R8)/spd
      oDateUB = mDateUB ! UB model date *before* reading new data ("old" UB)
       oSecUB =  mSecUB 
      if ( rDateLB <= rDateM .and. rDateM <= rDateUB ) then
         newData = .false. ! existing LB/UB data already brackets the model date
      else
         newData = .true.  ! existing LB/UB data does not bracket the model date
      end if
      if (debug > 0) then
         write(6,F01) 'Dates LB,model,UB: ',mDateLB,mSecLB,mDate,mSec,mDateUB,mSecUB
         if (      newData) write(6,F00) "reading in new LB/UB data is required"
         if (.not. newData) write(6,F00) "reading in new LB/UB data is not required"
      endif
   else
      !--- model just starting => no data read in yet => bundle dates not set ---
      newData = .true.  ! must read in some new (initial) data
      newFile = .true.  ! must open new (initial) data files
      if (debug > 0) then
         write(6,F00) 'LB/UB bundle dates have not yet been set (model starting up)'
         write(6,F00) "must read in initial LB/UB data"
      endif
      oDateUB = 0 ! UB model date *before* reading new data (undefined in this case)
       oSecUB = 0
   endif
   
   !----------------------------------------------------------------------------
   ! if new LB & UB data is needed, update the data bundles
   !----------------------------------------------------------------------------
   if (newData) then
      !--- determine what the bounding data is and where it is located ---
      call shr_stream_findBounds(stream,mDate,mSec,                 &
      &                          mDateLB,dDateLB,mSecLB,n_lb,fn_lb, &
      &                          mDateUB,dDateUB,mSecUB,n_ub,fn_ub  )
      
      newLBeqOldUB = (mDateLB == oDateUB .and. mSecLB == oSecUB) ! newLB =? oldUB

      !--- document input data files ---
      if (newFile) then
         !--- model just started up ---
         write(6,F00) "opened initial LB file: ",trim(fn_lb)
         write(6,F00) "opened initial UB file: ",trim(fn_ub)
      else if (trim(fn_lb) /= trim(fn_ub)) then
         !--- just opened next file in sequence ---
         !--- this logic assumes LB & UB are not in same file ---
         newFile = .true.
         write(6,F00) "opened new     UB file: ",trim(fn_ub)
      end if
      if (debug > 1) then 
         write(6,F00) "LB file = ",trim(fn_lb)
         write(6,F00) "UB file = ",trim(fn_ub)
      end if

      !--- acquire input data files ---
      call shr_stream_getFilePath(stream,path)
      inquire(file=trim(fn_lb),exist=fileExists)
      if (.not. fileExists) call shr_file_get(rCode,fn_lb,trim(path)//fn_lb)
      inquire(file=trim(fn_ub),exist=fileExists)
      if (.not. fileExists) call shr_file_get(rCode,fn_ub,trim(path)//fn_ub)

      !--- for each field, read data from file and put into bundle ---
      call dshr_bundle_getNumFields(bunLB,kFlds)
      do k=1,kFlds
         !--- get field name in data-stream & in bundle ---
         call shr_stream_getFileFieldName (stream,k,sFldName)
         call shr_stream_getModelFieldName(stream,k,bFldName)

         !--- debug info ---
         if (debug > 1) then
            write(6,F00) "read field: src, dest name = ",trim(sFldName),' ',trim(bFldName)
         end if

         !--- assign pointer to data array in bundle ---
         call dshr_bundle_assignPtr(bunLB,bFldName,dataLB)
         call dshr_bundle_assignPtr(bunUB,bFldName,dataUB)

         !--- read data directly into data array in bundle ---
         if (newLBeqOldUB) then !--- new UB == old LB ---
            if (debug > 1 .and. k==1) write(6,F00) "copy UB to LB, don't read ",trim(fn_ub)
            call dshr_bundle_getField(bunUB,dataLB,bFldName,rc=rCode)
         else
            call shr_ncread_tField(fn_lb, n_lb, sfldName, dataLB, rc=rCode)
            if (rCode /= 0) then
               write(str,'(a4)') "ERROR: reading field ",trim(sFldName)," from file ",trim(fn_lb)
               write(6,F00) trim(str)
               call shr_sys_abort(subName//trim(str))
            end if
         end if
         call shr_ncread_tField(fn_ub, n_ub, sfldName, dataUB, rc=rCode)
         if (rCode /= 0) then
            write(str,'(a4)') "ERROR: reading field ",trim(sFldName)," from file ",trim(fn_ub)
            write(6,F00) trim(str)
            call shr_sys_abort(subName//trim(str))
         end if
      end do

      !--- assign the appropriate model date to the new LB & UB bundles ---
      call dshr_bundle_putDate(bunLB,mDateLB,mSecLB)
      call dshr_bundle_putDate(bunUB,mDateUB,mSecUB)

      !--- determine previous & next data files in list of files ---
      call shr_stream_getPrevFileName(stream,fn_lb,fn_prev,path)
      call shr_stream_getNextFileName(stream,fn_ub,fn_next,path)

      !--- pre-fetch the next file ---
      if ( newFile ) then ! just opened a new file 
         inquire(file=trim(fn_next),exist=fileExists)
         if ( trim(fn_next) == "unknown" ) then
            write(6,F00) "not pre-fetching - couldn't determine next file"
         else if ( fileExists ) then
            write(6,F00) "exists/don't pre-fetch: ",trim(fn_next)
         else 
            call shr_file_get(rCode,fn_next,trim(path)//fn_next,async=.true.)
            write(6,F00) "pre-fetch         file: ",trim(fn_next)
         end if
      end if

      !---  remove the old file? ---
      if ( rmOldFile .and. newFile ) then ! just opened a new file
         if ( fn_prev/=fn_lb .and. fn_prev/=fn_ub .and. fn_prev /= fn_next ) then 
            !--- previous file is not in use and is not next in list ---
            inquire(file=trim(fn_prev),exist=fileExists)
            if ( fileExists ) then
               call shr_sys_system(" rm "//trim(fn_prev),rCode)
               write(6,F00) "remove previous   file: ",trim(fn_prev)
            end if
         end if
      end if

   end if

   call shr_sys_flush(6)

   if (present(rc)) rc=rCode

end subroutine dshr_dataIn_readLBUB

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_dataIn_setDebug -- Set local debug level
!
! !DESCRIPTION:
!    Set local debug level, 0 = production
!    \newline
!    General Usage: call dshr\_dataIn\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_dataIn_setDebug(level)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: level

!EOP

  !--- formats ---
  character(*),parameter :: subName =  "('dshr_dataIn_setDebug') "
  character(*),parameter :: F00     = "('(dshr_dataIn_setDebug) ',a) "
  character(*),parameter :: F01     = "('(dshr_dataIn_setDebug) ',a,i4) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  debug = level
  write(6,F01) "debug level reset to ",level

end subroutine dshr_dataIn_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_dataIn_getDebug -- get dshr_dataIn internal debug level
!
! !DESCRIPTION:
!    Get local debug level, 0 = production
!    \newline
!    General Usage: call dshr\_dataIn\_getDebug(level)
!
! !REVISION HISTORY:
!     2005-May-10  - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_dataIn_getDebug(level)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(out) :: level

!EOP

  !--- formats ---
  character(*),parameter :: subName =  "('dshr_dataIn_getDebug') "
  character(*),parameter :: F00     = "('(dshr_dataIn_getDebug) ',a) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  level = debug

end subroutine dshr_dataIn_getDebug

!===============================================================================

end module dshr_dataIn
