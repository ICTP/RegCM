!===============================================================================
! SVN $Id: cpl_iocdf_mod.F90 3380 2007-03-06 05:42:19Z robj $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_iocdf_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_iocdf_mod -- create and write-to netcdf data file.
!
! !DESCRIPTION:
!    This module creates and writes data in a {\tt bundle} to a file
!    in machine-independent netCDF format.  This module is intended to support 
!    the history writing sub-system of cpl6.

!
! !REMARKS:
!    Typically this module would be used to support the history file subsystem
!    of cpl6, but in fact it could be used to create netCDF data files for any
!    purpose.  The netCDF data file format might provide an alternative to 
!    binary data files.
!
! !REVISION HISTORY:
!     2002-Oct-22 - B. Kauffman - major upgrade, refactoring
!     2001-Dec-20 - B. Kauffman - first prototype
!
! !INTERFACE:  -----------------------------------------------------------------

module cpl_iocdf_mod

! !USES:

   use mct_mod       ! mct interface
   use cpl_comm_mod      ! mpi/mph communicator info
   use cpl_fields_mod    ! coupler/model data field indicies
   use cpl_bundle_mod    ! defines bundle
   use cpl_domain_mod    ! defines domain
   use cpl_kind_mod      ! defines F90 kinds
   use cpl_const_mod     ! defines constants (eg. spval)
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug
   use shr_sys_mod       ! share system routines
   use shr_date_mod      ! defines date data-type

   implicit none

#include <netcdf.inc>

   private ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_iocdf_create   ! create a new file (an empty file)
   public :: cpl_iocdf_open     ! open a named file
   public :: cpl_iocdf_close    ! close an open file
   public :: cpl_iocdf_set64bit ! select 32 or 64 bit real data in file
   public :: cpl_iocdf_append   ! add data to an existing file
 
! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   !----- module variables -----
   logical,save           :: flag64bit = .false. ! 32 or 64 bit netCDF reals?
   integer(IN),parameter  :: pid0 = 0            ! root process pid = zero
   character(*),parameter :: modName = "cpl_iocdf_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iocdf_create -- create a new file.
!
! !DESCRIPTION:
!    Create a new netCDF file with name {\tt fName} and no content other than
!    global attributes.  If optional argument {\tt desc} is present, it
!    will be placed in the ``description'' global attribute.
!
! !REVISION HISTORY:
!    2002-Oct-22 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iocdf_create(fName,desc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)          :: fName   ! file name
   character(*),intent(in),optional :: desc    ! description string

!EOP

   !----- local -----
   character(CL)       :: str      ! restart file ID string
   character(len= 8)   :: datestr  ! current date
   character(len=10)   :: timestr  ! current time
   integer(IN)         :: fid      ! file ID
   integer(IN)         :: vid      ! variable ID
   integer(IN)         :: vdid(3)  ! vector of dimension IDs
   integer(IN)         :: did      ! dimension ID
   integer(IN)         :: rcode    ! return code

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iocdf_create) '
   character(*),parameter :: F00 = "('(cpl_iocdf_create) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iocdf_create) ',3a,i4)"
   character(*),parameter :: F02 = "('(cpl_iocdf_create) ',a,i6,f10.2)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (dbug>1) write(6,F00) 'creating new file ',trim(fName)

   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe only
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0 ) return

   !----------------------------------------------------------------------------
   ! create the file, clobbering any existing file, add global attributes
   !----------------------------------------------------------------------------
   rcode = nf_create(trim(fName),NF_CLOBBER,fid)
   if (rcode /= NF_NOERR) then
      write(6,F00) "ERROR: creating file ",trim(fName) 
      write(6,F00) "ERROR: ",nf_strerror(rcode)
      call shr_sys_abort(subName // nf_strerror(rcode))
   end if
   call shr_sys_flush(6)

   !----------------------------------------------------------------------------
   ! add global attributes 
   !----------------------------------------------------------------------------
   str   = "cpl6 output netCDF data file"
   rcode = nf_put_att_text(fid,NF_GLOBAL,'title'      ,len_trim(str),str)

   str   = '(no description provided) '
   if (present(desc)) str = desc
   rcode = nf_put_att_text(fid,NF_GLOBAL,'description',len_trim(str),str)

   call date_and_time(datestr,timestr)
   str = 'File created: '                                            &
   &     //datestr(1:4)//'-'//datestr(5:6)//'-'//datestr(7:8)//' '   &
   &     //timestr(1:2)//':'//timestr(3:4)//':'//timestr(5:6)
   rcode = nf_put_att_text(fid,NF_GLOBAL,'history'    ,len_trim(str),str)

   str   ='CF-1.0 conventions'
   rcode = nf_put_att_text(fid,NF_GLOBAL,'Conventions',len_trim(str),str)

   str   ='SVN $Id: cpl_iocdf_mod.F90 3380 2007-03-06 05:42:19Z robj $'
   str   = trim(str) // ' ' // char(10)  ! ascii #10 is a new line
   str   = trim(str) // ' $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_iocdf_mod.F90 $'
   rcode = nf_put_att_text(fid,NF_GLOBAL,'RevisionControl',len_trim(str),str)

   !----------------------------------------------------------------------------
   ! special handling of time coords & related 
   !----------------------------------------------------------------------------
   rcode = nf_def_dim(fid,'time',NF_UNLIMITED,did)
   rcode = nf_def_var(fid,'time',NF_DOUBLE,1,did,vid)
   str   = 'time'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   str   = 'days since 0000-01-01 00:00:00'
   rcode = nf_put_att_text(fid,vid,'units'    ,len_trim(str),str)
   str   = 'noleap'
   rcode = nf_put_att_text(fid,vid,'calendar' ,len_trim(str),str)
   str   = 'time_bound'
   rcode = nf_put_att_text(fid,vid,"bounds"   ,len_trim(str),str)

   rcode = nf_def_dim(fid,'d2',2,did) ! upper & lower time bounds

   rcode = nf_inq_dimid(fid,'d2'  ,vdid(1))
   rcode = nf_inq_dimid(fid,'time',vdid(2))
   rcode = nf_def_var(fid,'time_bound',NF_DOUBLE,2,vdid,vid)
   str   = 'time average interval boundaries'
   rcode = nf_put_att_text(fid,vid,'long_name',len_trim(str),str)
   str   = 'days since 0000-01-01 00:00:00'
   rcode = nf_put_att_text(fid,vid,'units'    ,len_trim(str),str)
   str   = 'noleap'
   rcode = nf_put_att_text(fid,vid,'calendar' ,len_trim(str),str)

   !----------------------------------------------------------------------------
   ! close the file
   !----------------------------------------------------------------------------
   rcode = nf_close(fid)
   if (rcode /= NF_NOERR) then
      write(6,F00) "ERROR: closing file"
      write(6,F00) "ERROR: ",nf_strerror(rcode)
      call shr_sys_abort(subName // nf_strerror(rcode))
   end if

   call shr_sys_flush(6)

end subroutine cpl_iocdf_create

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iocdf_open -- open an existing file.
!
! !DESCRIPTION:
!    Open an existing file with name {\tt fName} and return the
!    NetCDF id in the output argument {\tt fid}.
!
! !REVISION HISTORY:
!    2002-Oct-22 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iocdf_open(fName,fid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: fName   ! file name
   integer(IN) ,intent(out) :: fid     ! file ID

!EOP

   !----- local -----
   integer(IN)         :: rCode    ! return code

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iocdf_open) '
   character(*),parameter :: F00 = "('(cpl_iocdf_open) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iocdf_open) ',3a,i4)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (dbug>1) write(6,F00) 'opening existing file ',trim(fName)

   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0) return

   !----------------------------------------------------------------------------
   ! open the file
   !----------------------------------------------------------------------------
   rcode = nf_open  (trim(fName),NF_WRITE  ,fid)
   if (rcode /= NF_NOERR) then
      write(6,F00) "ERROR: opening file ",trim(fName)
      write(6,F00) "ERROR: ",nf_strerror(rcode)
      call shr_sys_abort(subName // nf_strerror(rcode))
   end if

end subroutine cpl_iocdf_open

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iocdf_close -- close a file.
!
! !DESCRIPTION:
!    Close the netCDF file with netCDF id {\tt fid}.
!
! !REVISION HISTORY:
!    2002-Oct-22 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iocdf_close(fid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in) :: fid     ! file ID

!EOP

   !----- local -----
   integer(IN) :: rcode    ! return code

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iocdf_close) '
   character(*),parameter :: F00 = "('(cpl_iocdf_close) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iocdf_close) ',3a,i4)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (dbug>1) write(6,F00) 'closing data file'

   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0) return

   !----------------------------------------------------------------------------
   ! close the file
   !----------------------------------------------------------------------------
   rcode = nf_close(fid)

   call shr_sys_flush(6)

end subroutine cpl_iocdf_close

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iocdf_set64bit -- flags creation of 64 bit netCDF files.
!
! !DESCRIPTION:
!    Flags creation of 64 bit netCDF files, default is 32 bit.
!    If argument {\tt flag} is true, netCDF files with be 64 bit.
!
! !REMARKS:    
!    64 bit netCDF data was introduced for regression testing of coupled system.
!
! !REVISION HISTORY:
!    2004-Mar-31 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iocdf_set64bit(flag)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   logical,intent(in) :: flag    ! true <=> 64 bit

!EOP

   !----- formats -----
   character(*),parameter :: F01 = "('(cpl_iocdf_set64bit) ',a,L7)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   flag64bit = flag  ! set flag that is global within this module

   if (dbug>1) write(6,F01) 'reset 64 bit flag to ',flag64bit
   call shr_sys_flush(6)

end subroutine cpl_iocdf_set64bit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iocdf_append -- add data to an existing file.
!
! !DESCRIPTION:
!    Append all data in {\it bundle} {\tt bun} to an existing netCDF file 
!    with netCDF id {\tt fid}. If {\tt fName} is present, use 
!    {\tt  cpl\_iocdf\_open} to set {\tt fid}.
!
!    If the input time given by the range {\tt dateS, dateE} is not already 
!    in the time dimension of the file, it is appended to the end of the time
!    dimension.  This algorithm assumes that any input time value 
!    not already in the time dimension is greater than any value in the time 
!    dimension (so that time will be monotonically increasing).
!
!    If {\tt fName} is not present, {\tt fid} must be a valid netCDF file ID for an open
!    netCDF file.  This is presumably the normal mode of operation as it can
!    minimize the opening and closing of a file.
!
!    If {\tt fName} is present, the named file is opened, data is appended, and the
!    file is closed.  In this case, the input {\tt fid} is not used, although it's
!    value will be overwritten.  This mode is handy for ``one-off'' debug files.
!    In either case, the file being appended to must have previously been
!    created.
!
! !REVISION HISTORY:
!    2002-Oct-22 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iocdf_append(fid,date,bun,dateS,dateE,fName)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN)     ,intent(inout)       :: fid       ! file ID
   type(cpl_bundle),intent(in)          :: bun       ! bundle
   type(shr_date)  ,intent(in)          :: date      ! model date
   type(shr_date)  ,intent(in),optional :: dateS     ! time-bound start date
   type(shr_date)  ,intent(in),optional :: dateE     ! time-bound end   date
   character(*)    ,intent(in),optional :: fName     ! file name

!EOP

   !----- local -----
   character(CL)       :: str      ! generic text string
   character( 8)       :: datestr  ! current date
   character(10)       :: timestr  ! current time
   integer(IN)         :: start(3) ! vector of starting indicies
   integer(IN)         :: count(3) ! vector of edge lengths
   integer(IN)         :: vdid(3)  ! vector of dimension IDs
   integer(IN)         :: vid      ! variable ID
   integer(IN)         :: did      ! dimension ID
   type(mct_aVect) :: gData    ! global/gathered bundle data
   type(mct_aVect) :: gGrid    ! global/gathered bundle data
   integer(IN)         :: nflds    ! number of fields in a bundle
   integer(IN)         :: npts     ! number of points in a field
   integer(IN)         :: rcode    ! return code
   integer(IN)         :: n,k      ! generic index
   integer(IN)         :: kmax     ! maximum value for k index
   real(R8)            :: time     ! model time: eday + frac of day
   real(R8)            :: timeS    ! model time: time-bound start
   real(R8)            :: timeE    ! model time: time-bound end
   real(R8),parameter  :: time_eps = .01_R8/86400.0_R8  ! time epsilon = .01 seconds
   real(R8)            :: tmp      ! temporary real variable
   real(R8),allocatable:: XXX(:)   ! temporary real vector
   integer(IN)         :: eDay     ! model date: elapsed days
   integer(IN)         :: sec      ! model date: seconds

   type(mct_string) :: mctStr  ! temporary mct/mpeu string var
   character(80)        :: dSuffix ! domain suffix
   character(80)        :: dName   ! domain name
   character(80)        :: aName   ! attribute (field) name
   character(80)        :: bName   ! bundle name
   character(80)        :: nName   ! netCDF variable name
   character(80)        :: lName   ! netCDF long name
   character(80)        :: uName   ! netCDF units
   character(80)        ::  n_str  ! used to construct dimension var name
   character(80)        :: ni_str  ! used to construct dimension var name
   character(80)        :: nj_str  ! used to construct dimension var name
   character(80)        :: xc_str  ! used to coord array var name
   character(80)        :: yc_str  ! used to coord array var name

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iocdf_append) '
   character(*),parameter :: F00 = "('(cpl_iocdf_append) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iocdf_append) ',3a,i4)"
   character(*),parameter :: F02 = "('(cpl_iocdf_append) ',a,i6,3(f10.2,1x))"
   character(*),parameter :: F03 = "('(cpl_iocdf_append) ',a,3(f10.2,1x))"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (dbug>2) write(6,F00) &
   & 'adding data to file ',trim(fName),', bundle = ',trim(bun%name)

   !----------------------------------------------------------------------------
   ! do global gather to create non-distributed aVect (*part* of a bundle)
   !----------------------------------------------------------------------------
   call mct_aVect_gather(bun%data,gData,bun%dom%gsMap,pid0,cpl_comm_comp,rcode)
   call mct_aVect_gather(bun%dom%lGrid,gGrid,bun%dom%gsMap,pid0,cpl_comm_comp,rcode)
   
   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0) return

   !----------------------------------------------------------------------------
   ! if file name is present, open & close the file (note: fid is altered)
   !----------------------------------------------------------------------------
   if (present(fName)) call cpl_iocdf_open(str,fid) 

   !----------------------------------------------------------------------------
   ! add dimension & domain info (if it isn't there already)
   !----------------------------------------------------------------------------
   dSuffix = bun%dom%suffix
   dName   = bun%dom%name
    n_str =  "n_"//dSuffix
   ni_str = "ni_"//dSuffix
   nj_str = "nj_"//dSuffix

   rcode = nf_inq_dimid(fid,n_str,did)
   if ( rcode /= NF_NOERR ) then ! assume error means dim has not be defined yet
     
      !--- define dimension variables ---
      if (dbug>2) write(6,F00) 'adding domain info for domain = ',trim(dName)
      rcode = nf_redef(fid) ! (re) enter "define" mode

      rcode = nf_def_dim(fid,  n_str, bun%dom%n , did)
      rcode = nf_def_dim(fid, ni_str, bun%dom%ni, vdid(1))
      rcode = nf_def_dim(fid, nj_str, bun%dom%nj, vdid(2))

      !--- define domain variables ---
      do n=1,mct_aVect_nRAttr(gGrid)

         call mct_list_get(mctStr,n,gGrid%rList)
         aName = mct_string_toChar(mctStr)
         nName = "domain_" // trim(dSuffix) // "_" // trim(aName) 
         !--- work around for bug wrt "lat" meaning latitude vs latent ---
         str = aName
         if (trim(aName)=="lat") str = "lati" 
         call cpl_fields_getLongName(trim(str),lName,uName)

         rcode = nf_def_var  (fid,trim(nName),NF_FLOAT ,2,vdid,vid)
         rcode = nf_put_att_text(fid,vid,"long_name",len_trim(lName),lName)
         if ( trim(uName) /= "unitless" ) then
            rcode = nf_put_att_text(fid,vid,"units" ,len_trim(uName),uName)
         end if
         call mct_string_clean(mctStr)

      end do

      !--- put field data into nc file ---
      rcode = nf_enddef(fid) ! end "define" mode

      npts  = mct_aVect_lsize (gGrid)
      allocate(XXX(npts))

      do n=1,mct_aVect_nRAttr(gGrid)

         !--- construct netCDF var name ---
         call mct_list_get(mctStr,n,gGrid%rList)
         aName = mct_string_toChar(mctStr)
         nName = "domain_" // trim(dSuffix) // "_" // trim(aName) 
         call mct_string_clean(mctStr)

         !--- get variable's ID ---
         rcode = nf_inq_varid(fid, trim(nName) , vid)
         if (rcode /= NF_NOERR) then
            write(6,F00) "ERROR: getting vid for ",trim(nName)
            write(6,F00) "ERROR: ",nf_strerror(rcode)
            call shr_sys_abort(subName // nf_strerror(rcode))
         end if
   
         !--- put the data ---
         k = mct_aVect_indexRA(gGrid,trim(aName),perrWith=subName)
         XXX = gGrid%rAttr(k,:)
         rcode = nf_put_var_double(fid,vid,XXX)

      end do

      deallocate(XXX)

   else
      if (dbug>2) write(6,F00) 'domain info already exists for ',trim(dName)
   end if

   !----------------------------------------------------------------------------
   ! define fields (if they haven't been already)
   !----------------------------------------------------------------------------
   rcode = nf_redef(fid) ! re-enter "define" mode

   rcode = nf_inq_dimid(fid,ni_str,vdid(1))
   rcode = nf_inq_dimid(fid,nj_str,vdid(2))
   rcode = nf_inq_dimid(fid,"time",vdid(3))

   do n=1,mct_aVect_nRAttr(bun%data)

      !--- construct netCDF var name ---
      call mct_list_get(mctStr,n,bun%data%rList)
      aName = mct_string_toChar(mctStr)
      call mct_string_clean(mctStr)
      bName = bun%name
      nName = trim(bName) // "_" // trim(aName) 

      !--- test if var already exists (no error => var exists) ---
      rcode = nf_inq_varid(fid,nName,vid)

      if ( rcode == NF_NOERR ) then
         !--- no error => var is defined, assume all have been defined ---
         if ( n == 1 ) then
           if (dbug>2) write(6,F00) 'bundle variables have already been defined'
         end if
         exit
      else
         !--- assume an error means the var has not been defined ---
         if ( n == 1 ) then
           if (dbug>2) write(6,F00) 'defining bundle variables'
         end if

         !--- define units & long-name ---
         call cpl_fields_getLongName(trim(aName),lName,uName)
         if (flag64bit) then
            rcode = nf_def_var  (fid,trim(nName),NF_DOUBLE,3,vdid,vid)
         else
            rcode = nf_def_var  (fid,trim(nName),NF_FLOAT ,3,vdid,vid)
         endif
         rcode = nf_put_att_text(fid,vid,"long_name",len_trim(lName),lName)
         if ( trim(uName) /= "unitless" ) then
            rcode = nf_put_att_text(fid,vid,"units" ,len_trim(uName),uName)
         end if

         !--- define missing value ---
         if (flag64bit) then
            rcode = nf_put_att_double(fid,vid,"missing_value",NF_DOUBLE,1,cpl_const_spval)
            if (rcode /= NF_NOERR) write(6,F02) nName, nf_strerror(rcode)
            rcode = nf_put_att_double(fid,vid,"_FillValue"   ,NF_DOUBLE,1,cpl_const_spval)
            if (rcode /= NF_NOERR) write(6,F02) nName, nf_strerror(rcode)
         else
            rcode = nf_put_att_double(fid,vid,"missing_value",NF_FLOAT ,1,cpl_const_spval)
            if (rcode /= NF_NOERR) write(6,F02) nName, nf_strerror(rcode)
            rcode = nf_put_att_double(fid,vid,"_FillValue"   ,NF_FLOAT ,1,cpl_const_spval)
            if (rcode /= NF_NOERR) write(6,F02) nName, nf_strerror(rcode)
         end if

         !--- indicate time-mean data ---
         if ( present(dateS) .and. (present(dateE)) ) then
            call shr_date_getEDay(dateS,eDay,sec)
            timeS  = real(eDay,R8) + sec/86400.0_R8
            call shr_date_getEDay(dateE,eDay,sec)
            timeE  = real(eDay,R8) + sec/86400.0_R8
            if ( timeS /= timeE ) then
               str = "time: mean"
               rcode = nf_put_att_text(fid,vid,"cell_method",len_trim(str),str)
               if (rcode /= NF_NOERR) write(6,F02) nName, nf_strerror(rcode)
            end if
         end if
      end if

   end do

   !----------------------------------------------------------------------------
   ! find corresponding time index (lengthen time dim if necessary)
   !----------------------------------------------------------------------------
   rcode = nf_enddef(fid) ! end "define" mode

   !--- get current time & time-bounds, units are elapsed days ---
   call shr_date_getEDay(date,eDay,sec)
   time  = real(eDay,R8) + sec/86400.0_R8
   timeS = time
   timeE = time
   if (present(dateS)) then
      call shr_date_getEDay(dateS,eDay,sec)
      timeS  = real(eDay,R8) + sec/86400.0_R8
   end if
   if (present(dateE)) then
      call shr_date_getEDay(dateE,eDay,sec)
      timeE  = real(eDay,R8) + sec/86400.0_R8
   end if
   if ( time < timeS  .or.  timeE < time ) then
      write(6,F00) "WARNING: inconsistant time values"
      write(6,F03) "WARNING: ts,t,te = ",timeS,time,timeE
   end if

   !--- is current time already in time dimension? ---
   rcode = nf_inq_dimid (fid,"time",did)
   rcode = nf_inq_dimlen(fid,did,kmax)
   k   =   0 ! k==0 => corresponding time index not found
   tmp = 0.0_R8 
   do n=kmax,1,-1 ! desired time index is most likely n=kmax
      rcode = nf_get_var1_double(fid,did,n,tmp)
      if ( time >= (tmp + time_eps) ) exit ! time is monotonic, append new time
      if ( abs(time-tmp) < time_eps ) then
         k=n   ! k>0 => time index found
         exit
      end if
   end do

   !--- if corresponding time index not found, append to time dim ---
   if (k==0) then
      if (dbug>2) write(6,F00) 'lengthening time dimension'
      k=kmax+1

      !--- add time data ---
      start(1)=k
      count(1)=1
      rcode = nf_inq_varid(fid,'time',vid)
      if (rcode /= NF_NOERR) then
         write(6,F00) "ERROR: getting vid for time"
         write(6,F00) "ERROR: ",nf_strerror(rcode)
         call shr_sys_abort(subName // nf_strerror(rcode))
      end if
      rcode = nf_put_vara_double(fid, vid, start, count, time)
      if (rcode /= NF_NOERR) then
         write(6,F00) "ERROR: putting data into time"
         write(6,F00) "ERROR: ",nf_strerror(rcode)
         call shr_sys_abort(subName // nf_strerror(rcode))
      end if

      !--- add time_bounds data ---
      start(1)=1
      start(2)=k
      count(1)=1
      count(2)=1
      rcode = nf_inq_varid(fid,'time_bound',vid)
      if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
      rcode = nf_put_vara_double(fid, vid, start, count, timeS)
      if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
      start(1)=2
      rcode = nf_put_vara_double(fid, vid, start, count, timeE )
      if (rcode /= NF_NOERR) write(6,F00) nf_strerror(rcode)
   endif

   if (dbug>2) write(6,F02) 'time dimension index & time = ',k,time
   call shr_sys_flush(6)

   !----------------------------------------------------------------------------
   if (dbug>2) write(6,F00) 'put all bundle fields data into nc file'
   !----------------------------------------------------------------------------
   rcode = nf_enddef(fid) ! end "define" mode

   start(1) = 1
   start(2) = 1
   start(3) = k
   count(1) = bun%dom%ni
   count(2) = bun%dom%nj
   count(3) = 1

   npts  = mct_aVect_lsize (gGrid)
   allocate(XXX(npts))

   do n=1,mct_aVect_nRAttr(bun%data)

      !--- construct netCDF var name ---
      bName = bun%name
      call mct_list_get(mctStr,n,bun%data%rList)
      aName = mct_string_toChar(mctStr)
      call mct_string_clean(mctStr)
      nName = trim(bName) // "_" // trim(aName) 

      !--- get netCDF var ID ---
      rcode = nf_inq_varid(fid, trim(nName) , vid)
      if (rcode /= NF_NOERR) then
         write(6,F00) "ERROR: ",nf_strerror(rcode)
         write(6,F00) "ERROR: putting data into ",trim(nName)
         call shr_sys_abort(subName // nf_strerror(rcode))
      end if

      !--- put the data ---
      k = mct_aVect_indexRA(gData,trim(aName),perrWith=subName)
      XXX = gData%rAttr(k,:)
      rcode = nf_put_vara_double(fid,vid,start,count,XXX)
      if (rcode /= NF_NOERR) then
         write(6,F00) "ERROR: ",nf_strerror(rcode)
         write(6,F00) "ERROR: putting data into ",trim(nName)
         call shr_sys_abort(subName // nf_strerror(rcode))
      end if

   end do

   deallocate(XXX)

   call mct_aVect_clean(gData)
   call mct_aVect_clean(gGrid)
   call shr_sys_flush(6)

   !----------------------------------------------------------------------------
   ! if file name is present, open & close the file
   !----------------------------------------------------------------------------
   if (present(fName)) call cpl_iocdf_close(fid)

end subroutine cpl_iocdf_append

!===============================================================================
!===============================================================================

end module cpl_iocdf_mod
