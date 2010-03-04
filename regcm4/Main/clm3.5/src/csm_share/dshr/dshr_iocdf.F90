!===============================================================================
! SVN $Id: dshr_iocdf.F90 1051 2006-05-25 22:35:41Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_iocdf.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_iocdf -- create and write-to netcdf data file.
!
! !DESCRIPTION:
!    This module creates and writes data in a {\tt bundle} to a file
!    in machine-independent netCDF format.  This module is intended to support 
!    the history writing sub-system of dx7 models.  Assumes this is a 
!    single pe executable, no communication or pe management required.
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2005-Aug-01 - T. Craig - first version, from cpl_iocdf_mod.F90
!
! !INTERFACE:  -----------------------------------------------------------------

module dshr_iocdf

! !USES:

   use dshr_bundle        ! defines bundle
   use dshr_domain        ! defines domain
   use dshr_kind          ! defines F90 kinds
   use dshr_const         ! defines constants (eg. spval)
   use shr_sys_mod        ! share system routines
   use shr_date_mod       ! defines date data-type
   use shr_string_mod     ! string methods
   use shr_ncread_mod     ! ncread routines
   use netcdf             ! netcdf

   implicit none

   private ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: dshr_iocdf_create    ! create a new file (an empty file)
   public :: dshr_iocdf_open      ! open a named file
   public :: dshr_iocdf_close     ! close an open file
   public :: dshr_iocdf_set64bit  ! select 32 or 64 bit real data in file
   public :: dshr_iocdf_append    ! add data to an existing file
   public :: dshr_iocdf_appendAtt ! add a global attribute to an existing file
   public :: dshr_iocdf_read      ! read data from an existing file
   public :: dshr_iocdf_readAtt   ! read a global attribute from an existing file
 
! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   !----- module variables -----
   logical,save           :: flag64bit = .true.  ! 32 or 64 bit netCDF reals?
   logical,save           :: doabort = .true.    ! local abort flag
   integer(IN),save       :: debug = 4
   character(*),parameter :: modName = "dshr_iocdf"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_create -- create a new file.
!
! !DESCRIPTION:
!    Create a new netCDF file with name {\tt fName} and no content other than
!    global attributes.  If optional argument {\tt desc} is present, it
!    will be placed in the ``description'' global attribute.
!
! !REVISION HISTORY:
!    2005-Aug-01 - T. Craig, initial version
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_create(fName,title,desc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)          :: fName   ! file name
   character(*),intent(in),optional :: title   ! title       string
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
   character(*),parameter :: subName = '(dshr_iocdf_create) '
   character(*),parameter :: F00 = "('(dshr_iocdf_create) ',4a)"
   character(*),parameter :: F01 = "('(dshr_iocdf_create) ',3a,i4)"
   character(*),parameter :: F02 = "('(dshr_iocdf_create) ',a,i6,f10.2)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (debug>1) write(6,F00) 'creating new file ',trim(fName)

   !----------------------------------------------------------------------------
   ! create the file, clobbering any existing file, add global attributes
   !----------------------------------------------------------------------------
   rcode = nf90_create(trim(fName),NF90_CLOBBER,fid)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR create")

   !----------------------------------------------------------------------------
   ! add global attributes 
   !----------------------------------------------------------------------------
   str   = '(no title provided) '
   if (present(title)) str = title
   rcode = nf90_put_att(fid,NF90_GLOBAL,'title'      ,trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put global att")

   str   = '(no description provided) '
   if (present(desc)) str = desc
   rcode = nf90_put_att(fid,NF90_GLOBAL,'description',trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put global att")

   call date_and_time(datestr,timestr)
   str = 'File created: '                                            &
   &     //datestr(1:4)//'-'//datestr(5:6)//'-'//datestr(7:8)//' '   &
   &     //timestr(1:2)//':'//timestr(3:4)//':'//timestr(5:6)
   rcode = nf90_put_att(fid,NF90_GLOBAL,'history'    ,trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put global att")

   str   ='CF-1.0 conventions'
   rcode = nf90_put_att(fid,NF90_GLOBAL,'Conventions',trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put global att")

   str   = 'SVN $Id: dshr_iocdf.F90 1051 2006-05-25 22:35:41Z kauff $'
   str   = trim(str) // ' ' // char(10)  ! ascii #10 is a new line
   str   = trim(str) // ' $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_iocdf.F90 $'
   rcode = nf90_put_att(fid,NF90_GLOBAL,'RevisionControl',trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put global att")

   !----------------------------------------------------------------------------
   ! special handling of time coords & related 
   !----------------------------------------------------------------------------
   rcode = nf90_def_dim(fid,'time',NF90_UNLIMITED,did)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR def_dim time")
   rcode = nf90_def_var(fid,'time',NF90_DOUBLE,did,vid)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR def_var time")
   str   = 'time'
   rcode = nf90_put_att(fid,vid,'long_name',trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put_att long_name")
   str   = 'days since 0000-01-01 00:00:00'
   rcode = nf90_put_att(fid,vid,'units'    ,trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put_att units")
   str   = 'noleap'
   rcode = nf90_put_att(fid,vid,'calendar' ,trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put_att calendar")

   !----------------------------------------------------------------------------
   ! close the file
   !----------------------------------------------------------------------------
   rcode = nf90_close(fid)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR close")

   call shr_sys_flush(6)

end subroutine dshr_iocdf_create

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_open -- open an existing file.
!
! !DESCRIPTION:
!    Open an existing file with name {\tt fName} and return the
!    NetCDF id in the output argument {\tt fid}.
!
! !REVISION HISTORY:
!    2005-Aug-01 - T. Craig, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_open(fName,fid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: fName   ! file name
   integer(IN) ,intent(out) :: fid     ! file ID

!EOP

   !----- local -----
   integer(IN)         :: rCode    ! return code

   !----- formats -----
   character(*),parameter :: subName = '(dshr_iocdf_open) '
   character(*),parameter :: F00 = "('(dshr_iocdf_open) ',4a)"
   character(*),parameter :: F01 = "('(dshr_iocdf_open) ',3a,i4)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (debug>1) write(6,F00) 'opening existing file ',trim(fName)

   !----------------------------------------------------------------------------
   ! open the file
   !----------------------------------------------------------------------------
   rcode = nf90_open  (trim(fName),NF90_WRITE  ,fid)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR open "//trim(fName))

end subroutine dshr_iocdf_open

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_close -- close a file.
!
! !DESCRIPTION:
!    Close the netCDF file with netCDF id {\tt fid}.
!
! !REVISION HISTORY:
!    2005-Aug-01 - T. Craig, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_close(fid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in) :: fid     ! file ID

!EOP

   !----- local -----
   integer(IN) :: rcode    ! return code

   !----- formats -----
   character(*),parameter :: subName = '(dshr_iocdf_close) '
   character(*),parameter :: F00 = "('(dshr_iocdf_close) ',4a)"
   character(*),parameter :: F01 = "('(dshr_iocdf_close) ',3a,i4)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (debug>1) write(6,F00) 'closing data file'

   !----------------------------------------------------------------------------
   ! close the file
   !----------------------------------------------------------------------------
   rcode = nf90_close(fid)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR close")

   call shr_sys_flush(6)

end subroutine dshr_iocdf_close

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_set64bit -- flags creation of 64 bit netCDF files.
!
! !DESCRIPTION:
!    Flags creation of 64 bit netCDF files, default is 32 bit.
!    If argument {\tt flag} is true, netCDF files with be 64 bit.
!
! !REMARKS:    
!    64 bit netCDF data was introduced for regression testing of coupled system.
!
! !REVISION HISTORY:
!    2005-Aug-01 - T. Craig, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_set64bit(flag)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   logical,intent(in) :: flag    ! true <=> 64 bit

!EOP

   !----- formats -----
   character(*),parameter :: F01 = "('(dshr_iocdf_set64bit) ',a,L7)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   flag64bit = flag  ! set flag that is global within this module

   if (debug>1) write(6,F01) 'reset 64 bit flag to ',flag64bit
   call shr_sys_flush(6)

end subroutine dshr_iocdf_set64bit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_append -- add data to an existing file.
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
!    2005-Aug-01 - T. Craig, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_append(fid,date,bun,fName)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN)                 ,intent(inout)       :: fid    ! file ID
   type(dshr_bundle_bundleType),intent(in)          :: bun    ! bundle
   type(shr_date)              ,intent(in)          :: date   ! model date
   character(*)                ,intent(in),optional :: fName  ! file name

!EOP

   !----- local -----
   type(dshr_domain_domainType),pointer  :: dom    ! domain in bundle
   character(CX)       :: fldList  ! bundle field list
   character(CL)       :: str      ! generic text string
   integer(IN)         :: start(3) ! vector of starting indicies
   integer(IN)         :: count(3) ! vector of edge lengths
   integer(IN)         :: vdid2(2) ! vector of dimension IDs
   integer(IN)         :: vdid3(3) ! vector of dimension IDs
   integer(IN)         :: vid      ! variable ID
   integer(IN)         :: did      ! dimension ID
   integer(IN)         :: ni,nj    ! domain grid size
   integer(IN)         :: rcode    ! return code
   integer(IN)         :: n,k      ! generic index
   integer(IN)         :: kmax     ! maximum value for k index
   real(R8)            :: time     ! model time: eday + frac of day
   real(R8),parameter  :: time_eps = .01_R8/86400.0_R8  ! time epsilon = .01 seconds
   real(R8),allocatable:: tmp(:)   ! temporary real variable
   real(R8),pointer    :: XXX(:,:) ! temporary real array
   integer(IN)         :: eDay     ! model date: elapsed days
   integer(IN)         :: sec      ! model date: seconds

   character(CS)       :: dName    ! domain name
   character(CS)       :: aName    ! attribute (field) name
   character(CS)       :: bName    ! bundle name
   character(CS)       :: nName    ! netCDF variable name
   character(CS)       ::  n_str   ! used to construct dimension var name
   character(CS)       :: ni_str   ! used to construct dimension var name
   character(CS)       :: nj_str   ! used to construct dimension var name

   !----- formats -----
   character(*),parameter :: subName = '(dshr_iocdf_append) '
   character(*),parameter :: F00 = "('(dshr_iocdf_append) ',4a)"
   character(*),parameter :: F01 = "('(dshr_iocdf_append) ',3a,i4)"
   character(*),parameter :: F02 = "('(dshr_iocdf_append) ',a,i6,3(f10.2,1x))"
   character(*),parameter :: F03 = "('(dshr_iocdf_append) ',a,3(f10.2,1x))"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   bName = dshr_bundle_getName(bun)
   if (debug>0) write(6,F00) 'adding to file ',trim(fName),', bundle = ',trim(bName)
   call dshr_bundle_getFieldList(bun,fldList)
   call dshr_bundle_print(bun)
   call dshr_bundle_domainPtr(bun,dom)
   call dshr_domain_getDims(dom,ni,nj)

   !----------------------------------------------------------------------------
   ! if file name is present, open & close the file (note: fid is altered)
   !----------------------------------------------------------------------------
   if (present(fName)) call dshr_iocdf_open(fName,fid) 

   !----------------------------------------------------------------------------
   ! add dimension & domain info (if it isn't there already)
   !----------------------------------------------------------------------------
   call dshr_domain_getName(dom,dName)
   n_str  = "n_" //dName
   ni_str = "ni_"//dName
   nj_str = "nj_"//dName

   rcode = nf90_inq_dimid(fid,n_str,did)
   if ( rcode /= NF90_NOERR ) then ! assume error means dim has not be defined yet
     
      !--- define dimension variables ---
      if (debug>2) write(6,F00) 'adding domain info for domain = ',trim(dName)
      rcode = nf90_redef(fid) ! (re) enter "define" mode
      call dshr_iocdf_handleErr(rcode,subName//" ERROR redef")

      rcode = nf90_def_dim(fid,  n_str, ni*nj , did)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR def_dim n_str")
      rcode = nf90_def_dim(fid, ni_str, ni, vdid2(1))
      call dshr_iocdf_handleErr(rcode,subName//" ERROR def_dim ni_str")
      rcode = nf90_def_dim(fid, nj_str, nj, vdid2(2))
      call dshr_iocdf_handleErr(rcode,subName//" ERROR def_dim nj_str")

      !--- define domain variables ---
      aName = 'lat'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_def_var  (fid,trim(nName),NF90_DOUBLE ,vdid2,vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR def_var "//trim(nName))

      aName = 'lon'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_def_var  (fid,trim(nName),NF90_DOUBLE ,vdid2,vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR def_var "//trim(nName))

      aName = 'area'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_def_var  (fid,trim(nName),NF90_DOUBLE ,vdid2,vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR def_var "//trim(nName))

      aName = 'mask'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_def_var  (fid,trim(nName),NF90_DOUBLE ,vdid2,vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR def_var "//trim(nName))

      aName = 'index'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_def_var  (fid,trim(nName),NF90_DOUBLE ,vdid2,vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR def_var "//trim(nName))

      !--- put field data into nc file ---
      rcode = nf90_enddef(fid) ! end "define" mode
      call dshr_iocdf_handleErr(rcode,subName//" ERROR enddef")

      allocate(XXX(ni,nj))

      aName = 'lat'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_inq_varid(fid, trim(nName) , vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_varid "//trim(nName))
      call dshr_domain_getData(dom,XXX,'lat')
      rcode = nf90_put_var(fid,vid,XXX)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR put_var lat")

      aName = 'lon'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_inq_varid(fid, trim(nName) , vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_varid "//trim(nName))
      call dshr_domain_getData(dom,XXX,'lon')
      rcode = nf90_put_var(fid,vid,XXX)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR put_var lon")

      aName = 'area'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_inq_varid(fid, trim(nName) , vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_varid "//trim(nName))
      call dshr_domain_getData(dom,XXX,'area')
      rcode = nf90_put_var(fid,vid,XXX)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR put_var area")

      aName = 'mask'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_inq_varid(fid, trim(nName) , vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_varid "//trim(nName))
      call dshr_domain_getData(dom,XXX,'mask')
      rcode = nf90_put_var(fid,vid,XXX)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR put_var mask")

      aName = 'index'
      nName = trim(dName) // "_" // trim(aName) 
      rcode = nf90_inq_varid(fid, trim(nName) , vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_varid "//trim(nName))
      call dshr_domain_getData(dom,XXX,'index')
      rcode = nf90_put_var(fid,vid,XXX)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR put_var index")

      deallocate(XXX)

   else
      if (debug>2) write(6,F00) 'domain info already exists for ',trim(dName)
   end if

   !----------------------------------------------------------------------------
   ! define fields (if they haven't been already)
   !----------------------------------------------------------------------------
   rcode = nf90_redef(fid) ! re-enter "define" mode
   call dshr_iocdf_handleErr(rcode,subName//" ERROR redef")

   rcode = nf90_inq_dimid(fid,ni_str,vdid3(1))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_dimid ni_str")
   rcode = nf90_inq_dimid(fid,nj_str,vdid3(2))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_dimid nj_str")
   rcode = nf90_inq_dimid(fid,"time",vdid3(3))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_dimid time")

   do n=1,shr_string_listGetNum(fldList)

      !--- construct netCDF var name ---
      call shr_string_listGetName(fldList,n,aName)
      bName = dshr_bundle_getName(bun)
      nName = trim(bName) // "_" // trim(aName) 

      !--- test if var already exists (no error => var exists) ---
      rcode = nf90_inq_varid(fid,nName,vid)

      if ( rcode == NF90_NOERR ) then
         !--- no error => var is defined, assume all have been defined ---
         if ( n == 1 ) then
           if (debug>2) write(6,F00) 'bundle variables have already been defined'
         end if
         exit
      else
         !--- assume an error means the var has not been defined ---
         if ( n == 1 ) then
           if (debug>2) write(6,F00) 'defining bundle variables'
         end if

         !--- set type, float or double ---
         if (flag64bit) then
            rcode = nf90_def_var  (fid,trim(nName),NF90_DOUBLE,vdid3,vid)
            call dshr_iocdf_handleErr(rcode,subName//" ERROR def_var "//trim(nName))
         else
            rcode = nf90_def_var  (fid,trim(nName),NF90_FLOAT ,vdid3,vid)
            call dshr_iocdf_handleErr(rcode,subName//" ERROR def_var "//trim(nName))
         endif

         !--- define missing value ---
         if (flag64bit) then
            rcode = nf90_put_att(fid,vid,"missing_value",dshr_const_spval)
            call dshr_iocdf_handleErr(rcode,subName//" ERROR put_att missingvalue")
            rcode = nf90_put_att(fid,vid,"_FillValue"   ,dshr_const_spval)
            call dshr_iocdf_handleErr(rcode,subName//" ERROR put_att fillvalue")
         else
            rcode = nf90_put_att(fid,vid,"missing_value",dshr_const_spval)
            call dshr_iocdf_handleErr(rcode,subName//" ERROR put_att missingvalue")
            rcode = nf90_put_att(fid,vid,"_FillValue"   ,dshr_const_spval)
            call dshr_iocdf_handleErr(rcode,subName//" ERROR put_att fillvalue")
         end if
      endif
   end do

   !----------------------------------------------------------------------------
   ! find corresponding time index (lengthen time dim if necessary)
   !----------------------------------------------------------------------------
   rcode = nf90_enddef(fid) ! end "define" mode
   call dshr_iocdf_handleErr(rcode,subName//" ERROR enddef")

   !--- get current time, units are elapsed days ---
   call shr_date_getEDay(date,eDay,sec)
   time  = real(eDay,R8) + sec/86400.0_R8

   !--- is current time already in time dimension? ---
   rcode = nf90_inq_dimid (fid,"time",did)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_dimid time")
   rcode = nf90_inquire_dimension(fid,did,len=kmax)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR inquire_dimension")
   allocate(tmp(kmax))
   rcode = nf90_get_var(fid,did,tmp)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR get_var tmp")
   k = 0
   do n=kmax,1,-1 ! desired time index is most likely n=kmax
      if ( time >= (tmp(n) + time_eps) ) exit ! time is monotonic, append new time
      if ( abs(time-tmp(n)) < time_eps ) then
         k=n   ! k>0 => time index found
         exit
      end if
   end do
   deallocate(tmp)

   !--- if corresponding time index not found, append to time dim ---
   if (k==0) then
      if (debug>2) write(6,F00) 'lengthening time dimension'
      k=kmax+1

      !--- add time data ---
      rcode = nf90_inq_varid(fid,'time',vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR inq varid time")
      rcode = nf90_put_var(fid,vid,time,start=(/k/))
      call dshr_iocdf_handleErr(rcode,subName//" ERROR put var time")

   endif

   if (debug>2) write(6,F02) 'time dimension index & time = ',k,time
   call shr_sys_flush(6)

   !----------------------------------------------------------------------------
   if (debug>2) write(6,F00) 'put all bundle fields data into nc file'
   !----------------------------------------------------------------------------

   rcode = nf90_inq_dimid (fid,"time",did)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_dimid time")
   rcode = nf90_inquire_dimension(fid,did,len=kmax)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR inquire_dimension")
   allocate(tmp(kmax))
   rcode = nf90_get_var(fid,did,tmp)
   call dshr_iocdf_handleErr(rcode,subName//" ERROR get_var tmp")
   deallocate(tmp)

   start(1) = 1
   start(2) = 1
   start(3) = k
   count(1) = ni
   count(2) = nj
   count(3) = 1

   do n=1,shr_string_listGetNum(fldList)
      !--- construct netCDF var name ---
      call shr_string_listGetName(fldList,n,aName)
      bName = dshr_bundle_getName(bun)
      nName = trim(bName) // "_" // trim(aName) 

      !--- get netCDF var ID ---
      rcode = nf90_inq_varid(fid, trim(nName) , vid)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR inq_varid "//trim(nName))

      !--- put the data ---
      call dshr_bundle_assignPtr(bun,aName,XXX)
      rcode = nf90_put_var(fid,vid,XXX,start,count)
      call dshr_iocdf_handleErr(rcode,subName//" ERROR put_var XXX")

   end do

   !----------------------------------------------------------------------------
   ! if file name is present, open & close the file
   !----------------------------------------------------------------------------
   if (present(fName)) call dshr_iocdf_close(fid)

   call shr_sys_flush(6)

end subroutine dshr_iocdf_append

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_appendAtt -- add a global char attribute an existing file.
!
! !DESCRIPTION:
!    Add a global char attribute an existing file.
!
! !REVISION HISTORY:
!    2006-Feb-13 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_appendAtt(fid,attName,str,fName)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN) ,intent(inout)       :: fid     ! file ID
   character(*),intent(in)          :: attName ! attribute name
   character(*),intent(in)          :: str     ! attribute string
   character(*),intent(in),optional :: fName   ! file name

!EOP

   !----- local -----
   integer(IN) :: rCode     ! return Code

   !----- formats -----
   character(*),parameter :: subName = '(dshr_iocdf_appendAtt) '
   character(*),parameter :: F00 =   "('(dshr_iocdf_appendAtt) ',4a)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! if file name is present, open & close the file (note: fid is altered)
   !----------------------------------------------------------------------------
   if (present(fName)) call dshr_iocdf_open(fName,fid) 

   !----------------------------------------------------------------------------
   ! add global character attribute
   !----------------------------------------------------------------------------
   rcode = nf90_redef(fid) ! (re) enter "define" mode
   rcode = nf90_put_att(fid,nf90_global,trim(attName),trim(str))
   call dshr_iocdf_handleErr(rcode,subName//" ERROR put_att "//trim(attName))

   !----------------------------------------------------------------------------
   ! if file name is present, open & close the file
   !----------------------------------------------------------------------------
   if (present(fName)) call dshr_iocdf_close(fid)

   call shr_sys_flush(6)

end subroutine dshr_iocdf_appendAtt

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_read -- read data from an existing file.
!
! !DESCRIPTION:
!    Read all data in {\it bundle} {\tt bun} from an existing netCDF file
!    with filename {\tt fName}.   Use shr\_ncread extensively
!
!    Does not handles dates currently, uses the k=1 values
!
! !REVISION HISTORY:
!    2005-Aug-01 - T. Craig, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_read(fName,bun,date)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*)     ,intent(in)          :: fName     ! file name
   type(dshr_bundle_bundleType),intent(inout)  :: bun       ! bundle
   type(shr_date)   ,intent(in),optional :: date      ! model date

!EOP

   !----- local -----
   type(dshr_domain_domainType),pointer  :: dom    ! domain in bundle
   character(CX)       :: fldList  ! bundle field list
   integer(IN)         :: ni,nj    ! domain grid size from bundle
   integer(IN)         :: nd1,nd2  ! domain grid size from cdf file
   integer(IN)         :: tindex   ! time index of field
   integer(IN)         :: n,k      ! generic index
   real(R8),pointer    :: XXX(:,:) ! temporary real array

   character(80)        :: aName   ! attribute (field) name
   character(80)        :: bName   ! bundle name
   character(80)        :: nName   ! netCDF variable name

   !----- formats -----
   character(*),parameter :: subName = '(dshr_iocdf_read) '
   character(*),parameter :: F00 = "('(dshr_iocdf_read) ',4a)"
   character(*),parameter :: F01 = "('(dshr_iocdf_read) ',3a,i4)"
   character(*),parameter :: F02 = "('(dshr_iocdf_read) ',a,i6,3(f10.2,1x))"
   character(*),parameter :: F03 = "('(dshr_iocdf_read) ',a,3(f10.2,1x))"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (debug>2) write(6,F00) &
   & 'reading data from file ',trim(fName),', bundle = ',trim(dshr_bundle_getName(bun))
   call dshr_bundle_getFieldList(bun,fldList)
   call dshr_bundle_print(bun)
   call dshr_bundle_domainPtr(bun,dom)
   call dshr_domain_getDims(dom,ni,nj)

!   !--- get current time, units are elapsed days ---
   tindex = 1
!   call shr_date_getEDay(date,eDay,sec)
!   time  = real(eDay,R8) + sec/86400.0
   if (present(date)) then
     call shr_sys_abort(trim(subName)//' ERROR date option not supported')
   endif

   do n=1,shr_string_listGetNum(fldList)

      !--- construct netCDF var name ---
      call shr_string_listGetName(fldList,n,aName)
      bName = dshr_bundle_getName(bun)
      nName = trim(bName) // "_" // trim(aName) 

      call shr_ncread_varDimSize(fName,nName,1,nd1)
      call shr_ncread_varDimSize(fname,nName,2,nd2)
      if (nd1 /= ni .or. nd2 /= nj) then
         call dshr_iocdf_abort(trim(subName)//' ERROR ndim sizes of var '//trim(nName))
      endif

      call dshr_bundle_assignPtr(bun,aName,XXX)
      call shr_ncread_tField(fName,tindex,nName,XXX)

   end do

   call shr_sys_flush(6)

end subroutine dshr_iocdf_read

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_readAtt -- read a global char attribute from a file.
!
! !DESCRIPTION:
!    Read a global char attribute from a file.
!
! !REVISION HISTORY:
!    2006-Feb-13 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_readAtt(fName,attName,str)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: fName   ! file name
   character(*),intent(in)  :: attName ! attribute name
   character(*),intent(out) :: str     ! attribute string

!EOP

   !----- local -----
   integer(IN) :: fid   ! file ID
   integer(IN) :: n     ! attribute length
   integer(IN) :: rCode ! return code

   !----- formats -----
   character(*),parameter :: subName = '(dshr_iocdf_readAtt) '
   character(*),parameter :: F00 =   "('(dshr_iocdf_readAtt) ',8a)"
   character(*),parameter :: F01 =   "('(dshr_iocdf_readAtt) ',a,2i5)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   call dshr_iocdf_open(fName,fid) 

   rcode = nf90_inquire_attribute(fid,nf90_global,trim(attName),len=n)
   if ( rCode /= nf90_noerr ) then
      call dshr_iocdf_handleErr(rcode,subName//" ERROR inquire_att "//trim(attName))
   else if (len(str) < n) then
      write(6,F00) "ERROR: cannot read attribute  = ",trim(attName)
      write(6,F01) "ERROR: var & attribute length = ",len(str),n
      call dshr_iocdf_handleErr(rcode,subName//" ERROR att length "//trim(attName) )
   else
      str = " "
      rcode = nf90_get_att(fid,nf90_global,trim(attName),str)
      if ( rCode /= nf90_noerr )  &
         call dshr_iocdf_handleErr(rcode,subName//" ERROR get_att "//trim(attName) )
   end if
   if (debug>1) then
      write(6,F00) "file     : ",trim(fName)
      write(6,F00) "attribute: ",trim(attName)," = ",trim(str)
   end if

   call dshr_iocdf_close(fid)

   call shr_sys_flush(6)

end subroutine dshr_iocdf_readAtt

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_handleErr -- Print netCDF error message
!
! !DESCRIPTION:
!   Print the error message corresponding to the netCDF error status
!
! \newline
! General Usage:
!   call dshr_iocdf_handleErr(rCode,' check in xx call in subroutine yy ')
! \newline
! !REVISION HISTORY:
!     2005-Jan-31 - J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_handleErr(rCode, str)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN) ,intent (in) :: rCode
   character(*),intent (in) :: str

!EOP

   !----- formats -----
   character(*),parameter :: F00     = "('(dshr_iocdf_handleErr) ',4a)" 
   
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (rCode /= nf90_noerr) then
      write(6,F00) "netCDF error: ",trim(nf90_strerror(rCode))
      call dshr_iocdf_abort(str)
   end if

end subroutine dshr_iocdf_handleErr

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_setAbort -- Set local dshr_iocdf abort flag
!
! !DESCRIPTION:
!     Set local dshr_iocdf abort flag, true = abort, false = print and continue
! \newline
! General Usage:
!    call shr\_ncread\_setAbort(.false.)
! \newline
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_setAbort(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical,intent(in) :: flag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('dshr_iocdf_setAbort') "
  character(*),parameter :: F00     = "('(dshr_iocdf_setAbort) ',a) "

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

  doabort = flag

end subroutine dshr_iocdf_setAbort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_setDebug -- Set local dshr_iocdf debug level
!
! !DESCRIPTION:
!     Set local dshr_iocdf debug level, 0 = production
! \newline
! General Usage:
!    call shr\_ncread\_setDebug(2)
! \newline
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_setDebug(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: iflag

!EOP

  !--- local ---

  !--- formats ---
  character(*),parameter :: subName = "('dshr_iocdf_setDebug') "
  character(*),parameter :: F00     = "('(dshr_iocdf_setDebug) ',a) "

!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------

  debug = iflag

end subroutine dshr_iocdf_setDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_iocdf_abort -- local abort call
!
! !DESCRIPTION:
!     local abort call
! \newline
! General Usage:
!    call shr\_ncread\_abort(' ERROR in subroutine xyz ')
! \newline
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_iocdf_abort(string)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),optional,intent(IN) :: string

!EOP

  !--- local ---
  character(CL)           :: lstring
  character(*),parameter  :: subName =   "(dshr_iocdf_abort)"
  character(*),parameter  :: F00     = "('(dshr_iocdf_abort) ',a)"

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

end subroutine dshr_iocdf_abort

!===============================================================================
!===============================================================================

end module dshr_iocdf
