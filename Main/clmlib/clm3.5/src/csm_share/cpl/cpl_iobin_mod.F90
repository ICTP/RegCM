!===============================================================================
! SVN $Id: cpl_iobin_mod.F90 3380 2007-03-06 05:42:19Z robj $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_iobin_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_iobin_mod -- create, write-to, or read a binary data file.
!
! !DESCRIPTION:
!    This module creates, writes-to, and reads data in a {\it bundle} to/from a file
!    in machine-dependent binary format.  This module is intended to support 
!    the restart sub-system of cpl6.
!
! !REMARKS:
!    This module could be used to create binary data files for any
!    purpose.  The binary data file format is self-describing and extensible and
!    thus might provide an alternative to netCDF files.
!
! !REVISION HISTORY:
!     2002-nov-08 - B. Kauffman - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

module cpl_iobin_mod

! !USES:

   use mct_mod       ! mct interface
   use cpl_comm_mod      ! mpi/mph communicator info
   use cpl_fields_mod    ! coupler/model data field indices
   use cpl_bundle_mod    ! defines bundle
   use cpl_domain_mod    ! defines domain
   use cpl_kind_mod      ! defines F90 kinds
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug
   use shr_sys_mod       ! share system routines
   use shr_date_mod      ! defines date data-type
   use shr_mpi_mod       ! layer on MPI

   implicit none

   private ! except

! !PUBLIC TYPES:

   ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_iobin_create     ! COLLECTIVE: create a new file (an empty file)
   public :: cpl_iobin_open       ! COLLECTIVE: open a named file
   public :: cpl_iobin_close      ! COLLECTIVE: close an open file
   public :: cpl_iobin_appendBun  ! COLLECTIVE: add  bundle data to   a file
   public :: cpl_iobin_readBun    ! COLLECTIVE: read bundle data from a file
   public :: cpl_iobin_appendReal ! ROOT PID  : add  real array data   to a file
   public :: cpl_iobin_readReal   ! ROOT PID  : read real array data from a file
   public :: cpl_iobin_readDate   ! COLLECTIVE: read data date from a file
 
! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   !--- module variables ---
   integer(IN),parameter  :: pid0 = 0 ! root process pid = zero
   character(*),parameter :: modName = "cpl_iobin_mod"

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iobin_create -- create a new file.
!
! !DESCRIPTION:
!    Open a file with name {\tt fName} and write a small header
!    of 6 character strings each with length CL.  Use Fortran
!    unformatted write.  If optional argument {\tt desc} is present,
!    it will be included in the header.
!
! !REVISION HISTORY:
!    2002-Nov-08 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iobin_create(fName,desc)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)          :: fName   ! file name
   character(*),intent(in),optional :: desc    ! description string

!EOP

   !----- local -----
   integer(IN )  :: fid      ! file ID (unit number)
   character(CL) :: format   ! format string
   character(CL) :: str      ! generic string
   character(CL) :: dstr     ! F90 date string (wrt date_and_time)
   character(CL) :: tstr     ! F90 time string (wrt date_and_time)
   character(CL) :: comment  ! comment/description string
   character(CL),parameter :: svnID = "SVN $Id"
   logical       :: open     ! true => the unit is connected to a file

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_create) '
   character(*),parameter :: F00 = "('(cpl_iobin_create) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iobin_create) ',a,i2)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe only
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0 ) return

   !----------------------------------------------------------------------------
   ! find an unused unit number
   !----------------------------------------------------------------------------
   do fid=10,99
      inquire(fid,opened=open)
      if (.not. open) exit
   end do
   if (open) then
      write(6,F00) 'ERROR: couldn''t find unused unit number'
      call shr_sys_abort(subName//": all file unit numbers in use?")
   end if
   if (dbug>2) write(6,F01) 'using unit number ',fid

   !----------------------------------------------------------------------------
   ! create an empty file
   !----------------------------------------------------------------------------
   if (dbug>1) write(6,F00) 'creating new file ',trim(fName)

   format   = "format-str256:name:F90_date_and_time:caseDesc:svnId"
   str      = "Header for "// trim(fName)
   call date_and_time(dstr,tstr)
   comment  = 'description: (no description provided) '
   if (present(desc)) comment = desc

   open (fid,file=fName,form="UNFORMATTED",status="REPLACE")
   write(fid) format
   write(fid) str      
   write(fid) dstr,tstr
   write(fid) comment
   write(fid) svnId  
   close(fid)

   if (dbug>2) write(6,F00) 'created  new file ',trim(fName)

end subroutine cpl_iobin_create

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iobin_open -- open an existing file.
!
! !DESCRIPTION:
!    Open the pre-existing file with name {\tt fName} and assign it
!    the unit number {\tt fid}.  Also read the header information and
!    write it to stdout.
!
! !REVISION HISTORY:
!    2002-Nov-08 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iobin_open(fName,fid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)   :: fName   ! file name
   integer(IN) ,intent(out)  :: fid     ! file ID (file unit number)

!EOP

   !----- local -----
   character(CL) :: format   ! format    string
   character(CL) :: name     ! file name string
   character(CL) :: dstr     ! F90 date  string (wrt date_and_time)
   character(CL) :: tstr     ! F90 time  string (wrt date_and_time)
   character(CL) :: comment  ! comment   string
   character(CL) :: svnId    ! SVN Id    string
   character(CL) :: str      ! generic   string
   logical       :: open     ! true => unit number is in use

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_open) '
   character(*),parameter :: F00 = "('(cpl_iobin_open) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iobin_open) ',3a,i2)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe only
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0 ) return

   !----------------------------------------------------------------------------
   ! find an unused unit number
   !----------------------------------------------------------------------------
   do fid=10,99
      inquire(fid,opened=open)
      if (.not. open) exit
   end do
   if (open) then
      write(6,F00) 'ERROR: couldn''t find unused unit number'
      call shr_sys_abort(subName//": all units open?")
   end if

   !----------------------------------------------------------------------------
   ! open an existing file & print out the header info
   !----------------------------------------------------------------------------
   if (dbug>1) write(6,F01) 'open file ',trim(fName),', unit = ',fid

   !--- assume file was written with this format: --- 
!  format = "format-str256:name:F90_date_and_time:caseDesc:svnId"

   open(fid,file=fName,form="UNFORMATTED",status="OLD")
   read(fid) format     ! string len=256
   read(fid) name       ! string len=256
   read(fid) dstr,tstr  ! string len=256 * 2
   read(fid) comment    ! string len=256
   read(fid) svnID      ! string len=256

   write(6,F00) "format: ",trim(format)
   write(6,F00) "title: ",trim(name)
   str = 'File created: '                                   &
   &     //dstr(1:4)//'-'//dstr(5:6)//'-'//dstr(7:8)//' '   &
   &     //tstr(1:2)//':'//tstr(3:4)//':'//tstr(5:6)
   write(6,F00) trim(str)
   write(6,F00) "comment: ",trim(comment)
   write(6,F00) "SVN Id : ",trim(svnId)

   call shr_sys_flush(6)

end subroutine cpl_iobin_open

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iobin_close -- close an open file.
!
! !DESCRIPTION:
!    Call {\tt close} on file with unit number {\tt fid}.
!
! !REVISION HISTORY:
!    2002-Nov-08 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iobin_close(fid)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in)  :: fid    ! file ID (file unit number)

!EOP

   !----- local -----

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_close) '
   character(*),parameter :: F00 = "('(cpl_iobin_close) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iobin_close) ',a,i2)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe only
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0 ) return

   !----------------------------------------------------------------------------
   ! close it
   !----------------------------------------------------------------------------
   if (dbug>1) write(6,F01) 'close file, unit = ',fid
   close(fid)

end subroutine cpl_iobin_close

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iobin_appendBun -- add bundle data to an existing file.
!
! !DESCRIPTION:
!    Append the data in {\it bundle} {\tt bun} to the pre-existing file
!    with the unit number {\tt fid}.  Also write the date contained in
!    {\tt date} to the file.
!    {\tt fid} must be a valid fortran unit number for an open file.
!
!    All processors call this function and the root node will {\it MPI\_gather}
!    the data and write it to the file.
!
! !REMARKS:    
!    Domain data associated with the bundle is not written
!    but this functionality could be added, if desired.
!    The file format utilizes an extensible, self-describing data format.
!
!
! !REVISION HISTORY:
!    2002-Nov-08 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iobin_appendBun(fid,date,bun)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN)          ,intent(inout) :: fid     ! file ID
   type(shr_date)       ,intent(in)    :: date    ! model date
   type(cpl_bundle)     ,intent(in)    :: bun     ! bundle

!EOP

   !----- local -----
   character(CL)       :: format   ! format description text string
   character(CL)       :: name     ! field name text string
   type(mct_aVect) :: gData    ! global/gathered bundle data
   integer(IN)         :: nflds    ! number of fields in a bundle
   integer(IN)         :: npts     ! number of points in a field
   integer(IN)         :: rcode    ! return code
   integer(IN)         :: n,k      ! generic index
   real(R8),allocatable:: data(:)  ! temporary real vector
   integer(IN)         :: cDate    ! model date: yyyymmdd
   integer(IN)         :: sec      ! model date: seconds
   logical             :: open     ! true => the unit is connected to a file

   type(mct_string) :: mctStr  ! temporary mct/mpeu string var
   character(CL)        :: dSuffix ! domain suffix
   character(CL)        :: dName   ! domain name
   character(CL)        :: aName   ! attribute (field) name
   character(CL)        :: bName   ! bundle name


   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_appendBun) '
   character(*),parameter :: F00 = "('(cpl_iobin_appendBun) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iobin_appendBun) ',3a,i4)"
   character(*),parameter :: F02 = "('(cpl_iobin_appendBun) ',a,i6,f10.2)"

!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------

   if (dbug>0) write(6,F00) 'writing data for bundle = ',trim(bun%name)

   !----------------------------------------------------------------------------
   ! do global gather to create non-distributed aVect (*part* of a bundle)
   !----------------------------------------------------------------------------
   call mct_aVect_gather(bun%data,gData,bun%dom%gsMap,pid0,cpl_comm_comp,rcode)
   
   !----------------------------------------------------------------------------
   ! from here on, everything is done by the root pe
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid /= pid0) return

   !----------------------------------------------------------------------------
   ! append the data
   !----------------------------------------------------------------------------
   npts = mct_aVect_lsize (gData)
   allocate(data(npts))

   call shr_date_getCDate(date,cDate,sec) ! model time associate with data

   do n=1,mct_aVect_nRAttr(bun%data)

      !--- construct field var name ---
      if (dbug>2) write(6,F00) '* constructing name'
      call mct_list_get(mctStr,n,bun%data%rList)
      aName   = mct_string_toChar(mctStr) ! aVect var name
      call mct_string_clean(mctStr)
      bName   = bun%name                      ! bundle name
      dSuffix = bun%dom%suffix                ! domain suffix
!!!   name = trim(bName) // "_" // trim(aName) // "_" // trim(dSuffix) ! old
      name = trim(bName) // "_" // trim(aName) 

      !--- get data out of bundle ---
      if (dbug>2) write(6,F00) '* extracting data from aVect'
      k = mct_aVect_indexRA(gData,trim(aName),perrWith=subName)
      data = gData%rAttr(k,:)

      !--- append data to output file ---
      if (dbug>1) write(6,F00) '* writing data for variable = ',trim(name)
      call cpl_iobin_appendReal(fid,date,name,data,rcode)

   end do

   deallocate(data)
   call mct_aVect_clean(gData)

   if (dbug>2) write(6,F00) 'DONE adding data to file, EXIT routine'
   call shr_sys_flush(6)

end subroutine cpl_iobin_appendBun

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iobin_readBun -- read bundle data from a file.
!
! !DESCRIPTION:
!    Read data from file with unit number {\tt fid} and return in in
!    the {\it bundle} {\tt bun}.  Argument {\tt date} is currently
!    ignored.  
!    Data is read on node 0 and scattered using the information in
!    the {\it domain} associated with {\tt bun}.  On return, {\tt bun}
!    contains the data for points local to the calling processor.
!
! !REVISION HISTORY:
!    2002-Nov-08 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iobin_readBun(fid,date,bun)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN)     ,intent(in)    :: fid   ! file ID
   type(shr_date)  ,intent(inout) :: date  ! desired model date, file's date?
   type(cpl_bundle),intent(inout) :: bun   ! bundle to read

!EOP

   !----- local -----
   character(CL)       :: format   ! format description text string
   character(CL)       :: name     ! field name text string
   type(mct_aVect) :: gData    ! global/gathered bundle data
   integer(IN)         :: npts     ! number of points in a field
   integer(IN)         :: rcode    ! return code
   integer(IN)         :: n,k      ! generic index
   real(R8),allocatable:: data(:)  ! temporary real vector
   integer(IN)         :: cDate    ! model date: yyyymmdd
   integer(IN)         :: sec      ! model date: seconds

   type(mct_string) :: mctStr  ! temporary mct/mpeu string var
   character(CL)        :: dSuffix ! domain suffix
   character(CL)        :: dName   ! domain name
   character(CL)        :: aName   ! attribute (field) name
   character(CL)        :: bName   ! bundle name 
   character(CL)        :: vName   ! variable name in data file

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_readBun) '
   character(*),parameter :: F00 = "('(cpl_iobin_readBun) ',4a)"
   character(*),parameter :: F01 = "('(cpl_iobin_readBun) ',a,i10)"
   character(*),parameter :: F02 = "('(cpl_iobin_readBun) ',a,2es12.3)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (dbug>0) write(6,F00) 'reading data for bundle = ',trim(bun%name)

   !----------------------------------------------------------------------------
   ! read the data - master process only
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid == pid0) then

      if (dbug>2) write(6,F00) '* create temporary global data & aVect'
      npts = bun%dom%n                             ! global size of domain
      call mct_aVect_init(gData,bun%data,npts) ! it aVect, copy rList
      allocate(data(npts))                         ! read data into this array
      if (dbug>2) write(6,F01) '* allocated temp data array, size = ',npts

      do n=1,mct_aVect_nRAttr(bun%data)

         !--- construct field var name ---
         call mct_list_get(mctStr,n,bun%data%rList)
         aName   = mct_string_toChar(mctStr) ! aVect var name
         call mct_string_clean(mctStr)
         bName   = bun%name                      ! bundle name
         dSuffix = bun%dom%suffix                ! domain suffix
!!!      vName = trim(bName) // "_" // trim(aName) // "_" // trim(dSuffix) ! old
         vName = trim(bName) // "_" // trim(aName) 
   
         if (dbug>1) write(6,F00) '* read data for variable ',trim(aName)

         !--- locate and read a real array from data file ---
         call cpl_iobin_readReal(fid,date,vName,data,rcode)
         if (rcode /= 0) then
            write(6,F00) 'ERROR: reading from file, var = ',trim(vName)
            call shr_sys_abort(subName//": ERROR reading variable from file")
         end if
         if (dbug>2) then
            write(6,F02) '* min/max of data is ',minval(data),maxval(data)
         end if

         !--- add variable data to the gathered aVect ---
         k = mct_aVect_indexRA(gData,trim(aName),perrWith=subName)
         gData%rAttr(k,:) = data(:)
   
      end do
   end if

   !----------------------------------------------------------------------------
   ! scatter the global aVect 
   !----------------------------------------------------------------------------
   if (dbug>2) write(6,F00) '* scattering aVect'
   call mct_aVect_clean(bun%data)
   call mct_aVect_scatter(gData,bun%data,bun%dom%gsMap,pid0,cpl_comm_comp,rcode)

   if (cpl_comm_comp_pid == pid0) then
      if (dbug>2) write(6,F00) '* clean  temporary global data & aVect'
      deallocate(data)
      call mct_aVect_clean(gData)
   end if

   call shr_sys_flush(6)

end subroutine cpl_iobin_readBun

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iobin_appendReal -- Append real array data to file
!
! !DESCRIPTION:
!    Append data in real array {\tt data} to the alredy-open file
!    with unit number {\tt fid}.  Include in file the name of
!    the data {\tt vName} and the date {\tt date}.  {\tt rcode}
!    is 0 if successful and 1 if file is not open.
!
! !REMARKS:
!    This routine is run on root PE only.
!
! !REVISION HISTORY:
!    2003-Mar-07 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iobin_appendReal(fid,date,vName,data,rcode)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN)   ,intent(in)  :: fid     ! fortan unit number of open data file
   type(shr_date),intent(in)  :: date    ! desired model date
   character(*)  ,intent(in)  :: vName   ! name of var to read
   real(R8)      ,intent(in)  :: data(:) ! the data
   integer(IN)   ,intent(out) :: rcode   ! return code

!EOP

   !----- local -----
   character(CL) :: fmtStr    ! format text string
   character(CL) :: rcdStr    ! name-of-record text string
   integer(IN)   :: cDate,sec ! model date (yyyymmdd) & seconds
   logical       :: open      ! true iff file is opened

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_appendReal) '
   character(*),parameter :: F00 = "('(cpl_iobin_appendReal) ',4a)"
   character(*),parameter :: F02 = "('(cpl_iobin_appendReal) ',a,2es12.3)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rcode = 0       ! no errors yet

   !--- file is assumed to be open ---
   inquire(fid,OPENED=open)
   if (.not. open) then
      write(6,F00) 'ERROR: unit number must be assigned to an open file'
      rcode = 1
      return
   end if

   call shr_date_getCDate(date,cDate,sec) ! model time associate with data
   rcdStr = vName                         ! converts str length

   !----------------------------------------------------------------------------
   ! append data to output file
   !----------------------------------------------------------------------------
   if (dbug>2) write(6,F00) '* writing data for variable = ',trim(vName)

   fmtStr = "format-char256:name:model_time:data "
   write(fid) fmtStr          ! format description
   write(fid) rcdStr          ! name of data
   write(fid) cDate,sec       ! date/seconds of data
   write(fid) data            ! the data itself

end subroutine cpl_iobin_appendReal

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iobin_readReal -- read real array data from file
!
! !DESCRIPTION:
!    Read data with name {\tt vName} into real array {\tt data} 
!    from the alredy-open file with unit number {\tt fid}.
!    Argument {\tt date} is currently not used.  {\tt rcode}
!    is 0 if successful and 1 if file is not open or variable is not found.
!
! !REMARKS:
!    This routine is run on root PE only.
!
! !REVISION HISTORY:
!    2002-Nov-08 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iobin_readReal(fid,date,vName,data,rcode)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN)   ,intent(in)  :: fid     ! fortan unit number of open data file
   type(shr_date),intent(in)  :: date    ! desired data date
   character(*)  ,intent(in)  :: vName   ! name of var to read
   real(R8)      ,intent(out) :: data(:) ! array to put data into
   integer(IN)   ,intent(out) :: rcode   ! return code

!EOP

   !----- local -----
   character(CL) :: fmtStr    ! file format text string
   character(CL) :: rcdStr    ! name-of-record text string
   integer(IN)   :: cDate,sec ! model date (yyyymmdd) & seconds
   integer(IN)   :: n         ! generic index
   integer(IN)   :: nColon    ! number of colons in a format string
   logical       :: found     ! true iff variable was found in file
   logical       :: open      ! true iff file is opened

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_readReal) '
   character(*),parameter :: F00 = "('(cpl_iobin_readReal) ',4a)"
   character(*),parameter :: F02 = "('(cpl_iobin_readReal) ',a,2es12.3)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   rcode = 0  ! no errors yet

   !--- file is assumed to be open ---
   inquire(fid,OPENED=open)
   if (.not. open) then
      write(6,F00) 'ERROR: unit number must be assigned to an open file'
      rcode = 1
      return
   end if

   !----------------------------------------------------------------------------
   ! search for & read requested variable
   !----------------------------------------------------------------------------
   if (dbug>2) write(6,F00) 'searching for variable: ',trim(vName)

   rewind(fid)     ! start reading from beginning of file 
   found = .false. ! varable has not been found yet
   do while (.not. found) 
      read(fid,END=999) fmtStr  ! format description string
      read(fid,END=999) rcdStr  ! name-of-record string (variable name)
      if (trim(rcdStr) == trim(vName)) then
         !--- found requested variable => read it ---
         if (dbug>2) then
            write(6,F00) '  found/read  variable: ',trim(vName)
            call shr_sys_flush(6)
         end if
         found = .true.
         read(fid,END=999)      ! cDate,sec: date associated with data
         read(fid) data         ! the actual data
         if (dbug>2) then
            write(6,F02) '* min/max of data is ',minval(data),maxval(data)
         end if
      else
         !--- found some other variable, skip over it ---
         if (dbug>2) then
            write(6,F00) '  found/skip  variable: ',trim(rcdStr)
            call shr_sys_flush(6)
         end if
         found = .false.
         nColon = cpl_iobin_countChar(fmtStr,":")
         do n=2,nColon
            read(fid,END=999)
         end do
      end if
   end do

 999  continue

   if (.not. found) then
      write(6,F00) 'ERROR: couldn''t find variable in file, var = ',trim(vName)
      rcode = 1
   else
      if (dbug>2) write(6,F00) 'SUCCESS! exit routine'
   end if

end subroutine cpl_iobin_readReal

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_iobin_readDate - read data date from a file
!
! !DESCRIPTION:
!    Search the already open file with unit number {\tt fid} and find
!    the date string in the header.  If found, return date information
!    in output arguments {\tt cDate} and {\tt sec}.  {\tt rcode} is
!    0 if successful and 1 if file is not open.  This routine will
!    abort if the date string is not found.
!    Date is read on node 0 and broadcast to all processors calling this
!    routine.
!
! !REMARKS:
!    Assumes the date for all data is the same, hence returns 1st date found.
!
! !REVISION HISTORY:
!    2002-Dec-13 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_iobin_readDate(fid,cDate,sec,rcode)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in)  :: fid     ! fortan unit number of open data file
   integer(IN),intent(out) :: cDate   ! coded date of data in file
   integer(IN),intent(out) :: sec     ! seconds corresponding to cDate
   integer(IN),intent(out) :: rcode   ! return code

!EOP

   !----- local -----
   character(CL) :: fmtStr ! format text string
   integer(IN)   :: n      ! generic index
   integer(IN)   :: nColon ! number of colons in a format string
   logical       :: found  ! true iff variable was found in file
   logical       :: open   ! true iff file is opened

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_readDate) '
   character(*),parameter :: F00 = "('(cpl_iobin_readDate) ',4a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   found = .false. ! varable has not been found yet
   rcode = 0       ! no errors yet

   !----------------------------------------------------------------------------
   ! read the data - master process only
   !----------------------------------------------------------------------------
   if (cpl_comm_comp_pid == pid0) then

      !--- file is assumed to be open ---
      inquire(fid,OPENED=open)
      if (.not. open) then
         write(6,F00) 'ERROR: unit number must be assigned to an open file'
         rcode = 1
         return
      end if

      !--- start reading from beginning of file ---
      rewind(fid)

      if (dbug>2) write(6,F00) 'searching for cDate & sec'

      !--- search for a variable with an associated model date ---
      do while (.not. found) 
         read(fid,END=999) fmtStr  ! format string
         if (dbug>2) write(6,F00) '* fmtStr = ',trim(fmtStr)
         if ( index(fmtStr,"model_time") > 0 ) then
            !--- found it => read it ---
            found = .true.
            if (dbug>2) write(6,F00) '* found/read record with model_time'
            n = index(fmtStr,"model_time")
            nColon = cpl_iobin_countChar(fmtStr(1:n),":")
            do n=1,nColon-1
               read(fid,END=999)
            end do
            read(fid) cDate,sec     ! the actual data
         else
            found = .false.
            if (dbug>2) write(6,F00) '* found/skip record without model_time'
            nColon = cpl_iobin_countChar(fmtStr,":")
            do n=1,nColon
               read(fid,END=999)
            end do
         end if
      end do

 999  continue

      if (.not. found) then
         write(6,F00) 'ERROR: couldn''t find model date in file'
         call shr_sys_abort(subName//": couldn't find date in file")
      end if

   end if

   !----------------------------------------------------------------------------
   ! broadcast the model date
   !----------------------------------------------------------------------------
    call shr_mpi_bcast(cDate,cpl_comm_comp,subName//" MPI in cDate bcast")
    call shr_mpi_bcast(sec,  cpl_comm_comp,subName//" MPI in sec bcast")
  

end subroutine cpl_iobin_readDate

!===============================================================================
!=== private/internal module routine ===========================================
!===============================================================================

integer function cpl_iobin_countChar(str,char)

   implicit none

   character(*),intent(in) :: str   ! string to search
   character(*),intent(in) :: char  ! char to search for

   !----- local -----
   integer(IN) :: count    ! counts occurances of char
   integer(IN) :: n        ! generic index

   !----- formats -----
   character(*),parameter :: subName = '(cpl_iobin_countChar) '
   character(*),parameter :: F00 = "('(cpl_iobin_countChar) ',4a)"

!-------------------------------------------------------------------------------
!  count number of occurances of a character in a string
!-------------------------------------------------------------------------------

   count = 0
   do n = 1, len_trim(str)
      if (str(n:n) == char(1:1)) count = count + 1 
   end do
   cpl_iobin_countChar = count

end function cpl_iobin_countChar

!===============================================================================
!===============================================================================

end module cpl_iobin_mod
