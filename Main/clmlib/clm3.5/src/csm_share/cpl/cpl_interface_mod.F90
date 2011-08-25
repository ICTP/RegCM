!===============================================================================
! SVN $Id: cpl_interface_mod.F90 3380 2007-03-06 05:42:19Z robj $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/cpl/cpl_interface_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: cpl_interface_mod -- General model-coupler interaction.
!
! !DESCRIPTION:
!     This module represents a major subsystem of cpl6.
!
!     This module contains the highest level routines for communication
!     between models and the Coupler.
!
!     These routines present component models with an API using mostly native
!     Fortran 90 datatypes such as real an integer scalars and arrays.
!     Only one derived datatype, the {\it contract} appears in the arguments.
!     This satisifes a Coupler requirement to present a simple interface
!     to coupled model programmers.  Most of the routines in this module are
!     wrappers to other cpl6 subroutines.
!     
!     These routines are all that is  necessary for a component model to connect to,
!     and exchange data with, version 6 of the CCSM Coupler.
!
! !REVISION HISTORY:
!     2003-Jan-15 - T. Craig - change ibuf to infobuf datatype module
!     2002-Dec-05 - T. Craig - changed call from cpl_coupling to cpl_contract
!     2002-Sep-10 - T. Craig - abstracted functionality into cpl_coupling_mod.F90
!     2001-Aug-16 - B. Kauffman - reorganized code according to arch document.
!     2001-Mar-20 - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

module cpl_interface_mod

! !USES:

   use mct_mod         ! mct interface
   use cpl_comm_mod        ! mpi/mph communicator info
   use cpl_fields_mod      ! coupler/model data field indicies
   use cpl_bundle_mod      ! defines bundle
   use cpl_domain_mod      ! defines domain
   use cpl_infobuf_mod     ! defines infobuf
   use cpl_contract_mod    ! defines contract
   use cpl_kind_mod        ! defines cpl kinds
   use cpl_control_mod, only: dbug=>cpl_control_infoDBug
   use cpl_const_mod, only:  cpl_const_rearth2
   use shr_sys_mod         ! share system routines
   use shr_timer_mod       ! share timer routines
   use shr_mpi_mod         ! mpi layer

   implicit none

   private   ! except

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: cpl_interface_init
   public :: cpl_interface_finalize
   public :: cpl_interface_ibufSend
   public :: cpl_interface_ibufRecv
   public :: cpl_interface_infobufSend
   public :: cpl_interface_infobufRecv
   public :: cpl_interface_contractSend
   public :: cpl_interface_contractRecv
   public :: cpl_interface_contractInit
   public :: cpl_interface_contractIndex  ! return index of a field in contract
   public :: cpl_interface_contractField  ! return field of an index in contract
   public :: cpl_interface_contractNumatt ! return number of fields in contract
   public :: cpl_interface_dbugSet   ! set this module's internal dbug level

   interface cpl_interface_ibufSend; module procedure cpl_interface_infobufSend; end interface
   interface cpl_interface_ibufRecv; module procedure cpl_interface_infobufRecv; end interface

! !PUBLIC DATA MEMBERS:

  ! none

!EOP

   !--- module variables ---
   character(*),parameter :: modName = "cpl_interface_mod"

   integer(IN),save  :: timer01,timer02,timer03,timer04    ! timers
   integer(IN),save  :: timer11,timer12,timer13,timer14    ! timers

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_init -- initialize the coupling/mpi environment.
!
! !DESCRIPTION:
!    Wrapper routine for {\tt cpl\_comm\_init}.  Calls {\it MPI\_Init}
!    and reports model {\tt name} to the coupled system.  Returns
!    an {\it MPI\_Communicator} {\tt comm} for use in the calling model.
!
!    {\tt name} must be one of the component names in {\tt cpl\_fields\_mod}.
! 
! !REVISION HISTORY:
!     2001-Mar-20 - T. Craig, B. Kauffman, R. Jacob -- first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_init(name,comm)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: name ! name of component name
   integer(IN) ,intent(out) :: comm ! communicator group for component

!EOP

   integer(IN)          :: n    ! generic loop index

   character(*),parameter :: subName = '(cpl_interface_init) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call cpl_comm_init(name, comm)

end subroutine cpl_interface_init

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_finalize -- terminate the coupling/mpi environment.
!
! !DESCRIPTION:
!    Calls {\it MPI\_Finalize()} and disengages the model {\tt cname} from the CCSM.
!    {\tt cname} must be one of the component names in {\tt cpl\_fields\_mod}.
! 
! !REVISION HISTORY:
!    2001-mmm-dd -
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_finalize(cname)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in) :: cname  ! component name
   integer(IN)             :: rcode  ! return code

!EOP

   character(*),parameter :: subName = '(cpl_interface_finalize) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   call shr_mpi_finalize(subName//" MPI finalize")

end subroutine cpl_interface_finalize

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractInit -- Initialize contract
!
! !DESCRIPTION:
!     Initialize the {\tt contract} between model {\tt my\_name} and
!     model {\tt other\_name}.  This is a wrapper to {\tt cpl\_contract\_init}.
!     Only two sets of optional arguments are allowed.  A Model uses the ``send''
!     form and should provide {\tt buf}, which contains grid information
!     such as latitude and longitude values (see {\tt cpl\_contract\_initSend}).
!     The Coupler uses the ``recv'' form and does not provide {\tt buf}.
!     {\tt bunname} and {\tt fields} are used to initialize the {\it bundle} in
!     the output {\tt contract}.
!
! 
! !REVISION HISTORY:
!    2002-Jul-30 - T.Craig -- prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_contractInit(contract,my_name,other_name,fields,ibufi,buf,ibufr,bunname,decomp)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_contract)   ,intent(out)   :: contract   ! contract
   character(*)         ,intent(in)    :: my_name    ! my component name
   character(*)         ,intent(in)    :: other_name ! other component name
   character(*)         ,intent(in)    :: fields     ! fields char string for bun
   integer(IN) ,optional,intent(inout) :: ibufi(:)   ! info buffer ints
   real(R8)    ,optional,intent(in)    :: buf(:,:)   ! data buffer
   real(R8)    ,optional,intent(inout) :: ibufr(:)   ! info buffer reals
   character(*),optional,intent(in)    :: bunname
   integer(IN) ,optional,intent(in)    :: decomp     ! decomposition type

!EOP

   !--- local ---
   character(CL)          :: bn  ! bundle name
   character(*),parameter :: subName = '(cpl_interface_contractInit) '
   integer(IN)            :: decomp_type

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   decomp_type = 1
   if (present(decomp)) then
     decomp_type = decomp
   endif

   if (present(ibufi)) then
     contract%infobuf%ibuf = ibufi
   endif

   if (present(ibufr)) then
     contract%infobuf%rbuf = ibufr
   endif

   if (present(buf)) then
      call cpl_contract_Init('send',contract,my_name,other_name,buf,decomp=decomp_type)
   else
      call cpl_contract_Init('recv',contract,my_name,other_name,decomp=decomp_type)
   endif

   if (present(bunname)) then
     bn = bunname
   else
     bn = "undef"
   endif
   call cpl_bundle_init(contract%bundle,bn,fields,contract%domain)

end subroutine cpl_interface_contractInit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_infobufSend -- send an infobuffer using arrays
!
! !DESCRIPTION:
!     Load an {\it infobuffer} using the input arrays {\tt ibufi} and
!     {\tt ibufr} and send it to the root processor of the model with
!     name {\tt cname}.
! 
! !REVISION HISTORY:
!    2001-Aug-16 -
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_infobufSend(cname,ibufi,ibufr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)       :: cname      ! component name
   integer(IN),optional,intent(inout)    :: ibufi(:)   ! info buffer ints
   real(R8)   ,optional,intent(inout)    :: ibufr(:)   ! info buffer reals

!EOP

   integer(IN)             :: pid          ! mpi process id
   integer(IN),parameter   :: tag=1002     ! mpi msg tag
   type(cpl_infobuf)       :: infobuf      ! local info buffer

   character(*),parameter :: subName = '(cpl_interface_infobufSend) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   !--- identify cpl_comm_wrld pid for other component's pe0 ---
   if (cname == cpl_fields_atmname) then
      pid = cpl_comm_wrld_pe0_atm
   else if (cname == cpl_fields_icename) then
      pid = cpl_comm_wrld_pe0_ice
   else if (cname == cpl_fields_lndname) then
      pid = cpl_comm_wrld_pe0_lnd
   else if (cname == cpl_fields_ocnname) then
      pid = cpl_comm_wrld_pe0_ocn
   else if (cname == cpl_fields_cplname) then
      pid = cpl_comm_wrld_pe0_cpl
   else
      write(6,*) subName,'ERROR: this should never happen'
      write(6,*) subName,'unrecognized cname = ',cname
   endif

   if (present(ibufi)) then
     infobuf%ibuf = ibufi
   endif

   if (present(ibufr)) then
     infobuf%rbuf = ibufr
   endif

   !--- send ibuf (tell cpl how big local & global component data is) ---
   if (cpl_comm_comp_pid == 0) then
     call cpl_infobuf_send(infobuf,pid,tag,cpl_comm_wrld)
   endif


end subroutine cpl_interface_infobufSend

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_infobufRecv -- receive an infobuffer.
!
! !DESCRIPTION:
!     Receive an {\it infobuffer} from the model with name {\tt cname}
!     and return its integer contents in the {\tt ibufi} array and
!     the real contents in the {\tt ibufr} array.
! 
! !REVISION HISTORY:
!    2001-mmm-dd - 
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_infobufRecv(cname,ibufi,ibufr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)  :: cname  ! component name
   integer(IN),optional,intent(out) :: ibufi(cpl_infobuf_ibufSize) ! info-buffer
   real(R8)   ,optional,intent(out) :: ibufr(cpl_infobuf_rbufSize) ! info-buffer

!EOP

   integer(IN) :: pid                     ! mph process ID
   integer(IN),parameter :: pid0 = 0      ! component's root pid
   integer(IN),parameter :: tag = 1002    ! mpi msg tag
   type(cpl_infobuf)     :: infobuf       ! local infobuf

   character(*),parameter :: subName = '(cpl_interface_infobufRecv) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   if (cname == cpl_fields_atmname) then
      pid = cpl_comm_wrld_pe0_atm
   else if (cname == cpl_fields_icename) then
      pid = cpl_comm_wrld_pe0_ice
   else if (cname == cpl_fields_lndname) then
      pid = cpl_comm_wrld_pe0_lnd
   else if (cname == cpl_fields_ocnname) then
      pid = cpl_comm_wrld_pe0_ocn
   else if (cname == cpl_fields_cplname) then
      pid = cpl_comm_wrld_pe0_cpl
   else
      write(6,*) subName,'ERROR: this should never happen'
      write(6,*) subName,'unrecognized cname = ',cname
   endif

   !--- recv info-buffer data ---
   if (cpl_comm_comp_pid == pid0) then
     call cpl_infobuf_recv(infobuf,pid,tag,cpl_comm_wrld)
   endif
   call cpl_infobuf_bcast(infobuf,pid0,cpl_comm_comp)

   if (present(ibufi)) then
     ibufi = infobuf%ibuf
   endif

   if (present(ibufr)) then
     ibufr = infobuf%rbuf
   endif

end subroutine cpl_interface_infobufRecv

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractSend -- send information in a contract.
!
! !DESCRIPTION:
!     Send data in the {\it infobuffer} and {\it bundle} of the input {\tt contract}
!     to the component {\tt cname} (which could be a model or the Coupler).  If
!     optional arguments {\tt ibufi} and/or {\tt ibufr} are present, that information
!     will be placed in the {\tt contract}'s {\it infobuffer} and sent instead.
!     If optional argument {\tt buf} is present, then that data will be placed in
!     in the {\tt contract}'s {\it bundle} and sent to {\tt cname}.
!
!     In CCSM3, the Coupler calls this with only {\tt cname} and {\tt contract}
!     while Models use the simple {\tt ibufi} and {\tt buf} arguments. 
! 
! !REVISION HISTORY:
!    2001-mmm-dd - 
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_contractSend(cname,contract,ibufi,buf,ibufr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)    :: cname      ! component name
   type(cpl_contract)  ,intent(inout) :: contract   ! data buffer and domain
   integer(IN),optional,intent(inout) :: ibufi(:)   ! info buffer ints
   real(R8)   ,optional,intent(inout) :: ibufr(:)   ! info buffer reals
   real(R8)   ,optional,intent(inout) ::  buf(:,:)  ! data buffer

!EOP

   logical,save :: first_call = .true.    ! first time in subroutine
   integer(IN) :: AVsiz                   ! size of a field in the AttrVect
   integer(IN) :: AVnum                   ! number of fields in AttrVect
   integer(IN) :: i,j,n                   ! dummy variables
   integer(IN) :: pid                     ! mpi process ID
   integer(IN),parameter :: tag = 1003    ! msg tag

   character(*),parameter :: subName = '(cpl_interface_contractSend) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   !--- setup timers ---
   if (first_call) then
     first_call = .false.
     call shr_timer_get(timer01,'comp_send total')
     call shr_timer_get(timer02,'comp_send reorder')
     call shr_timer_get(timer03,'comp_send cpl_contract_Send')
     call shr_timer_get(timer04,'comp_send diagnostics')
   endif

   call shr_timer_start(timer01)

   !--- copy data into contract ---
   if (present(ibufi)) then
     contract%infobuf%ibuf = ibufi
   endif

   if (present(ibufr)) then
     contract%infobuf%rbuf = ibufr
   endif

   if (present(buf)) then
     AVnum = mct_aVect_nRAttr(contract%bundle%data)
     AVsiz = mct_aVect_lsize(contract%bundle%data)

     if (AVsiz /= size(buf,1) .or. AVnum /= size(buf,2)) then
       write(6,*) subName,' ERROR in buffer/contract buffer size:',  &
                  trim(contract%bundle%name)
       write(6,*) subName,' sending buffer, size:',size(buf,1),size(buf,2)
       write(6,*) subName,'contract buffer, size:',AVsiz,AVnum
       call shr_sys_flush(6)
     endif

     call shr_timer_start(timer02)
     !--- reorder bundle data as per mct data structure ---
#ifdef CPP_VECTOR
     do j=1,AVnum
     do i=1,AVsiz
#else
     do i=1,AVsiz
     do j=1,AVnum
#endif
       contract%bundle%data%rattr(j,i)=buf(i,j)
     enddo
     enddo
     call shr_timer_stop(timer02)
   endif

   !--- diagnostics ---
   if ( dbug >= 2) then
     call shr_timer_start(timer04)
     call cpl_bundle_gsum(contract%bundle,contract%bundle%dom%lGrid,'aream', &
                                          contract%bundle%dom%lGrid,'mask', &
                          scalar=cpl_const_rearth2,istr='send '//trim(cname))
     call shr_timer_stop(timer04)
   endif

   call shr_timer_start(timer03)
   
   !--- identify pid for component pe0 ---
   if (cname == cpl_fields_atmname) then
      pid = cpl_comm_wrld_pe0_atm
   else if (cname == cpl_fields_icename) then
      pid = cpl_comm_wrld_pe0_ice
   else if (cname == cpl_fields_lndname) then
      pid = cpl_comm_wrld_pe0_lnd
   else if (cname == cpl_fields_ocnname) then
      pid = cpl_comm_wrld_pe0_ocn
   else if (cname == cpl_fields_cplname) then
      pid = cpl_comm_wrld_pe0_cpl
   else
      write(6,*) subName,'ERROR: this should never happen'
      write(6,*) subName,'unrecognized cname = ',cname
   endif

   !--- send data ---
   call cpl_contract_Send(contract,cpl_comm_comp_pid,cpl_comm_comp,pid)

   call shr_timer_stop(timer03)
   call shr_timer_stop(timer01)

end subroutine cpl_interface_contractSend

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractRecv -- Receive information in a contract.
!
! !DESCRIPTION:
!     Recive data into the {\it infobuffer} and {\it bundle} of the input/output {\tt contract}
!     from the component {\tt cname} (which could be a model or the Coupler).  If
!     optional arguments {\tt ibufi} and/or {\tt ibufr} are present, that information
!     is copied out of the {\tt contract}'s {\it infobuffer} and placed into those arrays.
!     If optional argument {\tt buf} is present, then data  in the {\tt contract}'s
!     {\it bundle} after the receive will be returned in {\tt buf}.
!
!     In CCSM3, the Coupler calls this with only {\tt cname} and {\tt contract}
!     while Models use the simple {\tt ibufi} and {\tt buf} arguments. 
! 
! !REVISION HISTORY:
!    2001-mmm-dd - 
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_contractRecv(cname,contract,ibufi,buf,ibufr)

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   character(*)        ,intent(in)    :: cname      ! component name
   type(cpl_contract)  ,intent(inout) :: contract   ! data buffer and domain
   integer(IN),optional,intent(out)   :: ibufi(:)   ! info buffer ints
   real(R8)   ,optional,intent(out)   :: ibufr(:)   ! info buffer reals
   real(R8)   ,optional,intent(out)   ::  buf(:,:)  ! data buffer

!EOP

   logical,save:: first_call = .true.     ! first time in subroutine
   integer(IN) :: AVsiz                   ! size of a field in the AttrVect
   integer(IN) :: AVnum                   ! number of fields in AttrVect
   integer(IN) :: i,j,n                   ! dummy variables
   integer(IN) :: pid                     ! generic pid
   integer(IN),parameter :: pid0 = 0      ! component's root pid
   integer(IN),parameter :: tag = 1003    ! mpi msg tag

   character(*),parameter :: subName = '(cpl_interface_contractRecv) '

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!   write(6,*) subName,trim(cname)

   if (first_call) then
     call shr_timer_get(timer11,'comp_recv total')
     call shr_timer_get(timer12,'comp_recv reorder')
     call shr_timer_get(timer13,'comp_recv cpl_contract_Recv')
     call shr_timer_get(timer14,'comp_recv diagnostics')
     first_call = .false.
   endif

   if (present(buf)) then
     AVnum = mct_aVect_nRAttr(contract%bundle%data)
     AVsiz = mct_aVect_lsize(contract%bundle%data)
     if (AVsiz /= size(buf,1) .or. AVnum /= size(buf,2)) then
       write(6,*) subName,' ERROR in buffer/contract buffer size'
       write(6,*) subName,' sending buffer, size:',size(buf,1),size(buf,2)
       write(6,*) subName,'contract buffer, size:',AVsiz,AVnum
       call shr_sys_flush(6)
     endif
   endif

   call shr_timer_start(timer11)
   call shr_timer_start(timer13)

   !--- identify pid for component pe0 ---
   if (cname == cpl_fields_atmname) then
      pid = cpl_comm_wrld_pe0_atm
   else if (cname == cpl_fields_icename) then
      pid = cpl_comm_wrld_pe0_ice
   else if (cname == cpl_fields_lndname) then
      pid = cpl_comm_wrld_pe0_lnd
   else if (cname == cpl_fields_ocnname) then
      pid = cpl_comm_wrld_pe0_ocn
   else if (cname == cpl_fields_cplname) then
      pid = cpl_comm_wrld_pe0_cpl
   else
      write(6,*) subName,'ERROR: this should never happen'
      write(6,*) subName,'unrecognized cname = ',cname
   endif

   !--- receive data ---
   call cpl_contract_Recv(contract,cpl_comm_comp_pid,cpl_comm_comp,pid)
   call shr_timer_stop(timer13)

   !--- diagnostics ---
   if ( dbug >= 2) then
     call shr_timer_start(timer14)
     call cpl_bundle_gsum(contract%bundle,contract%bundle%dom%lGrid,'aream', &
                          scalar=cpl_const_rearth2,istr='recv '//trim(cname))
     call shr_timer_stop(timer14)
   endif

   !--- copy data out of contract ---
   if (present(ibufi)) then
     ibufi = contract%infobuf%ibuf
   endif

   if (present(ibufr)) then
     ibufr = contract%infobuf%rbuf
   endif

   if (present(buf)) then
     call shr_timer_start(timer12)
     !--- reorder bundle data as per mct data structure ---
#ifdef CPP_VECTOR
     do j=1,AVnum
     do i=1,AVsiz
#else
     do i=1,AVsiz
     do j=1,AVnum
#endif
       buf(i,j)=contract%bundle%data%rattr(j,i)
     enddo
     enddo
     call shr_timer_stop(timer12)
   endif

   call shr_timer_stop(timer11)

end subroutine cpl_interface_contractRecv

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractIndex -- return index in a contract
!
! !DESCRIPTION:
!     Return the index of an attribute in a contract
!     Facilitates field query logic in component models.
! 
! !REVISION HISTORY:
!    2003-01-31 - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

integer function cpl_interface_contractIndex(contract, item, perrWith, dieWith)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_contract)       ,intent(in) :: contract   ! data buffer and domain
   character(len=*)         ,intent(in) :: item       ! component name
   character(len=*),optional,intent(in) :: perrWith
   character(len=*),optional,intent(in) :: dieWith

!EOP

   integer(IN) :: indx

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(perrWith).and. (.not.present(dieWith))) then
      indx = mct_aVect_indexRA(contract%bundle%data,item,perrWith=perrWith)

   else if(.not.present(perrWith).and. present(dieWith)) then
      indx = mct_aVect_indexRA(contract%bundle%data,item,dieWith=dieWith)

   else if(present(perrWith).and. present(dieWith)) then
      indx = mct_aVect_indexRA(contract%bundle%data,item, &
                perrWith=perrWith,dieWith=dieWith)
   else
      indx = mct_aVect_indexRA(contract%bundle%data,item)
   endif

   cpl_interface_contractIndex =indx

end function cpl_interface_contractIndex

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractField -- return field in a contract
!
! !DESCRIPTION:
!     Return the field of an attribute in a contract
!     Facilitates field query logic in component models.
! 
! !REVISION HISTORY:
!    2005-09-20 - T. Craig
!
! !INTERFACE: ------------------------------------------------------------------

!character function cpl_interface_contractField(contract, index, perrWith, dieWith)
!character function cpl_interface_contractField(contract, index)
subroutine cpl_interface_contractField(contract, index, cfld)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_contract)       ,intent(in) :: contract   ! data buffer and domain
   integer(IN)              ,intent(in) :: index      ! field index
   character(len=*)         ,intent(out):: cfld       ! output field
!  character(len=*),optional,intent(in) :: perrWith
!  character(len=*),optional,intent(in) :: dieWith

!EOP

   type(mct_string) :: sfld

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call mct_aVect_getRList(sfld,index,contract%bundle%data)

   cfld = mct_string_toChar(sfld)
   call mct_string_clean(sfld)

end subroutine cpl_interface_contractField

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_contractNumAtt -- return number of attributes in 
!  a contract
!
! !DESCRIPTION:
!     Return the number of attributes in a contract. 
!     Facilitates dynamic array allocation in component models.
! 
! !REVISION HISTORY:
!    2003-01-31 - B. Kauffman 
!
! !INTERFACE: ------------------------------------------------------------------

integer function cpl_interface_contractNumatt(contract)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(cpl_contract),intent(in) :: contract   ! data buffer and domain

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   cpl_interface_contractNumatt = mct_aVect_nRAttr(contract%bundle%data)

end function cpl_interface_contractNumAtt

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: cpl_interface_dbugSet  -- set this module's internal debug level.
!
! !DESCRIPTION:
!    Set this module's internal debug level: 0,1,2,3 (lowest to highest). 
!    If debug level is 2 or greater, each call to {\tt cpl\_interface\_contractSend}
!    or {\tt cpl\_interface\_contractRecv} will output the global sum of
!    all data sent/received to stdout using {\tt cpl\_bundle\_gsum}.
!
! !REVISION HISTORY:
!    2003-Jan-21 - B. Kauffman, initial version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cpl_interface_dbugSet(level)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in) :: level  ! requested debug level

!EOP

   !----- local -----
   integer(IN):: newLevel  ! new debug level

   !----- formats -----
   character(*),parameter :: F00 = "('(cpl_interface_dbugSet) ',a,i1,a)"

!-------------------------------------------------------------------------------
!  Set module's internal debug level: 0,1,2, or 3 (lowest to highest)
!-------------------------------------------------------------------------------

   newLevel = max(0,min(3,level))

   !--- correct invalid level values ---
   if (newLevel /= level) then
      write(6,F00) 'WARNING: level ',level,' not in {0,1,2,3} '
      write(6,F00) 'WARNING: resetting level to ',newLevel
   end if

   if (dbug>0 .OR. dbug/=newLevel) write(6,F00) 'set debug level to ',newLevel

   dbug = newLevel

end subroutine cpl_interface_dbugSet

!===============================================================================
!===============================================================================

end module cpl_interface_mod
