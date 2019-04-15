!===============================================================================
! SVN $Id: dshr_hubInfo.F90 231 2006-02-06 23:33:22Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_hubInfo.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_hubInfo - defines scalar flags exchanged with a CCSM hub
!
! !DESCRIPTION:
!    In addition to 2d field data, CCSM components exchange a list of scalar data
!    that are control flags and other incidental scalar data.  This module 
!    defines that list of data.  This list of scalars is called an "info-buffer".
!
! !REVISION HISTORY:
!    2005 Jun 27 - B. Kauffman, adapted from cpl_fields_mod.F90
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_hubInfo

! !USES:

   use dshr_kind      ! kinds for strong typing

   private  ! default private

! !PUBLIC TYPES:

   public :: dshr_hubInfo_infoBuffType ! control flags to/from hub

! !PUBLIC MEMBER FUNCTIONS:

   public :: dshr_hubInfo_init         ! set initial info values

! !PUBLIC DATA MEMBERS:

   !----------------------------------------------------------------------------
   ! index into integer-valued "info-buffer" data
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: dshr_hubInfo_iBuf_total    = 100 ! size of info-buffer

   integer(IN),parameter,public :: dshr_hubInfo_iBuf_rCode    =   1 ! error code
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_cDate    =   2 ! current date: yymmdd
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_sec      =   3 ! elapsed sec on date
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_nCpl     =   4 ! cpl comm's per day
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_nFields  =  10 ! number of fields sent
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_gSize    =  11 ! global size of field
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_lSize    =  12 ! local  size of field
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_giSize   =  13 ! size of global i-dimension
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_gjSize   =  14 ! size of global j-dimension
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_liSize   =  15 ! size of local  i-dimension
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_ljSize   =  16 ! size of local  j-dimension
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_stopEod  =  19 ! stop at end-of-day flag
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_stopNow  =  20 ! stop now flag
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_restEod  =  21 ! restart file at EOD
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_restNow  =  22 ! restart file now
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_histEod  =  23 ! history file at EOD
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_histNow  =  24 ! history file now
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_histAvg  =  25 ! create monthly avg data
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_diagEod  =  26 ! diagnostics at EOD
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_diagNow  =  27 ! diagnostics now
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_infoTim  =  28 ! activate timing info
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_infoBug  =  29 ! debug level
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_precAdj  =  31 ! precip adjustment factor (* 1.0e+6)
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_aShift   =  32 ! albedo calculation time shift
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_nBasins  =  33 ! number of active runoff basins
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_xAlbic   =  34 ! request extra albedo solar init msg
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_dead     =  37 ! non-0 <=> dead model
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_domain   =  40
   integer(IN),parameter,public :: dshr_hubInfo_iBuf_useRest  =  41 ! non-0 <=> use restart data sent to cpl

   !----------------------------------------------------------------------------
   ! index into real-valued "info-buffer" data
   !----------------------------------------------------------------------------

   integer(IN),parameter,public :: dshr_hubInfo_rbuf_total    =  50 ! size of real info-buffer

   integer(IN),parameter,public :: dshr_hubInfo_rbuf_spval    =   1 ! special value (indicates invalid data)
   integer(IN),parameter,public :: dshr_hubInfo_rbuf_eccen    =  10 ! Earth's eccentricity
   integer(IN),parameter,public :: dshr_hubInfo_rbuf_obliqr   =  11 ! Earth's Obliquity
   integer(IN),parameter,public :: dshr_hubInfo_rbuf_lambm0   =  12 ! longitude of perihelion at v-equinox
   integer(IN),parameter,public :: dshr_hubInfo_rbuf_mvelpp   =  13 ! Earth's Moving vernal equinox of orbit +pi

!EOP

   type dshr_hubInfo_infoBuffType
      integer(IN) :: iRecv(dshr_hubInfo_iBuf_total) ! integers received from huecv'd
      integer(IN) :: iSend(dshr_hubInfo_iBuf_total) ! integers sent to hub
      real(R8)    :: rRecv(dshr_hubInfo_rbuf_total) ! reals received from hub
      real(R8)    :: rSend(dshr_hubInfo_rbuf_total) ! reals sent to hub
   end type dshr_hubInfo_infoBuffType

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_hubInfo_init -- initialize hubInfo data values
!
! !DESCRIPTION:
!     Initialize hubInfo data values -- default values are zero
!
! !REVISION HISTORY:
!     2005-Sep-30 - B. Kauffman - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_hubInfo_init(hubInfo,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS: 

   type(dshr_hubInfo_infoBuffType),intent(out) :: hubInfo  ! data to init
   integer(IN)           ,optional,intent(out) :: rc       ! return code

!EOP 

   !----- formats -----
   character(*),parameter :: subName = "(dshr_hubInfo_init) " 
   character(*),parameter ::   F00 = "('(dshr_hubInfo_init) ',8a)" 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   write(6,F00) "initializing all data to zero"

   hubInfo%iRecv(:) = 0      ! integers received from hub
   hubInfo%iSend(:) = 0      ! integers sent to hub
   hubInfo%rRecv(:) = 0.0_r8 ! reals received from hub
   hubInfo%rSend(:) = 0.0_r8 ! reals sent to hub

end subroutine dshr_hubInfo_init

!===============================================================================
!===============================================================================

end module dshr_hubInfo

