!===============================================================================
! SVN $Id: dshr_nml.F90 1049 2006-05-25 22:34:10Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_nml.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_nml -- defines and reads namelist common to all data models
!
! !DESCRIPTION:
!    Defines and reads a namelist common to all data models
!
! !REVISION HISTORY:
!    2005-Jun-27 - B. Kauffman - restricted to dshr namelist type and input
!    2004-Dec-15 - J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------


module dshr_nml 

   use shr_sys_mod   ! shared system calls
   use dshr_kind     ! kinds for strong typing

   implicit none

   private ! default private

! !PUBLIC TYPES:

   public :: dshr_nml_nmlType     ! namelist info

! !PUBLIC MEMBER FUNCTIONS:

   public :: dshr_nml_readNml     ! read dshr common namelist variables

! !PUBLIC DATA MEMBERS:

   integer(IN),parameter,public :: dshr_nml_nStreamsMax = 30 ! max number of streams
   integer(IN),parameter,public :: dshr_nml_nVectorsMax = 30 ! max number of vectors

!EOP

   type dshr_nml_nmlType            ! nml info common to all data models
      character(CS) :: restType     ! run type: initial, continue, branch
      character(CL) :: restPfile    ! IC data pointer file
      character(CL) :: restBfile    ! IC data file for branch runs, bundle data
      character(CL) :: restSfile    ! IC data file for branch runs, stream data
      character(CS) :: restFreq     ! restart file frequency
      integer(IN)   :: restDate     ! start date for initial or branch runs
      integer(IN)   :: restNday     ! restart file interval    (nday option
      integer(IN)   :: restOdate    ! restart file offset date (nday option)
      character(CS) :: caseName     ! case name
      character(CS) :: caseDesc     ! case description
      character(CL) :: dataMode     ! flags physics options wrt input data
      character(CS) :: domainFile   ! file   containing domain info
      character(CL) :: streamr      ! dlnd only: runoff stream description file name
      character(CL) :: streams(dshr_nml_nStreamsMax) ! stream description file names
      character(CL) :: vectors(dshr_nml_nVectorsMax) ! define vectors to vector map
      character(CS) :: mssDir       ! MSS output file directory
      character(CS) :: mssPass      ! MSS output file password
      character(CS) :: mssOpts      ! MSS output file mswrite options
      integer(IN)   :: mssRtpd      ! MSS output file MS retention period
      logical       :: mssRmlf      ! T => rm local copy when no longer needed
      integer(IN)   :: infoDbug     ! dbug level: 0=lowest, 3=highest
      logical       :: infoTime     ! T => print extra timing info
      real(R8)      :: infoSleep    ! simulates time to advance an active model
      integer(IN)   :: ncpl         ! number of msg pairs per day
   end type dshr_nml_nmlType

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: dshr_nml_readNml -- read control namelist
!
! !DESCRIPTION:
!     Reads namelist common to all data models
!
! !REVISION HISTORY:
!     2004-Dec-15 - J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine dshr_nml_readNml(nml,file,rc)

   implicit none

! !INPUT/OUTPUT PARAMETERS: 

   type(dshr_nml_nmlType),intent(out) :: nml  ! namelist data data-type
   character(*)          ,intent(in)  :: file ! file to read namelist from
   integer(IN),optional  ,intent(out) :: rc   ! return code

!EOP 

   integer(IN)    :: rCode         ! return code
   integer(IN)    :: nUnit         ! fortran i/o unit number

   !----- temporary/local namelist vars to read int -----
   character(CS)  :: restType      ! run type: initial, continue, branch
   character(CL)  :: restPfile     ! IC data pointer file
   character(CL)  :: restBfile     ! IC data file for branch runs, bundle data
   character(CL)  :: restSfile     ! IC data file for branch runs, stream data
   character(CS)  :: restFreq      ! restart file frequency
   integer(IN)    :: restDate      ! start date for initial or branch runs
   integer(IN)    :: restNday      ! restart file interval    (nday option
   integer(IN)    :: restOdate     ! restart file offset date (nday option)
   character(CL)  :: caseName      ! case name
   character(CL)  :: caseDesc      ! case description
   character(CL)  :: dataMode      ! flags physics options wrt input data
   character(CL)  :: domainFile    ! file   containing domain info
   character(CL)  :: streamr       ! dlnd only: roff stream description file name
   character(CL)  :: streams(dshr_nml_nStreamsMax) ! stream description file names
   character(CL)  :: vectors(dshr_nml_nvectorsMax) ! define vectors to vector map
   character(CL)  :: mssDir        ! MSS output file directory
   character(CS)  :: mssPass       ! MSS output file password
   character(CS)  :: mssOpts       ! MSS output file mswrite options
   integer(IN)    :: mssRtpd       ! MSS output file MS retention period
   logical        :: mssRmlf       ! T => rm local file after mswrite
   integer(IN)    :: infoDbug      ! dbug level: 0=lowest, 3=highest
   logical        :: infoTime      ! T => print extra timing info
   real(R8)       :: infoSleep     ! simulates time to advance an active model
   integer(IN)    :: ncpl          ! number of msg pairs per day

   !----- define namelist -----
   namelist / dshr_nml / &
        restType        &
      , restPfile       &
      , restBfile       &
      , restSfile       &
      , restFreq        &
      , restDate        &
      , restNday        &
      , restOdate       &
      , caseName        &
      , caseDesc        &
      , dataMode        &
      , domainFile      &
      , streamr         &
      , streams         &
      , vectors         &
      , mssDir          &
      , mssPass         &
      , mssOpts         &
      , mssRtpd         &
      , mssRmlf         &
      , infoDbug        &
      , infoTime        &
      , infoSleep       &
      , ncpl            

   !----- formats -----
   character(*),parameter :: subName = "(dshr_nml_readNml) " 
   character(*),parameter ::   F00 = "('(dshr_nml_readNml) ',8a)" 
   character(*),parameter ::   F01 = "('(dshr_nml_readNml) ',a,i6,a)" 
   character(*),parameter ::   F90 = "('(dshr_nml_readNml) ',58('-'))"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   !----------------------------------------------------------------------------
   ! set default values for namelist vars
   !----------------------------------------------------------------------------
   restType    = 'continue'  
   restPfile   = 'null'       
   restBfile   = 'null'       
   restSfile   = 'null'       
   restFreq    = 'never'        
   restDate    = 10101           
   restNday    = 365             
   restOdate   = 0101            
   caseName    = 'null'       
   caseDesc    = 'null'        
   dataMode    = 'COPYALL'        
   domainFile  = 'null'        
   streamr     = 'null'
   streams(:)  = 'null'
   vectors(:)  = 'null'
   mssDir      = 'null:'
   mssPass     = ' '
   mssOpts     = 'nowait,nomail'
   mssRtpd     = 180
   mssRmlf     = .false.
   infoDbug    = 1        
   infoTime    = .false.        
   infoSleep   = 0.0_R8
   ncpl        = 24             

   !----------------------------------------------------------------------------
   ! document default namelist values
   !----------------------------------------------------------------------------
   write(6,F90)
   write(6,F00) "Namelist values BEFORE reading namelist file (default values)..."
   write(6,dshr_nml)
   write(6,F90)

   !----------------------------------------------------------------------------
   ! read input namelist
   !----------------------------------------------------------------------------
   write(6,F00) 'reading input namelist file: ',trim(file)
   nUnit = shr_sys_ioUnit() ! get unused fortran i/o unit number
   open (nUnit,file=trim(file),status="old",action="read")
   read (nUnit,nml=dshr_nml,iostat=rCode)
   close(nUnit)
   if (rCode > 0) then
      write(6,F01) 'ERROR: reading input namelist, iostat=',rCode
      call shr_sys_abort(subName//": namelist read error")
   end if

   !----------------------------------------------------------------------------
   ! document namelist values
   !----------------------------------------------------------------------------
   write(6,F90)
   write(6,F00) "Namelist values AFTER reading namelist file..."
   write(6,dshr_nml)
   write(6,F90)

   !----------------------------------------------------------------------------
   ! fix variables as necessary
   !----------------------------------------------------------------------------
   restType = adjustL(restType)
   if      ( trim(restType) == "initial" ) then
      ! this is OK
   else if ( trim(restType) == "startup" ) then
      restType = "initial"  ! startup => initial
   else if ( trim(restType) == "continue") then
      ! this is OK
   else if ( trim(restType) == "branch"  ) then
      ! this is OK
   else 
      write(6,F00) "ERROR: invalid restType = ",trim(restType)
      call shr_sys_abort(subName//": invalid restType = "//trim(restType))
   end if

   !----------------------------------------------------------------------------
   ! copy temporary/local namelist vars into data structure
   !----------------------------------------------------------------------------
   nml%restType    = restType        
   nml%restPfile   = restPfile       
   nml%restBfile   = restBfile       
   nml%restSfile   = restSfile       
   nml%restFreq    = restFreq        
   nml%restDate    = restDate        
   nml%restNday    = restNday        
   nml%restOdate   = restOdate       
   nml%caseName    = caseName        
   nml%caseDesc    = caseDesc        
   nml%dataMode    = dataMode  
   nml%domainFile  = domainFile
   nml%streamr     = streamr
   nml%streams(:)  = streams(:)    
   nml%vectors(:)  = vectors(:)    
   nml%mssDir      = mssDir          
   nml%mssPass     = mssPass         
   nml%mssOpts     = mssOpts         
   nml%mssRtpd     = mssRtpd         
   nml%mssRmlf     = mssRmlf         
   nml%infoDbug    = infoDbug        
   nml%infoTime    = infoTime        
   nml%infoSleep   = infoSleep       
   nml%nCpl        = nCpl             

end subroutine dshr_nml_readNml

!===============================================================================
!===============================================================================

end module dshr_nml

