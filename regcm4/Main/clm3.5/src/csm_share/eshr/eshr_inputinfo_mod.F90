!===============================================================================
! SVN $Id$
! SVN $URL$
!===============================================================================

!BOP ===========================================================================
!
! !MODULE: eshr_inputinfo_mod --- Convert inputinfo object to/from ESMF State
!
! !DESCRIPTION:
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2006-Aug-24 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------

module eshr_inputinfo_mod

#ifdef SEQ_ESMF
! !USES:
   use shr_inputinfo_mod, only: shr_inputinfo_initType
   use shr_kind_mod,      only: SHR_KIND_CL, SHR_KIND_CS
   use eshr_rc_mod,       only: eshr_rc_check
   use eshr_estate_mod,   only: eshr_estate_printAttributes
   use ESMF_Mod

   implicit none

   private    ! default private

! ! PUBLIC TYPES:

  ! None

! ! PUBLIC MEMBER FUNCTIONS

  public :: eshr_inputinfo_info2EState ! Put CCSMInit object info on ESMF State
  public :: eshr_inputinfo_EState2Info ! Get information from ESMF State put into CCSMInit

! !PUBLIC DATA MEMBERS:



!EOP

   ! Private member functions:

   ! Private data members: 

   character(len=*), parameter :: CCSMInitEStateName = "CCSMInitEState"

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_inputinfo_info2EState -- put input CCSM-info values on ESMF State
!   
! !DESCRIPTION:
!   
!  Put the values on the input CCSM-info object onto the input ESMF State.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_inputinfo_info2EState( CCSMInit, eState, CCSMInitEState, print )
    use shr_inputinfo_mod, only: shr_inputinfo_initGetData
    use eshr_estate_mod,   only: eshr_estate_FLog2ELog
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(shr_inputinfo_initType), intent(IN)    :: CCSMInit       ! input object
    type(ESMF_State),             intent(INOUT) :: eState         ! Output ESMF State
    type(ESMF_State),             intent(INOUT) :: CCSMInitEState ! CCSMInit state object
    logical, optional,            intent(IN)    :: print          ! If print state at end
!EOP
    !----- local -----
    character(SHR_KIND_CL)  :: start_type     ! Type of startup
    character(SHR_KIND_CL)  :: case_desc      ! Long description of this case
    character(SHR_KIND_CS)  :: case_name      ! Short case identification
    character(SHR_KIND_CL)  :: restart_file   ! Full archive path to drv restart file
    character(SHR_KIND_CL)  :: archive_dir    ! Directory to archive to...
    character(SHR_KIND_CS)  :: logFilePostFix ! postfix for output log files
    character(SHR_KIND_CS)  :: outPathRoot    ! root for output log files
    character(SHR_KIND_CL)  :: mss_wpass      ! MSS write password
    logical                 :: atm_adiabatic  ! No surface models and atm adiabatic mode
    logical                 :: atm_ideal_phys ! No surface models and atm ideal-physics
    logical                 :: aqua_planet    ! No ice/lnd, analytic ocn, perpetual time
    logical                 :: brnch_retain_casename ! If branch should retain casename
    logical                 :: printIt        ! If should print state when done
    integer                 :: mss_irt        ! MSS retention period
    integer                 :: rc             ! ESMF error return code
!-------------------------------------------------------------------------------
! Notes: Just put the data that a sub-component would need. Not data needed at the
!        driver level only.
!-------------------------------------------------------------------------------
    
    if ( .not. present(print) )then
       printIt = .false.
    else
       printIt = print
    end if

    !------- Create a seperate state for the CCSMInit object -------------------

    CCSMInitEState = ESMF_StateCreate( CCSMInitEStateName, rc=rc )
    call eshr_rc_check( rc, "Error in creation of CCSMInit ESMF State" )

    !------  Add a long_name attribute to describe this state ------------------

    call ESMF_StateSetAttribute( CCSMInitEState, name="long_name", &
                                 value="shr_inputinfo_mod object ESMF State", rc=rc )
    call eshr_rc_check( rc, "Error adding long_name to CCSMInit State" )

    !------ Get the data from the CCSMInit object
    call shr_inputinfo_initGetData( CCSMInit, case_name=case_name,                &
                                    case_desc=case_desc, archive_dir=archive_dir, &
                                    mss_irt=mss_irt, mss_wpass=mss_wpass,         &
                                    start_type=start_type,                        &
                                    aqua_planet=aqua_planet,                      &
                                    atm_ideal_phys=atm_ideal_phys,                &
                                    brnch_retain_casename=brnch_retain_casename,  &
                                    atm_adiabatic=atm_adiabatic,                  &
                                    restart_file=restart_file,                    &
                                    logFilePostFix=logFilePostFix,                &
                                    outPathRoot=outPathRoot                       &
                                  )
    !
    !------ Add CCSMInit object data as attributes to this state
    !

    call ESMF_StateSetAttribute( CCSMInitEState, name="start_type", &
                                 value=start_type, rc=rc )
    call eshr_rc_check( rc, "Error adding start_type to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="case_name", &
                                 value=case_name, rc=rc )
    call eshr_rc_check( rc, "Error adding case_name to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="case_desc", &
                                 value=case_desc, rc=rc )
    call eshr_rc_check( rc, "Error adding case_desc to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="archive_dir", &
                                 value=archive_dir, rc=rc )
    call eshr_rc_check( rc, "Error adding archive_dir to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="mss_irt", &
                                 value=mss_irt, rc=rc )
    call eshr_rc_check( rc, "Error adding mss_irt to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="mss_wpass", &
                                 value=mss_wpass, rc=rc )
    call eshr_rc_check( rc, "Error adding mss_wpass to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="restart_file", &
                                 value=restart_file, rc=rc )
    call eshr_rc_check( rc, "Error adding restart_file to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="logFilePostFix", &
                                 value=logFilePostFix, rc=rc )
    call eshr_rc_check( rc, "Error adding logFilePostFix to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="outPathRoot", &
                                 value=outPathRoot, rc=rc )
    call eshr_rc_check( rc, "Error adding outPathRoot to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="aqua_planet", &
                                 value=eshr_estate_FLog2ELog(aqua_planet), rc=rc )
    call eshr_rc_check( rc, "Error adding aqua_planet to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="atm_ideal_phys", &
                                 value=eshr_estate_FLog2ELog(atm_ideal_phys), rc=rc )
    call eshr_rc_check( rc, "Error adding atm_ideal_phys to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="atm_adiabatic", &
                                 value=eshr_estate_FLog2ELog(atm_adiabatic), rc=rc )
    call eshr_rc_check( rc, "Error adding atm_adiabatic to CCSMInit State" )
    call ESMF_StateSetAttribute( CCSMInitEState, name="brnch_retain_casename", &
                                 value=eshr_estate_FLog2ELog(brnch_retain_casename), &
                                 rc=rc )
    call eshr_rc_check( rc, "Error adding brnch_retain_casename to CCSMInit State" )

    if ( printIt ) call eshr_estate_printAttributes( CCSMInitEState )
    !
    !------ Add the CCSMInit state to the input ESMF state ---------------------
    !

    call ESMF_StateAddState( eState, CCSMInitEState, rc=rc )
    call eshr_rc_check( rc, "Error in adding of CCSMInit ESMF State to input eState" )

END SUBROUTINE eshr_inputinfo_info2EState

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_inputinfo_EState2Info -- get input CCSM-info values from ESMF State
!   
! !DESCRIPTION:
!   
!  Get the values from the input ESMF State and put on the input CCSM-info object.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_inputinfo_eState2Info( eState, CCSMInit, print )
    use shr_inputinfo_mod, only: shr_inputinfo_initPutData, &
                                 shr_inputinfo_initSetDefaults
    use eshr_estate_mod,   only: eshr_estate_ELog2FLog
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),             intent(IN)  :: eState   ! Input ESMF State
    type(shr_inputinfo_initType), intent(OUT) :: CCSMInit ! Output object
    logical, optional,            intent(IN) :: print     ! If should print state when done
!EOP
    !----- local -----
    type(ESMF_State)        :: CCSMInitEState ! CCSMInit state object
    character(SHR_KIND_CL)  :: start_type     ! Type of startup
    character(SHR_KIND_CL)  :: case_desc      ! Long description of this case
    character(SHR_KIND_CS)  :: case_name      ! Short case identification
    character(SHR_KIND_CL)  :: restart_file   ! Full archive path to drv restart file
    character(SHR_KIND_CL)  :: archive_dir    ! Directory to archive to...
    character(SHR_KIND_CS)  :: logFilePostFix ! postfix for output log files
    character(SHR_KIND_CS)  :: outPathRoot    ! root for output log files
    character(SHR_KIND_CL)  :: mss_wpass      ! MSS write password
    character(SHR_KIND_CL)  :: string         ! Character string
    logical                 :: atm_adiabatic  ! No surface models and atm adiabatic mode
    logical                 :: atm_ideal_phys ! No surface models and atm ideal-physics
    logical                 :: aqua_planet    ! No ice/lnd, analytic ocn, perpetual time
    logical                 :: brnch_retain_casename ! If branch should retain casename
    logical                 :: printIt        ! If should print state when done
    integer                 :: mss_irt        ! MSS retention period
    integer                 :: rc             ! ESMF error return code
    type(ESMF_Logical)      :: ELogValue      ! ESMF logical value
!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( .not. present(print) )then
       printIt = .false.
    else
       printIt = print
    end if

    !------- Get the state for the CCSMInit object from the input EState -------

    call ESMF_StateGetState( eState, CCSMInitEStateName, CCSMInitEState, rc=rc)
    call eshr_rc_check( rc, "Error in getting CCSMInit ESMF State from input state" )

    if ( printIt ) call eshr_estate_printAttributes( CCSMInitEState )
    !
    !------ Get CCSMInit object data as attributes from this state
    !

    call ESMF_StateGetAttribute( CCSMInitEState, name="start_type", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting start_type from CCSMInit State" )
    start_type = string
    call ESMF_StateGetAttribute( CCSMInitEState, name="case_name", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting case_name from CCSMInit State" )
    case_name  = string
    call ESMF_StateGetAttribute( CCSMInitEState, name="case_desc", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting case_desc from CCSMInit State" )
    case_desc  = string
    call ESMF_StateGetAttribute( CCSMInitEState, name="archive_dir", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting archive_dir from CCSMInit State" )
    archive_dir  = string
    call ESMF_StateGetAttribute( CCSMInitEState, name="mss_irt", &
                                 value=mss_irt, rc=rc )
    call eshr_rc_check( rc, "Error getting mss_irt from CCSMInit State" )
    call ESMF_StateGetAttribute( CCSMInitEState, name="mss_wpass", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting mss_wpass from CCSMInit State" )
    mss_wpass  = string
    call ESMF_StateGetAttribute( CCSMInitEState, name="restart_file", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting restart_file from CCSMInit State" )
    restart_file  = string
    call ESMF_StateGetAttribute( CCSMInitEState, name="logFilePostFix", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting logFilePostFix from CCSMInit State" )
    logFilePostFix  = string
    call ESMF_StateGetAttribute( CCSMInitEState, name="outPathRoot", &
                                 value=string, rc=rc )
    call eshr_rc_check( rc, "Error getting outPathRoot from CCSMInit State" )
    outPathRoot  = string
    call ESMF_StateGetAttribute( CCSMInitEState, name="aqua_planet", &
                                 value=ELogValue, rc=rc )
    call eshr_rc_check( rc, "Error getting aqua_planet from CCSMInit State" )
    aqua_planet = eshr_estate_ELog2FLog(ELogValue)
    call ESMF_StateGetAttribute( CCSMInitEState, name="atm_ideal_phys", &
                                 value=ELogValue, rc=rc )
    call eshr_rc_check( rc, "Error getting atm_ideal_phys from CCSMInit State" )
    atm_ideal_phys = eshr_estate_ELog2FLog(ELogValue)
    call ESMF_StateGetAttribute( CCSMInitEState, name="atm_adiabatic", &
                                 value=ELogValue, rc=rc )
    call eshr_rc_check( rc, "Error getting atm_adiabatic from CCSMInit State" )
    atm_adiabatic = eshr_estate_ELog2FLog(ELogValue)
    call ESMF_StateGetAttribute( CCSMInitEState, name="brnch_retain_casename", &
                                 value=ELogValue, rc=rc )
    call eshr_rc_check( rc, "Error getting brnch_retain_casename from CCSMInit State" )
    brnch_retain_casename = eshr_estate_ELog2FLog(ELogValue)

    !------ Put the data on the CCSMInit object
    call shr_inputinfo_initSetDefaults( CCSMInit )  ! Set default values
    call shr_inputinfo_initPutData( CCSMInit, case_name=case_name,                &
                                    case_desc=case_desc, archive_dir=archive_dir, &
                                    mss_irt=mss_irt, mss_wpass=mss_wpass,         &
                                    start_type=start_type,                        &
                                    aqua_planet=aqua_planet,                      &
                                    atm_ideal_phys=atm_ideal_phys,                &
                                    brnch_retain_casename=brnch_retain_casename,  &
                                    atm_adiabatic=atm_adiabatic,                  &
                                    restart_file=restart_file,                    &
                                    logFilePostFix=logFilePostFix,                &
                                    outPathRoot=outPathRoot                       &
                                  )

END SUBROUTINE eshr_inputinfo_eState2Info

#endif

end module eshr_inputinfo_mod
