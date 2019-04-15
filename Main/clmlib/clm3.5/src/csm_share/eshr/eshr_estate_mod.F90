!===============================================================================
! SVN $Id$
! SVN $URL$
!===============================================================================

!BOP ===========================================================================
!
! !MODULE: eshr_estate_mod --- Functions to work on ESMF State objects
!
! !DESCRIPTION:
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2006-Aug-24 - E. Kluzek - creation
!
! !INTERFACE: ------------------------------------------------------------------

module eshr_estate_mod

! !USES:
   use shr_kind_mod,   only: SHR_KIND_CS, SHR_KIND_CL, SHR_KIND_R8
   use shr_string_mod, only: shr_string_listGetName
   use eshr_rc_mod,    only: eshr_rc_check
   use shr_sys_mod,    only: shr_sys_abort
   use ESMF_Mod

   implicit none

   private    ! default private

#ifdef SEQ_ESMF
! ! PUBLIC TYPES:

  ! None

! ! PUBLIC MEMBER FUNCTIONS

  public :: eshr_estate_init2DFromList   ! Create State of 2D Field from list of names
  public :: eshr_estate_initDomain       ! Create State of domain fields
  public :: eshr_estate_getDataPointer   ! Get the data pointer from a list index
  public :: eshr_estate_setDataPointer   ! Set the data pointer from a list index
  public :: eshr_estate_checkGridsMatch  ! Check that grids match
  public :: eshr_estate_getStats         ! Print basic statistics on input ESMF State
  public :: eshr_estate_fieldsAreSame    ! Return true if fields on two states are same
  public :: eshr_estate_matchValue       ! Return true if any of the data match a value
  public :: eshr_estate_copyFields       ! Copy fields from one ESMF state to another
  public :: eshr_estate_reDistInit       ! Initialize a Redist. from 1 E-state to another
  public :: eshr_estate_reDist           ! Redistribute fields from one E-state to another
  public :: eshr_estate_mergeInit        ! Initialize a Merge  from 3 E-states to another
  public :: eshr_estate_merge            ! Merge fields from 3 E-states to another
  public :: eshr_estate_destroy          ! Destroy a state and the fields that it contains
  public :: eshr_estate_getFilename      ! Get the given filename from the input estate
  public :: eshr_estate_putFilename      ! Put the given filename on the output estate
  public :: eshr_estate_printAttributes  ! Print out the list of attributes on the state
  public :: eshr_estate_printFieldAttr   ! Print out attributes on each field on the state
  public :: eshr_estate_ELog2FLog        ! Convert ESMF_Logical to FORTRAN Logical
  public :: eshr_estate_FLog2ELog        ! Convert FORTRAN Logical to ESMF_Logical
  public :: eshr_estate_Egrids2DAreSame  ! Check if two grids have same physical coords
  public :: eshr_estate_getSharedFNames  ! Get list of fieldnames shared on each state
  public :: eshr_estate_getFieldNames    ! Get list of fieldnames from state
  public :: eshr_estate_getFieldsNotNeed ! Get list of fieldnames that are not needed
  public :: eshr_estate_markRestNotNeed  ! Mark the other fields NOT needed

! !PUBLIC DATA MEMBERS:

  public :: domainFieldList              ! List of domain fieldnames
  public :: ICE_IDX, OCN_IDX, LND_IDX, XAO_IDX  ! Indices

!  INTERFACE BLOCKS


!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: eshr_estate_initDomain - Create a domain state for a model component.
!
! !INTERFACE:
   interface eshr_estate_initDomain
! !PRIVATE MEMBER FUNCTIONS:
      module procedure eshr_estate_init1DDomain
! !DESCRIPTION:
! This allows creation of a 1D (for arbitrary decomp) and 2D domain (for blocked)
! NOTE: Right now we only have the arbitrary decomposition version -- later the blocked
!       version will be added.
!EOPI
   end interface eshr_estate_initDomain

!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: eshr_estate_getDataPointer -- Get the pointer to data from a list and index
!
! !INTERFACE:
   interface eshr_estate_getDataPointer
! !PRIVATE MEMBER FUNCTIONS:
      module procedure eshr_estate_get1DDataPointer
! !DESCRIPTION:
! This interfaces to get methods that return data pointers based on a list and an index.
! NOTE: Right now we only have a 1D version for arbitrary decomp, later a 2D can be added.
!
!EOPI
   end interface eshr_estate_getDataPointer


!------------------------------------------------------------------------------
!BOPI
! !IROUTINE: eshr_estate_setDataPointer -- Set the pointer to data from a list and index
!
! !INTERFACE:
   interface eshr_estate_setDataPointer
! !PRIVATE MEMBER FUNCTIONS:
      module procedure eshr_estate_set1DDataPointer
! !DESCRIPTION:
! This interfaces to set methods that return data pointers based on a list and an index.
! NOTE: Right now we only have a 1D version for arbitrary decomp, later a 2D can be added.
!
!EOPI
   end interface eshr_estate_setDataPointer


! !PRIVATE MEMBER FUNCTIONS:

   private :: eshr_estate_egridGetLatLon    ! Get latitude and longitude from ESMF grid
   private :: eshr_estate_latLonAreSame     ! Check if latitude/longitudes are same

!EOP

   ! Private data members:
   real(SHR_KIND_R8),  parameter  :: eps = 1.0e-12    ! Acceptable difference
   character(len=*),   parameter  :: eshr_estate_dataIsFlux  = "flux"
   character(len=*),   parameter  :: eshr_estate_dataIsState = "state"
   integer,            parameter  :: OUT_IDX = 1      ! Index for output data
   integer,            parameter  :: ICE_IDX = 2      ! Index for sea-ice input data
   integer,            parameter  :: LND_IDX = 3      ! Index for land input data
   integer,            parameter  :: OCN_IDX = 4      ! Index for ocean input data
   integer,            parameter  :: XAO_IDX = 5      ! Index for atm/ocn flux input data
   !--- List of fields on domain states and their units and long names ----------
   character(len=*),   parameter  :: AreaName    = "area"
   character(len=*),   parameter  :: MaskName    = "mask"
   character(len=*),   parameter  :: MaxFracName = "maxfrac"
   character(len=*),   parameter  :: DomainFieldList   = AreaName//":"// &
                                                         MaskName//":"// &
                                                         MaxFracName
   character(len=*),   parameter  :: DomainFieldUnits  = "m^2:unitless:unitless"
   character(len=*),   parameter  :: DomainFieldLNames = &
                              "Area of grid cell:"// &
                              "If grid cell is active for this model:"// & 
                              "Maximum fraction that model can use in this grid cell"
   logical,            parameter :: OnlyUseNeededData = .false.

CONTAINS
!===============================================================================

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_init2DFromList -- Initialize State from a list of 2D fields
!   
! !DESCRIPTION:
!   
!  Initialize an ESMF State from a list of 2D fields.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_init2DFromList(  Grid, FieldList, FieldLNames,             &
                                        FieldUnits, notNeeded, initValue, eState, &
                                        usePTR )
    use shr_const_mod,   only: SHR_CONST_SPVAL
    use shr_string_mod,  only: shr_string_listGetNum, shr_string_listIsValid
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Grid)                        :: Grid        ! Grid fields are on
    character(len=*),           intent(IN) :: FieldList   ! Colen delimited list of fields
    character(len=*), optional, intent(IN) :: FieldLNames ! : delimited list of long_names
    character(len=*), optional, intent(IN) :: FieldUnits  ! : delimited list of units
    character(len=*), optional, intent(IN) :: NotNeeded   ! list of fields not needed
    real(SHR_KIND_R8),optional, intent(IN) :: initValue   ! Initial value
    type(ESMF_State),        intent(INOUT) :: eState      ! Output ESMF State
    logical,          optional, intent(IN) :: usePTR      ! Use model data pointer
!EOP
    !----- local -----
    character(len=*),parameter :: subName =   "(eshr_estate_init2DFromList) "
    integer                    :: i, nf        ! Indices
    integer                    :: rc           ! Return code
    character(len=SHR_KIND_CS) :: FieldName    ! Field name
    character(len=SHR_KIND_CL), allocatable :: FieldNames(:) ! Array of field names
    character(len=SHR_KIND_CS) :: long_name    ! Field long name
    character(len=SHR_KIND_CS) :: units        ! units of field
    character(len=SHR_KIND_CL) :: dataType     ! Type of data flux or state
    type(ESMF_ArraySpec)       :: arraySpec2d  ! 2D Array spec
    type(ESMF_ArraySpec)       :: arraySpec1d  ! 1D Array spec (for arbitrary dist lists)
    type(ESMF_GridStorage)     :: gridStore    ! Type of grid storage
    real(SHR_KIND_R8), pointer :: Data2D(:,:)  ! Pointer to 2D data
    real(SHR_KIND_R8), pointer :: Data1D(:)    ! Pointer to 1D data
    real(SHR_KIND_R8)          :: Value        ! Initial value to set data to
    logical                    :: ArbDecomp    ! If on arbitrary decomposition or not
    type(ESMF_Field), pointer  :: Fields(:)    ! Output fields in the state
    integer                    :: nFields      ! Output Number of fields
    type(ESMF_NeededFlag)      :: need         ! If field needed or not
    logical                    :: useDataPTR   ! Use model data pointer

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    if ( present(usePTR) )then
       useDataPTR = usePTR
    else
       useDataPTR = .false.
    end if
    !--- Check that field list is valid and get number of fields in list -------
    if ( .not. shr_string_listIsValid( FieldList ) )then
       call shr_sys_abort( subName//"ERROR: FieldList input is not a valid list" )
    end if
    nFields= shr_string_listGetNum( FieldList )
    !--- Check long_names ---
    if ( present(FieldLNames) )then
       if ( .not. shr_string_listIsValid( FieldLNames) )then
          call shr_sys_abort( subName//"ERROR: FieldLNames input is not a valid list" )
       end if
       nf = shr_string_listGetNum( FieldLNames)
       if ( nf /= nFields )then
          call shr_sys_abort( subName// &
                            "ERROR: number of FieldLNames do not match FieldList" )
       end if
    end if
    !--- Check units ---
    if ( present(FieldUnits) )then
       if ( .not. shr_string_listIsValid( FieldUnits) )then
          call shr_sys_abort( subName//"ERROR: FieldUnits input is not a valid list" )
       end if
       nf = shr_string_listGetNum( FieldUnits)
       if ( nf /= nFields )then
          call shr_sys_abort( subName// &
                              "ERROR: number of FieldLNames do not match FieldList" )
       end if
    end if
    !--- Set the initial value ---
    if ( present(initValue) )then
       Value = initValue
    else
       Value = SHR_CONST_SPVAL
    end if
    !--- Check if grid is arbitrary decomposition or block decomp ---
    call ESMF_GridGet( grid, gridStorage=gridStore, rc=rc )
    call eshr_rc_check( rc, "Error in get of grid storage" )
    if (      gridStore == ESMF_GRID_STORAGE_ARBITRARY )then
       ArbDecomp = .true.
    else
       ArbDecomp = .false.
    end if
    !--- Create array spec to base fields creation on ---
    if ( ArbDecomp )then
       call ESMF_ArraySpecSet( arraySpec1D, 1, ESMF_DATA_REAL, ESMF_R8, rc )
       call eshr_rc_check( rc, subName//'ERROR: setting 1D array spec' )
    else
       call ESMF_ArraySpecSet( arraySpec2D, 2, ESMF_DATA_REAL, ESMF_R8, rc )
       call eshr_rc_check( rc, subName//'ERROR: setting 2D array spec' )
    end if
    !--- Add attribute for fields that share the same grid ---
    call ESMF_StateSetAttribute( eState, name="fieldnames_with_shared_grid", &
                                 value=trim(FieldList), rc=rc )
    call eshr_rc_check( rc, subName// &
                        'ERROR: adding fieldnames_with_shared_grid to state' )
    !--- Check that field names aren't duplicated ---
    allocate( FieldNames(nFields) )
    do i = 1, nFields
       !--- Get the name of this field ---
       call shr_string_listGetName( FieldList, i, FieldName )
       FieldNames(i) = trim(FieldName)
       if ( i > 1 )then
          if ( any(FieldNames(1:i-1) == FieldName) )then
             call shr_sys_abort( subName//" ERROR: Input fieldname is duplicated" )
          end if
       end if
    end do
    !---------------------------------------------------------------------------
    ! Loop over each field
    !---------------------------------------------------------------------------
    allocate( Fields(nFields) )
    do i = 1, nFields
       !--- Get the name of this field ---
       call shr_string_listGetName( FieldList, i, FieldName )
       !--- Create the field ---
       if ( useDataPTR )then
          Fields(i) = ESMF_FieldCreateNoData(grid, arraySpec1D,      &
                                       horzRelloc=ESMF_CELL_CENTER,  &
                                       name=trim(FieldName), rc=rc)
          call eshr_rc_check( rc, subName//'ERROR: creating Field:'//trim(FieldName) )
       else if ( ArbDecomp )then
          Fields(i) = ESMF_FieldCreate(grid, arraySpec1D, allocFlag=ESMF_ALLOC, &
                                       horzRelloc=ESMF_CELL_CENTER,             &
                                       name=trim(FieldName), rc=rc)
          call eshr_rc_check( rc, subName//'ERROR: creating Field:'//trim(FieldName) )
       else
          Fields(i) = ESMF_FieldCreate(grid, arraySpec2D, allocFlag=ESMF_ALLOC, &
                                       horzRelloc=ESMF_CELL_CENTER,             &
                                       name=trim(FieldName), rc=rc)
          call eshr_rc_check( rc, subName//'ERROR: creating Field:'//trim(FieldName) )
       end if
       !
       !--- Add Attributes -----------------------------------------------------
       !
       !--- long_name ---
       if ( present(FieldLNames) )then
          call shr_string_listGetName( FieldLNames, i, long_name )
          call ESMF_FieldSetAttribute( Fields(i), "long_name", trim(long_name), rc )
          call eshr_rc_check( rc, subName//'ERROR: adding long_name to field' )
       end if
       !--- units ---
       if ( present(FieldUnits) )then
          call shr_string_listGetName( FieldUnits, i, units )
          call ESMF_FieldSetAttribute( Fields(i), "units", trim(units), rc )
          call eshr_rc_check( rc, subName//'ERROR: adding units to field' )
       end if
       !--- missing_value ---
       call ESMF_FieldSetAttribute( Fields(i), "missing_value", SHR_CONST_SPVAL, rc )
       call eshr_rc_check( rc, subName//'ERROR: adding missing_value to field' )
       !--- Type of data (Flux or State) ---
       if (      index(FieldName,"F") == 1 )then
          dataType = eshr_estate_dataIsFlux
       else if ( index(FieldName,"S") == 1 )then
          dataType = eshr_estate_dataIsState
       else if ( index(FieldName,"area") == 1 )then
          dataType = eshr_estate_dataIsState
       else if ( index(FieldName,"mask") == 1 )then
          dataType = eshr_estate_dataIsState
       else if ( index(FieldName,"maxfrac") == 1 )then
          dataType = eshr_estate_dataIsState
       else if ( index(FieldName,"lat") == 1 )then
          dataType = eshr_estate_dataIsState
       else if ( index(FieldName,"lon") == 1 )then
          dataType = eshr_estate_dataIsState
       else
          call shr_sys_abort( subName// &
               "ERROR: FieldName does NOT start with F or S, nor is it domain data : "// &
               trim(FieldName) )
       end if
       !
       !--- Set data to missing value ------------------------------------------
       !
       if ( .not. useDataPTR )then
          if ( ArbDecomp )then
             call ESMF_FieldGetDataPointer( Fields(i), Data1D, &
                                            copyflag=ESMF_DATA_REF, rc=rc)
             call eshr_rc_check( rc, subName//'ERROR: Getting data pointer from field' )
             Data1D(:) = Value
             nullify( data1D )
          else
             call ESMF_FieldGetDataPointer( Fields(i), Data2D, &
                                            copyflag=ESMF_DATA_REF, rc=rc)
             call eshr_rc_check( rc, subName//'ERROR: Getting data pointer from field' )
             Data2D(:,:) = Value
             nullify( data2D )
          end if
       end if
    end do
    !--- Add the fields to the state --------------------------------------------
    call ESMF_StateAddField( EState, nFields, Fields, rc )
    call eshr_rc_check( rc, "Error, adding fields to state" )
    !--- Turn some fields to not needed if a list is provided -------------------
    if ( present(NotNeeded) )then
       nf = shr_string_listGetNum( NotNeeded )
       do i = 1, nf
          call shr_string_listGetName( NotNeeded, i, Fieldname )
          if ( .not. any(Fieldname == Fieldnames(:) ) ) &
               call shr_sys_abort( subName// &
               "ERROR: Fieldname in notneeded list NOT in list of fields" )
          call ESMF_StateSetNeeded( EState, itemname=Fieldname, &
                                    neededflag=ESMF_NOTNEEDED,  &
                                    rc=rc )
          call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
       end do
    end if
    deallocate( FieldNames )

END SUBROUTINE eshr_estate_init2DFromList

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_init1DDomain -- Initialize domain
!   
! !DESCRIPTION:
!   
!  Initialize a domain for a model with an arbitrary decomposition.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_init1DDomain(  modeltype, Grid, area, rmask, maxfrac, eSDomain, &
                                      AllowMiss )
    use shr_const_mod,   only: SHR_CONST_SPVAL
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Grid)                      :: Grid       ! Grid fields are on
    character(len=*),        intent(IN)  :: modelType  ! Type of model (atm,lnd,ice,ocn)
    real(SHR_KIND_R8),       intent(IN)  :: area(:)    ! Area of grid cells
    real(SHR_KIND_R8),       intent(IN)  :: rmask(:)   ! Mask if grid cells active
    real(SHR_KIND_R8),       intent(IN)  :: maxfrac(:) ! Maximum areal fraction
    type(ESMF_State),        intent(OUT) :: eSDomain   ! Output ESMF Domain State
    logical,       optional, intent(IN)  :: AllowMiss  ! Allow missing values in input
!EOP
    !----- local -----
    character(len=*),parameter :: subName =   "(eshr_estate_init1DDomain) "
    integer                    :: rc           ! Return code
    real(SHR_KIND_R8), pointer :: Data1D(:)    ! Pointer to 1D data
    type(ESMF_GridStorage)     :: gridStore    ! Type of grid storage
    logical                    :: Allow        ! Allow missing values

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    !
    !----- Make sure grid is arbitrary decomposition ---------------------------
    !
    call ESMF_GridGet( grid, gridStorage=gridStore, rc=rc )
    call eshr_rc_check( rc, "Error in get of grid storage" )
    if (      gridStore /= ESMF_GRID_STORAGE_ARBITRARY )then
       call shr_sys_abort( subname//" : Error, input grid is not arbitrary decomposition")
    end if
    !---- Make sure modeltype is recognized -----
    if ( trim(modeltype) /= "atm" .and. trim(modeltype) /= "lnd" .and.  &
         trim(modeltype) /= "ice" .and. trim(modeltype) /= "ocn" )then
       call shr_sys_abort( subname//" : Error, model type not recognized: "// &
                           trim(modeltype) )
    end if
    !
    !----- Create the state ----------------------------------------------------
    !
    eSDomain = ESMF_StateCreate( "domain_"//trim(modeltype), &
                                 statetype=ESMF_STATE_EXPORT, rc=rc )
    call eshr_rc_check( rc, "Error, creating domain state for:"//trim(modeltype) )
    call eshr_estate_init2DFromList( Grid, DomainFieldList, DomainFieldLNames,        &
                                     DomainFieldUnits, estate=eSDomain )
    !
    !--- Fill state with data --------------------------------------------------
    !
    call ESMF_StateGetDataPointer( eSDomain, itemname="area",   &
                                   dataPointer=Data1D,       &
                                   copyflag=ESMF_DATA_REF, rc=rc )
    call eshr_rc_check( rc, subName//" : Error, in getting area from domain state:"// &
                        trim(modeltype) )
    if ( size(Data1D) /= size(area) )then
       call shr_sys_abort( subname//" : Error, size of area different than input area" )
    end if
    Data1D(:) = area(:)
    nullify( Data1D )
    call ESMF_StateGetDataPointer( eSDomain, itemname="mask",   &
                                   dataPointer=Data1D,       &
                                   copyflag=ESMF_DATA_REF, rc=rc )
    call eshr_rc_check( rc, subName//" : Error, in getting mask from domain state:"// &
                        trim(modeltype) )
    if ( size(Data1D) /= size(rmask) )then
       call shr_sys_abort( subname//" : Error, size of rmask different than input rmask" )
    end if
    Data1D(:) = rmask(:)
    nullify( Data1D )
    call ESMF_StateGetDataPointer( eSDomain, itemname="maxfrac", &
                                   dataPointer=Data1D,           &
                                   copyflag=ESMF_DATA_REF, rc=rc )
    call eshr_rc_check( rc, subName//" : Error, in getting maxfrac from domain state:"// &
                        trim(modeltype) )
    if ( size(Data1D) /= size(maxfrac) )then
       call shr_sys_abort( subname//" : Error, size of maxfrac different than input maxfrac" )
    end if
    Data1D(:) = maxfrac(:)
    nullify( Data1D )
    !--- Check ----
    Allow = .false.
    if ( present( AllowMiss) )then
       if ( AllowMiss ) Allow = .true.
    end if
    if ( .not. Allow )then
       if ( eshr_estate_matchValue( eSDomain, SHR_CONST_SPVAL) ) &
          call shr_sys_abort( subname//" : Error, domain has missing values" )
    end if

END SUBROUTINE eshr_estate_init1DDomain

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_get1DDataPointer -- Get pointer to data from list index
!   
! !DESCRIPTION:
!   
!  Get the pointer to the data from a list and a index into that list.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_get1DDataPointer( eState, list, Listindex, edata1D, expectSize )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State), intent(IN) :: eState     ! Input ESMF state to get data pointer
    character(len=*), intent(IN) :: list       ! Colen delimited list of fieldnames
    integer,          intent(IN) :: ListIndex  ! Integer index into list of field to get
    real(SHR_KIND_R8), pointer   :: edata1D(:) ! Output pointer to field data
    integer, optional,intent(IN) :: expectSize ! Expected size of data
!EOP
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_get1DDataPointer) "
    integer                     :: rc            ! Return code
    integer                     :: ierr          ! Error code
    character(len=SHR_KIND_CL)  :: Fieldn        ! Fieldname

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call shr_string_listGetName( list, ListIndex, Fieldn, rc=ierr)
    call ESMF_StateGetDataPointer( eState, itemname=Fieldn, dataPointer=edata1D,   &
                                   copyflag=ESMF_DATA_REF, rc=rc )
    call eshr_rc_check( rc, subName// &
                            " : Error, in getting pointer to field data on state: "// &
                            trim(Fieldn) )
    if ( present(expectSize) )then
       if ( size(edata1D) /= expectSize ) &
          call shr_sys_abort( subName//" : Data != expected size" )
    end if

END SUBROUTINE eshr_estate_get1DDataPointer

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_set1DDataPointer -- Set pointer to data from list index
!   
! !DESCRIPTION:
!   
!  Get the pointer to the data from a list and a index into that list.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_set1DDataPointer( eState, list, Listindex, data1D )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State), intent(IN) :: eState     ! Input ESMF state to get data pointer
    character(len=*), intent(IN) :: list       ! Colen delimited list of fieldnames
    integer,          intent(IN) :: ListIndex  ! Integer index into list of field to get
    real(SHR_KIND_R8), pointer   :: data1D(:)  ! Input pointer to field data
!EOP
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_set1DDataPointer) "
    integer                     :: rc            ! Return code
    integer                     :: ierr          ! Error code
    character(len=SHR_KIND_CL)  :: Fieldn        ! Fieldname
    type(ESMF_Field)            :: eField        ! ESMF field

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call shr_string_listGetName( list, ListIndex, Fieldn, rc=ierr)
    call ESMF_StateGetField( eState, FieldName=Fieldn, Field=eField, rc=rc )
    call eshr_rc_check( rc, subName// &
                        " : Error, in setting pointer to field data on state: "// &
                        trim(Fieldn) )
    call ESMF_FieldSetDataPointer( eField, dataPointer=data1D,   &
                                   copyFlag=ESMF_DATA_REF, rc=rc )
    call eshr_rc_check( rc, subName// &
                        " : Error, in setting pointer to field data on state: "// &
                        trim(Fieldn) )

END SUBROUTINE eshr_estate_set1DDataPointer

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_checkGrid -- check that input grids match
!   
! !DESCRIPTION:
!   
!  Check that the four input grids from each component match
!      
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION eshr_estate_checkGridsMatch( x2i_i, x2a_a, x2o_o, x2l_l, gridsExact, &
                                              Debug)
    use shr_string_mod, only: shr_string_listGetIndexF
    use shr_const_mod,  only: SHR_CONST_SPVAL
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State), intent(IN) :: x2i_i          ! Ice import state
    type(ESMF_State), intent(IN) :: x2a_a          ! Atmosphere import state
    type(ESMF_State), intent(IN) :: x2o_o          ! Ocean import state
    type(ESMF_State), intent(IN) :: x2l_l          ! Land import state
    logical,          intent(IN) :: gridsExact     ! If all grids should be identical
    logical, optional,intent(IN) :: Debug          ! If should print info on grids
!EOP
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_checkGridsMatch) "
    integer                             :: rc             ! Return code
    integer                             :: i              ! Index
    integer                             :: nsize_atm      ! Number of atm grid points
    integer                             :: nFields        ! Number of fields in state
    integer                             :: indexMaxFrac   ! Index to maximum fraction
    integer                             :: indexRMask     ! Index to real mask
    integer,                    pointer :: Mask_lnd(:)    ! Land mask on atm decomp
    character(len=SHR_KIND_CL), pointer :: FieldNames(:)  ! Names of fields  in state
    character(len=SHR_KIND_CL)          :: AtmFieldName   ! Field name on atmosphere
    character(len=SHR_KIND_CL)          :: IceFieldName   ! Field name on sea-ice 
    character(len=SHR_KIND_CL)          :: LndFieldName   ! Field name on land
    character(len=SHR_KIND_CL)          :: OcnFieldName   ! Field name on ocean
    type(ESMF_Field)                    :: AtmField       ! Atmosphere field
    type(ESMF_Field)                    :: IceField       ! Sea-ice field
    type(ESMF_Field)                    :: LndField       ! Land field
    type(ESMF_Field)                    :: OcnField       ! Ocean field
    type(ESMF_Grid)                     :: AtmGrid        ! Atmosphere grid
    type(ESMF_Grid)                     :: IceGrid        ! Sea-ice grid
    type(ESMF_Grid)                     :: LndGrid        ! Land grid
    type(ESMF_Grid)                     :: OcnGrid        ! Ocean grid
    type(ESMF_State)                    :: domain_ice     ! Ice domain
    type(ESMF_State)                    :: domain_atm     ! Atmosphere domain
    type(ESMF_State)                    :: domain_ocn     ! Ocean domain
    type(ESMF_State)                    :: domain_lnd     ! Land domain
    type(ESMF_State)                    :: LndDomain_atm  ! Land domain on atm grid
    type(ESMF_Bundle)                   :: EBundleIn      ! Input bundle
    type(ESMF_Bundle)                   :: EBundleOut     ! Output bundle
    type(ESMF_RouteHandle)              :: eRoute         ! Route between land and atm
    real(SHR_KIND_R8),          pointer :: maxfrac_atm(:) ! Maximum atm fraction
    real(SHR_KIND_R8),          pointer :: maxfrac_ocn(:) ! Maximum ocean fraction
    real(SHR_KIND_R8),          pointer :: maxfrac_ice(:) ! Maximum sea-ice fraction
    real(SHR_KIND_R8),          pointer :: maxfrac_lnd(:) ! Maximum land fraction
    real(SHR_KIND_R8),          pointer :: area(:)        ! Area
    real(SHR_KIND_R8),          pointer :: rmask(:)       ! Mask
    real(SHR_KIND_R8),          pointer :: maxfrac(:)     ! Maximum fraction
    logical                             :: DebugInfo      ! If should print info on grids

!-------------------------------------------------------------------------------
! Notes: Later need to allow grids to be different between atmosphere and ocean.
!-------------------------------------------------------------------------------

   if ( present(Debug) )then
     DebugInfo = Debug
   else
     DebugInfo = .false.
   end if
   eshr_estate_checkGridsMatch = .true.
   !
   !------  Get first field name from each state -------------------------------
   !
   call eshr_estate_getFieldNames( x2a_a, FieldNames, nFields )
   if ( nFields == 0 )then
      eshr_estate_checkGridsMatch = .false.
      return
   end if
   atmFieldName = FieldNames(1)
   deallocate( FieldNames )
   call eshr_estate_getFieldNames( x2i_i, FieldNames, nFields )
   if ( nFields == 0 )then
      eshr_estate_checkGridsMatch = .false.
      return
   end if
   iceFieldName = FieldNames(1)
   deallocate( FieldNames )
   call eshr_estate_getFieldNames( x2l_l, FieldNames, nFields )
   if ( nFields == 0 )then
      eshr_estate_checkGridsMatch = .false.
      return
   end if
   lndFieldName = FieldNames(1)
   deallocate( FieldNames )
   call eshr_estate_getFieldNames( x2o_o, FieldNames, nFields )
   if ( nFields == 0 )then
      eshr_estate_checkGridsMatch = .false.
      return
   end if
   ocnFieldName = FieldNames(1)
   deallocate( FieldNames )
   !
   !------  Get field and grid from each state ---------------------------------
   ! 
   call ESMF_StateGetField( x2a_a, fieldname=atmFieldName, field=atmField, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting field from atm import state" )
   call ESMF_FieldGet( atmField, grid=atmGrid, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting grid from atm import state field" )
   call ESMF_StateGetField( x2i_i, fieldname=iceFieldName, field=iceField, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting field from ice import state" )
   call ESMF_FieldGet( iceField, grid=iceGrid, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting grid from ice import state field" )
   call ESMF_StateGetField( x2l_l, fieldname=lndFieldName, field=lndField, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting field from lnd import state" )
   call ESMF_FieldGet( lndField, grid=lndGrid, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting grid from lnd import state field" )
   call ESMF_StateGetField( x2o_o, fieldname=ocnFieldName, field=ocnField, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting field from ocn import state" )
   call ESMF_FieldGet( ocnField, grid=ocnGrid, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting grid from ocn import state field" )
   !
   !------  Check that atm, ice, ocn grids coordinates are the same ------------
   ! 
   if ( DebugInfo ) write(6,*) "Compare atm and ice grid"
   if ( .not. eshr_estate_Egrids2DAreSame( atmGrid, iceGrid, GridsExact, &
                                    SameDecomp=.true., Debug=DebugInfo ) )then
      if ( DebugInfo ) write(6,*) "atm and ice grids or decomposition are different"
      eshr_estate_checkGridsMatch = .false.
      return
   end if
   if ( DebugInfo ) write(6,*) "Compare ice and ocn grid"
   if ( .not. eshr_estate_Egrids2DAreSame( iceGrid, ocnGrid, GridsExact, &
                                    SameDecomp=.true., Debug=DebugInfo ) )then
      if ( DebugInfo ) write(6,*) "ice and ocn grids or decomposition are different"
      eshr_estate_checkGridsMatch = .false.
      return
   end if
   !--- Compare basic attributes of atm and lnd grid ----
   if ( DebugInfo ) write(6,*) "Compare basic attributes of atm and lnd grid"
   if ( .not. eshr_estate_Egrids2DAreSame( atmGrid, lndGrid, GridsExact, &
                                    SameDecomp=.true., NotData=.true.,   &
                                    Debug=DebugInfo ) )then
      if ( DebugInfo ) write(6,*) "atm and lnd grids or decomposition are different"
      eshr_estate_checkGridsMatch = .false.
      return
   end if
   !
   !--- Get domains for each model ---------------------------------------------
   !
   call ESMF_StateGetState( x2a_a, "domain_atm", domain_atm, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting domain_atm from atm import state" )
   call ESMF_StateGetState( x2i_i, "domain_ice", domain_ice, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting domain_ice from ice import state" )
   call ESMF_StateGetState( x2l_l, "domain_lnd", domain_lnd, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting domain_lnd from lnd import state" )
   call ESMF_StateGetState( x2o_o, "domain_ocn", domain_ocn, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting domain_ocn from ocn import state" )
   !
   !------ Check that maxfrac on domains are consistent ------------------------
   !
   indexMaxFrac = shr_string_listGetIndexF( DomainFieldList, MaxFracName )
   call eshr_estate_getDataPointer( domain_atm, DomainFieldList, indexMaxFrac, &
                                    maxfrac_atm )
   call eshr_estate_getDataPointer( domain_ice, DomainFieldList, indexMaxFrac, &
                                    maxfrac_ice )
   call eshr_estate_getDataPointer( domain_ocn, DomainFieldList, indexMaxFrac, &
                                    maxfrac_ocn )
   nsize_atm = size(maxfrac_atm)
   if ( any(maxfrac_atm(:) /= 1.0_SHR_KIND_R8) )then
      if ( DebugInfo ) write(6,*) "atmosphere max fraction not equal to 1."
      eshr_estate_checkGridsMatch = .false.
   end if
   if ( any(maxfrac_ice(:) /= maxfrac_ocn) )then
      if ( DebugInfo ) write(6,*) "ice and ocean max fractions are different"
      eshr_estate_checkGridsMatch = .false.
   end if
   !
   !--- Compare the land and atmosphere grid --------------------------------
   !
   if ( DebugInfo ) write(6,*) "Compare atm and lnd grid"
   !--- Regrid land domain to atmosphere grid -----
   allocate( area(nsize_atm), rmask(nsize_atm), maxfrac(nsize_atm) )
   area(:)    = SHR_CONST_SPVAL
   rmask(:)   = SHR_CONST_SPVAL
   maxfrac(:) = 0.0_SHR_KIND_R8
   call eshr_estate_initDomain(  "lnd", atmGrid, area, rmask, maxfrac, LndDomain_atm, &
                                 AllowMiss=.true. )
   deallocate( area, rmask, maxfrac )
   call eshr_estate_reDistInit( domain_lnd, LndDomain_atm, EBundleIn, EBundleOut, &
                                eRoute, noGridCheck=.true. )
   call eshr_estate_reDist( eBundleIn, eBundleOut, eRoute )
   call eshr_estate_getDataPointer( LndDomain_atm, DomainFieldList, indexMaxFrac, &
                                    maxfrac_lnd, expectSize=nsize_atm )
   indexRMask = shr_string_listGetIndexF( DomainFieldList, MaskName )
   call eshr_estate_getDataPointer( LndDomain_atm, DomainFieldList, indexRMask,   &
                                    rmask, expectSize=nsize_atm )
   allocate( Mask_lnd(nsize_atm) )
   do i = 1, nsize_atm
      if ( (rmask(i) == SHR_CONST_SPVAL) .or. (rmask(i) == 0.0_SHR_KIND_R8) )then
         Mask_lnd(i) = 0
      else
         Mask_lnd(i) = 1
      end if
   end do
   
   if ( .not. eshr_estate_Egrids2DAreSame( atmGrid, lndGrid, GridsExact, &
                                           Grid1Mask=Mask_lnd, Debug=DebugInfo ) )then
      eshr_estate_checkGridsMatch = .false.
   else
      if ( (size(maxfrac_ocn) /= nsize_atm) .or. (size(maxfrac_ocn) /= nsize_atm) )then
         eshr_estate_checkGridsMatch = .false.
      else
         do i = 1, nsize_atm
            if ( Mask_lnd(i) == 1 )then
               if ( abs(maxfrac_lnd(i)+maxfrac_ocn(i)-1.0_SHR_KIND_R8) > eps ) &
                   eshr_estate_checkGridsMatch = .false.
            else
               if( abs(maxfrac_ocn(i)-1.0_SHR_KIND_R8) > eps ) &
                  eshr_estate_checkGridsMatch = .false.
            end if
         end do
      end if
   end if
   deallocate( Mask_lnd )
   nullify( maxfrac_lnd, maxfrac_ocn, maxfrac_atm )
   call ESMF_StateDestroy( LndDomain_atm, rc=rc )
   call eshr_rc_check( rc, subName//" : Error destroying land domain on atm grid" )

END FUNCTION eshr_estate_checkGridsMatch

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_getStats -- Get the basic stats on the input ESMF State
!   
! !DESCRIPTION:
!   
!  Get basic statistics on all of the fields in the input ESMF state.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_getStats( eState, mpicom, rootid )
    use shr_kind_mod,    only: SHR_KIND_CL
    use shr_const_mod,   only: SHR_CONST_SPVAL
    use shr_mpi_mod,     only: shr_mpi_gatherv, shr_mpi_gathScatVInit, &
                               shr_mpi_commrank, shr_mpi_barrier
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),  intent(IN) :: EState      ! Input state to get the stats on
    integer, optional, intent(IN) :: mpicom      ! Mpi communicator to gather to
    integer, optional, intent(IN) :: Rootid      ! Rootid to gather on
!EOP
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_getStats) "
    integer                             :: rc            ! Return code
    integer                             :: j, n          ! Indices
    integer                             :: nFields       ! Number of fields in state
    integer                             :: rank          ! Rank of field
    integer                             :: mpiRank       ! MPI rank
    character(len=SHR_KIND_CL), pointer :: FieldNames(:) ! Names of fields  in state
    character(len=SHR_KIND_CL)          :: EStateName    ! Name of state
    logical                             :: printIt       ! If printing
    type(ESMF_Field)                    :: Field         ! Field to query
    type(ESMF_Array)                    :: array         ! Pointer to array in field
    type(ESMF_DataType)                 :: dataType      ! Type of data in array
    type(ESMF_DataKind)                 :: kind          ! Kind of data in array
    real(ESMF_KIND_R8),         pointer :: data2D(:,:)   ! Pointer 2D array in field
    logical,                    pointer :: mask2D(:,:)   ! 2D mask
    real(ESMF_KIND_R8),         pointer :: edata1D(:)    ! Pointer 1D array in field
    real(ESMF_KIND_R8),         pointer :: data1D(:)     ! Pointer to global 1D data
    logical,                    pointer :: mask1D(:)     ! 1D mask
    integer,                    pointer :: displs(:)     ! Displacement array
    integer,                    pointer :: globSize(:)   ! Global size
    real(ESMF_KIND_R8)                  :: dataMin       ! Data minimum value
    real(ESMF_KIND_R8)                  :: dataMax       ! Data maximum value
    real(ESMF_KIND_R8)                  :: dataAvg       ! Data average value
    real(ESMF_KIND_R8)                  :: dataStDev     ! Data standard-deviation 
    type(ESMF_NeededFlag)               :: need          ! If field needed or not
    integer                             :: lbnd(2)       ! Lower bound of data arrays
    integer                             :: ubnd(2)       ! Upper bound of data arrays
    real(ESMF_KIND_R8)                  :: missing_val   ! Missing value

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( present(mpicom) .and. present(rootid) )then
      call shr_mpi_commrank( mpicom, mpiRank )
      if ( mpiRank /= rootid )then
         printIt = .false.
      else
         printIt = .true.
      end if
   else
      printIt = .true.
   end if
   call ESMF_StateGet( state=EState, name=EStateName, rc=rc )
   call eshr_rc_check( rc, subName//' : Error getting state name' )
   if ( printIt ) write(6,*) 'Get statistics on fields on state: ', trim(EStateName)
   !
   !------  Get number and list of field names in state ------------------------
   !
   call eshr_estate_getFieldNames( eState, FieldNames, nFields )
   !
   !------  Loop through all fields in state and extract it's info -------------
   !
   do j = 1, nFields
      call ESMF_StateGetField( state=EState, field=Field, fieldname=FieldNames(j), &
                               rc=rc )
      call eshr_rc_check( rc, subName//' : Error getting field '//FieldNames(j)// &
                              ' on input state' )

      !------ Check if field is needed and if not exit early -------------------
      call ESMF_StateGetNeeded( EState, itemname=Fieldnames(j), neededflag=need, &
                                rc=rc )
      call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
      if ( need == ESMF_NOTNEEDED) cycle    !--- Exit early if this field NOT needed ------

#ifdef ESMF_3
      !call ESMF_FieldGetInternArray( field=Field, array=array, rc=rc )
      call shr_sys_abort( "function not correct for version 3 of ESMF" )
#else
      call ESMF_FieldGetArray( field=Field, array=array, rc=rc )
#endif
      call eshr_rc_check( rc, subName//' : Error getting array from field '// &
                          FieldNames(j)// ' on input state' )
      call ESMF_ArrayGet( array=array, rank=rank, type=dataType, kind=kind, rc=rc )
      call eshr_rc_check( rc, subName//': Error in getting array type info from field' )
      !----- Get missing value ----
      call ESMF_FieldGetAttributeInfo( Field, "missing_value", count=n, rc=rc )
      call eshr_rc_check( rc, subName//': Error in getting attribute info from state' )
      if ( n == 1 )then
         call ESMF_FieldGetAttribute( Field, "missing_value", value=missing_val, &
                                      rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting missing value from state' )
      else
         missing_val = SHR_CONST_SPVAL
      end if
      if ( (dataType == ESMF_DATA_REAL) .and. (kind == ESMF_R8) .and. (rank == 2) )then
         call ESMF_FieldGetDataPointer( field=field,   &
                                        ptr=data2D,    &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to field data' )
         lbnd(1:2) = lbound(data2D)
         ubnd(1:2) = ubound(data2D)
         allocate( mask2D(lbnd(1):ubnd(1),lbnd(2):ubnd(2))  )
         mask2D = (data2D /= missing_val)
         if ( count( mask2D ) > 0 )then
            dataMin   = minval(data2D, mask=mask2D )
            dataMax   = maxval(data2D, mask=mask2D )
            dataAvg   =    sum(data2D, mask=mask2D ) / count( mask2D )
            dataStDev =    sum((data2D-dataAvg)**2, mask=mask2D ) / count( mask2D )
            dataStDev = sqrt( dataStDev )
         else
            dataMin   = missing_val
            dataMax   = missing_val
            dataAvg   = missing_val
            dataStDev = missing_val
            dataStDev = missing_val
            if ( printIt ) write(6,*) subName//" : Warning " &
                           //trim(FieldNames(j))// " all missing_val"
         end if
         nullify( data2D )
         deallocate( mask2D )
      else if ( (dataType == ESMF_DATA_REAL) .and. (kind == ESMF_R8) .and. (rank == 1) &
      )then
         call ESMF_FieldGetDataPointer( field=field,   &
                                        ptr=edata1D,   &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to field data' )
         if ( present(mpicom) .and. present(rootid) )then
            call shr_mpi_gathScatvInit( mpicom, rootid, edata1D, data1D, globSize, &
                                        displs )
            call shr_mpi_barrier( mpicom )
            call shr_mpi_gatherv( edata1D, size(edata1D), data1D, globSize, displs,  &
                                  rootid, mpicom )
            call shr_mpi_barrier( mpicom )
         else
            data1D => edata1D
         end if
         lbnd(1:1) = lbound(data1D)
         ubnd(1:1) = ubound(data1D)
         allocate( mask1D(lbnd(1):ubnd(1)) )
         mask1D = (data1D /= missing_val)
         if ( count( mask1D ) > 0 )then
            dataMin   = minval(data1D, mask=mask1D )
            dataMax   = maxval(data1D, mask=mask1D )
            dataAvg   =    sum(data1D, mask=mask1D ) / count( mask1D )
            dataStDev =    sum((data1D-dataAvg)**2, mask=mask1D ) / count( mask1D )
            dataStDev = sqrt( dataStDev )
         else
            dataMin   = missing_val
            dataMax   = missing_val
            dataAvg   = missing_val
            dataStDev = missing_val
            dataStDev = missing_val
            if ( printIt ) write(6,*) subName//" : Warning " &
                           //trim(FieldNames(j))// " all missing_val"
         end if
         if ( present(mpicom) .and. present(rootid) )then
            deallocate( data1D, globSize, displs )
         else
            nullify( data1D )
         end if
         nullify( edata1D )
         deallocate( mask1D )
      else
         call shr_sys_abort( subName// &
                             " : Error rank not equal to [1 or 2] or not real-r8")
      end if
      if ( dataMin /= missing_val .and. printIt ) &
         write(6,*) trim(FieldNames(j)),            &
                          ' std-dev: ', dataStDev,  &
                          ' min:     ', dataMin,    &
                          ' avg:     ', dataAvg,    &
                          ' max:     ', dataMax
   end do
   if ( nFields > 0 ) deallocate( FieldNames )

END SUBROUTINE eshr_estate_getStats

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_fieldsAreSame -- Return true if fields are identical
!   
! !DESCRIPTION:
!   
!  Check that all the fields on the two input states have the same names
!  type of data, and sizes, and contain identical values. If the field lists
!  are different or, the types different or sizes different -- will return false.
!  Optionally can enter a list of fields, in this case both states MUST have all
!  of the fields on the input list -- if NOT return false.
!      
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION eshr_estate_fieldsAreSame( eState1, eState2, notData, namesOnly, &
                                            onlyNeeded, FieldNames, Mask1D )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),  intent(IN)    :: EState1    ! First state
    type(ESMF_State),  intent(IN)    :: EState2    ! State to compare to
    logical, optional, intent(IN)    :: notData    ! Only compare names, sizes and types
    logical, optional, intent(IN)    :: namesOnly  ! Only compare names not sizes
    logical, optional, intent(IN)    :: onlyNeeded ! Only compare fields that are needed
    character(len=SHR_KIND_CL), optional, intent(IN) :: FieldNames(:)     ! Field names 1
    integer, optional, intent(IN)    :: Mask1D(:)  ! 1D mask of grid points to compare
!EOP
    !----- local -----
    character(len=*),          parameter :: subName = " (eshr_estate_fieldsAreSame) "
    integer                              :: rc                ! Return code
    integer                              :: i, j              ! Indices
    integer                              :: arank             ! Array rank
    integer                              :: nFields           ! No. fields
    integer                              :: nFields2          ! No. fields 2
    character(len=SHR_KIND_CL), pointer  :: FieldNames1(:)    ! Field names 1
    character(len=SHR_KIND_CL), pointer  :: FieldNames2(:)    ! Field names 2
    real(SHR_KIND_R8),          pointer  :: Data1D1(:)        ! Pointer to data 1
    real(SHR_KIND_R8),          pointer  :: Data1D2(:)        ! Pointer to data 2
    real(SHR_KIND_R8),          pointer  :: Data2D1(:,:)      ! Pointer to data 1
    real(SHR_KIND_R8),          pointer  :: Data2D2(:,:)      ! Pointer to data 2
    type(ESMF_DataType)                  :: dataType          ! Type of data
    type(ESMF_DataKind)                  :: dataKind          ! Kind of data
    type(ESMF_Field)                     :: Field1            ! field 1
    type(ESMF_Field)                     :: Field2            ! field 2
    type(ESMF_Array)                     :: array1            ! array 1
    type(ESMF_Array)                     :: array2            ! array 2
    type(ESMF_NeededFlag)                :: need              ! If field needed or not
    logical                              :: different         ! If difference found
    logical                              :: neededOnly        ! Only needed

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( present(OnlyNeeded) )then
      neededOnly = OnlyNeeded
   else
      neededOnly = .false.
   end if
   arank     = 1
   different = .false.
   !----- Make sure list of fields on both states are the same ----
   call eshr_estate_getFieldNames( eState1,  FieldNames1,  nFields  )
   call eshr_estate_getFieldNames( eState2,  FieldNames2,  nFields2 )
   if ( (nFields  == 0) .and. (nFields2 == 0) )then
      different = .false.
      return
   else if ( nFields  == 0 )then
      different = .true.
      deallocate( FieldNames2 )
      return
   else if ( nFields2 == 0 )then
      different = .true.
      deallocate( FieldNames1 )
      return
   end if
   if ( .not. present(FieldNames) )then
      if ( nFields /= nFields2 ) different = .true.
      do i = 1, nFields
         if ( different ) exit
         if ( trim(FieldNames1(i)) /= trim(FieldNames2(i)) ) different = .true.
      end do
   !------ If list of fieldnames entered make sure both states have everything on list ---
   else
      do i = 1, size(FieldNames)
         if ( .not. any(FieldNames1(:) == FieldNames(i) ) )then
            different = .true.
            exit
         end if
      end do
      do i = 1, size(FieldNames)
         if ( .not. any(FieldNames2(:) == FieldNames(i) ) )then
            different = .true.
            exit
         end if
      end do
      !---- Copy fieldnames to fieldnames1 array ----
      nFields = size(FieldNames)
      do i = 1, nFields
         if ( different ) exit
         FieldNames1(i) = FieldNames(i)
      end do
   end if
   if ( nFields2 > 0 ) deallocate(FieldNames2)
   ! ---- Make sure dimension sizes same, and 1D or 2D list of SHR_KIND_R8 fields ----
   do i = 1, nFields
      !----- Check reasons to end loop early -----
      if ( different ) exit
      if ( present(namesOnly) )then
         if ( namesOnly ) exit
      end if
      if ( neededOnly )then
         call ESMF_StateGetNeeded( EState1, itemname=Fieldnames1(i), neededflag=need, &
                                   rc=rc )
         call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
         if ( need==ESMF_NOTNEEDED ) cycle
         call ESMF_StateGetNeeded( EState2, itemname=Fieldnames1(i), neededflag=need, &
                                   rc=rc )
         call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
         if ( need==ESMF_NOTNEEDED ) cycle
      end if
      !----- Get fields and compare type, rank, and size information
      call ESMF_StateGetField( state=EState1, field=Field1, &
                               fieldname=FieldNames1(i), rc=rc )
      call eshr_rc_check( rc, subName//' : Error getting field '//FieldNames1(i)// &
                              ' on input state' )
      call ESMF_StateGetField( state=EState2, field=Field2,  &
                               fieldname=FieldNames1(i), rc=rc )
      call eshr_rc_check( rc, subName//' : Error getting field '//FieldNames1(i)// &
                              ' on input state' )
#ifdef ESMF_3
      !call ESMF_FieldGetInternArray( field=Field1, array=array1, rc=rc )
      call shr_sys_abort( "function not correct for version 3 of ESMF" )
#else
      call ESMF_FieldGetArray( field=Field1, array=array1, rc=rc )
#endif
      call eshr_rc_check( rc, subName//' : Error getting array from field '// &
                          FieldNames1(i)// ' on input state' )
      call ESMF_ArrayGet( array=array1, rank=arank, type=dataType, kind=dataKind, &
                          rc=rc )
      call eshr_rc_check( rc, subName// &
                          ': Error in getting array type info from field' )
      if ( (dataType /= ESMF_DATA_REAL) .and. (dataKind /= ESMF_R8) )then
         call shr_sys_abort( subName//" : Error field not real-r8" )
      end if
      if ( arank /= 1 .and. arank /= 2 )then
         call shr_sys_abort( subName//" : Error field not 1 or 2D" )
      end if
      if ( arank == 1 )then
         call ESMF_FieldGetDataPointer( field=field1, &
                                        ptr=data1D1,  &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting ptr to fieldIn data' )
         call ESMF_FieldGetDataPointer( field=field2, &
                                        ptr=data1D2,  &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting ptr to fieldOut data' )
         if ( size(data1D1) /= size(data1D2) ) different = .true.
         nullify( data1D1 )
         nullify( data1D2 )
      else
         call ESMF_FieldGetDataPointer( field=field1, &
                                        ptr=data2D1,  &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting ptr to fieldIn data' )
         call ESMF_FieldGetDataPointer( field=field2, &
                                        ptr=data2D2,  &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting ptr to fieldOut data' )
         if ( size(data2D1,1) /= size(data2D2,1) .or. &
              size(data2D1,2) /= size(data2D2,2) ) different = .true.
         nullify( data2D1 )
         nullify( data2D2 )
      end if
   end do
   ! ---- Actually compare the data ------
   if ( arank == 1 )then
      do i = 1, nFields
         !----- Check reasons to end loop early ------
         if ( different ) exit
         if ( present(notData) )then
            if ( notData ) exit
         end if
         if ( present(namesOnly) )then
            if ( namesOnly ) exit
         end if
         if ( neededOnly )then
            call ESMF_StateGetNeeded( EState1, itemname=Fieldnames1(i), neededflag=need, &
                                      rc=rc )
            call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
            if ( need==ESMF_NOTNEEDED ) cycle
            call ESMF_StateGetNeeded( EState2, itemname=Fieldnames1(i), neededflag=need, &
                                      rc=rc )
            call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
            if ( need==ESMF_NOTNEEDED ) cycle
         end if
         !----- Get pointers to data and compare the data ------
         call ESMF_StateGetDataPointer( state=eState1,            &
                                        itemName=FieldNames1(i),  &
                                        dataPointer=data1D1,      &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
         call ESMF_StateGetDataPointer( state=eState2,            &
                                        itemName=FieldNames1(i),  &
                                        dataPointer=data1D2,      &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
         if ( .not. present(Mask1D) )then
            if ( any( data1D1(:) /= data1D2(:) ) ) different = .true.
         else
            do j = 1, size(data1D1)
               if ( (Mask1D(j) == 1) .and. (data1D1(j) /= data1D2(j)) )then
                  different = .true.
                  exit
               end if
            end do
         end if
         nullify( data1D1 )
         nullify( data1D2 )
      end do
   else
      do i = 1, nFields
         !----- Check reasons to end loop early ------
         if ( different ) exit
         if ( present(notData) )then
            if ( notData ) exit
         end if
         if ( present(namesOnly) )then
            if ( namesOnly ) exit
         end if
         if ( neededOnly )then
            call ESMF_StateGetNeeded( EState1, itemname=Fieldnames1(i), neededflag=need, &
                                      rc=rc )
            call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
            if ( need==ESMF_NOTNEEDED ) cycle
            call ESMF_StateGetNeeded( EState2, itemname=Fieldnames1(i), neededflag=need, &
                                      rc=rc )
            call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
            if ( need==ESMF_NOTNEEDED ) cycle
         end if
         !----- Get pointers to data and compare the data ------
         call ESMF_StateGetDataPointer( state=eState1,            &
                                        itemName=FieldNames1(i),  &
                                        dataPointer=data2D1,      &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
         call ESMF_StateGetDataPointer( state=eState2,            &
                                        itemName=FieldNames1(i),  &
                                        dataPointer=data2D2,      &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
         if ( any( data2D1(:,:) /= data2D2(:,:) ) ) different = .true.
         nullify( data2D1 )
         nullify( data2D2 )
      end do
   end if
   if ( nFields > 0 ) deallocate( FieldNames1 )
   eshr_estate_fieldsAreSame = .not. different
   return

END FUNCTION  eshr_estate_fieldsAreSame

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_matchValue -- Return true if matches a given value.
!   
! !DESCRIPTION:
!   
!      
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION eshr_estate_matchValue( eState, Value, PrintNames )
    use shr_sys_mod, only: shr_sys_flush
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),  intent(IN)    :: EState     ! State to query
    real(SHR_KIND_R8), intent(IN)    :: Value      ! Value to data compare to
    logical, optional, intent(IN)    :: PrintNames ! Print names of data matching
!EOP
    !----- local -----
    character(len=*),          parameter :: subName = " (eshr_estate_matchValue) "
    character(len=SHR_KIND_CL), pointer :: FieldNames(:) ! Fieldnames to run through
    integer                             :: nFields       ! No. fields on estate
    integer                             :: i             ! Indices
    integer                             :: rc            ! Return code
    real(SHR_KIND_R8),         pointer  :: Data1D(:)     ! Pointer to estate data
    type(ESMF_NeededFlag)               :: need          ! If needed or not

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    eshr_estate_matchValue = .false.
    call eshr_estate_getFieldNames( eState,  FieldNames,   nFields )
    do i = 1, nFields
       call ESMF_StateGetDataPointer( state=eState,             &
                                      itemName=FieldNames(i),   &
                                      dataPointer=data1D,       &
                                      copyflag=ESMF_DATA_REF, rc=rc )
       call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
       if ( rc /= ESMF_SUCCESS ) cycle
       call ESMF_StateGetNeeded( eState, itemname=Fieldnames(i), neededflag=need, &
                                 rc=rc )
       call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
       if ( (need == ESMF_NEEDED) .and. any(data1D(:) == Value) )then
          eshr_estate_matchValue = .true.
          if ( present( PrintNames) )then
             if ( PrintNames )then
                write(6,*) "Matching field = ", trim(FieldNames(i))
                call shr_sys_flush(6)
             end if
          end if
       end if
       nullify( data1D )
       if ( eshr_estate_matchValue .and. .not. present(PrintNames) ) exit
    end do
    deallocate( Fieldnames )

END FUNCTION eshr_estate_matchValue

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_copyFields -- Copy fields from one ESMF state to another.
!   
! !DESCRIPTION:
!   
!  Copy the fields that are shared on two ESMF States from one to the other.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_copyFields( Init, rank, eStateIn, eStateOut, FieldNames, &
                                   nFields )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    logical,              intent(IN)    :: Init          ! If doing initialization phase
    integer,              intent(IN)    :: rank          ! Expected rank of fields to copy
    type(ESMF_State),     intent(IN)    :: EStateIn      ! Input state to copy
    type(ESMF_State),     intent(INOUT) :: EStateOut     ! Output state to copy to
    character(len=SHR_KIND_CL), pointer :: FieldNames(:) ! Fieldnames to copy
    integer,              intent(INOUT) :: nFields       ! No. input fields
!EOP
    !----- local -----
    character(len=*),          parameter :: subName = " (eshr_estate_copyFields) "
    integer                              :: rc                ! Return code
    integer                              :: i                 ! Indices
    integer                              :: arank             ! Array rank
    real(SHR_KIND_R8),          pointer  :: Data1DIn(:)       ! Pointer to input data
    real(SHR_KIND_R8),          pointer  :: Data1DOut(:)      ! Pointer to output data
    real(SHR_KIND_R8),          pointer  :: Data2DIn(:,:)     ! Pointer to input data
    real(SHR_KIND_R8),          pointer  :: Data2DOut(:,:)    ! Pointer to output data
    type(ESMF_DataType)                  :: dataType          ! Type of data
    type(ESMF_DataKind)                  :: dataKind          ! Kind of data
    type(ESMF_Field)                     :: EFieldIn          ! Input field
    type(ESMF_Grid)                      :: EGridIn           ! Input grid
    type(ESMF_Field)                     :: EFieldOut         ! Output field
    type(ESMF_Grid)                      :: EGridOut          ! Output grid
    character(len=SHR_KIND_CL), pointer  :: FieldNamesOut(:)  ! Output Fieldnames
    integer                              :: nFieldsOut        ! No. output fields

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   !
   !--- If on initialization phase do checking to make sure states are consistent
   !
   if ( Init )then
      !--- Find names that are shared between the two states ---
      call eshr_estate_getSharedFNames( eStateIn, eStateOut, FieldNames, nFields, &
                                        onlyNeeded=OnlyUseNeededData )
      if ( nFields == 0 )then
         write(6,*) subName// " : The two input states do NOT share any field names"
         return
      end if
      !--- Now check that for shared names data sizes are correct ----
      if ( .not. eshr_estate_fieldsAreSame( eStateIn, eStateOut, notData=.true., &
      FieldNames=FieldNames ) ) &
         call shr_sys_abort( subName// &
                      " : Dimensions of fields on input states are different" )
      call ESMF_StateGetField( eStateIn, FieldName=FieldNames(1), field=EFieldIn, &
                               rc=rc )
      call eshr_rc_check( rc, subName// &
                       ': Error in getting in-field' )
      call ESMF_StateGetField( eStateOut, FieldName=FieldNames(1), field=EFieldOut, &
                               rc=rc )
      call eshr_rc_check( rc, subName// &
                       ': Error in getting out-field' )
      call ESMF_FieldGet( EFieldIn, grid=EGridIn, rc=rc )
      call eshr_rc_check( rc, subName// &
                       ': Error in getting in-field grid' )
      call ESMF_FieldGet( EFieldOut, grid=EGridOut, rc=rc )
      call eshr_rc_check( rc, subName// &
                       ': Error in getting out-field grid' )
      !--- Check that grids are same and on same decomposition -----
      if ( .not. eshr_estate_Egrids2DAreSame( EGridIn, EGridOut, Exact=.true., &
                                              SameDecomp=.true. ) ) &
         call shr_sys_abort( subName// &
                      " : Input grids of states to copy are different" )
   end if
   !---- Actually copy the data ------
   if ( rank == 1 )then
      do i = 1, nFields
         call ESMF_StateGetDataPointer( state=eStateIn,           &
                                        itemName=FieldNames(i),   &
                                        dataPointer=data1DIn,     &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
         call ESMF_StateGetDataPointer( state=eStateOut,          &
                                        itemName=FieldNames(i),   &
                                        dataPointer=data1DOut,    &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
         data1DOut(:) = data1dIn(:)
         nullify( data1DIn  )
         nullify( data1DOut )
      end do
   else
      do i = 1, nFields
         call ESMF_StateGetDataPointer( state=eStateIn,           &
                                        itemName=FieldNames(i),   &
                                        dataPointer=data2DIn,     &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
         call ESMF_StateGetDataPointer( state=eStateOut,          &
                                        itemName=FieldNames(i),   &
                                        dataPointer=data2DOut,    &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state' )
         data2DOut(:,:) = data2DIn(:,:)
         nullify( data2DIn  )
         nullify( data2DOut )
      end do
   end if

END SUBROUTINE eshr_estate_copyFields

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_reDistInit -- Initialize a ReDist. from 1 ESMF state to another.
!   
! !DESCRIPTION:
!   
!  Initialize the re-distribute the fields that are shared on two ESMF States 
!  from one to the other.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_reDistInit( eStateIn, eStateOut, EBundleIn, EBundleOut, eRoute, &
                                   noGridCheck  )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),       intent(INOUT) :: EStateIn     ! Input state to redist 
    type(ESMF_State),       intent(INOUT) :: EStateOut    ! Output state to redist to
    type(ESMF_Bundle),      intent(INOUT) :: EBundleIn    ! In bundle of fields to redist
    type(ESMF_Bundle),      intent(INOUT) :: EBundleOut   ! Out bundle of fields to redist
    type(ESMF_RouteHandle), intent(OUT)   :: ERoute       ! Comm. handle to do the redist
    logical, optional,      intent(IN)    :: noGridCheck  ! Do NOT check grids
!EOP
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_redistInit) "
    integer                              :: rc                ! Return code
    integer                              :: i                 ! Indices
    integer                              :: nFields           ! Number of input fields
    character(len=SHR_KIND_CL), pointer  :: FieldNames(:)     ! Input field names
    type(ESMF_VM)                        :: currentVM         ! This Virtual machine
    type(ESMF_Grid)                      :: EGridIn           ! Input grid
    type(ESMF_Grid)                      :: EGridOut          ! Output grid
    type(ESMF_Field),           pointer  :: EFieldIn(:)       ! Input fields to redist 
    type(ESMF_Field),           pointer  :: EFieldOut(:)      ! Output fields to redist to

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   !----- Make sure list of fields on both states are the same ----
   if ( .not. eshr_estate_fieldsAreSame( eStateIn, eStateOut, notData=.true., &
   namesOnly=.true. ) ) &
      call shr_sys_abort( subName// &
                   " : names of fields on input states are different" )
   call eshr_estate_getSharedFNames( eStateIn, eStateOut, FieldNames, nFields, &
                                     onlyNeeded=OnlyUseNeededData )
   if ( nFields == 0 )then
      call shr_sys_abort( subName// " : no fields on input state" )
   end if
   ! ---- Get list of fields from both states ---
   allocate( EFieldIn (nFields) )
   allocate( EFieldOut(nFields) )
   do i = 1, nFields
      call ESMF_StateGetField( eStateIn, FieldName=FieldNames(i), field=EFieldIn(i), &
                               rc=rc )
      call eshr_rc_check( rc, subName// &
                       ': Error in getting in-field for redist between two states' )
      call ESMF_StateGetField( eStateOut, FieldName=FieldNames(i), field=EFieldOut(i), &
                               rc=rc )
      call eshr_rc_check( rc, subName// &
                       ': Error in getting out-field for redist between two states' )
   end do
   !------- Check that grids are at least reasonably close to the same ---
   if ( .not. present(noGridCheck) )then
      call ESMF_FieldGet( EFieldIn(1), grid=EGridIn, rc=rc )
      call eshr_rc_check( rc, subName// &
                          ': Error in getting in-field grid' )
      call ESMF_FieldGet( EFieldOut(1), grid=EGridOut, rc=rc )
      call eshr_rc_check( rc, subName// &
                          ': Error in getting out-field grid' )
      if ( .not. eshr_estate_Egrids2DAreSame( EGridIn, EGridOut, Exact=.false., &
                               NotData=.true. )  ) call shr_sys_abort( subName// &
                         " : Input grids of states to copy are different" )
   else if ( noGridCheck .eqv. .false. )then
      call shr_sys_abort( subName//" : If noGridCheck given it should be set to true" )
   end if
   !------ Create bundles from both sets of fields ----
   EBundleIn = ESMF_BundleCreate( nFields, EFieldIn, &
                                  name="Input field bundle", rc=rc )
   call eshr_rc_check( rc, subName// &
                       ': Error in creating in field bundle' )
   EBundleOut = ESMF_BundleCreate( nFields, EFieldOut, &
                                   name="Output field bundle", rc=rc )
   call eshr_rc_check( rc, subName// &
                       ': Error in creating out field bundle' )
   !------ Get the route to redist between the two bundles ---
   call ESMF_VMGetCurrent( currentVM, rc=rc)
   call eshr_rc_check( rc, subName// ': Error getting VM for redist operation' )
   call ESMF_BundleRedistStore( EBundleIn, EBundleOut, currentVM, &
                                eRoute, rc=rc)
   call eshr_rc_check( rc, subName// &
                       ': Error in getting route for redist between two states' )
   !------ Deallocate fieldnames ----
   deallocate( FieldNames )
   deallocate( EFieldIn   )
   deallocate( EFieldOut  )

END SUBROUTINE eshr_estate_redistInit
!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_reDist-- ReDistribute fields from one ESMF state to another.
!   
! !DESCRIPTION:
!   
!  Re-distribute the fields that are shared on two ESMF States from one to the other.
!  Must first do a eshr_estate_redistInit and use the input from that for this call.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_reDist( eBundleIn, eBundleOut, eRoute )
    use shr_sys_mod, only: shr_sys_flush
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Bundle),      intent(INOUT) :: EBundleIn  ! In bundle of fields to redist
    type(ESMF_Bundle),      intent(INOUT) :: EBundleOut ! Out bundle of fields to redist
    type(ESMF_RouteHandle), intent(INOUT) :: ERoute     ! Handle to route to do the redist
!EOP
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_redist) "
    type(ESMF_VM)               :: VM         ! This Virtual machine
    integer                     :: rc         ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   call ESMF_BundleRedist( srcBundle=EBundleIn, dstBundle=EBundleOut, &
                           routeHandle=eRoute, rc=rc)
   call eshr_rc_check( rc, subName// &
                       ': Error in doing redist between two states' )

END SUBROUTINE eshr_estate_redist
!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_mergeInit -- Initialize a Merge from 3 ESMF states to another.
!   
! !DESCRIPTION:
!   
!  Initialize the merge from ice, ocn, and land fields to output merged fields.
!  Will only merge fields with the correct all-surface name "x" on the output
!  state and that have component specific counterparts on ice, ocean, and land
!  states.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_mergeInit( eStateIce, eStateLnd, eStateOcn, eStateOut,       &
                                  eStateXao, FieldNames, IceFrac, LndFrac, OcnFrac, &
                                  nFields, nGrPts )
    use shr_string_mod,  only: shr_string_listGetNum, shr_string_listGetName
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),       intent(INOUT) :: EStateLnd       ! Input land state to merge
    type(ESMF_State),       intent(INOUT) :: EStateIce       ! Input ice state to merge
    type(ESMF_State),       intent(INOUT) :: EStateOcn       ! Input ocean state to merge
    type(ESMF_State),       intent(INOUT) :: EStateOut       ! Output state merging
    type(ESMF_State),       intent(INOUT) :: EStateXao       ! Input atm/ocean flux to merge
    character(len=SHR_KIND_CL), pointer   :: FieldNames(:,:) ! Merge Fieldnames
    real(SHR_KIND_R8),          pointer   :: IceFrac(:)      ! Ice area fraction
    real(SHR_KIND_R8),          pointer   :: LndFrac(:)      ! Land area fraction
    real(SHR_KIND_R8),          pointer   :: OcnFrac(:)      ! Ocean area fraction
    integer,                intent(OUT)   :: nFields         ! Number of fields to merge
    integer,                intent(OUT)   :: nGrPts          ! Number of grid points
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_mergeInit) "
    character(len=SHR_KIND_CS), parameter :: LndPrefixList = & ! Land prefix to replace
              "Sl_:Fall_"
    character(len=SHR_KIND_CS), parameter :: OcnPrefixList = & ! Ocean prefix to replace
              "So_:Faoo_"
    character(len=SHR_KIND_CS), parameter :: XaoPrefixList = & ! A/O flx prefix to replace
              "So_:Faox_"
    character(len=SHR_KIND_CS), parameter :: IcePrefixList = & ! Ice prefix to replace
              "Si_:Faii_"
    character(len=SHR_KIND_CS), parameter :: OutPrefixList = & ! Output prefix to replace
              "Sx_:Faxx_"
    character(len=SHR_KIND_CL), parameter :: IgnoreList    = & ! List of fields to ignore
              "Sx_snowh:Sx_fv:Sx_ram1:Faxx_flxdst1:Faxx_flxdst2:Faxx_flxdst3:Faxx_flxdst4"
    character(len=SHR_KIND_CL), parameter :: iceFName      = & ! Fractional ice field name
              "Sx_ifrac"
    character(len=SHR_KIND_CL), parameter :: ocnFName      = & ! Fractional ocn field name
              "Sx_ofrac"
    character(len=SHR_KIND_CL), parameter :: lndFName      = & ! Fractional lnd field name
              "Sx_lfrac"
    integer                               :: rc                ! Return code
    integer                               :: ifld, i, j, n     ! Indices
    integer                               :: nOutFields        ! Number of output fields
    character(len=SHR_KIND_CL), pointer   :: OutFieldNames(:)  ! Output field names
    character(len=SHR_KIND_CL)            :: Prefix            ! Prefix testing
    character(len=SHR_KIND_CL)            :: OutPrefix         ! Output prefix
    character(len=SHR_KIND_CL)            :: Fieldn            ! Fieldname to ignore
    character(len=SHR_KIND_CL)            :: IceFieldn         ! Ice Fieldname
    character(len=SHR_KIND_CL)            :: LndFieldn         ! Land Fieldname
    character(len=SHR_KIND_CL)            :: OcnFieldn         ! Ocean Fieldname
    character(len=SHR_KIND_CL), pointer   :: IceFieldNames(:)  ! List of ice Fieldnames
    character(len=SHR_KIND_CL), pointer   :: OcnFieldNames(:)  ! List of ocean Fieldnames
    character(len=SHR_KIND_CL), pointer   :: LndFieldNames(:)  ! List of land Fieldnames
    character(len=SHR_KIND_CL), pointer   :: OceanType(:)      ! If field in ocean or xao
    real(SHR_KIND_R8),          pointer   :: edataIn(:)        ! Pointer to data in state
    integer                               :: itype             ! Type index
    type(ESMF_StateItemType)              :: eItemType         ! ESMF Item type

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
          call eshr_estate_getFieldNames( eStateOut,  OutFieldNames,  nOutFields, &
                                          onlyNeeded=OnlyUseNeededData  )
          if ( nOutFields > 0 ) allocate( IceFieldNames(nOutFields) )
          if ( nOutFields > 0 ) allocate( OcnFieldNames(nOutFields) )
          if ( nOutFields > 0 ) allocate( LndFieldNames(nOutFields) )
          if ( nOutFields > 0 ) allocate( OceanType    (nOutFields) )
          nFields = 0
          nGrPts  = 0
fldloop:  do ifld = 1, nOutFields
              do i = 1, shr_string_listGetNum( OutPrefixList )
                 call shr_string_listGetName(OutPrefixList, i, Prefix, rc=rc )
                 !--------------------------------------------------------------
                 ! If begining part of fieldname matches one of the prefixes -- process it
                 !--------------------------------------------------------------
                 if ( index(OutFieldNames(ifld),trim(Prefix)) == 1 )then
                    !--- Make sure doesn't match one of the fieldnames to ignore ---
                    do j = 1, shr_string_listGetNum( IgnoreList )
                        call shr_string_listGetName( IgnoreList, j, Fieldn, rc=rc )
                        if ( trim(OutFieldNames(ifld)) == trim(Fieldn) ) cycle fldloop
                    end do
                    !--- Make sure doesn't match fractional names ---
                    if ( trim(OutFieldNames(ifld)) == trim(IceFname) ) cycle fldloop
                    if ( trim(OutFieldNames(ifld)) == trim(OcnFname) ) cycle fldloop
                    if ( trim(OutFieldNames(ifld)) == trim(LndFname) ) cycle fldloop
                    !--- If make it to here -- increment number of merged fields ---
                    nFields = nFields + 1
                    OutFieldnames(nFields) = OutFieldnames(ifld)
                    !--- Replace output prefix with ice prefix ---
                    call shr_string_listGetName(IcePrefixList, i, OutPrefix, rc=rc )
                    n = len_trim(Prefix)
                    if ( len_trim(OutPrefix) /= n )then
                       call shr_sys_abort( SubName// &
                                " : Error, length of prefixes different" )
                    end if
                    IceFieldnames(nFields) = trim(OutFieldNames(ifld))
                    IceFieldnames(nFields)(1:n) = trim(OutPrefix)
                    !--- Now make sure fields exist on ice state ---
                    call ESMF_StateGetItemInfo( eStateIce, name=IceFieldnames(nFields), &
                      stateItemType=eItemType, rc=rc )
                    if ( eItemType == ESMF_STATEITEM_NOTFOUND )then
                       nFields = nFields - 1
                       cycle fldloop
                    end if
                    call eshr_rc_check( rc,subName// &
                                  ': Error in getting info on field from ice-state')
                    !--- Replace output prefix with ocn prefix ---
                    call shr_string_listGetName(OcnPrefixList, i, OutPrefix, rc=rc )
                    n = len_trim(Prefix)
                    if ( len_trim(OutPrefix) /= n )then
                       call shr_sys_abort( SubName// &
                                 " : Error, length of prefixes different" )
                    end if
                    OcnFieldnames(nFields) = trim(OutFieldNames(ifld))
                    OcnFieldnames(nFields)(1:n) = trim(OutPrefix)
                    !---Now make sure fields exist on ocn state ---
                    call ESMF_StateGetItemInfo( eStateOcn, name=OcnFieldnames(nFields), &
                      stateItemType=eItemType, rc=rc )
                    !--- Check xao if not found ---
                    if ( eItemType == ESMF_STATEITEM_NOTFOUND )then
                       !--- Replace output prefix with A/O flux prefix ---
                       call shr_string_listGetName(XaoPrefixList, i, OutPrefix, rc=rc )
                       n = len_trim(Prefix)
                       if ( len_trim(OutPrefix) /= n )then
                          call shr_sys_abort( SubName//" : Error, length of prefixes different" )
                       end if
                       OcnFieldnames(nFields) = trim(OutFieldNames(ifld))
                       OcnFieldnames(nFields)(1:n) = trim(OutPrefix)
                       call ESMF_StateGetItemInfo( eStateXao, &
                            name=OcnFieldnames(nFields), stateItemType=eItemType, rc=rc )
                       if ( eItemType == ESMF_STATEITEM_NOTFOUND )then
                          nFields = nFields - 1
                          cycle fldloop
                       else
                          OceanType(nFields) = "xao"
                       end if
                    else
                       OceanType(nFields) = "ocn"
                    end if
                    call eshr_rc_check( rc,subName// &
                          ': Error in getting info on field from ocn-state')
                    !--- Replace output prefix with lnd prefix ---
                    call shr_string_listGetName(LndPrefixList, i, OutPrefix, rc=rc )
                    n = len_trim(Prefix)
                    if ( len_trim(OutPrefix) /= n )then
                       call shr_sys_abort( SubName// &
                                           " : Error, length of prefixes different" )
                    end if
                    LndFieldnames(nFields) = trim(OutFieldNames(ifld))
                    LndFieldnames(nFields)(1:n) = trim(OutPrefix)
                    !--- Now make sure fields exist on lnd state ---
                    call ESMF_StateGetItemInfo( eStateLnd, name=LndFieldnames(nFields), &
                                                stateItemType=eItemType, rc=rc )
                    if ( eItemType == ESMF_STATEITEM_NOTFOUND )then
                       nFields = nFields - 1
                       cycle fldloop
                    end if
                    call eshr_rc_check( rc,subName// &
                           ': Error in getting info on field from lnd-state')
                 end if
              end do
          end do fldloop
          !---Now loop through list of fields to get the merge list of fields ---
          if ( nFields > 0 ) allocate(FieldNames(5,nFields))
          do ifld = 1, nFields
             FieldNames(OUT_IDX,ifld) = trim(OutFieldNames(ifld))
             FieldNames(OCN_IDX,ifld) = trim(OcnFieldNames(ifld))
             FieldNames(ICE_IDX,ifld) = trim(IceFieldNames(ifld))
             FieldNames(LND_IDX,ifld) = trim(LndFieldNames(ifld))
             FieldNames(XAO_IDX,ifld) = trim(OceanType(ifld))
          end do
          !---------------------------------------------------------------------
          ! Do error checking
          !---------------------------------------------------------------------
          do ifld = 1, nFields
             do itype = 1, 4
                if ( itype == OUT_IDX )then
                   call ESMF_StateGetDataPointer( state=eStateOut,               &
                                                  itemName=OutFieldNames(ifld),  &
                                                  dataPointer=edataIn,           &
                                                  copyflag=ESMF_DATA_REF, rc=rc )
                   if ( ifld == 1 ) nGrPts = size(edataIn)
                else if ( itype == ICE_IDX )then
                   call ESMF_StateGetDataPointer( state=eStateIce,                 &
                                                  itemName=FieldNames(itype,ifld), &
                                                  dataPointer=edataIn,             &
                                                  copyflag=ESMF_DATA_REF, rc=rc )
                else if ( itype == LND_IDX )then
                   call ESMF_StateGetDataPointer( state=eStateLnd,                 &
                                                  itemName=FieldNames(itype,ifld), &
                                                  dataPointer=edataIn,             &
                                                  copyflag=ESMF_DATA_REF, rc=rc )
                else if ( itype == OCN_IDX )then
                   if ( trim(FieldNames(XAO_IDX,ifld)) == "ocn" )then
                       call ESMF_StateGetDataPointer( state=eStateOcn,             &
                                                  itemName=FieldNames(itype,ifld), &
                                                  dataPointer=edataIn,             &
                                                  copyflag=ESMF_DATA_REF, rc=rc )
                   else
                       call ESMF_StateGetDataPointer( state=eStateXao,             &
                                                  itemName=FieldNames(itype,ifld), &
                                                  dataPointer=edataIn,             &
                                                  copyflag=ESMF_DATA_REF, rc=rc )
                   end if
                end if
                call eshr_rc_check( rc, subName//': Error in getting F90 ptr to state '// &
                                 trim(FieldNames(itype,ifld)) )
                n = size(edataIn)
                if ( n /= nGrPts )then
                   call shr_sys_abort( SubName//" : Error, inconsistent grid size" )
                end if
                nullify( edataIn )
             end do
          end do
          !---- Get data pointers for areal fractions -----
          call ESMF_StateGetDataPointer( state=eStateOut,                 &
                                         itemName=IceFname,               &
                                         dataPointer=IceFrac,             &
                                         copyflag=ESMF_DATA_REF, rc=rc )
          call eshr_rc_check( rc, subName//': Error in getting icefrac F90 ptr to state' )
          call ESMF_StateGetDataPointer( state=eStateOut,                 &
                                         itemName=LndFname,               &
                                         dataPointer=LndFrac,             &
                                         copyflag=ESMF_DATA_REF, rc=rc )
          call eshr_rc_check( rc, subName//': Error in getting lndfrac F90 ptr to state' )
          call ESMF_StateGetDataPointer( state=eStateOut,                 &
                                         itemName=OcnFname,               &
                                         dataPointer=OcnFrac,             &
                                         copyflag=ESMF_DATA_REF, rc=rc )
          call eshr_rc_check( rc, subName//': Error in getting ocnfrac F90 ptr to state' )

          if ( nFields > 0 ) deallocate( OutFieldNames )
          if ( nFields > 0 ) deallocate( IceFieldNames )
          if ( nFields > 0 ) deallocate( OcnFieldNames )
          if ( nFields > 0 ) deallocate( LndFieldNames )
          if ( nFields > 0 ) deallocate( OceanType     )

END SUBROUTINE eshr_estate_mergeInit

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_merge -- Merge fields from input ESMF states to another.
!   
! !DESCRIPTION:
!   
!  Merge the fields that are shared on the input land, ice, and ocean ESMF States 
!  to the output ESMF state. Must first do a eshr_estate_mergeInit and use the 
!  input from that for this call.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_merge( EStateIce, EStateLnd, EStateOcn, EStateXao, nFields,      &
                              nGrPts, FieldNames, IceFrac, LndFrac, OcnFrac, EStateOut, &
                              check )
    use shr_const_mod, only: SHR_CONST_SPVAL
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),           intent(IN) :: EStateLnd       ! Input land state to merge
    type(ESMF_State),           intent(IN) :: EStateIce       ! Input ice state to merge
    type(ESMF_State),           intent(IN) :: EStateOcn       ! Input ocean state to merge
    type(ESMF_State),           intent(IN) :: EStateXao       ! Input atm/ocean flux to merge
    integer,                    intent(IN) :: nFields         ! No. of fields to merge
    integer,                    intent(IN) :: nGrPts          ! No. of grid points
    character(len=SHR_KIND_CL), intent(IN) :: FieldNames(5,nFields) ! Merge Fieldnames
    real(SHR_KIND_R8),          intent(IN) :: IceFrac(nGrPts) ! Ice area fraction
    real(SHR_KIND_R8),          intent(IN) :: LndFrac(nGrPts) ! Land area fraction
    real(SHR_KIND_R8),          intent(IN) :: OcnFrac(nGrPts) ! Ocean area fraction
    type(ESMF_State),        intent(INOUT) :: EStateOut       ! Output state merging
    logical,          optional, intent(IN) :: check           ! Extra checking
!EOP
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_merge) "
    integer                     :: rc            ! Return code
    integer                     :: ifld, ig      ! Indices
    logical                     :: checking      ! If doing extra checking
    real(SHR_KIND_R8)           :: AreaFrac      ! Fractional surface
    real(SHR_KIND_R8), pointer  :: edataIce(:)   ! Pointer to data in ice input state
    real(SHR_KIND_R8), pointer  :: edataOcn(:)   ! Pointer to data in ocean input state
    real(SHR_KIND_R8), pointer  :: edataLnd(:)   ! Pointer to data in land input state
    real(SHR_KIND_R8), pointer  :: edataOut(:)   ! Pointer to data in output state
    real(SHR_KIND_R8)           :: iceMissVal    ! Sea-Ice missing value
    real(SHR_KIND_R8)           :: ocnMissVal    ! Ocean missing value
    real(SHR_KIND_R8)           :: lndMissVal    ! Land missing value

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
    if ( present(check) )then
       checking = check
    else
       checking = .false.
    end if
    iceMissVal = SHR_CONST_SPVAL
    lndMissVal = SHR_CONST_SPVAL
    ocnMissVal = SHR_CONST_SPVAL
    ! This is a fairly complex loop -- run it through with Open-MP
    !$OMP PARALLEL DO PRIVATE (ifld,edataOut,edataIce,edataLnd,edataOcn,ig,AreaFrac,rc)
    do ifld = 1, nFields
       !--- Output data pointers ---
       call ESMF_StateGetDataPointer( state=eStateOut,               &
                                      itemName=FieldNames(OUT_IDX,ifld),  &
                                      dataPointer=edataOut,          &
                                      copyflag=ESMF_DATA_REF, rc=rc )
       call eshr_rc_check( rc, subName//': Error in getting F90 ptr to output state' )
       !--- Ice data pointers ---
       call ESMF_StateGetDataPointer( state=eStateIce,               &
                                      itemName=FieldNames(ICE_IDX,ifld),  &
                                      dataPointer=edataIce,          &
                                      copyflag=ESMF_DATA_REF, rc=rc )
       call eshr_rc_check( rc, subName//': Error in getting F90 ptr to ice state' )
       !--- Land data pointers ---
       call ESMF_StateGetDataPointer( state=eStateLnd,               &
                                      itemName=FieldNames(LND_IDX,ifld),  &
                                      dataPointer=edataLnd,          &
                                      copyflag=ESMF_DATA_REF, rc=rc )
       call eshr_rc_check( rc, subName//': Error in getting F90 ptr to land state' )
       if ( trim(FieldNames(XAO_IDX,ifld)) == "ocn" )then
          !--- Ocean data pointers ---
          call ESMF_StateGetDataPointer( state=eStateOcn,               &
                                         itemName=FieldNames(OCN_IDX,ifld),  &
                                         dataPointer=edataOcn,          &
                                         copyflag=ESMF_DATA_REF, rc=rc )
          call eshr_rc_check( rc, subName//': Error in getting F90 ptr to ocean state' )
       else
          !--- Xao data pointers ---
          call ESMF_StateGetDataPointer( state=eStateXao,               &
                                         itemName=FieldNames(OCN_IDX,ifld),  &
                                         dataPointer=edataOcn,          &
                                         copyflag=ESMF_DATA_REF, rc=rc )
          call eshr_rc_check( rc, subName//': Error in getting F90 ptr to xao state' )
       end if
       !--- Initialize to zero ---
       eDataOut(:) = 0.0_SHR_KIND_R8
       do ig = 1, nGrPts
          !--- Add in land data ---
          AreaFrac = LndFrac(ig)
          if ( checking .and. (AreaFrac < 0.0_SHR_KIND_R8) .or. &
                              (AreaFrac > 1.0_SHR_KIND_R8) )then
             call shr_sys_abort( subName//" : Error land fraction out of bounds" )
          end if
          if ( AreaFrac > 0._SHR_KIND_R8 )then
             eDataOut(ig) = eDataOut(ig) + eDataLnd(ig) * AreaFrac
             if ( checking .and. (eDataLnd(ig) == lndMissVal) )then
                call shr_sys_abort( subName//" : Error, using land missing value" )
             end if
          end if
          !--- Add in ice data ---
          AreaFrac = IceFrac(ig)
          if ( checking .and. (AreaFrac < 0.0_SHR_KIND_R8) .or. &
                              (AreaFrac > 1.0_SHR_KIND_R8) )then
             call shr_sys_abort( subName//" : Error ice fraction out of bounds" )
          end if
          if ( AreaFrac > 0._SHR_KIND_R8 )then
             eDataOut(ig) = eDataOut(ig) + eDataIce(ig) * AreaFrac
             if ( checking .and. (eDataIce(ig) == iceMissVal) )then
                call shr_sys_abort( subName//" : Error, using ice missing value" )
             end if
          end if
          !--- Add in ocean data ---
          AreaFrac = OcnFrac(ig)
          if ( checking .and. (AreaFrac < 0.0_SHR_KIND_R8) .or. &
                              (AreaFrac > 1.0_SHR_KIND_R8) )then
             call shr_sys_abort( subName//" : Error ocean fraction out of bounds" )
          end if
          if ( AreaFrac > 0._SHR_KIND_R8 )then
             eDataOut(ig) = eDataOut(ig) + eDataOcn(ig) * AreaFrac
             if ( checking .and. (eDataOcn(ig) == ocnMissVal) )then
                call shr_sys_abort( subName//" : Error, using ocean missing value" )
             end if
          end if
       end do
       nullify(eDataIce)
       nullify(eDataLnd)
       nullify(eDataOcn)
       nullify(eDataOut)
    end do


END SUBROUTINE eshr_estate_merge

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_destroy -- Destroy a state and the fields it contains
!   
! !DESCRIPTION:
!   
!  Destroy a state and the fields that are contained in it.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_destroy( EState )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),           intent(INOUT) :: EState     ! State to destroy
!EOP
    !----- local -----
    character(len=*), parameter            :: subName = " (eshr_estate_destroy) "
    integer                                :: rc            ! Return code
    integer                                :: ifld          ! Indices
    integer                                :: nFields       ! Number of fields
    character(len=SHR_KIND_CL), pointer    :: FieldNames(:) ! Fieldsnames on this state
    type(ESMF_Field),           pointer    :: Fields(:)     ! Fields on state

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
    call eshr_estate_getFieldNames( eState,  FieldNames,  nFields  )

    if ( NFields > 0 )then
       allocate( Fields(nFields) )
       do ifld = 1, nFields
          call ESMF_StateGetField( eState, FieldNames(ifld), Fields(ifld), rc=rc )
          call eshr_rc_check( rc, subName//': Error in getting field from state' )
       end do
       call ESMF_StateDestroy( eState )
       do ifld = 1, nFields
          call ESMF_FieldDestroy( Fields(ifld) )
       end do
       deallocate( Fields )
    else
       call ESMF_StateDestroy( eState )
    end if

END SUBROUTINE eshr_estate_destroy

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_getFilename -- get the given filename from the ESMF state
!   
! !DESCRIPTION:
!   
!  Get the filename from the ESMF state.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_getFilename( eState, FileType, Filename )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State), intent(IN)  :: eState     ! Input ESMF State
    character(len=*), intent(IN)  :: FileType   ! Type of file looking for
    character(len=*), intent(OUT) :: Filename   ! File-path of file looking for
!EOP
    !----- local -----
    character(len=*), parameter  :: subName = " (eshr_estate_getFilename) "
    integer                      :: rc           ! ESMF Return code
    character(len=len(Filename)) :: FilenameIn   ! File-path of file looking for
!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
    if ( len_trim(FileType) == 0 ) &
      call shr_sys_abort( subName//" : Error filetype is empty" )
    call ESMF_StateGetAttribute( eState, name=FileType, value=FilenameIn, rc=rc )
    call eshr_rc_check( rc, "Error, on getting filename from ESMF State: "//FileType )
    Filename = trim(FilenameIn)

END SUBROUTINE eshr_estate_getFilename

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_putFilename -- put the given filename on the ESMF state
!   
! !DESCRIPTION:
!   
!  Put the given filename on the ESMF_State
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_putFilename( FileType, Filename, eState )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(IN)    :: FileType   ! Type of file looking for
    character(len=*), intent(IN)    :: Filename   ! File-path of file looking for
    type(ESMF_State), intent(INOUT) :: eState     ! Output ESMF State
!EOP
    !----- local -----
    character(len=*), parameter            :: subName = " (eshr_estate_putFilename) "
    integer :: rc   ! ESMF Return code
    character(len=len(Filename)+25) :: value     ! File-path of file looking for

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------
    value=trim(Filename)
    if ( len_trim(FileType) == 0 ) &
      call shr_sys_abort( subName//" : Error filetype is empty" )
    call ESMF_StateSetAttribute( eState, name=FileType, value=value, rc=rc )
    call eshr_rc_check( rc, "Error, on setting filename to ESMF State: "//FileType )

END SUBROUTINE eshr_estate_putFilename

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_printAttributes -- print the list of attributes on the state
!   
! !DESCRIPTION:
!   
!  Print out the list of attributes on the state
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_printAttributes( eState )
    use shr_kind_mod,   only: SHR_KIND_R8
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State), intent(IN)  :: eState     ! Input ESMF State
!EOP
    !----- local -----
    integer             :: rc     ! ESMF Return code
    integer             :: count  ! Count of attributes
    integer             :: i      ! Indices
    character(len=256)  :: name   ! Name of state or attribute
    character(len=256)  :: value  ! Value of attribute
    type(ESMF_DataType) :: dtype  ! Data-type
    type(ESMF_Logical)  :: lvalue ! ESMF Logical value
    integer             :: ivalue ! integer attribute value
    real(SHR_KIND_R8)   :: rvalue ! real attribute value

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call ESMF_StateGet( eState, name=name, rc=rc )
    call eshr_rc_check( rc, "Error in Get estate name in printAttributes" )
    write(6,'(A,A)') "Name of state to print attributes on: ", trim(name)
    call ESMF_StateGetAttributeCount( eState, count=count, rc=rc )
    call eshr_rc_check( rc, "Error in Get Attribute count in printAttributes" )
    write(6,'(A,I3)') "Number of attributes = ", count
    !---- Loop through all attributes on state ---------------------------------
    do i = 1, count
       call ESMF_StateGetAttributeInfo( eState, attributeIndex=i, name=name, &
                                        datatype=dtype, rc=rc )
       call eshr_rc_check( rc, "Error in Get Attribute info in printAttributes" )
       write(6,'(A,I3,A,A,$)') "Attribute #: ", i, " = ", trim(name)
       !--- Character data ----
       if (      dtype == ESMF_DATA_CHARACTER )then
           call ESMF_StateGetAttribute( eState, name=name, value=value, rc=rc )
           call eshr_rc_check( rc, "Error in Get Attribute in printAttributes" )
           write(6,'(A,A)') ' = ', trim(value)
       !--- Integer data ----
       else if ( dtype == ESMF_DATA_INTEGER  )then
           call ESMF_StateGetAttribute( eState, name=name, value=ivalue, rc=rc )
           call eshr_rc_check( rc, "Error in Get Attribute in printAttributes" )
           write(6,'(A,I9)') ' = ', ivalue
       !--- Logical data ----
       else if ( dtype == ESMF_DATA_LOGICAL )then
           call ESMF_StateGetAttribute( eState, name=name, value=lvalue, rc=rc )
           call eshr_rc_check( rc, "Error in Get Attribute in printAttributes" )
           if ( lvalue == ESMF_TRUE )then
              write(6,'(A,L1)') ' = ', .true.
           else
              write(6,'(A,L1)') ' = ', .false.
           end if
       !--- Real data ----
       else if ( dtype == ESMF_DATA_REAL )then
           call ESMF_StateGetAttribute( eState, name=name, value=rvalue, rc=rc )
           call eshr_rc_check( rc, "Error in Get Attribute in printAttributes" )
           write(6,'(A,g15.5)') ' = ', rvalue
       !--- Bad data ----
       else
           write(6,'(A)') 'Bad-type'
       end if
    end do

END SUBROUTINE eshr_estate_printAttributes

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_printFieldAttr -- print attributes on the fields on this state
!   
! !DESCRIPTION:
!   
!  Print out the list of attributes on each field that is on this state.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_printFieldAttr( eState )
    use shr_kind_mod,   only: SHR_KIND_R8, SHR_KIND_CL
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State), intent(IN)  :: eState     ! Input ESMF State
!EOP
    !----- local -----
    integer                                 :: rc            ! ESMF Return code
    integer                                 :: count         ! Count of attributes
    integer                                 :: nFields       ! Number of fields on state
    integer                                 :: i, j          ! Indices
    character(len=SHR_KIND_CL)              :: name          ! Name of state or attribute
    character(len=SHR_KIND_CL)              :: value         ! Value of attribute
    type(ESMF_DataType)                     :: dtype         ! Data-type
    type(ESMF_Logical)                      :: lvalue        ! ESMF Logical value
    integer                                 :: ivalue        ! integer attribute value
    real(SHR_KIND_R8)                       :: rvalue        ! real attribute value
    character(len=SHR_KIND_CL), pointer     :: FieldNames(:) ! Names of fields in state
    type(ESMF_Field)                        :: Field         ! Field to query
    character(len=*),parameter :: subName =   "(eshr_estate_printFieldAttr) "

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call ESMF_StateGet( eState, name=name, rc=rc )
    call eshr_rc_check( rc, subName//"Error in Get estate name" )
    write(6,'(A,A)') "Name of state to print field attributes on: ", trim(name)
    !
    !------  Get number and list of field names in state ------------------------
    !
    call eshr_estate_getFieldNames( eState, FieldNames, nFields )
    do j = 1, nFields
       call ESMF_StateGetField( eState, Fieldnames(j), field, rc=rc)
       call eshr_rc_check( rc, subName//' : Error getting field on input state' )
       call ESMF_FieldGetAttributeCount( Field, count=count, rc=rc )
       call eshr_rc_check( rc, subName//"Error in Get Attribute count" )
       write(6,'(T5,A, A)') "Field: ", trim(FieldNames(j))
       do i = 1, count
          call ESMF_FieldGetAttributeInfo( Field, attributeIndex=i, name=name, &
                                           datatype=dtype, rc=rc )
          call eshr_rc_check( rc, subName//"Error in Get Attribute info" )
          write(6,'(T10,A,$)') trim(name)
          if (      dtype == ESMF_DATA_CHARACTER )then
              call ESMF_FieldGetAttribute( Field, name=name, value=value, rc=rc )
              call eshr_rc_check( rc, subName//"Error in Get Attribute" )
              write(6,'(A,A)') ' = ', trim(value)
          else if ( dtype == ESMF_DATA_INTEGER  )then
              call ESMF_FieldGetAttribute( Field, name=name, value=ivalue, rc=rc )
              call eshr_rc_check( rc, subName//"Error in Get Attribute" )
              write(6,'(A,I9)') ' = ', ivalue
          else if ( dtype == ESMF_DATA_LOGICAL )then
              call ESMF_FieldGetAttribute( Field, name=name, value=lvalue, rc=rc )
              call eshr_rc_check( rc, subName//"Error in Get Attribute" )
              if ( lvalue == ESMF_TRUE )then
                 write(6,'(A,L1)') ' = ', .true.
              else
                 write(6,'(A,L1)') ' = ', .false.
              end if
          else if ( dtype == ESMF_DATA_REAL )then
              call ESMF_FieldGetAttribute( Field, name=name, value=rvalue, rc=rc )
              call eshr_rc_check( rc, subName//"Error in Get Attribute" )
              write(6,'(A,g15.5)') ' = ', rvalue
          else
              write(6,'(A)') 'Bad-type'
          end if
       end do
    end do
    if ( nFields > 0 ) deallocate( FieldNames )

END SUBROUTINE eshr_estate_printFieldAttr

!===============================================================================
! !IROUTINE: eshr_estate_ELog2FLog -- convert ESMF_Logical data into FORTRAN Logical
!   
! !DESCRIPTION:
!   
!  Convert the input ESMF Logical data value into FORTRAN Logical
!      
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION eshr_estate_ELog2FLog( ELogValue )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Logical),  intent(IN)  :: ELogValue ! Input ESMF Logical data
!EOP
    !----- local -----
!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( ELogValue == ESMF_TRUE )then
      eshr_estate_ELog2FLog = .true.
   else
      eshr_estate_ELog2FLog = .false.
   end if
END FUNCTION eshr_estate_ELog2FLog

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_FLog2ELog -- convert FORTRAN Logical data into ESMF Logical
!   
! !DESCRIPTION:
!   
!  Convert the input FORTRAN Logical data value into ESMF Logical
!      
! !INTERFACE: ------------------------------------------------------------------
type(ESMF_Logical) FUNCTION eshr_estate_FLog2ELog( FLogValue )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    logical,  intent(IN)  :: FLogValue ! Input FORTRAN  Logical data
!EOP
    !----- local -----
!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( FLogValue )then
      eshr_estate_FLog2ELog = ESMF_TRUE
   else
      eshr_estate_FLog2ELog = ESMF_FALSE
   end if
END FUNCTION eshr_estate_FLog2ELog

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_getFieldNames -- Get the field names on the state
!   
! !DESCRIPTION:
!   
!  Get the list of field names on the input state.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_getFieldNames( eState, FieldNames, nFields, onlyNeeded )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),       intent(IN)  :: eState        ! Input ESMF State
    character(len=SHR_KIND_CL), pointer :: Fieldnames(:) ! Output fieldnames
    integer,                intent(OUT) :: nFields       ! Number of output fields
    logical, optional,      intent(IN)  :: onlyNeeded    ! Only get names of needed fields
!EOP
    !----- local ----
    character(len=*), parameter :: subName = " (eshr_estate_getFieldNames) "
    integer                                 :: rc            ! Return code
    integer                                 :: i, j          ! Indices
    integer                                 :: nItems        ! Number of items in state
    character(len=SHR_KIND_CL), allocatable :: names(:)      ! Names of items in state
    type(ESMF_StateItemType),   allocatable :: itemType(:)   ! Type of items in state
    type(ESMF_DataType)                     :: dataType      ! Type of data in array
    type(ESMF_NeededFlag)                   :: need          ! If needed or not
    logical                                 :: neededOnly    ! If only extract needed
    !
    if ( present(onlyNeeded) )then
       neededOnly = onlyNeeded
    else
       neededOnly = .false.
    end if
    !
    !------  Get number and list of field names in state ------------------------
    !
    call ESMF_StateGet( state=EState, itemCount=nItems, rc=rc )
    call eshr_rc_check( rc, subName//' : Error getting no. items on input state' )
    if ( nItems == 0 )then
       nFields = 0
    else
       allocate( itemType(nItems) )
       allocate( names(nItems) )
       call ESMF_StateGet( state=EState, itemNameList=names, &
                           stateitemtypeList=itemType, rc=rc )
       call eshr_rc_check( rc, subName//' : Error getting information on input state' )
       nFields = 0
       do j = 1, nItems
          if ( itemType(j) == ESMF_STATEITEM_FIELD )then
             nFields = nFields + 1
             !----  Remove field if only needed and this field is not needed ------
             if ( neededOnly )then
                 call ESMF_StateGetNeeded( EState, itemname=names(j), neededflag=need, &
                                           rc=rc )
                 call eshr_rc_check( rc, subname// &
                                     " : Error in getting needed flag on state" )
                if ( need == ESMF_NOTNEEDED) nFields = nFields - 1
             endif
          end if
       end do
       if ( nFields > 0 )then
          allocate( FieldNames(nFields) )
          i = 0
          do j = 1, nItems
             call ESMF_StateGetNeeded( EState, itemname=names(j), neededflag=need, &
                                       rc=rc )
             call eshr_rc_check( rc, subname// &
                                 " : Error in getting needed flag on state" )
             if ( itemType(j) == ESMF_STATEITEM_FIELD )then
                if ( .not. neededOnly .or. (need == ESMF_NEEDED) )then
                   i = i + 1
                   FieldNames(i) = names(j)
                end if
             end if
          end do
          deallocate( itemType )
       end if
       deallocate( names )
    end if
END SUBROUTINE eshr_estate_getFieldNames

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_getFieldsNotNeed -- Get list of fields not needed
!   
! !DESCRIPTION:
!   
!  Get the list of field names that are not needed.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_getFieldsNotNeed( eState, FieldsNotNeeded, nFieldsNotNeed )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State), intent(IN)  :: eState          ! Input ESMF State
    character(len=*), intent(OUT) :: FieldsNotNeeded ! Output fields not needed
    integer,          intent(OUT) :: nFieldsNotNeed  ! Number of fields not needed
!EOP
    !----- local ----
    character(len=*), parameter :: subName = " (eshr_estate_getFieldsNotNeed) "
    integer                                 :: rc            ! Return code
    integer                                 :: i, j, n       ! Indices
    integer                                 :: nFields       ! Number of fields in state
    character(len=SHR_KIND_CL), pointer :: Fieldnames(:)     ! Output fieldnames
    type(ESMF_NeededFlag)                   :: need          ! If needed or not
    !
    call eshr_estate_getFieldNames( eState, FieldNames, nFields )
    j = 0
    FieldsNotNeeded = " "
    do i = 1, nFields
       call ESMF_StateGetNeeded( EState, itemname=FieldNames(i), neededflag=need, &
                                 rc=rc )
       call eshr_rc_check( rc, subname//" : Error in getting needed flag on state" )
       if ( need == ESMF_NOTNEEDED)then
          j = j + 1
          if ( j > 1 ) FieldsNotNeeded = trim(FieldsNotNeeded)//":"
          !--- Test that length of input string is long enough ----
          n = len_trim(FieldsNotNeeded) + len_trim(FieldNames(i)) + 1
          if ( n > len(FieldsNotNeeded) ) &
             call shr_sys_abort( subName// &
                 " : Error,  input field length too short for non-needed fields" )
          !---- Add next NOT needed field to the end of the list ---
          FieldsNotNeeded = trim(FieldsNotNeeded)//trim(FieldNames(i))
       end if
    end do
    nFieldsNotNeed = j
    if ( nFields /= 0 ) deallocate( FieldNames )
END SUBROUTINE eshr_estate_getFieldsNotNeed

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_markRestNotNeed  -- Mark other fields entered as not needed
!   
! !DESCRIPTION:
!   
!  Mark the other fields on state as not needed.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_markRestNotNeed( eState, FieldsNeeded, nFieldsNeeded )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State), intent(INOUT) :: eState          ! Input ESMF State
    character(len=*), intent(IN)    :: FieldsNeeded(:) ! Fields needed
    integer,          intent(IN)    :: nFieldsNeeded   ! Number of fields needed
!EOP
    !----- local ----
    character(len=*), parameter :: subName = " (eshr_estate_markRestNotNeed) "
    integer                                 :: rc            ! Return code
    integer                                 :: i             ! Indices
    integer                                 :: nFields       ! Number of fields in state
    character(len=SHR_KIND_CL), pointer :: Fieldnames(:)     ! Output fieldnames
    !
    if ( nFieldsNeeded > 0 )then
       if ( size(FieldsNeeded) < nFieldsNeeded ) &
          call shr_sys_abort( subName//" : Error size of input fieldsnotneed too small" )
       call eshr_estate_getFieldNames( eState, FieldNames, nFields )
       do i = 1, nFields
          if ( .not. any(trim(FieldNames(i)) == FieldsNeeded) )then
             call ESMF_StateSetNeeded( EState, itemname=FieldNames(i), &
                                       neededflag=ESMF_NOTNEEDED, rc=rc )
             call eshr_rc_check( rc, subname//" : Error in setting needed flag on state" )
          end if
       end do
       if ( nFields /= 0 ) deallocate( FieldNames )
    end if
END SUBROUTINE eshr_estate_markRestNotNeed
!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_getSharedFNames -- Get the field names shared on states.
!   
! !DESCRIPTION:
!   
!  Get the list of field names shared between the two input states.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_getSharedFNames( eStateIn, eStateOut, FieldNames, nFields, &
                                        onlyNeeded )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_State),        intent(IN)  :: EStateIn      ! Input state to copy
    type(ESMF_State),        intent(IN)  :: EStateOut     ! Output state to copy to
    character(len=SHR_KIND_CL), pointer  :: FieldNames(:) ! Shared Fieldnames
    integer,                 intent(OUT) :: nFields       ! No. of shared fields
    logical, optional,       intent(IN)  :: onlyNeeded    ! Only get needed fields
!EOP
    !----- local ----
    character(len=SHR_KIND_CL), pointer  :: FieldNamesIn(:)   ! Input Fieldnames
    character(len=SHR_KIND_CL), pointer  :: FieldNamesTmp(:)  ! Temporary fieldnames
    integer                              :: nFieldsIn         ! No. input fields
    character(len=SHR_KIND_CL), pointer  :: FieldNamesOut(:)  ! Output Fieldnames
    integer                              :: nFieldsOut        ! No. output fields
    integer                              :: i, j              ! Indices
    logical,                    pointer  :: taken(:)          ! If outname already taken
    logical                              :: neededOnly        ! If only extract needed
    !
    if ( present(onlyNeeded) )then
       neededOnly = onlyNeeded
    else
       neededOnly = .false.
    end if

    call eshr_estate_getFieldNames( eStateIn,   FieldNamesIn,  nFieldsIn,  &
                                    onlyNeeded=neededOnly  )
    call eshr_estate_getFieldNames( eStateOut,  FieldNamesOut, nFieldsOut, &
                                    onlyNeeded=neededOnly  )
    if ( (nFieldsIn == 0) .or. (nFieldsOut == 0) )then
       nFields = 0
       if ( nFieldsIn  > 0 ) deallocate( FieldNamesIn )
       if ( nFieldsOut > 0 ) deallocate( FieldNamesOut )
       return
    end if
    allocate( FieldNamesTmp(min(nFieldsIn,nFieldsOut)) )
    allocate( taken(nFieldsOut) )
    taken(:) = .false.
    nFields = 0
    do i = 1, nFieldsIn
       do j = 1, nFieldsOut
          if ( taken(j) )then
             cycle
          else if ( trim(FieldNamesOut(j)) == trim(FieldNamesIn(i)) )then
             nFields                = nFields + 1
             FieldNamesTmp(nFields) = trim(FieldNamesOut(j))
             taken(j)               = .true.
          end if
       end do
    end do
    allocate( FieldNames(nFields) )
    FieldNames(:) = FieldNamesTmp(1:nFields)

    deallocate( FieldNamesTmp )
    deallocate( taken         )
    deallocate( FieldNamesOut )
    deallocate( FieldNamesIn  )
END SUBROUTINE eshr_estate_getSharedFNames

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_Egrids2DAreSame -- compare two grids to see if identical
!   
! !DESCRIPTION:
!   
!  Compare to see if the two input grids are identical or not.
!      
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION eshr_estate_Egrids2DAreSame( grid1, grid2, Exact, Grid1Mask, &
                                              Debug, SameDecomp, NotData )
    implicit none

! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Grid),   intent(IN) :: grid1        ! First grid
    type(ESMF_Grid),   intent(IN) :: grid2        ! Grid to compare to
    logical,           intent(IN) :: Exact        ! If grids should be bit-for-bit same
    integer, optional, intent(IN) :: Grid1Mask(:) ! Mask on grid1 of grid2 data
    logical, optional, intent(IN) :: Debug        ! If should print info on grids
    logical, optional, intent(IN) :: SameDecomp   ! If decomposition needs to be identical
    logical, optional, intent(IN) :: NotData      ! Do NOT check exact latitude/longitudes
!EOP
    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_Egrids2DAreSame) "
    character(len=*), parameter :: gridList = "lat:lon"  ! List of grid fields
    integer                     :: i                ! Indices
    integer                     :: rc               ! Return code
    integer                     :: dc               ! Dimensions count
    type(ESMF_GridStorage)      :: GridStore1       ! Type of grid storage
    type(ESMF_GridStorage)      :: GridStore2       ! Type of grid storage
    real(SHR_KIND_R8)           :: MinExtents_1(3)  ! Minimum extents of 1st grid
    real(SHR_KIND_R8)           :: MinExtents_2(3)  ! Minimum extents of 2nd grid
    real(SHR_KIND_R8)           :: MaxExtents_1(3)  ! Minimum extents of 1st grid
    real(SHR_KIND_R8)           :: MaxExtents_2(3)  ! Minimum extents of 2nd grid
    real(SHR_KIND_R8),  pointer :: Lat_1(:)         ! 1st grid latitude
    real(SHR_KIND_R8),  pointer :: Lon_1(:)         ! 1st grid longitude
    real(SHR_KIND_R8),  pointer :: Lat_2(:)         ! 2nd grid latitude
    real(SHR_KIND_R8),  pointer :: Lon_2(:)         ! 2nd grid longitude
    real(SHR_KIND_R8),  pointer :: edata(:)         ! Pointer to state data
    integer                     :: dimCounts_1(3)   ! Dimension sizes, grid 1
    integer                     :: dimCounts_2(3)   ! Dimension sizes, grid 2
    logical                     :: ArbDecomp1       ! If grid1 Arbitrary decomposition
    logical                     :: ArbDecomp2       ! If grid2 Arbitrary decomposition
    integer                     :: nlon1, nlat1     ! grid 1 dimensions for long/latitude
    integer                     :: nlon2, nlat2     ! grid 2 dimensions for long/latitude
    logical,            pointer :: Mask1(:)         ! Mask of grid2 on grid 1 decomp.
    logical                     :: PrintInfo        ! If should print info on grids
    logical                     :: DecompIdentical  ! If decomposition has to be identical
    logical                     :: Same             ! If two grids are same
    type(ESMF_Bundle)           :: EBundle1         ! In bundle of fields to redist
    type(ESMF_Bundle)           :: EBundle2         ! Out bundle of fields to redist
    type(ESMF_RouteHandle)      :: ERoute           ! Comm. handle to do the redist
    type(ESMF_State)            :: grid1State       ! Grid values on grid1 state
    type(ESMF_State)            :: grid2State       ! Grid values on grid2 state

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

   if ( .not. present(Debug) )then
      printInfo = .false.
   else
      printInfo = Debug
   end if
   if ( .not. present(SameDecomp) )then
      decompIdentical = .false.
   else
      decompIdentical = SameDecomp
   end if
   eshr_estate_Egrids2DAreSame = .true.
   !
   !------  First make sure coordinates line up --------------------------------
   ! 
   call ESMF_GridGet( Grid1, horzrelloc=ESMF_CELL_CENTER,       &
                      minGlobalCoordPerDim=MinExtents_1(:),     &
                      maxGlobalCoordPerDim=MaxExtents_1(:),     &
                      globalCellCountPerDim=dimCounts_1(:),     &
                      gridStorage=gridStore1,                   &
                      dimCount=dc, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting grid extents from grid1" )
   if ( dc < 2 .or. dc > 3 )then
      call shr_sys_abort( subName//" : Error no. of dimensions in from grid1" )
   end if
   call ESMF_GridGet( Grid2, horzrelloc=ESMF_CELL_CENTER,       &
                      minGlobalCoordPerDim=MinExtents_2(:),     &
                      maxGlobalCoordPerDim=MaxExtents_2(:),     &
                      globalCellCountPerDim=dimCounts_2(:),     &
                      gridStorage=gridStore2,                   &
                      dimCount=dc, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting grid extents from grid2" )
   if ( dc < 2 .or. dc > 3 )then
      call shr_sys_abort( subName//" : Error no. of dimensions in grid2" )
   end if
   !----- Check if grid is arbitrary decomposition or block decomp ----
   if (      gridStore1 == ESMF_GRID_STORAGE_ARBITRARY )then
      ArbDecomp1 = .true.
   else
      ArbDecomp1 = .false.
   end if
   if (      gridStore2 == ESMF_GRID_STORAGE_ARBITRARY )then
      ArbDecomp2 = .true.
   else
      ArbDecomp2 = .false.
   end if
   
   if ( PrintInfo ) write(6,*) subName, " First check that grid extents are the same"
   if ( Exact )then
      if ( any(MinExtents_1(:2) /= MinExtents_2(:2)) ) &
                               eshr_estate_Egrids2DAreSame = .false.
      if ( any(MaxExtents_1(:2) /= MaxExtents_2(:2)) ) &
                               eshr_estate_Egrids2DAreSame = .false.
   else
      if ( any(abs(MinExtents_1(:2) - MinExtents_2(:2))> eps) ) &
               eshr_estate_Egrids2DAreSame = .false.
      if ( any(abs(MaxExtents_1(:2) - MaxExtents_2(:2))> eps) ) &
               eshr_estate_Egrids2DAreSame = .false.
   end if
   if ( any(dimCounts_1(:2)  /= dimCounts_2(:2))  ) eshr_estate_Egrids2DAreSame = .false.
   !
   !------ If want same decomposition and and arbitrary decomposition different then FALSE
   !
   if ( (ArbDecomp1 .neqv. ArbDecomp2) .and. DecompIdentical ) &
                                       eshr_estate_Egrids2DAreSame = .false.

   if ( present(NotData) )then
      if ( NotData ) return
   end if
   !
   !------  Now get coordinates from each grid ---------------------------------
   ! 
   if ( eshr_estate_Egrids2DAreSame )then
      call eshr_estate_egridGetLatLon( Grid1, Lat_1, Lon_1 )
      call eshr_estate_egridGetLatLon( Grid2, Lat_2, Lon_2 )
      nlon1 = size(Lon_1)
      nlat1 = size(Lat_1)
      nlon2 = size(Lon_2)
      nlat2 = size(Lat_2)
      if ( PrintInfo ) write(6,*) subName, " Now compare coordinates"
      !------ First compute mask if it wasn't sent in
      allocate( mask1(nlon1) )
      if ( present(Grid1Mask) )then
         mask1(:) = (Grid1Mask(:) == 1)
      else
         mask1(:) = .true.
      end if
      !
      !------ If same decomposition required -- check if same -----
      !
      if ( DecompIdentical )then
          eshr_estate_Egrids2DAreSame = eshr_estate_latlonAreSame( Lat_1, nlat1,     &
                             Lon_1, nlon1, Lat_2, nlat2, Lon_2, nlon2, Exact, Mask1, &
                             PrintInfo )
      else
          !
          !------ If same decomposition NOT required check if is the same anyway -----
          !
          same = eshr_estate_latlonAreSame( Lat_1, nlat1, Lon_1, nlon1, Lat_2, &
                                            nlat2, Lon_2, nlon2, Exact, Mask1, & 
                                            PrintInfo )
          if ( .not. same ) then
             !
             !------ Decomposition is different so redistribute to same grid -----
             !
             if ( PrintInfo ) write(6,*) subName, " Decomps are different redist to same"
             !--- First create states to redist between -----
             grid1state = ESMF_StateCreate( "grid1", rc=rc )
             grid2state = ESMF_StateCreate( "grid2", rc=rc )
             call eshr_estate_init2DFromList( Grid1, FieldList=gridList,  &
                                              eState=grid1state )
             call eshr_estate_init2DFromList( Grid2, FieldList=gridList,  &
                                              eState=grid2state )
             !--- Put latitude and longitude on 2nd state ----
             call ESMF_StateGetDataPointer( state=grid2State,                &
                                            itemName="lat",                  &
                                            dataPointer=edata,               &
                                            copyflag=ESMF_DATA_REF, rc=rc )
             call eshr_rc_check( rc, subName//': Error in getting lat F90 ptr to state' )
             edata(:) = Lat_2(:)
             nullify( edata )
             call ESMF_StateGetDataPointer( state=grid2State,                &
                                            itemName="lon",                  &
                                            dataPointer=edata,               &
                                            copyflag=ESMF_DATA_REF, rc=rc )
             call eshr_rc_check( rc, subName//': Error in getting lon F90 ptr to state' )
             edata(:) = Lon_2(:)
             nullify( edata )
             ! Now do redist from state-2 to state-1 so can compare directly
             ! First create states to redist between
             if ( PrintInfo ) write(6,*) subName, " Do the redistInit"
             call eshr_estate_reDistInit( grid2State, grid1State, EBundle1, EBundle2, &
                                          eRoute, noGridCheck=.true. )
             if ( PrintInfo ) write(6,*) subName, " Do the redist"
             call eshr_estate_reDist( EBundle1, EBundle2, eRoute )
             !----- And compare now ----
             nullify( Lat_2 )
             nullify( Lon_2 )
             if ( PrintInfo ) write(6,*) subName, " Compare again"
             call ESMF_StateGetDataPointer( state=grid1State,                &
                                            itemName="lat",                  &
                                            dataPointer=Lat_2,               &
                                            copyflag=ESMF_DATA_REF, rc=rc )
             call eshr_rc_check( rc, subName//': Error in getting lat F90 ptr to state' )
             call ESMF_StateGetDataPointer( state=grid1State,                &
                                            itemName="lon",                  &
                                            dataPointer=Lon_2,               &
                                            copyflag=ESMF_DATA_REF, rc=rc )
             call eshr_rc_check( rc, subName//': Error in getting lon F90 ptr to state' )
             nlon2 = size(Lon_2)
             nlat2 = size(Lat_2)
             eshr_estate_Egrids2DAreSame = eshr_estate_latlonAreSame( Lat_1, nlat1, &
                         Lon_1, nlon1, Lat_2, nlat2, Lon_2, nlon2, Exact, Mask1,    &
                         PrintInfo )
         end if
      end if
      deallocate( mask1 )
      nullify( Lon_1 )
      nullify( Lat_1 )
      nullify( Lon_2 )
      nullify( Lat_2 )
   end if
   !
   return

END FUNCTION eshr_estate_Egrids2DAreSame

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_egridGetLatLon -- Get the latitudes and longitude from the grid
!   
! !DESCRIPTION:
!   
!  Get the local latitude and longitudes from the grid.
!      
! !INTERFACE: ------------------------------------------------------------------
SUBROUTINE eshr_estate_egridGetLatLon( Grid, Lat, Lon )
    implicit none
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Grid),   intent(IN) :: grid       ! Input grid
    real(SHR_KIND_R8),  pointer      :: Lat(:)  ! Grid latitude
    real(SHR_KIND_R8),  pointer      :: Lon(:)  ! Grid longitude
!EOP

    !----- local -----
    character(len=*), parameter :: subName = " (eshr_estate_egridGetLatLon) "
    type(ESMF_Array)                 :: CoordsG(3)      ! Grid coordinates
    integer                          :: rc               ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

#ifdef ESMF_3
   call ESMF_GridGetCoord( Grid, 1, horzRelloc=ESMF_CELL_CENTER, &
                           centerCoord=Lon, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting coords from grid" )
   call ESMF_GridGetCoord( Grid, 2, horzRelloc=ESMF_CELL_CENTER, &
                           centerCoord=Lat, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting coords from grid" )
#else
   call ESMF_GridGetCoord( Grid, horzRelloc=ESMF_CELL_CENTER, &
                           centerCoord=CoordsG, rc=rc )
   call eshr_rc_check( rc, subName//" : Error getting coords from grid 1" )
   call ESMF_ArrayGetData( CoordsG(1), fptr=Lon, docopy=ESMF_DATA_REF, rc=rc)
   call eshr_rc_check( rc, subName//" : Error getting ptr to coords from grid" )
   call ESMF_ArrayGetData( CoordsG(2), fptr=Lat, docopy=ESMF_DATA_REF, rc=rc)
   call eshr_rc_check( rc, subName//" : Error getting ptr to coords from grid" )
#endif
END SUBROUTINE eshr_estate_egridGetLatLon

!===============================================================================

!===============================================================================
! !IROUTINE: eshr_estate_latlonAreSame -- Check to see if the lat/longs are the same
!   
! !DESCRIPTION:
!   
!  Check to see if the input latitudes and longitudes are the same
!      
! !INTERFACE: ------------------------------------------------------------------
logical FUNCTION eshr_estate_latlonAreSame( Lat_1, nlat1, Lon_1, nlon1, Lat_2, &
                                            nlat2, Lon_2, nlon2, Exact, Mask1, &
                                            PrintInfo )
   implicit none
! !INPUT/OUTPUT PARAMETERS:
   integer,           intent(in) :: nlat1          ! Size of latitude 1 array
   integer,           intent(in) :: nlon1          ! Size of longitude 1 array
   integer,           intent(in) :: nlat2          ! Size of latitude 2 array
   integer,           intent(in) :: nlon2          ! Size of longitude 2 array
   real(SHR_KIND_R8), intent(in) :: Lat_1(nlat1)   ! Latitude 1 array
   real(SHR_KIND_R8), intent(in) :: Lon_1(nlon1)   ! Longitude 1 array
   real(SHR_KIND_R8), intent(in) :: Lat_2(nlat2)   ! Latitude 2 array
   real(SHR_KIND_R8), intent(in) :: Lon_2(nlon2)   ! Longitude 2 array
   logical,           intent(in) :: Exact          ! If comparision should be exact
   logical,           intent(in) :: Mask1(:)       ! Mask of grid2 points on grid1
   logical,           intent(in) :: PrintInfo      ! If should print extra debugging info
   !----- local -----
   character(len=*), parameter :: subName = " (eshr_estate_latlonAreSame) "
   integer                       :: i               ! Indices
!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------


   eshr_estate_latlonAreSame = .true.
   !
   !-----  Make sure dimensions are same ------
   !
   if ( nlon1 /= nlon2 ) eshr_estate_latlonAreSame = .false.
   if ( nlat1 /= nlat2 ) eshr_estate_latlonAreSame = .false.
   if ( eshr_estate_latlonAreSame )then
      !-------------------------------------------------------------------------
      ! If require exact comparision check that grid points are identical
      !-------------------------------------------------------------------------
      if ( Exact )then
         do i = 1, nlon1
            if ( Mask1(i) )then
                if ( Lon_1(i) /= Lon_2(i) ) eshr_estate_latlonAreSame = .false.
            end if
            if ( PrintInfo .and. .not. eshr_estate_latlonAreSame ) then
                write(6,*) "Lon_1 = ", Lon_1(i), ' i = ', i
                write(6,*) "Lon_2 = ", Lon_2(i), ' i = ', i
            end if
            if ( .not. eshr_estate_latlonAreSame ) exit
         end do
         if ( PrintInfo .and. eshr_estate_latlonAreSame )  &
                                         write(6,*) subName, " Longitude identical: "
         do i = 1, nlat1
            if ( Mask1(i) )then
               if ( Lat_1(i) /= Lat_2(i) ) eshr_estate_latlonAreSame = .false.
               if ( PrintInfo .and. .not. eshr_estate_latlonAreSame ) then
                  write(6,*) "Lat_1 = ", Lat_1(i), ' i = ', i
                  write(6,*) "Lat_2 = ", Lat_2(i), ' i = ', i
               end if
            end if
            if ( .not. eshr_estate_latlonAreSame ) exit
         end do
         if ( PrintInfo .and. eshr_estate_latlonAreSame ) &
                                         write(6,*) subName, " Latitude identical: "
      !-------------------------------------------------------------------------
      ! If not check that grid points are identical to within acceptable value
      !-------------------------------------------------------------------------
      else
         do i = 1, nlon1
           if ( Mask1(i) )then
              if ( abs(Lon_1(i) - Lon_2(i)) > eps ) eshr_estate_latlonAreSame = .false.
              if ( PrintInfo .and. .not. eshr_estate_latlonAreSame )then
                 write(6,*) "Lon_1 = ", Lon_1(i), ' i = ', i
                 write(6,*) "Lon_2 = ", Lon_2(i), ' i = ', i
              end if
           end if
           if ( .not. eshr_estate_latlonAreSame ) exit
         end do
         if ( PrintInfo .and. eshr_estate_latlonAreSame ) &
                                      write(6,*) subName, " Longitude identical: "
         do i = 1, nlat1
           if ( Mask1(i) )then
              if ( abs(Lat_1(i) - Lat_2(i)) > eps) eshr_estate_latlonAreSame = .false.
              if ( PrintInfo .and. .not. eshr_estate_latlonAreSame )then
                 write(6,*) "Lat_1 = ", Lat_1(i), ' i = ', i
                 write(6,*) "Lat_2 = ", Lat_2(i), ' i = ', i
              end if
           end if
           if ( .not. eshr_estate_latlonAreSame ) exit
         end do
         if ( PrintInfo .and. eshr_estate_latlonAreSame ) &
                                      write(6,*) subName, " Latitude identical: "
      end if
   end if
END FUNCTION eshr_estate_latlonAreSame

!===============================================================================

#endif

end module eshr_estate_mod
