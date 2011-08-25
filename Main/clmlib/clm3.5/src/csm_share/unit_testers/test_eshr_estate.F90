module test_eshr_inputinfo_mod
  use eshr_inputinfo_mod, only: eshr_inputinfo_eState2Info, eshr_inputinfo_Info2EState
  use shr_inputinfo_mod,  only: shr_inputinfo_initType
  use eshr_rc_mod,        only: eshr_rc_check
  use eshr_estate_mod,    only: eshr_estate_printAttributes
  use ESMF_Mod
  implicit none
contains

logical function check_ccsminitConversion( CCSMInit )
  use shr_inputinfo_mod, only: shr_inputinfo_initIsSame
  use ESMF_Mod
  implicit none
  type(shr_inputinfo_initType), intent(IN) :: CCSMInit

  type(ESMF_State) :: EState
  type(ESMF_State) :: CCSMInitEState
  type(shr_inputinfo_initType) :: CCSMInit2
  integer :: rc

  eState = ESMF_StateCreate( "Temp state", rc=rc )
  call eshr_rc_check( rc, "Error in creation of ESMF State" )
  call eshr_inputinfo_info2EState( CCSMInit, Estate, CCSMInitEState, print=.true. )
  call eshr_estate_printAttributes( Estate )
  call ESMF_StatePrint( EState )
  call eshr_inputinfo_EState2Info( Estate, CCSMInit2, print=.true. )
  if ( shr_inputinfo_initisSame( CCSMInit, CCSMInit2 ) ) then
     check_ccsminitConversion = .true.
  else
     check_ccsminitConversion = .false.
  end if
  call ESMF_StateDestroy( eState )
end function check_ccsminitConversion

end module test_eshr_inputinfo_mod

module test_eshr_timemgr_mod
  use eshr_timemgr_mod, only: eshr_timemgr_eState2Info, eshr_timemgr_Info2EState
  !                            eshr_timemgr_clockInfoInit
  use eshr_rc_mod,      only: eshr_rc_check
  use eshr_estate_mod,  only: eshr_estate_printAttributes
  use ESMF_Mod
  implicit none
contains

logical function check_clockInfoConversion( clockInfo  )
  use eshr_timemgr_mod, only: eshr_timemgr_clockInfoIsSame, &
                              eshr_timemgr_clockInfoType
  use ESMF_Mod
  implicit none
  type(eshr_timemgr_clockInfoType), intent(IN) :: clockInfo

  type(ESMF_State) :: EState
  type(ESMF_State) :: ClockInfoEState
  type(eshr_timemgr_clockInfoType) :: clockInfo2
  integer :: rc

  eState = ESMF_StateCreate( "Temp state", rc=rc )
  call eshr_rc_check( rc, "Error in creation of ESMF State" )
  call eshr_timemgr_info2EState( ClockInfo, Estate, ClockInfoEState, print=.true. )
  call eshr_estate_printAttributes( Estate )
  call ESMF_StatePrint( EState )
  call eshr_timemgr_EState2Info( Estate, ClockInfo2, print=.true. )
  if ( eshr_timemgr_clockInfoisSame( ClockInfo, ClockInfo2 ) ) then
     check_clockInfoConversion = .true.
  else
     check_clockInfoConversion = .false.
  end if
  call ESMF_StateDestroy( eState )
end function check_clockInfoConversion

end module test_eshr_timemgr_mod

module test_eshr_estate_mod
  use shr_const_mod,     only: SHR_CONST_SPVAL
  use eshr_estate_mod,   only: eshr_estate_getFilename,         &
                               eshr_estate_putFilename,         &
                               eshr_estate_init2DFromList,      &
                               eshr_estate_getStats,            &
                               eshr_estate_printAttributes,     &
                               eshr_estate_printFieldAttr,      &
                               eshr_estate_getfieldnames,       &
                               eshr_estate_getfieldsNotNeed,    &
                               eshr_estate_getSharedFNames,     &
                               eshr_estate_getDataPointer,      &
                               eshr_estate_setDataPointer,      &
                               eshr_estate_markRestNotNeed,     &
                               eshr_estate_copyFields,          &
                               eshr_estate_redistInit,          &
                               eshr_estate_redist,              &
                               eshr_estate_fieldsAreSame,       &
                               eshr_estate_matchValue,          &
                               eshr_estate_mergeInit,           &
                               eshr_estate_merge,               &
                               eshr_estate_checkGridsMatch,     &
                               eshr_estate_Egrids2DAreSame,     &
                               eshr_estate_destroy,             &
                               eshr_estate_initDomain
  use eshr_rc_mod,       only: eshr_rc_check
  use shr_sys_mod,       only: shr_sys_abort
  use shr_kind_mod,      only: SHR_KIND_CS, SHR_KIND_CL, SHR_KIND_R8, &
                               SHR_KIND_CX
  use ESMF_Mod
  implicit none

  integer, parameter :: RANDOM_FILL   = 1
  integer, parameter :: CONSTANT_FILL = 2
  integer, parameter :: SLOPE_FILL    = 3
  integer, parameter :: STEP_FILL     = 4

  real(SHR_KIND_R8), parameter :: icefrac = 0.2_SHR_KIND_R8
  real(SHR_KIND_R8), parameter :: ocnfrac = 0.5_SHR_KIND_R8
  real(SHR_KIND_R8), parameter :: lndfrac = 0.3_SHR_KIND_R8
contains

logical function check_FilenameConversion( FileType, Filename )
  use ESMF_Mod
  implicit none
  character(len=*), intent(in) :: FileType
  character(len=*), intent(in) :: Filename

  type(ESMF_State) :: EState
  character(len=256) :: Filename2
  integer :: rc

  eState = ESMF_StateCreate( "Temp state", rc=rc )
  call eshr_rc_check( rc, "Error in creation of ESMF State" )
  if ( len_trim(FileType) == 0 )then
     check_FilenameConversion = .false.
  else
     call eshr_estate_putFilename( FileType, Filename, eState )
     call eshr_estate_getFilename( estate, FileType, Filename2 )
     if ( trim(FileName) == trim(FileName2) )then
        check_FilenameConversion = .true.
     else
        check_FilenameConversion = .false.
     end if
  end if
  call ESMF_StateDestroy( eState )

end function check_FilenameConversion

logical function check_EstateInit( grid, FieldList, FieldLNames, FieldUnits, arbitrary )
    use shr_string_mod,          only: shr_string_listGetName, &
                                       shr_string_listGetNum
    type(ESMF_Grid)               :: Grid        ! Grid fields are on
    character(len=*),  intent(IN) :: FieldList   ! Colen delimited list of fields
    character(len=*),  intent(IN) :: FieldLNames ! : delimited list of long_names
    character(len=*),  intent(IN) :: FieldUnits  ! : delimited list of units
    logical,           intent(IN) :: arbitrary   ! If distributed on an arbitrary decomposition

    type(ESMF_State)            :: eState          ! Output ESMF State
    integer                     :: rc              ! ESMF return code
    integer                     :: i               ! Index
    integer                     :: nfields         ! Number of fields
    integer                     :: itemCount       ! Items in state
    character(len=256)          :: FieldName       ! Field name looking at
    character(len=256), pointer :: itemNameList(:) ! List of field names
    real(SHR_KIND_R8),  pointer :: Data1D(:)       ! Pointer to 1D data
    real(SHR_KIND_R8),  pointer :: Data2D(:,:)     ! Pointer to 2D data
    real(SHR_KIND_R8)           :: missval         ! Missing value
    type(ESMF_GridStorage)      :: gridStore       ! Grid storage
!
    ! Make sure grid is correctly labeled according to arbitrary flag sent in
    call ESMF_GridGet( grid, gridStorage=gridStore, rc=rc )
    call eshr_rc_check( rc, "Error in get of grid storage" )
    if ( arbitrary .and. (gridStore /= ESMF_GRID_STORAGE_ARBITRARY) )then
       call shr_sys_abort( "Input arbitrary flag inconsistent with grid sent in" )
    end if
    if ( (.not. arbitrary) .and. (gridStore == ESMF_GRID_STORAGE_ARBITRARY) )then
       call shr_sys_abort( "Input arbitrary flag inconsistent with grid sent in" )
    end if
    missval = SHR_CONST_SPVAL
    check_EStateInit  = .false.
    eState = ESMF_StateCreate( "Temp state", rc=rc )
    call eshr_rc_check( rc, "Error in creation of ESMF State" )
    if (      len_trim( FieldLNames ) > 0 .and. len_trim( FieldUnits) > 0 )then
       call eshr_estate_init2DFromList(  Grid, FieldList, FieldLNames=FieldLNames, &
                                         FieldUnits=FieldUnits, eState=eState )
    else if ( len_trim( FieldLNames ) > 0 )then
       call eshr_estate_init2DFromList(  Grid, FieldList, FieldLNames=FieldLNames,      &
                                         eState=eState )
    else if ( len_trim( FieldUnits ) > 0 )then
       call eshr_estate_init2DFromList(  Grid, FieldList, FieldUnits=FieldUnits,        &
                                         eState=eState )
    else
       call eshr_estate_init2DFromList(  Grid, FieldList, eState=eState )
    end if
    if ( .not. eshr_estate_matchValue( eState, missval, PrintNames=.true. ) )then
       call shr_sys_abort( "NOT finding missing value on freshly initialized state" )
    end if
    ! ------ Check the state out ------
    call ESMF_StateValidate( EState, rc=rc )
    call eshr_rc_check( rc, "Error, on state validate" )
    if ( rc /= ESMF_SUCCESS ) return
    nFields = shr_string_listGetNum( FieldList )
    allocate( itemNameList(nFields) )
    call ESMF_StateGet( EState, itemCount=itemCount, itemNameList=itemNameList, rc=rc )
    call eshr_rc_check( rc, "Error, getting items and names from state" )
    if ( rc /= ESMF_SUCCESS ) return
    if ( itemCount /= nFields )then
       write(6,*) 'Error, number of items on state is different from expected'
       return
    end if
    do i = 1, nfields
       call shr_string_listGetName( FieldList, i, FieldName )
       if ( trim(itemNameList(i)) /= trim(FieldName) )then
          write(6,*) 'Error, field name on state different than expected'
          write(6,*) 'FieldName = ', trim(FieldName)
          write(6,*) 'itemName  = ', trim(itemNameList(i))
          return
       end if
       if ( arbitrary )then
          call ESMF_StateGetDataPointer( eState, FieldName, Data1D, &
                                         copyflag=ESMF_DATA_REF, rc=rc)
       else
          call ESMF_StateGetDataPointer( eState, FieldName, Data2D, &
                                         copyflag=ESMF_DATA_REF, rc=rc)
       end if
       call eshr_rc_check( rc, 'ERROR: Getting data pointer from state' )
       if ( rc /= ESMF_SUCCESS ) return
       if ( .not. ESMF_StateIsNeeded(estate, FieldName, rc) )then
          write(6,*) 'Error, field name on state NOT needed -- terminate with error'
          return
       end if
       if ( arbitrary )then
          if ( any(Data1D(:) /= missval) ) then
             write(6,*) 'Error, field name on state different than expected'
             return
          end if
          nullify( Data1D )
       else
          if ( any(Data2D(:,:) /= missval) ) then
             write(6,*) 'Error, field name on state different than expected'
             return
          end if
          nullify( Data2D )
       end if
    end do
    deallocate( itemNameList )
    ! ------ Print info on the state and contained fields ----
    call ESMF_StatePrint( eState, rc=rc )
    call eshr_rc_check( rc, "Error, printing state" )
    if ( rc /= ESMF_SUCCESS ) return
    call eshr_estate_PrintAttributes( eState )
    call eshr_estate_PrintFieldAttr( eState )

    ! ------- Destroy the state and then contained fields ----
    call ESMF_StateDestroy( eState, rc )
    call eshr_rc_check( rc, "Error, destroying state" )
    if ( rc /= ESMF_SUCCESS ) return

    check_EStateInit  = .true.
end function check_EstateInit

subroutine create_domains( agrid, igrid, lgrid, ogrid, lMask_a, domain_a, domain_i, &
                           domain_l, domain_o )
  implicit none
  type(ESMF_Grid),  intent(IN)  :: agrid
  type(ESMF_Grid),  intent(IN)  :: igrid
  type(ESMF_Grid),  intent(IN)  :: lgrid
  type(ESMF_Grid),  intent(IN)  :: ogrid
  integer,          intent(IN)  :: lMask_a(:)
  type(ESMF_State), intent(OUT) :: domain_a, domain_i, domain_l, domain_o
  real(SHR_KIND_R8),         pointer     :: area(:)
  real(SHR_KIND_R8),         pointer     :: rmask(:)
  real(SHR_KIND_R8),         pointer     :: maxfrac(:)
  real(SHR_KIND_R8),         pointer     :: data1D(:)
  real(SHR_KIND_R8),         pointer     :: lat(:)
  real(SHR_KIND_R8)                      :: x
  integer                                 :: ngrid, nl, i
  integer                                 :: rc
  type(ESMF_Array) :: CoordsG(2)
  ! Create domains
  ! Atm
  call ESMF_GridGetCoord( aGrid, horzRelloc=ESMF_CELL_CENTER, &
                          centerCoord=CoordsG, rc=rc )
  call eshr_rc_check( rc, " : Error getting coords from grid" )
  call ESMF_ArrayGetData( CoordsG(2), fptr=Lat, docopy=ESMF_DATA_REF, rc=rc)
  call eshr_rc_check( rc, " : Error getting ptr to coords from grid" )
  ngrid = size(Lat)
  nullify( Lat )
  allocate( area(ngrid), rmask(ngrid), maxfrac(ngrid) )
  area(:)    = 1.0_SHR_KIND_R8
  rmask(:)   = 1.0_SHR_KIND_R8
  maxfrac(:) = 1.0_SHR_KIND_R8
  call eshr_estate_initDomain( "atm", agrid, area, rmask, maxfrac, domain_a )
  deallocate( area, rmask, maxfrac )
  if ( eshr_estate_matchValue( domain_a, SHR_CONST_SPVAL) )then
     call shr_sys_abort( "Found missing value on atm domain" )
  end if
  ! Lnd
  call ESMF_GridGetCoord( lGrid, horzRelloc=ESMF_CELL_CENTER, &
                          centerCoord=CoordsG, rc=rc )
  call eshr_rc_check( rc, " : Error getting coords from grid" )
  call ESMF_ArrayGetData( CoordsG(2), fptr=Lat, docopy=ESMF_DATA_REF, rc=rc)
  call eshr_rc_check( rc, " : Error getting ptr to coords from grid" )
  ngrid = size(lat)
  nullify( lat )
  allocate( area(ngrid), rmask(ngrid), maxfrac(ngrid) )
  area(:) = 1.0_SHR_KIND_R8
  if ( eshr_estate_Egrids2DAreSame( agrid, lgrid, Exact=.true., Grid1Mask=lMask_a, &
                                    SameDecomp=.true. ) )then
     rmask(:) = lMask_a(:)
     where( lMask_a > 0 )
        maxfrac = lndfrac
     elsewhere
        maxfrac = 1.0_SHR_KIND_R8
     end where
  else
     maxfrac(:) = lndfrac
     rmask(:) = 1.0_SHR_KIND_R8
  end if
  call eshr_estate_initDomain( "lnd", lgrid, area, rmask, maxfrac, domain_l )
  if ( eshr_estate_matchValue( domain_l, SHR_CONST_SPVAL) )then
     call shr_sys_abort( "Found missing value on lnd domain" )
  end if
  deallocate( area, rmask, maxfrac )
  ! Ocn
  call ESMF_GridGetCoord( oGrid, horzRelloc=ESMF_CELL_CENTER, &
                          centerCoord=CoordsG, rc=rc )
  call eshr_rc_check( rc, " : Error getting coords from grid" )
  call ESMF_ArrayGetData( CoordsG(2), fptr=Lat, docopy=ESMF_DATA_REF, rc=rc)
  call eshr_rc_check( rc, " : Error getting ptr to coords from grid" )
  ngrid = size(lat)
  nullify( lat )
  allocate( area(ngrid), rmask(ngrid), maxfrac(ngrid) )
  area(:) = 1.0_SHR_KIND_R8
  ! Calc maxfrac and mask
  rmask(:) = lMask_a
  where( lMask_a > 0 )
    maxfrac = 1.0_SHR_KIND_R8 - lndfrac
  elsewhere
    maxfrac = 1.0_SHR_KIND_R8
  end where
  call eshr_estate_initDomain( "ocn", ogrid, area, rmask, maxfrac, domain_o )
  if ( eshr_estate_matchValue( domain_o, SHR_CONST_SPVAL) )then
     call shr_sys_abort( "Found missing value on ocn domain" )
  end if
  ! Ice
  call ESMF_GridGetCoord( iGrid, horzRelloc=ESMF_CELL_CENTER, &
                          centerCoord=CoordsG, rc=rc )
  call eshr_rc_check( rc, " : Error getting coords from grid" )
  call ESMF_ArrayGetData( CoordsG(2), fptr=Lat, docopy=ESMF_DATA_REF, rc=rc)
  call eshr_rc_check( rc, " : Error getting ptr to coords from grid" )
  ngrid = size(lat)
  nullify( lat )
  call eshr_estate_initDomain( "ice", igrid, area, rmask, maxfrac, domain_i )
  deallocate( area, rmask, maxfrac )
  if ( eshr_estate_matchValue( domain_i, SHR_CONST_SPVAL) )then
     call shr_sys_abort( "Found missing value on ice domain" )
  end if
end subroutine create_domains

logical function check_EstateImportCreate( a2x_fields, i2x_fields, l2x_fields,     &
                                           o2x_fields, agrid, igrid, lgrid, ogrid, &
                                           lMask_a )
  use shr_string_mod, only: shr_string_listGetNum
  use shr_kind_mod,   only: SHR_KIND_CL
  implicit none
  character(len=*), intent(IN) :: a2x_fields
  character(len=*), intent(IN) :: i2x_fields
  character(len=*), intent(IN) :: l2x_fields
  character(len=*), intent(IN) :: o2x_fields
  type(ESMF_Grid),  intent(IN) :: agrid
  type(ESMF_Grid),  intent(IN) :: igrid
  type(ESMF_Grid),  intent(INOUT) :: lgrid
  type(ESMF_Grid),  intent(IN) :: ogrid
  integer,          intent(IN) :: lMask_a(:)

  type(ESMF_State)             :: a2x_a
  type(ESMF_State)             :: i2x_i
  type(ESMF_State)             :: l2x_l
  type(ESMF_State)             :: o2x_o
  type(ESMF_State)             :: domain_a, domain_i, domain_l, domain_o

  type(ESMF_State)             :: a2x_l
  type(ESMF_State)             :: a2x_i
  type(ESMF_State)             :: a2x_o
  type(ESMF_Array)                        :: CoordsG(3)    ! Grid coordinates
  integer                                 :: rc
  integer                                 :: ntag
  integer                                 :: nlnd          ! Size of land arrays
  integer                                 :: i, j, n       ! Indices
  integer                                 :: nFields       ! Number of fields in state
  integer                                 :: rank          ! Rank of field
  character(len=SHR_KIND_CL), pointer     :: FieldNames(:) ! Names of fields  in state
  character(len=SHR_KIND_CL), pointer     :: fname         ! First fieldname
  type(ESMF_Field)                        :: Field         ! Field to query
  type(ESMF_Array)                        :: array         ! Pointer to array in field
  type(ESMF_DataType)                     :: dataType      ! Type of data in array
  type(ESMF_DataKind)                     :: kind          ! Kind of data in array
  real(ESMF_KIND_R8)                      :: x             ! Random number
  real(ESMF_KIND_R8)                      :: xsum, xmin, xmax, xsq ! Statistics
  real(ESMF_KIND_R8),         pointer     :: data2D(:,:)   ! Pointer 2D array in field
  real(ESMF_KIND_R8),         pointer     :: data1D(:)     ! Pointer 1D array in field
  real(ESMF_KIND_R8),         pointer     :: ptr1(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr2(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr3(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr4(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr5(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr6(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr7(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr8(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr9(:)       ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr10(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr11(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr12(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr13(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr14(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr15(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr16(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr17(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr18(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr19(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: ptr20(:)      ! FORTRAN ptr for land data
  real(ESMF_KIND_R8),         pointer     :: Lon(:)        ! FORTRAN ptr to longitudes

  check_EstateImportCreate = .false.
  a2x_a = ESMF_StateCreate( statename="atm import", statetype=ESMF_STATE_IMPORT, rc=rc )
  call eshr_rc_check( rc, "Error, creating state" )
  if ( rc /= ESMF_SUCCESS ) return
  i2x_i = ESMF_StateCreate( statename="ice import", statetype=ESMF_STATE_IMPORT, rc=rc )
  call eshr_rc_check( rc, "Error, creating state" )
  if ( rc /= ESMF_SUCCESS ) return
  l2x_l = ESMF_StateCreate( statename="lnd import", statetype=ESMF_STATE_IMPORT, rc=rc )
  call eshr_rc_check( rc, "Error, creating state" )
  if ( rc /= ESMF_SUCCESS ) return
  o2x_o = ESMF_StateCreate( statename="ocn import", statetype=ESMF_STATE_IMPORT, rc=rc )
  call eshr_rc_check( rc, "Error, creating state" )
  if ( rc /= ESMF_SUCCESS ) return
  x = 0.0_ESMF_KIND_R8
  call eshr_estate_init2DFromList(  aGrid, a2x_fields, initValue=x, eState=a2x_a )
  call eshr_estate_init2DFromList(  iGrid, i2x_fields, initValue=x, eState=i2x_i )
  call eshr_estate_init2DFromList(  lGrid, l2x_fields, initValue=x, eState=l2x_l, &
                                    UsePTR=.true. )
  call eshr_estate_init2DFromList(  oGrid, o2x_fields, initValue=x, eState=o2x_o )

  check_EstateImportCreate = .true.

  if ( .not. eshr_estate_Egrids2DAreSame( aGrid, lGrid, Exact=.false., &
             Grid1Mask=lMask_a ) )then
     write(6,*) "Grids do NOT match"
     check_EstateImportCreate = .false.
  end if
  if ( .not. eshr_estate_Egrids2DAreSame( aGrid, iGrid, Exact=.true., Grid1Mask=lMask_a, &
             SameDecomp=.true. ) )then
     write(6,*) "Grids do NOT match"
     check_EstateImportCreate = .false.
  end if
  if ( .not. eshr_estate_Egrids2DAreSame( aGrid, oGrid, Exact=.true., Grid1Mask=lMask_a, &
             SameDecomp=.true. ) )then
     write(6,*) "Grids do NOT match"
     check_EstateImportCreate = .false.
  end if
  if ( check_EstateImportCreate )then
     call create_domains( agrid, igrid, lgrid, ogrid, lMask_a, domain_a, domain_i, &
                          domain_l, domain_o )

     call ESMF_StateAddState( i2x_i, domain_i, rc=rc )
     call eshr_rc_check( rc, " : Error adding domain_i to i2x_i" )
     call ESMF_StateAddState( a2x_a, domain_a, rc=rc )
     call eshr_rc_check( rc, " : Error adding domain_a to a2x_a" )
     call ESMF_StateAddState( l2x_l, domain_l, rc=rc )
     call eshr_rc_check( rc, " : Error adding domain_l to l2x_l" )
     call ESMF_StateAddState( o2x_o, domain_o, rc=rc )
     call eshr_rc_check( rc, " : Error adding domain_o to o2x_o" )

     if ( .not. eshr_estate_checkGridsMatch( i2x_i, a2x_a, o2x_o, l2x_l, &
            gridsExact=.true.) )then
        write(6,*) "Grids do NOT match"
        check_EstateImportCreate = .false.
     end if
  else
     domain_a = ESMF_StateCreate( statename="atm domain", rc=rc )
     call eshr_rc_check( rc, "Error, creating state" )
     domain_l = ESMF_StateCreate( statename="lnd domain", rc=rc )
     call eshr_rc_check( rc, "Error, creating state" )
     domain_i = ESMF_StateCreate( statename="ice domain", rc=rc )
     call eshr_rc_check( rc, "Error, creating state" )
     domain_o = ESMF_StateCreate( statename="ocn domain", rc=rc )
     call eshr_rc_check( rc, "Error, creating state" )
     call ESMF_StateAddState( i2x_i, domain_i, rc=rc )
     call eshr_rc_check( rc, " : Error adding domain_i to i2x_i" )
     call ESMF_StateAddState( a2x_a, domain_a, rc=rc )
     call eshr_rc_check( rc, " : Error adding domain_a to a2x_a" )
     call ESMF_StateAddState( l2x_l, domain_l, rc=rc )
     call eshr_rc_check( rc, " : Error adding domain_l to l2x_l" )
     call ESMF_StateAddState( o2x_o, domain_o, rc=rc )
     call eshr_rc_check( rc, " : Error adding domain_o to o2x_o" )
  end if
  !
  ! Set land data pointers
  !
  call ESMF_GridGetCoord( lGrid, horzRelloc=ESMF_CELL_CENTER, &
                          centerCoord=CoordsG, rc=rc )
  call eshr_rc_check( rc, " : Error getting coords from land grid" )
  call ESMF_ArrayGetData( CoordsG(1), fptr=Lon, docopy=ESMF_DATA_REF, rc=rc )
  call eshr_rc_check( rc, " : Error getting long coords from land grid" )
  nlnd = size(Lon)
  nullify( Lon )
  write(6,*) "nlnd=", nlnd
  allocate( ptr1(nlnd),  ptr2(nlnd),  ptr3(nlnd),  ptr4(nlnd),  ptr5(nlnd)  )
  allocate( ptr6(nlnd),  ptr7(nlnd),  ptr8(nlnd),  ptr9(nlnd),  ptr10(nlnd) )
  allocate( ptr11(nlnd), ptr12(nlnd), ptr13(nlnd), ptr14(nlnd), ptr15(nlnd)  )
  allocate( ptr16(nlnd), ptr17(nlnd), ptr18(nlnd), ptr19(nlnd), ptr20(nlnd) )
  if ( shr_string_listGetNum( l2x_fields) /= 20 ) &
     call shr_sys_abort( "Not the right number of elements in l2x_fields" )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  1, ptr1  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  2, ptr2  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  3, ptr3  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  4, ptr4  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  5, ptr5  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  6, ptr6  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  7, ptr7  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  8, ptr8  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields,  9, ptr9  )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 10, ptr10 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 11, ptr11 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 12, ptr12 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 13, ptr13 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 14, ptr14 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 15, ptr15 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 16, ptr16 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 17, ptr17 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 18, ptr18 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 19, ptr19 )
  call eshr_estate_setDataPointer( l2x_l, l2x_fields, 20, ptr20 )
  ptr1(:) = 1.0
  ptr2(:) = 2.0
  ptr3(:) = 3.0
  call eshr_estate_getDataPointer( l2x_l, l2x_fields, 1, data1D, expectSize=nlnd )
  if ( .not. all(data1D == ptr1) ) &
     call shr_sys_abort( "Pointers not handled correctly" )
  nullify( data1D )
  call eshr_estate_getDataPointer( l2x_l, l2x_fields, 2, data1D, expectSize=nlnd )
  if ( .not. all(data1D == ptr2) ) &
     call shr_sys_abort( "Pointers not handled correctly" )
  nullify( data1D )
  call eshr_estate_getDataPointer( l2x_l, l2x_fields, 3, data1D, expectSize=nlnd )
  if ( .not. all(data1D == ptr3) ) &
     call shr_sys_abort( "Pointers not handled correctly" )
  nullify( data1D )
  deallocate( ptr1,  ptr2,  ptr3,  ptr4,  ptr5  )
  deallocate( ptr6,  ptr7,  ptr8,  ptr9,  ptr10 )
  deallocate( ptr11, ptr12, ptr13, ptr14, ptr15  )
  deallocate( ptr16, ptr17, ptr18, ptr19, ptr20 )
  ! 
  ! Fill in data to test get stats
  !
  call eshr_estate_getFieldNames( a2x_a, FieldNames, nFields )
  xmin = 100.0
  xmax = 0.0
  xsum = 0.0
  xsq  = 0.0
  ntag = 0
  do j = 1, nFields
     if ( .not. check_EstateImportCreate ) exit
     call ESMF_StateGetField( state=a2x_a, field=Field, fieldname=FieldNames(j), &
                              rc=rc )
     call eshr_rc_check( rc, ' : Error getting field '//FieldNames(j)// &
                             ' on input state' )
#ifdef ESMF_3
     !call ESMF_FieldGetInternArray( field=Field, array=array, rc=rc )
     call shr_sys_abort( "Do not have get array function" )
#else
     call ESMF_FieldGetArray( field=Field, array=array, rc=rc )
#endif
     call eshr_rc_check( rc, ' : Error getting array from field '// &
                         FieldNames(j)// ' on input state' )
     call ESMF_ArrayGet( array=array, rank=rank, type=dataType, kind=kind, rc=rc )
     call eshr_rc_check( rc, ': Error in getting array type info from field' )
     if ( (dataType == ESMF_DATA_REAL) .and. (kind == ESMF_R8) .and. (rank == 2) )then
         call ESMF_FieldGetDataPointer( field=field,   &
                                        ptr=data2D,    &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, ': Error in getting F90 ptr to field data' )
         do n = 1, size(Data2D,2)
             do i = 1, size(Data2D,1)
                call random_number( x )
                x           = x*100.0
                data2D(i,n) = x
                if ( x < xmin ) xmin = x
                if ( x > xmax ) xmax = x
                xsum = xsum + x
                ntag = ntag + 1
             end do
         end do
         do n = 1, size(Data2D,2)
             do i = 1, size(Data2D,1)
                xsq = xsq + (data2D(i,n) - (xsum/ntag))**2
             end do
         end do
         nullify( Data2D )
     else if ( (dataType == ESMF_DATA_REAL) .and. (kind == ESMF_R8) .and. (rank == 1) &
     )then
         call ESMF_FieldGetDataPointer( field=field,   &
                                        ptr=data1D,    &
                                        copyflag=ESMF_DATA_REF, rc=rc )
         call eshr_rc_check( rc, ': Error in getting F90 ptr to field data' )
         do i = 1, size(Data1D)
            call random_number( x )
            x         = x*100.0
            data1D(i) = x
            if ( x < xmin ) xmin = x
            if ( x > xmax ) xmax = x
            xsum = xsum + x
            ntag = ntag + 1
         end do
         do i = 1, size(Data1D)
            xsq = xsq + (data1D(i) - (xsum/ntag))**2
         end do
         nullify( Data1D )
     end if

  end do
  if ( eshr_estate_matchValue( a2x_a, SHR_CONST_SPVAL) )then
     call shr_sys_abort( "Found missing value on atm state" )
  end if
  if ( ntag > 0 )then
    write(6,1010) xmin, xsum/ntag, xmax, sqrt(xsq/ntag)
 1010 format( "Statistics on state should be near: min=", f4.1, ", avg=", f5.1, &
               ", max=", f5.1, ", std-dev=", f5.1 )
    call eshr_estate_getStats( a2x_a )
  end if
  ! Destroy states and then fields
  call eshr_estate_destroy( a2x_a )
  call eshr_estate_destroy( i2x_i )
  call eshr_estate_destroy( l2x_l )
  call eshr_estate_destroy( o2x_o )
  call eshr_estate_destroy( domain_a )
  call eshr_estate_destroy( domain_i )
  call eshr_estate_destroy( domain_l )
  call eshr_estate_destroy( domain_o )
  if ( nfields > 0 ) deallocate( FieldNames )
  return
end function check_EstateImportCreate

logical function check_EstateCopy( estateIn, eStateOut, nf, arb )
  type(ESMF_State), intent(inout)    :: estateIn  ! EState to copy
  type(ESMF_State), intent(inOUT) :: estateOut ! Output EState
  integer, intent(in) :: nf               ! Number of output fields expected
  logical, intent(in) :: arb              ! If arbitrary grid or not

  integer :: rank, j
  character(len=SHR_KIND_CL), pointer       :: FieldNames(:) ! Fieldnames to copy

  logical :: check
  integer :: nfields

  check_EstateCopy = .true.
  rank = 2
  if ( arb ) rank = 1
  ! ----- If list of names different
  if ( .not. eshr_estate_fieldsAreSame( eStateIn, eStateOut, notData=.true. ) )then
     call eshr_estate_getSharedFNames( eStateIn, eStateOut, Fieldnames, nFields )
     if ( nFields == 0 ) check_EstateCopy = .false.
  else
     call eshr_estate_getFieldNames( eStateIn, Fieldnames, nFields )
     if ( nFields == 0 ) check_EstateCopy = .false.
  end if
  ! ----- Check that list of fields is correct ------
  if ( .not. eshr_estate_fieldsAreSame( eStateIn, eStateOut, notData=.true., &
                                        FieldNames=FieldNames ) ) &
     check_EstateCopy = .false.
  ! -- Copy the fields ----
  if ( check_EstateCopy ) &
      call eshr_estate_copyFields( Init=.true., rank=rank, eStateIn=eStateIn,  &
                               eStateOut=eStateOut, FieldNames=FieldNames, &
                               nFields=nfields )
  if ( .not. eshr_estate_fieldsAreSame( eStateIn, eStateOut, FieldNames=FieldNames ) ) &
     check_EstateCopy = .false.
  if ( nf /= nfields ) check_EstateCopy = .false.
  if ( check_EstateCopy )then
     ! Loop over a bunch of times so that can test if Open-MP is faster
     do j = 1, 1000
        ! ---- Copy fields without initialization
        if ( check_EstateCopy ) &
        call eshr_estate_copyFields( Init=.false., rank=rank, eStateIn=eStateIn, &
                                     eStateOut=eStateOut, FieldNames=FieldNames, &
                                     nFields=nfields )
        if ( .not. eshr_estate_fieldsAreSame( eStateIn, eStateOut, FieldNames=FieldNames ) ) &
           check_EstateCopy = .false.
     end do
  end if
  ! Clean up
  if ( associated( FieldNames ) ) deallocate( FieldNames )
  return
end function check_EstateCopy

logical function check_EstateMerge( estateIce, estateLnd, estateOcn, eStateOut, &
                                    eStateXao, eStateExpect )
  use shr_string_mod,  only: shr_string_listGetName
  use eshr_estate_mod, only: ICE_IDX, OCN_IDX, LND_IDX, ICE_IDX, XAO_IDX
  implicit none
  type(ESMF_State), intent(inout)    :: estateIce ! Ice EState
  type(ESMF_State), intent(inout)    :: estateLnd ! Land  EState
  type(ESMF_State), intent(inout)    :: estateOcn ! Ocean EState
  type(ESMF_State), intent(inout)    :: estateXao ! Atm/Ocean EState
  type(ESMF_State), intent(inOUT) :: estateOut    ! Output EState
  type(ESMF_State), intent(inOUT) :: estateExpect ! What expect the output estate to look like

  integer :: rank, j, n, n1, i
  character(len=SHR_KIND_CL), pointer :: FieldNames(:,:) ! Fieldnames to merge
  character(len=SHR_KIND_CL), pointer :: FieldNamesE(:) ! Fieldnames on expect
  character(len=SHR_KIND_CL), pointer :: FieldNeedO(:)   ! OCN Fields needed
  character(len=SHR_KIND_CL), pointer :: FieldNeedXao(:)   ! XAO Fields needed
  character(len=SHR_KIND_CX)          :: FieldsNot          ! Fields not needed
  character(len=SHR_KIND_CL)          :: FieldName ! Fieldname

  integer :: nfields, nGrPts, nfieldsE, nFieldsNot
  real(SHR_KIND_R8), pointer :: IceFrac(:), LndFrac(:), OcnFrac(:)

  ! if grids don't match -- check is false
  check_EstateMerge = eshr_estate_checkGridsMatch( eStateIce, eStateOut, estateOcn, &
                                                   eStateLnd, gridsExact=.false. )
  if ( check_EstateMerge ) &
  check_EstateMerge = eshr_estate_checkGridsMatch( eStateIce, eStateOut, estateXao, &
                                                   eStateLnd, gridsExact=.false. )
  ! Initialize the merge
  if ( check_EstateMerge )then
     call eshr_estate_mergeInit( eStateIce, eStateLnd, eStateOcn, eStateOut,       &
                                 eStateXao, Fieldnames, IceFrac, LndFrac, OcnFrac, &
                                 nFields, nGrPts )
     ! Check if the initial merge failed
     if ( nFields == 0 ) check_EstateMerge = .false.
     if ( nGrPts  == 0 ) check_EstateMerge = .false.
     if ( .not. associated(IceFrac) ) check_EstateMerge = .false.
     if ( .not. associated(OcnFrac) ) check_EstateMerge = .false.
     if ( .not. associated(LndFrac) ) check_EstateMerge = .false.
     if ( .not. check_EstateMerge ) write(6,*) "MergeInit did not run successfully"
     ! If successful mark extra fields as not needed
     if ( check_EstateMerge )then
        call eshr_estate_markRestNotNeed( eStateIce, FieldNames(ICE_IDX,:), nFields )
        call eshr_estate_getFieldsNotNeed(  eStateIce, FieldsNot, nFieldsNot )
        do i = 1, nFields
           do j = 1, nFieldsNot
              call shr_string_listGetName( FieldsNot, j, FieldName )
              if ( trim(FieldNames(ICE_IDX,i)) == trim(FieldName) ) &
                 call shr_sys_abort( " : Error, not needed fields NOT correct" )
           end do
        end do
        call eshr_estate_markRestNotNeed( eStateLnd, FieldNames(LND_IDX,:), nFields )
        call eshr_estate_getFieldsNotNeed(  eStateLnd, FieldsNot, nFieldsNot )
        do i = 1, nFields
           do j = 1, nFieldsNot
              call shr_string_listGetName( FieldsNot, j, FieldName )
              if ( trim(FieldNames(LND_IDX,i)) == trim(FieldName) ) &
                 call shr_sys_abort( " : Error, not needed fields NOT correct" )
           end do
        end do
        allocate( FieldNeedO(nFields), FieldNeedXao(nFields) )
        n  = 0
        n1 = 0
        do i = 1, nFields
          if ( trim(FieldNames(XAO_IDX,i)) == "xao" )then
            n = n + 1
            FieldNeedXao(n) = FieldNames(OCN_IDX,i)
          else
            n1 = n1 + 1
            FieldNeedO(n1) = FieldNames(OCN_IDX,i)
          end if
        end do
        if ( n+n1 /= nFields ) &
           call shr_sys_abort( " : Error, ocean fields not split correctly" )
        call eshr_estate_markRestNotNeed( eStateOcn, FieldNeedO(:n1), n1 )
        call eshr_estate_getFieldsNotNeed(  eStateOcn, FieldsNot, nFieldsNot )
        do i = 1, nFields
           do j = 1, nFieldsNot
              call shr_string_listGetName( FieldsNot, j, FieldName )
              if ( trim(FieldNeedO(i)) == trim(FieldName) ) &
                 call shr_sys_abort( " : Error, not needed fields NOT correct" )
           end do
        end do
        call eshr_estate_markRestNotNeed( eStateXao, FieldNeedXao(:n), n )
        call eshr_estate_getFieldsNotNeed(  eStateXao, FieldsNot, nFieldsNot )
        do i = 1, nFields
           do j = 1, nFieldsNot
              call shr_string_listGetName( FieldsNot, j, FieldName )
              if ( trim(FieldNeedXao(i)) == trim(FieldName) ) &
                 call shr_sys_abort( " : Error, not needed fields NOT correct" )
           end do
        end do
        deallocate( FieldNeedO, FieldNeedXao )
     end if
  end if
  ! Actually do the merge
  if ( check_EstateMerge )then
     call eshr_estate_getFieldNames( eStateExpect,  FieldNamesE,  nFieldsE )
     if ( nFieldsE == 0 ) check_EstateMerge = .false.
     ! Loop over a bunch of times so that can test if Open-MP is faster
     do j = 1, 1000
        if ( check_EstateMerge ) &
           call eshr_estate_merge( EStateIce, EStateLnd, EStateOcn, EStateXao, nFields, &
                                   nGrPts, FieldNames, IceFrac, LndFrac, OcnFrac, EStateOut )
        ! Now check the merge
        if ( check_EstateMerge )then
           check_EstateMerge = eshr_estate_fieldsAreSame( eStateOut, estateExpect, &
                                                    Fieldnames=FieldnamesE,        &
                                                    onlyNeeded=.true. )
        end if
     end do
  end if
  ! Deallocate data
  if ( associated(Fieldnames) ) deallocate( Fieldnames )
  if ( associated(IceFrac) ) nullify( IceFrac )
  if ( associated(LndFrac) ) nullify( LndFrac )
  if ( associated(OcnFrac) ) nullify( OcnFrac )
  return
end function check_EstateMerge


subroutine prepare_MergeExpect( lMask_a, eStateIce, eStateOcn, eStateLnd, eStateOut, eStateXao, eStateExpect )

  implicit none
  integer,          intent(in)    :: lMask_a(:)
  type(ESMF_State), intent(inout) :: eStateIce
  type(ESMF_State), intent(inout) :: eStateOcn
  type(ESMF_State), intent(inout) :: eStateLnd
  type(ESMF_State), intent(inout) :: eStateOut
  type(ESMF_State), intent(inout) :: eStateXao
  type(ESMF_State), intent(inout) :: eStateExpect

  integer :: i, rc, n, j, m, l
  character(len=SHR_KIND_CL) :: fname, fname2
  real(SHR_KIND_R8), pointer :: data1D(:), data1DOut(:)
  integer                       :: nFieldsIce, nFieldsOcn, nFieldsLnd, nFieldsExpect
  integer                       :: nFieldsOut, nFieldsXao
  character(len=SHR_KIND_CS), pointer :: FieldNamesIce(:)
  character(len=SHR_KIND_CS), pointer :: FieldNamesOcn(:)
  character(len=SHR_KIND_CS), pointer :: FieldNamesOut(:)
  character(len=SHR_KIND_CS), pointer :: FieldNamesLnd(:)
  character(len=SHR_KIND_CS), pointer :: FieldNamesXao(:)
  character(len=SHR_KIND_CS), pointer :: FieldNamesExpect(:)
  type(ESMF_State)                    :: domain_i, domain_a, domain_l, domain_o
  type(ESMF_Grid)                     :: agrid, igrid, lgrid, ogrid
  type(ESMF_Field)                    :: Field

  !
  ! Get fields from states
  !
  ! Get names and  number of fields on states
  call eshr_estate_getFieldNames( eStateIce,    FieldNamesIce,    nFieldsIce )
  call eshr_estate_getFieldNames( eStateOcn,    FieldNamesOcn,    nFieldsOcn )
  call eshr_estate_getFieldNames( eStateOut,    FieldNamesOut,    nFieldsOut, &
                                  onlyNeeded=.true. )
  call eshr_estate_getFieldNames( eStateLnd,    FieldNamesLnd,    nFieldsLnd )
  call eshr_estate_getFieldNames( eStateXao,    FieldNamesXao,    nFieldsXao )
  call eshr_estate_getFieldNames( eStateExpect, FieldNamesExpect, nFieldsExpect )
  ! Set input fractions to constants
  do i = 1, nFieldsOut
     fname = FieldNamesOut(i)
     if ( index(fname,"ifrac") /= 0 )then
       call ESMF_StateGetDataPointer( estateOut,      &
                                      fname,          &
                                      data1D,    &
                                      copyflag=ESMF_DATA_REF, rc=rc )
       call eshr_rc_check( rc, ': Error in getting data ptr' )
       data1D(:) = icefrac
     else if ( index(fname,"lfrac") /= 0 )then
       call ESMF_StateGetDataPointer( estateOut,     &
                                      fname,          &
                                      data1D,    &
                                      copyflag=ESMF_DATA_REF, rc=rc )
       call eshr_rc_check( rc, ': Error in getting data ptr' )
       data1D(:) = lndfrac
     else if ( index(fname,"ofrac") /= 0 )then
       call ESMF_StateGetDataPointer( estateOut,   &
                                      fname,          &
                                      data1D,    &
                                      copyflag=ESMF_DATA_REF, rc=rc )
       call eshr_rc_check( rc, ': Error in getting data ptr' )
       data1D(:) = ocnfrac
     end if
     nullify(data1D)
  end do
  ! Calculate average
  do i = 1, min(nFieldsExpect,nFieldsOut)
     fname = FieldNamesExpect(i)
     call ESMF_StateGetDataPointer( estateExpect,   &
                                    fname,          &
                                    data1DOut,    &
                                    copyflag=ESMF_DATA_REF, rc=rc )
     call eshr_rc_check( rc, ': Error in getting data ptr' )
     n = index(fname,"_")
     n = n + 1
     l = index(fname," ")
     data1DOut(:) = 0.0_SHR_KIND_R8
     ! Add in land
     do j = 1, nFieldsLnd
        fname2 = FieldNamesLnd(j)
        m = index(fname2,"_")
        m = m + 1
        if ( trim(fname2(m:l)) == trim(fname(n:l)) )then
           call ESMF_StateGetDataPointer( eStateLnd,   &
                                          fname2,          &
                                          data1D,    &
                                          copyflag=ESMF_DATA_REF, rc=rc )
           call eshr_rc_check( rc, ': Error in getting data ptr' )
           data1DOut(:) = data1dOut(:) + lndfrac*data1D(:size(data1DOut))
           nullify(data1D)
        end if
     end do
     ! Add in ice
     do j = 1, nFieldsIce
        fname2 = FieldNamesIce(j)
        m = index(fname2,"_")
        m = m + 1
        if ( trim(fname2(m:l)) == trim(fname(n:l)) )then
           call ESMF_StateGetDataPointer( eStateIce,   &
                                          fname2,          &
                                          data1D,    &
                                          copyflag=ESMF_DATA_REF, rc=rc )
           call eshr_rc_check( rc, ': Error in getting data ptr' )
           data1DOut(:) = data1dOut(:) + icefrac*data1D(:size(data1DOut))
           nullify(data1D)
        end if
     end do
     ! Add in ocean
     do j = 1, nFieldsOcn
        fname2 = FieldNamesOcn(j)
        m = index(fname2,"_")
        m = m + 1
        if ( trim(fname2(m:l)) == trim(fname(n:l)) )then
           call ESMF_StateGetDataPointer( eStateOcn,  &
                                          fname2,          &
                                          data1D,    &
                                          copyflag=ESMF_DATA_REF, rc=rc )
           call eshr_rc_check( rc, ': Error in getting data ptr' )
           data1DOut(:) = data1dOut(:) + ocnfrac*data1D(:size(data1DOut))
           nullify(data1D)
        end if
     end do
     ! Add in xao
     do j = 1, nFieldsXao
        fname2 = FieldNamesXao(j)
        m = index(fname2,"_")
        m = m + 1
        if ( trim(fname2(m:l)) == trim(fname(n:l)) )then
           call ESMF_StateGetDataPointer( eStateXao,  &
                                          fname2,          &
                                          data1D,    &
                                          copyflag=ESMF_DATA_REF, rc=rc )
           call eshr_rc_check( rc, ': Error in getting data ptr' )
           data1DOut(:) = data1dOut(:) + ocnfrac*data1D(:size(data1DOut))
           nullify(data1D)
        end if
     end do
     nullify(data1DOut)
  end do
  ! Create domains and put on states
  call ESMF_StateGetField( eStateExpect, FieldName=FieldNamesExpect(1), Field=Field, rc=rc )
  call eshr_rc_check( rc, ': Error in getting field' )
  call ESMF_FieldGet( Field, grid=agrid, rc=rc )
  call eshr_rc_check( rc, ': Error in getting agrid from field' )
  call ESMF_StateGetField( eStateLnd, FieldName=FieldNamesLnd(1), Field=Field, rc=rc )
  call eshr_rc_check( rc, ': Error in getting field' )
  call ESMF_FieldGet( Field, grid=lgrid, rc=rc )
  call eshr_rc_check( rc, ': Error in getting lgrid from field' )
  call ESMF_StateGetField( eStateIce, FieldName=FieldNamesIce(1), Field=Field, rc=rc )
  call eshr_rc_check( rc, ': Error in getting field' )
  call ESMF_FieldGet( Field, grid=igrid, rc=rc )
  call eshr_rc_check( rc, ': Error in getting igrid from field' )
  call ESMF_StateGetField( eStateOcn, FieldName=FieldNamesOcn(1), Field=Field, rc=rc )
  call eshr_rc_check( rc, ': Error in getting field' )
  call ESMF_FieldGet( Field, grid=ogrid, rc=rc )
  call eshr_rc_check( rc, ': Error in getting ogrid from field' )
  call create_domains( agrid, igrid, lgrid, ogrid, lMask_a, domain_a, domain_i, &
                       domain_l, domain_o )

  call ESMF_StateAddState( eStateIce, domain_i, rc=rc )
  call eshr_rc_check( rc, " : Error adding domain_i to eStateIce" )
  call ESMF_StateAddState( eStateOut, domain_a, rc=rc )
  call eshr_rc_check( rc, " : Error adding domain_a to eSTateOut" )
  call ESMF_StateAddState( eStateExpect, domain_a, rc=rc )
  call eshr_rc_check( rc, " : Error adding domain_a to eSTateExpect" )
  call ESMF_StateAddState( eStateLnd, domain_l, rc=rc )
  call eshr_rc_check( rc, " : Error adding domain_l to EStateLnd" )
  call ESMF_StateAddState( eStateOcn, domain_o, rc=rc )
  call eshr_rc_check( rc, " : Error adding domain_o to EStateOcn" )
  call ESMF_StateAddState( eStateXao, domain_o, rc=rc )
  call eshr_rc_check( rc, " : Error adding domain_o to eStateXao" )
  ! de-Allocate
  if ( nFieldsIce    > 0 ) deallocate( FieldNamesIce    )
  if ( nFieldsOcn    > 0 ) deallocate( FieldNamesOcn    )
  if ( nFieldsLnd    > 0 ) deallocate( FieldNamesLnd    )
  if ( nFieldsXao    > 0 ) deallocate( FieldNamesXao    )
  if ( nFieldsExpect > 0 ) deallocate( FieldNamesExpect )

end subroutine prepare_MergeExpect


logical function check_EstateRedist( estateIn, estateOut, Mask )
  use shr_string_mod,  only: shr_string_listGetNum, shr_string_listMerge
  type(ESMF_State), intent(INout)  :: estateIn  ! EState to redist
  type(ESMF_State), intent(inOUT)  :: estateOut ! Output EState
  integer,          intent(IN)     :: Mask(:)

  type(ESMF_Bundle)           :: EBundleIn     ! In bundle of fields to redist
  type(ESMF_Bundle)           :: EBundleOut    ! Out bundle of fields to redist
  type(ESMF_Bundle)           :: EBundle1      ! In bundle of fields to redist backward
  type(ESMF_Bundle)           :: EBundle2      ! Out bundle of fields to redist backward
  type(ESMF_RouteHandle)      :: ERoute        ! Comm. handle to do the redist
  type(ESMF_RouteHandle)      :: ERouteBack    ! Comm. handle to do the redist backwards
  type(ESMF_State)            :: estateBack    ! Output EState from do redist backwards
  type(ESMF_Field)            :: Field         ! Field to get off input state
  type(ESMF_Grid)             :: Grid          ! Grid of input estate
  type(ESMF_Grid)             :: GridOut       ! Grid of output estate
  integer                     :: nfields       ! Number of fields on back estate
  integer                     :: nf            ! Number of fields on back estate
  integer                     :: nf1, nf2      ! Number of fields not needed
  integer                     :: i             ! Index
  character(len=SHR_KIND_CS), pointer   :: FieldNames(:) ! Fieldnames to copy
  character(len=SHR_KIND_CX)  :: FieldList     ! Fieldnames to copy
  character(len=SHR_KIND_CX)  :: NotNeededList ! Fieldnames not needed
  character(len=SHR_KIND_CX)  :: NotNeededList1! Fieldnames not needed
  character(len=SHR_KIND_CX)  :: NotNeededList2! Fieldnames not needed
  integer :: rc
!
  check_EstateRedist = .true.
  ! ----- Check that list of fields is correct ------
  if ( .not. eshr_estate_fieldsAreSame( eStateIn, eStateOut, notData=.true., &
  namesOnly=.true. ) )then
     check_EstateRedist = .false.
     return
  end if
  ! ----- Create and Put list of fields on Back estate ------
  call eshr_estate_getSharedFNames( eStateIn, eStateOut, FieldNames, nFields )
  if ( nfields == 0 )then
     deallocate( FieldNames )
     check_EstateRedist = .false.
     return
  end if
  FieldList = trim(FieldNames(1))
  do i = 2, nfields
     FieldList = trim(FieldList)//":"//trim(FieldNames(i))
  end do
  nf = shr_string_listGetNum( FieldList )
  if ( nf /= nfields ) call shr_sys_abort( "Field names constructed wrong" )
  call ESMF_StateGetField( eStateIn, FieldName=FieldNames(1), Field=Field, rc=rc )
  call eshr_rc_check( rc, ': Error in getting field from estate' )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_FieldGet( Field, grid=Grid, rc=rc )
  call eshr_rc_check( rc, ': Error in getting grid from field' )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  eStateBack = ESMF_StateCreate( "Temp state to compare to Input state", rc=rc )
  call eshr_rc_check( rc, "Error in creation of ESMF State" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  ! Get output grid
  call ESMF_StateGetField( eStateOut, FieldName=FieldNames(1), Field=Field, rc=rc )
  call eshr_rc_check( rc, ': Error in getting field from estate' )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_FieldGet( Field, grid=GridOut, rc=rc )
  call eshr_rc_check( rc, ': Error in getting grid from field' )
  ! Compare grids -- make sure same
  if( .not. eshr_estate_Egrids2DAreSame( Grid, GridOut, Exact=.true., Grid1Mask=Mask ) &
  )then
     check_EstateRedist = .false.
     deallocate( Fieldnames )
     return
  end if
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call eshr_estate_getFieldsNotNeed( eStateIn,   notNeededList1, nf1 )
  call eshr_estate_getFieldsNotNeed( eStateOut,  notNeededList2, nf2 )
  call shr_string_listMerge( notNeededList1, notNeededList2, notNeededList )
  write(6,*) "FieldList = ", trim(FieldList)
  if ( nf1+nf2 > 0 )then
     write(6,*) "notNeededList = ", trim(notNeededList)
     call eshr_estate_init2DFromList(  Grid, FieldList,  &
                                       eState=eStateBack, notNeeded=notNeededList )
     nf  = shr_string_listGetNum( notNeededList )
     nf1 = shr_string_listGetNum( FieldList )
     if ( nf == nf1 )then
        check_EstateRedist = .false.
        return
     end if
  else
     call eshr_estate_init2DFromList(  Grid, FieldList, eState=eStateBack )
  end if
  ! ------ Initialize the redistribution -------
  call eshr_estate_reDistInit( eStateIn,  eStateOut,  &
                               EBundleIn, EBundleOut, eRoute, noGridCheck=.true. )
  call eshr_estate_reDistInit( eStateOut, eStateBack, &
                               EBundle1,  EBundle2,   eRouteBack, noGridCheck=.true. )
  ! ------- Run the redistribution --------
  call eshr_estate_reDist( eBundleIn, eBundleOut, eRoute     )
  call eshr_estate_reDist( eBundle1,  eBundle2,   eRouteBack )
  ! ------- Check that resultant back state is same as input state ----
  if ( .not. eshr_estate_fieldsAreSame( eStateIn, eStateBack, onlyNeeded=.true., &
       Mask1D=Mask ) ) &
        check_EstateRedist = .false.
  ! ----- Clean up
  call ESMF_StateDestroy(  EStateBack )
  call ESMF_BundleDestroy( EBundleIn  )
  call ESMF_BundleDestroy( EBundleOut )
  call ESMF_BundleDestroy( EBundle1   )
  call ESMF_BundleDestroy( EBundle2   )
  deallocate( Fieldnames )
  call ESMF_BundleRedistRelease( eRoute     )
  call ESMF_BundleRedistRelease( eRouteBack )
  return
end function check_EstateRedist

subroutine fill_EstateField( eState, arb, type )
   use shr_const_mod,           only: SHR_CONST_SPVAL
   implicit none
   type(ESMF_State), intent(inout) :: eState
   logical, intent(in) :: arb
   integer, intent(in), optional :: type

   integer :: nfields
   real(SHR_KIND_R8), pointer :: data1d(:)
   real(SHR_KIND_R8) :: x
   character(len=SHR_KIND_CL), pointer :: FieldNames(:)
   integer :: i, j, n, FillType
   integer :: rc

   if ( .not. present(type) )then
      FillType = RANDOM_FILL
   else
      FillType = type
   end if
   call eshr_estate_getFieldNames( eState, FieldNames, nFields, onlyNeeded=.true. )
   if ( nfields == 0 ) return
   if ( .not. arb ) call shr_sys_abort( "Fill only works for Arbitrary decomp fields" )
   do n = 1, nfields
      call ESMF_StateGetDataPointer( eState, FieldNames(n), data1d, copyflag=ESMF_DATA_REF, &
                                     rc=rc )
      if ( FillType == RANDOM_FILL )then
         do i = 1, size(Data1D)
            call random_number( x )
            data1D(i) = x*1000.0
         end do
         ! Some random missing data points
         do j = 1, 10
            call random_number( x )
            i = nint(x*(size(Data1D)-1)) + 1
            data1D(i) = SHR_CONST_SPVAL
         end do
      else if ( FillType == CONSTANT_FILL )then
         call random_number( x )
         x = x * 1.e5
         data1D(:) = x
      else
         call shr_sys_abort( "Invalid fill type" )
      end if
      nullify( Data1D )
   end do
   if ( FillType == RANDOM_FILL )then
       if ( .not. eshr_estate_matchValue( estate, SHR_CONST_SPVAL) ) &
           call shr_sys_abort( "missing value NOT found on a2x_a state" )
   else if ( FillType == CONSTANT_FILL )then
       if ( eshr_estate_matchValue( estate, SHR_CONST_SPVAL) ) &
           call shr_sys_abort( "missing value found on a2x_a state" )
   end if
   deallocate( FieldNames )

end subroutine fill_EstateField

end module test_eshr_estate_mod

program test_eshr_estate
!
! Program to test the eshr EState modules
!
  use shr_kind_mod,            only: r8=>SHR_KIND_R8, SHR_KIND_CX
  use shr_inputinfo_mod,       only: shr_inputinfo_initType,        &
                                     shr_inputinfo_initSetDefaults, &
                                     shr_inputinfo_initPutData
  use test_eshr_inputinfo_mod, only: check_ccsminitConversion
  use test_eshr_timemgr_mod,   only: check_clockInfoConversion
  use test_eshr_estate_mod,    only: check_FilenameConversion,      &
                                     check_estateInit,              &
                                     check_EstateImportCreate,      &
                                     check_EstateCopy,              &
                                     check_EstateRedist,            &
                                     fill_EstateField,              &
                                     CONSTANT_FILL,                 &
                                     check_EstateMerge,             &
                                     prepare_MergeExpect
  use eshr_estate_mod,         only: eshr_estate_Egrids2DAreSame,   &
                                     eshr_estate_init2DFromList,    &
                                     eshr_estate_getStats,          &
                                     eshr_estate_destroy,           &
                                     eshr_estate_getFieldnames,     &
                                     eshr_estate_markRestNotNeed      
  use shr_sys_mod,             only: shr_sys_abort, shr_sys_flush
  use eshr_rc_mod,             only: eshr_rc_check, eshr_rc_set
  use eshr_timemgr_mod,        only: eshr_timemgr_clockInfoType,     &
                                     eshr_timemgr_NMLinfoSetDefault, &
                                     eshr_timemgr_NMLinfoPutData,    &
                                     eshr_timemgr_clockInitNMLinfo,  &
                                     eshr_timemgr_clockGet,          &
                                     eshr_timemgr_NMLinfoType,       &
                                     eshr_timemgr_clockType
  use shr_orb_mod,             only: SHR_ORB_UNDEF_INT
  use GridnDecomp_mod,         only: GridnDecomp_CenterToVertices, &
                                     GridnDecomp_CalcDecomp,       &
                                     GridnDecomp_GetLocalGIndices, &
                                     setRandomSeed
  use shr_string_mod,          only: shr_string_listGetNum
  use ESMF_Mod
  use shr_mpi_mod,             only: shr_mpi_init, shr_mpi_commrank, shr_mpi_finalize, &
                                     shr_mpi_barrier
  implicit none
  type(shr_inputinfo_initType) :: CCSMInit
  type(eshr_timemgr_clockInfoType) :: ClockInfo
  type(eshr_timemgr_clockType) :: Clock
  type(eshr_timemgr_NMLinfoType) :: clocknmlinfo
  character(len=80) :: Filetype, Filename
  character(len=80), pointer :: Fieldnames(:)
  integer :: nfields
  character(len=*), parameter :: no_rtnFMT = '(A20," ", A40," ...... ", T70, $)'
  character(len=80) :: name
   !----------------------------------------------------------------------------
   character(*), parameter :: seq_flds_i2x_states = &
         'Si_t'        &    ! temperature                     DEF
      //':Si_tref'     &    ! 2m reference temperature        DEF
      //':Si_qref'     &    ! 2m reference specific humidity  DEF
      //':Si_avsdr'    &    ! albedo: visible, direct         DEF
      //':Si_anidr'    &    ! albedo: near ir, direct         DEF
      //':Si_avsdf'    &    ! albedo: visible, diffuse        DEF
      //':Si_anidf'    &    ! albedo: near ir, diffuse        DEF
      //':Si_ifrac'    &    ! fractional ice cov wrt gridcell DEF
      //':Si_aice'     &    ! fractional ice cov wrt ocean    DEF
      //':Si_sicthk'        ! sea ice thickness (m)           DEF
  ! Fluxes
   character(*), parameter :: seq_flds_i2x_fluxes = &
         'Faii_taux'   &    ! wind stress, zonal              DEF
      //':Faii_tauy'   &    ! wind stress, meridional         DEF
      //':Faii_lat'    &    ! latent          heat flux       DEF
      //':Faii_sen'    &    ! sensible        heat flux       DEF
      //':Faii_lwup'   &    ! upward longwave heat flux       DEF
      //':Faii_evap'   &    ! evaporation    water flux       DEF
      //':Faii_swnet'  &    ! shortwave: net absorbed         DEF
      //':Fioi_melth'       ! heat  flux from melting ice     DEF

   ! States
   character(*), parameter :: seq_flds_l2x_states = &
         'Sl_t'        &    ! temperature                     DEF
      //':Sl_tref'     &    ! 2m reference temperature        DEF
      //':Sl_qref'     &    ! 2m reference specific humidity  DEF
      //':Sl_avsdr'    &    ! albedo: direct , visible        DEF
      //':Sl_anidr'    &    ! albedo: direct , near-ir        DEF
      //':Sl_avsdf'    &    ! albedo: diffuse, visible        DEF
      //':Sl_anidf'    &    ! albedo: diffuse, near-ir        DEF
      //':Sl_snowh'    &    ! snow height                     DEF
      //':Sl_landfrac'       ! fractional land                 DEF
  ! Fluxes
   character(*), parameter :: seq_flds_l2x_fluxes = &
         'Fall_taux'    &   ! wind stress, zonal              DEF
      //':Fall_tauy'    &   ! wind stress, meridional         DEF
      //':Fall_lat'     &   ! latent          heat flux       DEF
      //':Fall_sen'     &   ! sensible        heat flux       DEF
      //':Fall_lwup'    &   ! upward longwave heat flux       DEF
      //':Fall_evap'    &   ! evaporation    water flux       DEF
      //':Fall_swnet'   &   ! shortwave: net absorbed         DEF
      //':Fall_flxdst1' &   ! dust flux bin 1
      //':Fall_flxdst2' &   ! dust flux bin 2
      //':Fall_flxdst3' &   ! dust flux bin 3
      //':Fall_flxdst4'     ! dust flux bin 4

   character(*), parameter :: seq_flds_o2x_states = &
         'So_tsocn'    &    ! ocean layer temperature         DEF
      //':So_t'             ! temperature                     DEF
   ! Fluxes
   character(*), parameter :: seq_flds_o2x_fluxes = &
         'Fioo_q'            ! heat of fusion (q>0) melt pot (q<0)  DEF 
   ! States
   character(*), parameter :: seq_flds_xao_states = &
         'So_tref'     &    ! 2m reference temperature        DEF
      //':So_qref'     &    ! 2m reference specific humidity  DEF
      //':So_avsdr'    &    ! albedo: visible, direct         DEF
      //':So_anidr'    &    ! albedo: near ir, direct         DEF
      //':So_avsdf'    &    ! albedo: visible, diffuse        DEF
      //':So_anidf'    &    ! albedo: near ir, diffuse        DEF
      //':So_ustar'    &    ! ustar                           DEF
      //':So_ssq'      &    ! surface saturation spec. hum.   DEF
      //':So_re'            ! sqrt of exch. coeff (tracers)   DEF

   ! Fluxes
   character(*), parameter :: seq_flds_xao_fluxes = &
         'Fioo_q'       &    ! heat of fusion (q>0) melt pot (q<0)  DEF 
      //':Faox_taux'    &   ! wind stress, zonal              DEF
      //':Faox_tauy'    &   ! wind stress, meridional         DEF
      //':Faox_lat'     &   ! latent          heat flux       DEF
      //':Faox_sen'     &   ! sensible        heat flux       DEF
      //':Faox_lwup'    &   ! upward longwave heat flux       DEF
      //':Faox_evap'    &   ! evaporation    water flux       DEF
      //':Faox_swnet'       ! shortwave: net absorbed         DEF
      
   ! States
   character(*), parameter :: seq_flds_a2x_states = &
         'Sa_z'        &    ! bottom atm level height         DEF
      //':Sa_u'        &    ! bottom atm level zon wind       DEF
      //':Sa_v'        &    ! bottom atm level mer wind       DEF
      //':Sa_tbot'     &    ! bottom atm level temp           DEF
      //':Sa_ptem'     &    ! bottom atm level pot temp       DEF
      //':Sa_shum'     &    ! bottom atm level spec hum       DEF
      //':Sa_dens'     &    ! bottom atm level air den        DEF
      //':Sa_pbot'     &    ! bottom atm level pressurea      DEF
      //':Sa_pslv'          ! sea level atm pressure          DEF
   character(*), parameter :: seq_flds_a2x_state_units = &
         'm     '     &    ! bottom atm level height         DEF
      //':m/s   '     &    ! bottom atm level zon wind       DEF
      //':m/s   '     &    ! bottom atm level mer wind       DEF
      //':K     '     &    ! bottom atm level temp           DEF
      //':K     '     &    ! bottom atm level pot temp       DEF
      //':kg/kg '     &    ! bottom atm level spec hum       DEF
      //':kg/m^3'     &    ! bottom atm level air den        DEF
      //':Pa    '     &    ! bottom atm level pressurea      DEF
      //':Pa    '          ! sea level atm pressure          DEF
   character(*), parameter :: seq_flds_a2x_state_lnames = &
         'bottom atm level height    ' &
      //':bottom atm level zon wind  ' &
      //':bottom atm level mer wind  ' &
      //':bottom atm level temp      ' &
      //':bottom atm level pot temp  ' &
      //':bottom atm level spec hum  ' &
      //':bottom atm level air den   ' &
      //':bottom atm level pressurea ' &
      //':sea level atm pressure     '

   ! Fluxes
   character(*), parameter :: seq_flds_a2x_fluxes = &
         'Faxa_lwdn'   &    ! downward lw heat flux           DEF
      //':Faxa_rainc'  &    ! prec: liquid "convective"       DEF
      //':Faxa_rainl'  &    ! prec: liquid "large scale"      DEF
      //':Faxa_snowc'  &    ! prec: frozen "convective"       DEF
      //':Faxa_snowl'  &    ! prec: frozen "large scale"      DEF
      //':Faxa_swndr'  &    ! sw: nir direct  downward        DEF
      //':Faxa_swvdr'  &    ! sw: vis direct  downward        DEF
      //':Faxa_swndf'  &    ! nir diffuse downward            DEF
      //':Faxa_swvdf'  &    ! sw: vis diffuse downward        DEF
      //':Faxa_swnet'       ! sw: net                         DEF
   character(*), parameter :: seq_flds_a2x_flux_units = &
         'W/m^2 '  &    ! downward lw heat flux           DEF
      //':kg/m^2'  &    ! prec: liquid "convective"       DEF
      //':kg/m^2'  &    ! prec: liquid "large scale"      DEF
      //':kg/m^2'  &    ! prec: frozen "convective"       DEF
      //':kg/m^2'  &    ! prec: frozen "large scale"      DEF
      //':W/m^2 '  &    ! sw: nir direct  downward        DEF
      //':W/m^2 '  &    ! sw: vis direct  downward        DEF
      //':W/m^2 '  &    ! nir diffuse downward            DEF
      //':W/m^2 '  &    ! sw: vis diffuse downward        DEF
      //':W/m^2 '       ! sw: net                         DEF
   ! States
   character(*), parameter :: seq_flds_x2a_states = &
         'Sx_tref'     &    ! 2m reference temperature        DEF
      //':Sx_qref'     &    ! 2m reference specific humidity  DEF
      //':Sx_avsdr'    &    ! albedo, visible, direct         DEF
      //':Sx_anidr'    &    ! albedo, near-ir, direct         DEF
      //':Sx_avsdf'    &    ! albedo, visible, diffuse        DEF
      //':Sx_anidf'    &    ! albedo, near-ir, diffuse        DEF
      //':Sx_t'        &    ! surface temperature             DEF
      //':So_t'        &    ! sea surface temperature         DEF
      //':Sx_snowh'    &    ! surface snow depth              DEF
      //':Sx_lfrac'    &    ! surface land fraction           DEF
      //':Sx_ifrac'    &    ! surface ice fraction            DEF
      //':Sx_ofrac'    &    ! surface ocn fraction            DEF
      //':So_ustar'    &    ! needed for isoptope calc        DEF
      //':So_re'       &    ! needed for isoptope calc        DEF
      //':So_ssq'      &    ! needed for isoptope calc        DEF
      //':Sx_fv'       &    !
      //':Sx_ram1'          !

   ! Fluxes
   character(*), parameter :: seq_flds_x2a_fluxes = &
         'Faxx_taux'    &   ! wind stress, zonal              DEF
      //':Faxx_tauy'    &   ! wind stress, meridional         DEF
      //':Faxx_lat'     &   ! latent          heat flux       DEF
      //':Faxx_sen'     &   ! sensible        heat flux       DEF
      //':Faxx_lwup'    &   ! upward longwave heat flux       DEF
      //':Faxx_evap'    &   ! evaporation    water flux       DEF
      //':Faxx_flxdst1' &   ! dust flux bin 1
      //':Faxx_flxdst2' &   ! dust flux bin 2
      //':Faxx_flxdst3' &   ! dust flux bin 3
      //':Faxx_flxdst4'     ! dust flux bin 4
   character(*), parameter :: seq_flds_x2a_mrg_fields = &
         'Sx_tref'     &    ! 2m reference temperature        DEF
      //':Sx_qref'     &    ! 2m reference specific humidity  DEF
      //':Sx_avsdr'    &    ! albedo, visible, direct         DEF
      //':Sx_anidr'    &    ! albedo, near-ir, direct         DEF
      //':Sx_avsdf'    &    ! albedo, visible, diffuse        DEF
      //':Sx_anidf'    &    ! albedo, near-ir, diffuse        DEF
      //':Sx_t'        &    ! surface temperature             DEF
      //':Faxx_taux'    &   ! wind stress, zonal              DEF
      //':Faxx_tauy'    &   ! wind stress, meridional         DEF
      //':Faxx_lat'     &   ! latent          heat flux       DEF
      //':Faxx_sen'     &   ! sensible        heat flux       DEF
      //':Faxx_lwup'    &   ! upward longwave heat flux       DEF
      //':Faxx_evap'        ! evaporation    water flux       DEF

   character(*), parameter :: seq_flds_a2x_flux_lnames = &
         'downward lw heat flux   ' &
      //':prec liquid "convective ' &
      //':prec liquid "large scale' &
      //':prec frozen "convective ' &
      //':prec frozen "large scale' &
      //':sw nir direct  downward ' &
      //':sw vis direct  downward ' &
      //':nir diffuse downward    ' &
      //':sw vis diffuse downward ' &
      //':sw net                  '
   ! States
   character(*), parameter :: seq_flds_x2l_states = &
         'Sa_z'        &    ! bottom atm level height         DEF
      //':Sa_u'        &    ! bottom atm level zon wind       DEF
      //':Sa_v'        &    ! bottom atm level mer wind       DEF
      //':Sa_tbot'     &    ! bottom atm level temp           DEF
      //':Sa_ptem'     &    ! bottom atm level pot temp       DEF
      //':Sa_shum'     &    ! bottom atm level spec hum       DEF
      //':Sa_pbot'          ! bottom atm level pressure       DEF

   ! Fluxes
   character(*), parameter :: seq_flds_x2l_fluxes = &
         'Faxa_lwdn'   &    ! downward longwave heat flux     DEF
      //':Faxa_rainc'  &    ! precip: liquid, convective      DEF
      //':Faxa_rainl'  &    ! precip: liquid, large-scale     DEF
      //':Faxa_snowc'  &    ! precip: frozen, convective      DEF
      //':Faxa_snowl'  &    ! precip: frozen, large-scale     DEF
      //':Faxa_swndr'  &    ! shortwave: nir direct  down     DEF
      //':Faxa_swvdr'  &    ! shortwave: vis direct  down     DEF
      //':Faxa_swndf'  &    ! shortwave: nir diffuse down     DEF
      //':Faxa_swvdf'       ! shortwave: vis diffuse down     DEF
   ! States
   character(*), parameter :: seq_flds_x2i_states = &
         'So_t'        &    ! ocean layer temperature         DEF
      //':Sa_z'        &    ! atm bottom layer height         DEF
      //':Sa_u'        &    ! atm u velocity                  DEF
      //':Sa_v'        &    ! atm v velocity                  DEF
      //':Sa_ptem'     &    ! atm potential temp              DEF
      //':Sa_tbot'     &    ! atm bottom temp                 DEF
      //':Sa_pbot'     &    ! atm bottom pressure             DEF
      //':Sa_shum'     &    ! atm specfic humidity            DEF
      //':Sa_dens'          ! atm air density                 DEF
      
   ! Fluxes
   character(*), parameter :: seq_flds_x2i_fluxes = &
         'Fioo_q'      &    ! ocn freeze or melt heat         DEF
      //':Faxa_swndr'  &    ! atm sw near-ir, direct          DEF
      //':Faxa_swvdr'  &    ! atm sw visable, direct          DEF
      //':Faxa_swndf'  &    ! atm sw near-ir, diffuse         DEF
      //':Faxa_swvdf'  &    ! atm sw visable, diffuse         DEF
      //':Faxa_lwdn'   &    ! long-wave down                  DEF
      //':Faxa_rain'   &    ! prec: liquid                    DEF
      //':Faxa_snow'        ! prec: frozen                    DEF

   character(*), parameter :: seq_flds_a2x_fields = &
      trim(seq_flds_a2x_states)//":"//trim(seq_flds_a2x_fluxes)
   character(*), parameter :: seq_flds_x2a_fields = &
      trim(seq_flds_x2a_states)//":"//trim(seq_flds_x2a_fluxes)
   character(*), parameter :: seq_flds_x2l_fields = &
      trim(seq_flds_x2l_states)//":"//trim(seq_flds_x2l_fluxes)
   character(*), parameter :: seq_flds_x2i_fields = &
      trim(seq_flds_x2i_states)//":"//trim(seq_flds_x2i_fluxes)
   character(*), parameter :: seq_flds_i2x_fields = &
      trim(seq_flds_i2x_states)//":"//trim(seq_flds_i2x_fluxes)
   character(*), parameter :: seq_flds_l2x_fields = &
      trim(seq_flds_l2x_states)//":"//trim(seq_flds_l2x_fluxes)
   character(*), parameter :: seq_flds_o2x_fields = &
      trim(seq_flds_o2x_states)//":"//trim(seq_flds_o2x_fluxes)
   character(*), parameter :: seq_flds_xao_fields = &
      trim(seq_flds_xao_states)//":"//trim(seq_flds_xao_fluxes)
   integer, parameter :: lname_len = (SHR_KIND_CX * 5)
   character(len=lname_len) :: seq_flds_a2x_field_lnames = &
      trim(seq_flds_a2x_state_lnames)//":"//trim(seq_flds_a2x_flux_lnames)
   character(len=lname_len) :: seq_flds_a2x_field_units =  &
      trim(seq_flds_a2x_state_units)//":"//trim(seq_flds_a2x_flux_units)
   character(len=lname_len) :: blankString = " "
   integer, parameter :: nInitTests = 3
  character(len=80) :: FieldType(nInitTests)  = (/   &
                                          "simple      ",       &
                                          "lnames      ",       &
                                          "lnames&units" /)
  character(len=lname_len) :: FieldList(nInitTests)  = (/ &
                                     seq_flds_a2x_fields,       &
                                     seq_flds_a2x_fields,       &
                                     seq_flds_a2x_fields        &
                                   /)
  character(len=lname_len) :: FieldUnits(nInitTests)
  character(len=lname_len) :: FieldLNames(nInitTests)
  integer :: rc, nf
  integer :: i, MPIRank, j
  logical :: MasterTask
  type(ESMF_DELayout) :: deLayout
  type(ESMF_VM) :: vm
  type(ESMF_Grid) :: gridArbDecomp
  type(ESMF_Grid) :: gridArb2Decomp
  type(ESMF_Grid) :: FVgridArbDecomp
  type(ESMF_Grid) :: FVgridArb2Decomp
  integer, parameter :: nlats = 8
  integer, parameter :: nlons = 16
  real(ESMF_KIND_R8) :: coordLat(nlats+1)
  real(ESMF_KIND_R8) :: coordLon(nlons+1)
  real(r8), parameter :: centerLon(nlons) = &
      (/ 0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., &
         247.5, 270., 292.5, 315., 337.5 /)
  real(r8), parameter :: centerLat(nlats) = &
    (/ -73.7992136285632, -52.8129431899943, -31.704091745008, -10.5698823125761, &
       10.5698823125761, 31.704091745008, 52.8129431899943, 73.7992136285632 /)
  integer, parameter :: FVnlats = 19
  integer, parameter :: FVnlons = 24
  real(ESMF_KIND_R8) :: FVcoordLat(FVnlats+1)
  real(ESMF_KIND_R8) :: FVcoordLon(FVnlons+1)
  real(r8), parameter :: FVcenterLon(FVnlons) = &
      (/ 0., 15., 30., 45., 60., 75., 90., 105., 120., 135., 150., 165.,  &
         180., 195., 210., 225., 240., 255., 270., 285., 300., 315., 330., 345./)
  real(r8), parameter :: FVcenterLat(FVnlats) = &
    (/ -90., -80., -70., -60., -50., -40., -30., -20., -10., 0., 10., 20., 30.,   &
        40., 50., 60., 70., 80., 90. /)
  integer,  pointer    :: GlobalMask(:,:), FVGlobalMask(:,:)
  integer,  pointer    :: Mask(:), FVMask(:)
  integer              :: myCount, myCount1
  integer              :: DeComp(2)
  integer, allocatable :: myIndices(:,:)
  integer              :: npes, id
  integer              :: rank
  integer              :: MasterRank
  integer              :: mpicom
  type(ESMF_ArraySpec) :: arraySpec1D  ! 1D Array spec for arbitrary decomp
  type(ESMF_ArraySpec) :: arraySpec2D  ! 2D Array spec
  type(ESMF_Field)     :: FieldArb     ! field in the state
  type(ESMF_State)     :: eState       ! ESMF State
  type(ESMF_State)     :: eStateIn     ! ESMF State input state for copy or redist
  type(ESMF_State)     :: eStateIce, eStateLnd, eStateOcn, eStateXao
  type(ESMF_State)     :: eStateOut    ! ESMF State output state for copy or redist
  type(ESMF_State)     :: eStateExpect ! What expect merge state to look like

#include <mpif.h>


  call shr_mpi_init( )
  call shr_mpi_commrank( MPI_COMM_WORLD, MPIRank )
  MasterTask = (MPIRank == 0)
  MasterRank = 0
  rank = MPIRank
  mpicom = MPI_COMM_WORLD
  FieldUnits(1:nInitTests-1)           = blankString
  FieldUnits(nInitTests)               = seq_flds_a2x_field_units
  FieldLNames(1:nInitTests-2)          = blankString
  FieldLNames(nInitTests-1:nInitTests) = seq_flds_a2x_field_lnames

  if ( MasterTask ) write(6,*) "Initialize ESMF....."
  call eshr_rc_set( continueOnFailure = .true. )
  call ESMF_Initialize( rc=rc )
  call eshr_rc_check( rc, "Error on ESMF_Init" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )

  ! ----- Get the global VM and create a deLayout and Grid that can be used
  call ESMF_VMGetGlobal(vm, rc)
  call eshr_rc_check( rc, "Error on vm get global" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_VMGet(vm, petcount=npes, localpet=id, rc=rc)
  call eshr_rc_check( rc, "Error on vm get" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  deLayout = ESMF_DELayoutCreate( vm, rc=rc )
  call eshr_rc_check( rc, "Error on delayout create" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_DELayoutPrint( delayout, rc=rc )
  call eshr_rc_check( rc, "Error on delayout print" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  ! Spectral grids
  if ( MasterTask ) write(6,*) "center to vertex....."
  call GridnDecomp_CenterToVertices( centerLon, centerLat, coordLon, coordLat )
  call setRandomSeed( 1 )
  if ( MasterTask ) write(6,*) "calc decomp....."
  call GridnDecomp_CalcDecomp( MasterRank, npes, rank, mpicom, MyCount )
  if ( MasterTask ) write(6,*) "create grid....."
  gridArbDecomp = ESMF_GridCreateHorzLatLon( coordLon, coordLat,                     &
                                    horzstagger=ESMF_GRID_HORZ_STAGGER_A,   &
                                    dimNames=(/"longitude", "latitude "/),  &
                                    dimUnits=(/"degrees",   "degrees"  /),  &
                                    coordorder=ESMF_COORD_ORDER_XYZ,        &
                                    periodic=(/ESMF_TRUE, ESMF_FALSE/),     &
                                    name="Test grid arb decomp", rc=rc)
  call eshr_rc_check( rc, "Error on grid create" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  gridArb2Decomp = ESMF_GridCreateHorzLatLon( coordLon, coordLat,                     &
                                    horzstagger=ESMF_GRID_HORZ_STAGGER_A,   &
                                    dimNames=(/"longitude", "latitude "/),  &
                                    dimUnits=(/"degrees",   "degrees"  /),  &
                                    coordorder=ESMF_COORD_ORDER_XYZ,        &
                                    periodic=(/ESMF_TRUE, ESMF_FALSE/),     &
                                    name="Test grid arb decomp", rc=rc)
  call eshr_rc_check( rc, "Error on grid create" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  allocate( MyIndices(MyCount,2) )
  call GridnDecomp_GetLocalGIndices( MyCount, rank, MyIndices )
  call shr_mpi_barrier(mpicom)
  if ( MasterTask ) write(6,*) "Grid distribute....."
!  write(6,*) "id = ", id, " MyCount = ", MyCount
!  write(6,1080) "id = ", id, ( MyIndices(i,1), MyIndices(i,2), i = 1, MyCount )
!1080 format( A, i3, 128( "(",i3,",",i3,") " ) )
!  call shr_sys_flush(6)
  call ESMF_GridDistribute( gridArbDecomp, delayout, myCount=MyCount, &
                            myIndices=myIndices, rc=rc)
  call shr_mpi_barrier(mpicom)
  allocate( Mask(MyCount) )
  MyCount1 = MyCount
  call GridnDecomp_CenterToVertices( centerLon, centerLat, coordLon, coordLat )
  call setRandomSeed( 5 )
  call GridnDecomp_CalcDecomp( MasterRank, npes, rank, mpicom, MyCount, &
                               landGrid=.true., GlobalMask=GlobalMask )
  do i = 1, MyCount1
     Mask(i) = GlobalMask(MyIndices(i,1),MyIndices(i,2))
  end do
  deallocate( MyIndices )
  allocate( MyIndices(MyCount,2) )
  call GridnDecomp_GetLocalGIndices( MyCount, rank, MyIndices )
  call shr_mpi_barrier(mpicom)
  if ( MasterTask ) write(6,*) "Grid distribute....."
  call ESMF_GridDistribute( gridArb2Decomp, delayout, myCount=MyCount, &
                            myIndices=myIndices, rc=rc)
  deallocate( MyIndices )
  call ESMF_GridValidate( gridArbDecomp, rc=rc )
  call eshr_rc_check( rc, "Error on grid validate" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_GridPrint( gridArbDecomp, rc=rc )
  call eshr_rc_check( rc, "Error on grid print" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_GridValidate( gridArb2Decomp, rc=rc )
  call eshr_rc_check( rc, "Error on grid validate" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_GridPrint( gridArb2Decomp, rc=rc )
  call eshr_rc_check( rc, "Error on grid print" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  ! FV grid
  call shr_mpi_barrier(mpicom)
  call GridnDecomp_CenterToVertices( FVcenterLon, FVcenterLat, FVcoordLon, FVcoordLat )
  call setRandomSeed( 15 )
  call GridnDecomp_CalcDecomp( MasterRank, npes, rank, mpicom, MyCount )
  call eshr_rc_check( rc, "Error on grid create" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  FVgridArbDecomp = ESMF_GridCreateHorzLatLon( FVcoordLon, FVcoordLat,        &
                                    horzstagger=ESMF_GRID_HORZ_STAGGER_A,   &
                                    dimNames=(/"longitude", "latitude "/),  &
                                    dimUnits=(/"degrees",   "degrees"  /),  &
                                    coordorder=ESMF_COORD_ORDER_XYZ,        &
                                    periodic=(/ESMF_TRUE, ESMF_FALSE/),     &
                                    name="Test FV grid arb decomp", rc=rc)
  call eshr_rc_check( rc, "Error on grid create" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  FVgridArb2Decomp = ESMF_GridCreateHorzLatLon( FVcoordLon, FVcoordLat,        &
                                    horzstagger=ESMF_GRID_HORZ_STAGGER_A,   &
                                    dimNames=(/"longitude", "latitude "/),  &
                                    dimUnits=(/"degrees",   "degrees"  /),  &
                                    coordorder=ESMF_COORD_ORDER_XYZ,        &
                                    periodic=(/ESMF_TRUE, ESMF_FALSE/),     &
                                    name="Test FV grid arb decomp", rc=rc)
  call eshr_rc_check( rc, "Error on grid create" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  allocate( MyIndices(MyCount,2) )
  call GridnDecomp_GetLocalGIndices( MyCount, rank, MyIndices )
  call ESMF_GridDistribute( FVgridArbDecomp, delayout, myCount=MyCount, &
                            myIndices=myIndices, rc=rc)
  call shr_mpi_barrier(mpicom)
  call GridnDecomp_CenterToVertices( FVcenterLon, FVcenterLat, FVcoordLon, FVcoordLat )
  call setRandomSeed( 50 )
  MyCount1 = MyCount
  call GridnDecomp_CalcDecomp( MasterRank, npes, rank, mpicom, MyCount, &
                               landGrid=.true., GlobalMask=FVGlobalMask )
  allocate( FVMask(MyCount1) )
  do i = 1, MyCount1
     FVMask(i) = FVGlobalMask(MyIndices(i,1),MyIndices(i,2))
  end do
  deallocate( MyIndices )
  allocate( MyIndices(MyCount,2) )
  call GridnDecomp_GetLocalGIndices( MyCount, rank, MyIndices )
  if ( MasterTask ) write(6,*) "FV Grid distribute....."
  call ESMF_GridDistribute( FVgridArb2Decomp, delayout, myCount=MyCount, &
                            myIndices=myIndices, rc=rc)
  call shr_mpi_barrier(mpicom)
  call ESMF_GridValidate( FVgridArbDecomp, rc=rc )
  call eshr_rc_check( rc, "Error on grid validate" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_GridPrint( FVgridArbDecomp, rc=rc )
  call eshr_rc_check( rc, "Error on grid print" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  !
  ! Compare the two grids
  !
  if ( MasterTask ) write(6,*) "Compare spectral to FV grid"
  if ( eshr_estate_Egrids2DAreSame( gridArbDecomp, FVgridArbDecomp, exact=.true., &
       Debug=MasterTask ) )then
     write(6,*) "FV and Spectral grids compare same and should NOT"
     call shr_sys_abort( "FAIL" )
  end if
  if ( MasterTask ) write(6,*) "Compare spectral grids with different decomp"
  if ( .not. eshr_estate_Egrids2DAreSame( gridArbDecomp, gridArb2Decomp, exact=.true., &
  Debug=MasterTask, Grid1Mask=Mask ) )then
     write(6,*) "Arbitrary decomp 1 and 2 Spectral grids compare different"
     call shr_sys_abort( "FAIL" )
  end if
  if ( MasterTask ) write(6,*) "Compare FV grids with different decomp"
  if ( .not. eshr_estate_Egrids2DAreSame( FVgridArbDecomp, FVgridArb2Decomp, &
  exact=.true., Debug=MasterTask, Grid1Mask=FVMask ) )then
     write(6,*) "Arbitrary decomp FV grids 1 and 2 compare different"
     call shr_sys_abort( "FAIL" )
  end if
  if ( MasterTask ) write(6,*) "Compare spectral grids with different decomp, non-exact"
  if ( .not. eshr_estate_Egrids2DAreSame( gridArbDecomp, gridArb2Decomp, exact=.false., &
  Debug=MasterTask, Grid1Mask=Mask ) )then
     write(6,*) "Arbitrary decomp 1 and 2 Spectral grids compare different non exact"
     call shr_sys_abort( "FAIL" )
  end if
  if ( MasterTask ) write(6,*) "Compare FV grids with different decomp, non-exact"
  if ( .not. eshr_estate_Egrids2DAreSame( FVgridArbDecomp, FVgridArb2Decomp, &
  exact=.false., Debug=MasterTask, Grid1Mask=FVMask ) )then
     write(6,*) "Arbitrary decomp FV grids 1 and 2 compare different non exact"
     call shr_sys_abort( "FAIL" )
  end if
  !
  ! Simple field and state create with above grid
  !
  eState = ESMF_StateCreate( "Temp state", statetype=ESMF_STATE_IMPORT, rc=rc )
  call eshr_rc_check( rc, "Error in creation of ESMF State" )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_ArraySpecSet( arraySpec1D, rank=1, type=ESMF_DATA_REAL, kind=ESMF_R8, &
                          rc=rc )
  call eshr_rc_check( rc, 'ERROR: setting 1D arb array spec' )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  FieldArb = ESMF_FieldCreate(gridArbDecomp, arraySpec1D, allocFlag=ESMF_ALLOC, &
                           horzRelloc=ESMF_CELL_CENTER,                &
                           name="2D_FieldTest_arb_decomp", rc=rc)
  call ESMF_ArraySpecSet( arraySpec2D, rank=2, type=ESMF_DATA_REAL, kind=ESMF_R8, &
                          rc=rc )
  call eshr_rc_check( rc, 'ERROR: setting 2D array spec' )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  call ESMF_StateAddField( eState, fieldArb, rc=rc )
  call eshr_rc_check( rc, 'ERROR: add field to state' )
  if ( rc /= ESMF_SUCCESS ) call shr_sys_abort( "FAIL" )
  !
  ! Loop through estateInit tests that should work for arbitrary decomp
  !
  do i = 1, nInitTests
     if ( MasterTask ) write(6,no_rtnFMT) "EstateInit arbitrary decomp test ", &
                                           trim(FieldType(i))
     if ( .not. check_EstateInit( gridArbDecomp, FieldList(i), FieldLNames(i), &
                                  FieldUnits(i), arbitrary=.true. ) )then
        write(6,*) "EstateInit NOT successfully initialized on list"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
     if ( MasterTask ) write(6,no_rtnFMT) "EstateInit FV arbitrary decomp test ", &
                                           trim(FieldType(i))
     if ( .not. check_EstateInit( FVgridArbDecomp, FieldList(i), FieldLNames(i), &
                                  FieldUnits(i), arbitrary=.true. ) )then
        write(6,*) "EstateInit NOT successfully initialized on list"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
  end do
  if ( .not. check_EstateImportCreate( seq_flds_a2x_fields, seq_flds_i2x_fields, &
                          seq_flds_l2x_fields, seq_flds_o2x_fields, &
                          gridArbDecomp, gridArbDecomp, gridArbDecomp, gridArbDecomp, &
                          Mask) )then
     write(6,*) "EstateImportCreate NOT successfully ran on identical grids"
     call shr_sys_abort( "FAIL" )
  end if
  if (       check_EstateImportCreate( seq_flds_a2x_fields, seq_flds_i2x_fields, &
                      seq_flds_l2x_fields, seq_flds_o2x_fields, &
                      gridArbDecomp, FVgridArbDecomp, gridArbDecomp, FVgridArbDecomp, &
                      Mask) )then
     write(6,*) "EstateImportCreate ran successfully on different grids"
     call shr_sys_abort( "FAIL" )
  end if
  if ( .not. check_EstateImportCreate( seq_flds_a2x_fields, seq_flds_i2x_fields, &
                      seq_flds_l2x_fields, seq_flds_o2x_fields, &
                      gridArbDecomp, gridArbDecomp, gridArb2Decomp, gridArbDecomp, &
                      Mask) )then
     write(6,*) "EstateImportCreate NOT successfully ran on identical grids with different land decomp"
     call shr_sys_abort( "FAIL" )
  end if
  !
  ! Loop through Estate-copy tests that should work
  !
  do i = 1, 3
     eStateIn = ESMF_StateCreate( "Temp input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateOut = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_a2x_fields, &
                                       eState=eStateIn )
     nf = shr_string_listGetNum( seq_flds_a2x_fields )
     if ( i < 3 )then
       call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_a2x_fields, &
                                         eState=eStateOut )
     end if
     if ( i == 1 )then
        name = "simple-random fields"
        call fill_EstateField( eStateIn, arb=.true. )
     else if ( i == 2 )then
        name = "simple-constant fields"
        call fill_EstateField( eStateIn, arb=.true., type=CONSTANT_FILL )
     else if ( i == 3 )then
        name = "different-lists random fields"
        call fill_EstateField( eStateIn, arb=.true. )
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_x2l_fields, &
                                          eState=eStateOut )
        call fill_EstateField( eStateOut, arb=.true. )
        nf = shr_string_listGetNum( seq_flds_x2l_fields )
     else if ( i == 4 )then
        name = "constant fluxes not needed"
        call fill_EstateField( eStateIn, arb=.true., type=CONSTANT_FILL )
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_a2x_fields, &
                                          notNeeded=seq_flds_a2x_fluxes,                &
                                          eState=eStateOut )
        call fill_EstateField( eStateOut, arb=.true., type=CONSTANT_FILL )
        nf = shr_string_listGetNum( seq_flds_a2x_states )
     else if ( i == 5 )then
        name = "different-lists that share random fields"
        call fill_EstateField( eStateIn, arb=.true. )
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_x2l_fields, &
                                          eState=eStateOut )
        call fill_EstateField( eStateOut, arb=.true. )
        nf = shr_string_listGetNum( seq_flds_x2l_fields )
     else
        write(6,*) "Index out of bounds"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,no_rtnFMT) "Estate-copy tests ", trim(name)
     call eshr_estate_getStats( eStateIn, mpicom=mpicom, rootid=0  )
     if ( .not. check_EstateCopy( estateIn, eStateOut, nf, arb=.true. ) )then
        write(6,*) "EstateImportCopy NOT successfully ran"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
     call eshr_estate_getStats( eStateOut )
     ! Clean-up 
     call eshr_estate_destroy( eStateIn )
     call eshr_estate_destroy( eStateOut )
  end do
  !
  ! Loop through Estate-copy tests that should NOT work
  !
  do i = 1, 3
     eStateIn = ESMF_StateCreate( "Temp input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateOut = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     if ( i == 1 )then
        name = "mismatched field names"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_a2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_i2x_fields, &
                                          eState=eStateOut )
        nf = 0
     else if ( i == 2 )then
        cycle     ! This aborts on error when run
        if ( npes == 1 ) cycle
        name = "mismatched decomposition"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_i2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  GridArb2Decomp, FieldList=seq_flds_i2x_fields, &
                                          eState=eStateOut )
        nf = shr_string_listGetNum( seq_flds_i2x_fields )
     else if ( i == 3 )then
        name = "mismatched grids"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_l2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  FVGridArbDecomp, FieldList=seq_flds_l2x_fields, &
                                          eState=eStateOut )
        nf = shr_string_listGetNum( seq_flds_i2x_fields )
     else if ( i == 4 )then
        name = "constant no fields needed"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_l2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_l2x_fields, &
                                          notNeeded=seq_flds_l2x_fields,                &
                                          eState=eStateOut )
        nf = 0
     else
        write(6,*) "Index out of bounds"
        call shr_sys_abort( "FAIL" )
     end if
     call fill_EstateField( eStateIn, arb=.true. )
     if ( MasterTask ) write(6,no_rtnFMT) "Estate-copy fail tests ", trim(name)
     if ( check_EstateCopy( estateIn, eStateOut, nf, arb=.true. ) )then
        write(6,*) "EstateImportCopy successfully ran and should NOT have been"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
     ! Clean-up 
     call eshr_estate_destroy( eStateIn  )
     call eshr_estate_destroy( eStateOut )
  end do
  !
  ! Loop through Estate-redist tests that should work
  !
  do i = 1, 2
     eStateIn = ESMF_StateCreate( "Temp input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateOut = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     if ( i == 1)then
        name = "different grid decomp"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_a2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  GridArb2Decomp, FieldList=seq_flds_a2x_fields, &
                                          eState=eStateOut )
     else if ( i == 2 )then
        name = "same grids"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_o2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_o2x_fields, &
                                          eState=eStateOut )
     else if ( i == 3)then
        name = "different decomp some not needed"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_a2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  GridArb2Decomp, FieldList=seq_flds_a2x_fields, &
                                          eState=eStateOut, notNeeded=seq_flds_a2x_states  )
     else
        write(6,*) "Index out of bounds"
        call shr_sys_abort( "FAIL" )
     end if
     call fill_EstateField( eStateIn, arb=.true. )
     if ( MasterTask ) write(6,no_rtnFMT) "Estate-redist tests ", trim(name)
     if ( .not. check_EstateRedist( estateIn, eStateOut, Mask ) )then
        write(6,*) "EstateImportRedist NOT successfully ran"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
     ! Clean-up 
     call eshr_estate_destroy( eStateIn  )
     call eshr_estate_destroy( eStateOut )
  end do
  !
  ! Loop through Estate-redist tests that should NOT work
  !
  do i = 1, 2
     eStateIn = ESMF_StateCreate( "Temp input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateOut = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     if ( i == 1 )then
        name = "mismatched field names"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_l2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  GridArb2Decomp, FieldList=seq_flds_i2x_fields, &
                                          eState=eStateOut )
     else if ( i == 2 )then
        name = "mismatched grids"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_o2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  FVGridArbDecomp, FieldList=seq_flds_o2x_fields, &
                                          eState=eStateOut )
     else if ( i == 3 )then
        name = "all fields not needed"
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_o2x_fields, &
                                          eState=eStateIn )
        call eshr_estate_init2DFromList(  GridArb2Decomp, FieldList=seq_flds_o2x_fields, &
                                          notNeeded=seq_flds_o2x_fields, estate=eStateOut )
     else
        write(6,*) "Index out of bounds"
        call shr_sys_abort( "FAIL" )
     end if
     call fill_EstateField( eStateIn, arb=.true. )
     if ( MasterTask ) write(6,no_rtnFMT) "Estate-redist fail tests ", trim(name)
     if ( check_EstateRedist( estateIn, eStateOut, Mask ) )then
        write(6,*) "EstateImportRedist successfully ran and should NOT have been"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
     ! Clean-up 
     call eshr_estate_destroy( eStateIn  )
     call eshr_estate_destroy( eStateOut )
  end do
  !
  ! Loop through Estate-merge tests that should work
  !
  do i = 1, 2
     eStateIce = ESMF_StateCreate( "Temp ice input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateLnd = ESMF_StateCreate( "Temp land input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateOcn = ESMF_StateCreate( "Temp ocean input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateOut = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateXao = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateExpect = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_x2a_fields, &
                                       eState=eStateOut )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_x2a_mrg_fields, &
                                       eState=eStateExpect )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_i2x_fields, &
                                       eState=eStateIce )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_o2x_fields, &
                                       eState=eStateOcn )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_xao_fields, &
                                       eState=eStateXao )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_l2x_fields, &
                                       eState=eStateLnd )
     if ( i == 1 )then
        name = "simple avg"
        call fill_EstateField( eStateIce, arb=.true. )
        call fill_EstateField( eStateOcn, arb=.true. )
        call fill_EstateField( eStateLnd, arb=.true. )
        call fill_EstateField( eStateOut, arb=.true. )
        call fill_EstateField( eStateXao, arb=.true. )
     else if ( i == 2 )then
        name = "constant avg"
        call fill_EstateField( eStateIce, arb=.true., type=CONSTANT_FILL )
        call fill_EstateField( eStateOcn, arb=.true., type=CONSTANT_FILL )
        call fill_EstateField( eStateLnd, arb=.true., type=CONSTANT_FILL )
        call fill_EstateField( eStateOut, arb=.true., type=CONSTANT_FILL )
        call fill_EstateField( eStateXao, arb=.true., type=CONSTANT_FILL )
     else if ( i == 3 )then
        name = "simple avg -- some not needed"
        call fill_EstateField( eStateIce, arb=.true. )
        call fill_EstateField( eStateOcn, arb=.true. )
        call fill_EstateField( eStateLnd, arb=.true. )
        call fill_EstateField( eStateOut, arb=.true. )
        call fill_EstateField( eStateXao, arb=.true. )
        call eshr_estate_getFieldNames( estateOut, Fieldnames, nfields )
        nfields = nfields/2
        call eshr_estate_markRestNotNeed( estateOut, Fieldnames(:nfields), nFields )
        deallocate( Fieldnames )
     else
        write(6,*) "Index out of bounds"
        call shr_sys_abort( "FAIL" )
     end if
     call prepare_MergeExpect( Mask, eStateIce, eStateOcn, eStateLnd, eStateOut, &
                               eStateXao, eStateExpect )
     if ( MasterTask ) write(6,no_rtnFMT) "Estate-merge tests ", trim(name)
     if ( .not. check_EstateMerge( estateIce, estateLnd, estateOcn, eStateOut, &
                                   estateXao, estateExpect ) )then
        write(6,*) "EstateImportMerge NOT successfully ran"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
     call eshr_estate_getStats( eStateOut )
     ! Clean-up 
     call eshr_estate_destroy( eStateIce  )
     call eshr_estate_destroy( eStateLnd  )
     call eshr_estate_destroy( eStateOcn  )
     call eshr_estate_destroy( eStateOut )
     call eshr_estate_destroy( eStateXao )
  end do
  !
  ! Loop through Estate-merge tests that should NOT work
  !
  do i = 1, 2
     eStateIce = ESMF_StateCreate( "Temp ice input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateLnd = ESMF_StateCreate( "Temp land input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateOcn = ESMF_StateCreate( "Temp ocean input state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateOut = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateXao = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     eStateExpect = ESMF_StateCreate( "Temp output state", rc=rc )
     call eshr_rc_check( rc, "Error in creation of ESMF State" )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_x2a_fields, &
                                       eState=eStateOut )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_x2a_mrg_fields, &
                                       eState=eStateExpect )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_i2x_fields, &
                                       eState=eStateIce )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_o2x_fields, &
                                       eState=eStateOcn )
     call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_xao_fields, &
                                       eState=eStateXao )
     if ( i == 1 )then
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_a2x_fields, &
                                          eState=eStateLnd )
        name = "no lnd flds"
     else if ( i == 2 )then
        call eshr_estate_init2DFromList(  FVGridArb2Decomp, FieldList=seq_flds_l2x_fields, &
                                          eState=eStateLnd )
        name = "mismatched grids"
     else if ( i == 3 )then
        call eshr_estate_init2DFromList(  GridArb2Decomp, FieldList=seq_flds_l2x_fields, &
                                          eState=eStateLnd )
        call eshr_estate_destroy( eStateOut )
        eStateOut = ESMF_StateCreate( "Temp output state", rc=rc )
        call eshr_rc_check( rc, "Error in creation of ESMF State" )
        call eshr_estate_init2DFromList(  GridArbDecomp, FieldList=seq_flds_x2a_fields, &
                                          notNeeded=seq_flds_x2a_fields, eState=eStateOut )
        name = "All out flds not needed"
     else
        write(6,*) "Index out of bounds"
        call shr_sys_abort( "FAIL" )
     end if
     call fill_EstateField( eStateIce, arb=.true. )
     call fill_EstateField( eStateOcn, arb=.true. )
     call fill_EstateField( eStateLnd, arb=.true. )
     call fill_EstateField( eStateOut, arb=.true. )
     call fill_EstateField( eStateXao, arb=.true. )
     call prepare_MergeExpect( Mask, eStateIce, eStateOcn, eStateLnd, eStateOut, &
                               eStateXao, eStateExpect )
     if ( MasterTask ) write(6,no_rtnFMT) "Estate-merge fail tests ", trim(name)
     if ( check_EstateMerge( estateIce, estateLnd, estateOcn, eStateOut, &
                             estateXao, estateExpect ) )then
        write(6,*) "EstateMerge successfully ran and should NOT have been"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
     call eshr_estate_getStats( eStateOut )
     ! Clean-up 
     call eshr_estate_destroy( eStateIce  )
     call eshr_estate_destroy( eStateLnd  )
     call eshr_estate_destroy( eStateOcn  )
     call eshr_estate_destroy( eStateOut )
     call eshr_estate_destroy( eStateXao )
  end do
  !
  ! Loop through inputInfo tests that should work
  !
  do i = 1, 6
     call shr_inputInfo_initSetDefaults( CCSMInit )
     if (      i == 1 )then
        name = "simple"
        call shr_inputinfo_initPutData( CCSMInit, start_type="startup",      &
                                  case_name="case", case_desc="Description", &
                                  mss_wpass="password",                      &
                                  archive_dir="/CCSM/csm/case/drv" )
     else if ( i == 2 )then
        name = "default"
     else if ( i == 3 )then
        name = "long_names"
        call shr_inputinfo_initPutData( CCSMInit, start_type="branch",             &
            case_name="long_casename123456789",                                    &
            case_desc="Description of a very long case name and long description", &
            mss_wpass="p", mss_irt=0, atm_ideal_phys=.true.,                       &
            brnch_retain_casename=.true.,                                          &
            archive_dir="mss:/CCSM_LONGPATHCASETESTING/csm/long_casename123456789/drv" )
     else if ( i == 4 )then
        name = "aqua"
        call shr_inputinfo_initPutData( CCSMInit, start_type="continue",           &
            mss_irt=4292, aqua_planet=.true.,                                      &
            archive_dir="null:" )
     else if ( i == 5 )then
        name = "ideal"
        call shr_inputinfo_initPutData( CCSMInit, start_type="continue",           &
            mss_irt=1, atm_ideal_phys=.true., mss_wpass=" ",                       &
            archive_dir="file:/home/user/test/drv" )
     else if ( i == 6 )then
        name = "restart_file_misc"
        call shr_inputinfo_initPutData( CCSMInit, start_type="branch",             &
        restart_file="restart_file_drv", logFilePostFix=".ccsm.log",               &
        outPathRoot="/ptmp/user/ccsmtest/ccsm",                                    &
        archive_dir="mss:/CCSM_LONGPATHCASETESTING/csm/long_casename123456789/drv" )
     else
        write(6,*) "Index beyond list using: ", i
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,no_rtnFMT) "CCSMInit test ", trim(name)
     if ( .not. check_ccsminitConversion( CCSMInit ) )then
        write(6,*) "CCSMInit NOT successfully changed to EState and back"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
  end do
  !
  ! Error test -- these variables won't be duplicated, as such these tests will fail
  !
  do i = 1, 2
     call shr_inputInfo_initSetDefaults( CCSMInit )
     if ( i == 1 )then
        name = "restart_pfile"
        call shr_inputinfo_initPutData( CCSMInit, start_type="branch",         &
        restart_pfile = "thing",                                               &
        archive_dir=" " )
     else if ( i == 2 )then
        name = "restart_file_override"
        call shr_inputinfo_initPutData( CCSMInit, start_type="branch",         &
        mss_wpass="passw", mss_irt=32,                                         &
        restart_file_override="mss_wpass:mss_irt",                             &
        archive_dir="mss:/CCSM_LONGPATHCASETESTING/csm/long_casename123456789/drv" )
     else
        write(6,*) "index out of bounds", i
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,no_rtnFMT) "CCSMInit Error test", trim(name)
     if ( check_ccsminitConversion( CCSMInit ) )then
        write(6,*) "CCSMInit successfully changed to EState and back -- on variables not expected"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
  end do
  !
  ! Filename tests
  !
  do i = 1, 3
     if ( i == 1 )then
        name = "simple"
        Filetype="NLFilename"
        Filename="namelist"
     else if ( i == 2 )then
        name = "long"
        Filetype="NLFilename"
        Filename="/home/user/longpath/namelist"
     else if ( i == 3 )then
        name = "blank"
        Filetype="configuration_file_for_sequential_ccsm"
        Filename=" "
     else
        write(6,*) "index out of bounds", i
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,no_rtnFMT) "Filename test", trim(name)
     if ( .not. check_FilenameConversion( FileType, Filename ) ) then
        write(6,*) "Filename NOT successfully changed to EState and back"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
  end do
  !
  ! Filename failure tests
  !
  do i = 1, 1
     if ( i == 1 )then
        name = "blank_filetype"
        Filetype=" "
        Filename="namelist"
     else
        write(6,*) "index out of bounds", i
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,no_rtnFMT) "Filename error test", trim(name)
     if ( check_FilenameConversion( FileType, Filename ) ) then
        write(6,*) "Filename successfully changed to EState and back -- and should not be"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,no_rtnFMT) "Filename error test", trim(name)
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
  end do
  !
  ! Clockinfo time-manager tests
  !
  do i = 1, 5
     call eshr_timemgr_NMLinfoSetDefault( clockNMLinfo=clockNMLinfo )
     call eshr_timemgr_NMLinfoPutData( clockNMLinfo, start_ymd=20060901,      &
                                       stop_option="date", stop_ymd=20060910, &
                                       restart_option="end", dtime=1200,      &
                                       orb_iyear_AD=2006 )
     if (      i == 1 )then
        name = "default"
     else if ( i == 2 )then
        name = "simple"
        call eshr_timemgr_NMLinfoPutData( clockNMLinfo, desc="My clock", &
                                          restart_option="none" )
     else if ( i == 3 )then
        name = "orbit"
        call eshr_timemgr_NMLinfoPutData( clockNMLinfo, desc="Change Orb",     &
                                          orb_iyear_AD=SHR_ORB_UNDEF_INT,      &
                                          orb_eccen=0.1_r8, orb_mvelp=278._r8, &
                                          orb_obliq=27.9_r8 )
     else if ( i == 4 )then
        name = "perpetual"
        call eshr_timemgr_NMLinfoPutData( clockNMLinfo, desc="Change Perpetual", &
                                          perpetual_ymd=20060321,                &
                                          perpetual_run=.true. )
     else if ( i == 5 )then
        name = "calendar"
        call eshr_timemgr_NMLinfoPutData( clockNMLinfo, desc=" ", &
                                          calendar="GREGORIAN" )
     else
        write(6,*) "Index beyond list using: ", i
        call shr_sys_abort( "FAIL" )
     end if
     call eshr_timemgr_clockInitNMLinfo( clockNMLinfo, LogPrint=.false., clockOut=clock )
     if ( MasterTask ) write(6,no_rtnFMT) "ClockInfo test ", trim(name)
     call eshr_timemgr_clockGet( clock, Info=clockInfo )
     if ( .not. check_clockInfoConversion( ClockInfo ) )then
        write(6,*) "ClockInfo NOT successfully changed to EState and back"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )
  end do
  !
  ! Clockinfo time-manager failure tests
  !
  do i = 1, 1
     call eshr_timemgr_NMLinfoSetDefault( clockNMLinfo=clockNMLinfo )
     call eshr_timemgr_NMLinfoPutData( clockNMLinfo, start_ymd=20060901,      &
                                       stop_option="date", stop_ymd=20060910, &
                                       restart_option="end", dtime=1200,      &
                                       orb_iyear_AD=2006 )
     if (      i == 1 )then
        name = "MasterSyncClock"
        call eshr_timemgr_NMLinfoPutData( clockNMLinfo, desc="My clock", &
                                          MasterSyncClock=.false. )
     else
        write(6,*) "Index beyond list using: ", i
        call shr_sys_abort( "FAIL" )
     end if
     call eshr_timemgr_clockInitNMLinfo( clockNMLinfo, LogPrint=.false., clockOut=clock )
     if ( MasterTask ) write(6,no_rtnFMT) "ClockInfo failure test ", trim(name)
     call eshr_timemgr_clockGet( clock, Info=clockInfo )
     if ( check_clockInfoConversion( ClockInfo ) )then
        write(6,*) "ClockInfo successfully changed to EState and back when should NOT have"
        call shr_sys_abort( "FAIL" )
     end if
     if ( MasterTask ) write(6,'(A)') "PASS"
     if ( MasterTask ) call shr_sys_flush( 6 )

  end do

  call ESMF_Finalize( terminationflag=ESMF_KEEPMPI )
  if ( MasterTask ) write(6,'(///,A)') "PASS"

  call shr_mpi_finalize( )

end program test_eshr_estate
