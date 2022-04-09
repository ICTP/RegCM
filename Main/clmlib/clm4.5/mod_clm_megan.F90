!================================================================================
! Handles MEGAN VOC emissions metadata for CLM produced chemical emissions
! MEGAN = Model of Emissions of Gases and Aerosols from Nature
!
! This reads the megan_emis_nl namelist in drv_flds_in and makes the relavent
! information available to CAM, CLM, and driver. The driver sets up CLM to CAM
! communication for the  VOC flux fields. CLM needs to know what specific VOC
! fluxes need to be passed to the coupler and how to assimble the fluxes.
! CAM needs to know what specific VOC fluxes to expect from CLM.
!
! Francis Vitt -- 26 Oct 2011
!================================================================================
module mod_clm_megan

  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_dynparam
  use mod_mppparam
  use mod_dynparam
  use mod_stdio

  implicit none

  private

  save

  public :: shr_megan_readnl      ! reads megan_emis_nl namelist
  ! points to an array of chemical compounds (in CAM-Chem mechanism)
  ! that have MEGAN emissions
  public :: shr_megan_mechcomps
  ! number of unique compounds in the CAM chemical mechanism
  ! that have MEGAN emissions
  public :: shr_megan_mechcomps_n
  public :: shr_megan_megcomps_n   ! number of unique MEGAN compounds
  public :: shr_megan_megcomp_t    ! MEGAN compound data type
  ! data type for chemical compound in CAM mechanism than has MEGAN emissions
  public :: shr_megan_mechcomp_t
  ! points to linked list of shr_megan_comp_t objects
  public :: shr_megan_linkedlist
  ! switch to use mapped emission factors
  public :: shr_megan_mapped_emisfctrs
  public :: shr_megan_comp_ptr

  ! First drydep fields token
  character(len=80), public :: shr_megan_fields_token = ''
  character(len=256), public :: shr_megan_factors_file = ''

  integer(ik4) , parameter :: max_specifier_len = 1024

  ! MEGAN compound data structure (or user defined type)
  type shr_megan_megcomp_t
    ! MEGAN compound name (in MEGAN input table)
    character(len=16) :: name
    integer(ik4) :: index
    ! function of plant-function-type (PFT)
    real(rk8), pointer :: emis_factors(:)
    integer(ik4) :: class_number    ! MEGAN class number
    ! molecular weight of the MEGAN compound (g/mole)
    real(rk8) :: molec_weight
    ! points to next member in the linked list
    type(shr_megan_megcomp_t) , pointer :: next_megcomp
  endtype shr_megan_megcomp_t

  type shr_megan_comp_ptr
    type(shr_megan_megcomp_t), pointer :: ptr
  endtype shr_megan_comp_ptr

  ! chemical compound in CAM mechanism than has MEGAN emissions
  type shr_megan_mechcomp_t
    character(len=16) :: name           ! compound name
    ! an array of pointers to megan emis compounds
    type(shr_megan_comp_ptr), pointer :: megan_comps(:)
    ! number of megan emis compounds than make up the emissions for
    ! this mechanis compound
    integer(ik4) :: n_megan_comps
  end type shr_megan_mechcomp_t

  ! array of chemical compounds (in CAM mechanism) than have MEGAN emissions
  type(shr_megan_mechcomp_t), pointer :: shr_megan_mechcomps(:)
  ! points to linked list top
  type(shr_megan_megcomp_t),  pointer :: shr_megan_linkedlist

  ! number of unique megan compounds
  integer(ik4) :: shr_megan_megcomps_n  = 0
  ! number of unique compounds in the CAM chemical mechanism than have
  ! MEGAN emissions
  integer(ik4) :: shr_megan_mechcomps_n = 0

  ! switch to use mapped emission factors
  logical :: shr_megan_mapped_emisfctrs = .false.

  ! private data
  type parser_items_t
    character(len=16),pointer :: megan_comp_names(:)
    character(len=16) :: mech_comp_name
    integer(ik4) :: n_megan_comps
  end type parser_items_t

  contains

  !-------------------------------------------------------------------------
  !
  ! This reads the megan_emis_nl namelist group in drv_flds_in and parses the
  ! namelist information for the driver, CLM, and CAM.
  !
  ! Namelist variables:
  !   megan_specifier, shr_megan_mapped_emisfctrs, megan_factors_file
  !
  ! megan_specifier is a series of strings where each string contains one
  !  CAM chemistry constituent name (left of = sign) and one or more MEGAN
  !  compounds (seperated by + sign if more than one).  The specification of
  !  the MEGAN compounds to the right of the = signs tells the MEGAN VOC
  !  model within CLM how to construct the VOC fluxes using the factors in
  !  megan_factors_file and land surface state.
  !
  ! megan_factors_file read by CLM contains valid MEGAN compound names,
  !  MEGAN class groupings and scalar emission factors
  !
  ! shr_megan_mapped_emisfctrs switch is used to tell the MEGAN model to use
  !  mapped emission factors read in from the CLM surface data input file
  !  rather than the scalar factors from megan_factors_file
  !
  ! Example:
  ! &megan_emis_nl
  !  megan_specifier = 'ISOP = isoprene',
  !     'C10H16 = myrcene + sabinene + limonene + carene_3 +
  !               ocimene_t_b + pinene_b + ...',
  !     'CH3OH = methanol',
  !     'C2H5OH = ethanol',
  !     'CH2O = formaldehyde',
  !     'CH3CHO = acetaldehyde',
  ! ...
  !  megan_factors_file = '$datapath/megan_emis_factors.nc'
  ! /
  !-------------------------------------------------------------------------
!  subroutine shr_megan_readnl( NLFileName, megan_fields )
  subroutine shr_megan_readnl( NLFileName )
    implicit none
    character(len=*), intent(in)  :: NLFileName
!    character(len=*), intent(out) :: megan_fields
    character(len=256) ::       megan_fields
    integer(ik4) :: unitn            ! namelist unit number
    integer(ik4) :: i                ! Loop index
    integer(ik4) :: ierr             ! error code
    logical :: exists           ! if file exists or not
    integer(ik4), parameter :: maxspc = 150
    character(len=max_specifier_len) :: megan_specifier(maxspc) = ' '
    character(len=256) :: megan_factors_file = ' '
    character(*),parameter :: F00   = "('(seq_drydep_read) ',2a)"

    namelist /megan_emis_nl/ megan_specifier , megan_factors_file ,  &
                     shr_megan_mapped_emisfctrs , shr_megan_megcomps_n , &
                     shr_megan_mechcomps_n

    if ( myid == iocpu ) then
      inquire( file=trim(NLFileName), exist=exists)
      if ( exists ) then
        unitn = file_getUnit()
        open( unitn, file=trim(NLFilename), status='old' )
        if ( debug_level > 0 ) write(stdout,F00) &
             'Read in megan_emis_readnl namelist from: ', trim(NLFilename)
        read(unitn, megan_emis_nl, iostat=ierr)
        if (ierr > 0) then
          write(stderr,*) 'megan_emis_nl namelist read error'
          call fatal(__FILE__,__LINE__, &
             'problem on read of megan_emis_nl namelist in shr_megan_readnl' )
        end if
      end if
      call file_freeUnit( unitn )
    end if

    call bcast(shr_megan_megcomps_n)
    call bcast(shr_megan_mechcomps_n)
    call bcast(megan_factors_file,256)
    call bcast(shr_megan_mapped_emisfctrs)
    do i = 1 , maxspc
      call bcast(megan_specifier(i),max_specifier_len)
    end do

    shr_megan_factors_file = megan_factors_file

    if ( myid == iocpu ) then
      write(stdout,*) 'MEGAN FACTOR FILE : ',trim(megan_factors_file)
    end if
    ! parse the namelist info and initialize the module data
    call shr_megan_init( megan_specifier, megan_fields )

  end subroutine shr_megan_readnl

  !-------------------------------------------------------------------------
  ! module data initializer
  !-------------------------------------------------------------------------
  subroutine shr_megan_init( specifier, megan_fields )
    implicit none
    character(len=*), intent(in) :: specifier(:)
    character(len=*), intent(out) :: megan_fields
    integer(ik4) :: n_entries
    integer(ik4) :: i, j, k
    integer(ik4) :: spc_len
    type(parser_items_t), pointer :: items
    character(len=12) :: token   ! megan field name to add

    nullify(shr_megan_linkedlist)

    n_entries = size(specifier)
    allocate(shr_megan_mechcomps(n_entries))
    shr_megan_mechcomps(:)%n_megan_comps = 0
    megan_fields = ' '
    do i = 1,n_entries
      spc_len=len_trim(specifier(i))
      if ( spc_len > 0 ) then
        items => get_parser_items( specifier(i) )

        do k = 1 , shr_megan_mechcomps_n
          if ( trim(shr_megan_mechcomps(k)%name) == &
               trim(items%mech_comp_name) ) then
            call fatal(__FILE__,__LINE__, &
                'shr_megan_init : duplicate compound names : '// &
                trim(items%mech_comp_name))
          end if
        end do

        shr_megan_mechcomps(i)%name = items%mech_comp_name
        shr_megan_mechcomps(i)%n_megan_comps = items%n_megan_comps
        allocate(shr_megan_mechcomps(i)%megan_comps(items%n_megan_comps))

        do j = 1 , items%n_megan_comps
          shr_megan_mechcomps(i)%megan_comps(j)%ptr => &
                  add_megan_comp( items%megan_comp_names(j) )
        end do
        shr_megan_mechcomps_n = shr_megan_mechcomps_n

        call destroy_parser_items( items )

        ! Need to explicitly add Fl_ based on naming convention
        write(token,"('Fall_voc',i3.3)") shr_megan_mechcomps_n

        if ( shr_megan_mechcomps_n == 1 ) then
          ! no not prepend ":" to the string for the first token
          megan_fields = trim(token)
          shr_megan_fields_token = token
        else
          megan_fields = trim(megan_fields)//':'//trim(token)
        end if
      end if
    end do

  end subroutine shr_megan_init

  !-------------------------------------------------------------------------
  ! private methods...
  !-------------------------------------------------------------------------
  function get_parser_items( spec_entry ) result(items)
    implicit none
    character(len=*), intent(in) :: spec_entry
    type(parser_items_t), pointer :: items ! items returned
    integer(ik4) :: ndxs(512)
    integer(ik4) :: nelem, j, i
    character(len=max_specifier_len) :: tmp_str

    j = scan( spec_entry, '=' )

    nelem = 1
    ndxs(nelem) = j

    tmp_str = trim( spec_entry(j+1:) )
    j = scan( tmp_str, '+' )

    do while( j > 0 )
      nelem = nelem+1
      ndxs(nelem) = ndxs(nelem-1) + j
      tmp_str = tmp_str(j+1:)
      j = scan( tmp_str, '+' )
    end do
    ndxs(nelem+1) = len(spec_entry)+1
    allocate(items)
    allocate(items%megan_comp_names(nelem))
    items%mech_comp_name = trim(adjustl( spec_entry(:ndxs(1)-1)))
    items%n_megan_comps = nelem
    do i = 1 , nelem
       items%megan_comp_names(i) = &
               trim(adjustl( spec_entry(ndxs(i)+1:ndxs(i+1)-1)))
    end do
  end function get_parser_items

  subroutine destroy_parser_items( items )
    implicit none
    type(parser_items_t), pointer :: items

    deallocate( items%megan_comp_names )
    deallocate( items )
    nullify( items )
  end subroutine destroy_parser_items

  function add_megan_comp( name ) result(megan_comp)
    implicit none
    character(len=16), intent(in) :: name
    type(shr_megan_megcomp_t), pointer :: megan_comp

    megan_comp => get_megan_comp_by_name(shr_megan_linkedlist, name)
    if(associated(megan_comp)) then
      ! already in the list so return...
      return
    end if

    ! create new megan compound and add it to the list
    allocate(megan_comp)

    ! element%index = lookup_element( name )
    ! element%emis_factors = get_factors( list_elem%index )

    megan_comp%index = shr_megan_megcomps_n+1

    megan_comp%name = trim(name)
    nullify(megan_comp%next_megcomp)

    call add_megan_comp_to_list(megan_comp)

  end function add_megan_comp

  recursive function get_megan_comp_by_name(list_comp, name) result(megan_comp)
    implicit none
    type(shr_megan_megcomp_t), pointer  :: list_comp
    character(len=*), intent(in) :: name  ! variable name
    type(shr_megan_megcomp_t), pointer  :: megan_comp ! returned object

    if ( associated(list_comp) ) then
      if ( list_comp%name == name ) then
        megan_comp => list_comp
      else
        megan_comp => get_megan_comp_by_name(list_comp%next_megcomp, name)
      end if
    else
      nullify(megan_comp)
    end if
  end function get_megan_comp_by_name

  subroutine add_megan_comp_to_list( new_megan_comp )
    implicit none
    type(shr_megan_megcomp_t), target, intent(in) :: new_megan_comp
    type(shr_megan_megcomp_t), pointer :: list_comp

    if ( associated(shr_megan_linkedlist) ) then
      list_comp => shr_megan_linkedlist
      do while ( associated(list_comp%next_megcomp) )
        list_comp => list_comp%next_megcomp
      end do
      list_comp%next_megcomp => new_megan_comp
    else
      shr_megan_linkedlist => new_megan_comp
    end if
    shr_megan_megcomps_n = shr_megan_megcomps_n + 1
  end subroutine add_megan_comp_to_list

end module mod_clm_megan
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
