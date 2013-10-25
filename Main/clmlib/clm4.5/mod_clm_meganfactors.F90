module mod_clm_meganfactors
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: MEGANFactorsMod
!
! !DESCRIPTION:
! Manages input of MEGAN emissions factors from netCDF file
!
! !USES:
  use mod_realkinds
  use mod_stdio
  use mod_clm_nchelper
!
  implicit none
  private
  save
!
! !PUBLIC MEMBERS:
  public :: megan_factors_init
  public :: megan_factors_get
  public :: comp_names
!
! !PUBLIC DATA:
  real(rk8), public, allocatable :: LDF(:)  ! light dependent fraction
  real(rk8), public, allocatable :: Agro(:) ! growing leaf age factor
  real(rk8), public, allocatable :: Amat(:) ! mature leaf age factor
  real(rk8), public, allocatable :: Anew(:) ! new leaf age factor
  real(rk8), public, allocatable :: Aold(:) ! old leaf age factor
  real(rk8), public, allocatable :: betaT(:)! temperature factor
  real(rk8), public, allocatable :: ct1(:)  ! temperature coefficient 1
  real(rk8), public, allocatable :: ct2(:)  ! temperature coefficient 2
  real(rk8), public, allocatable :: Ceo(:)  ! Eopt coefficient
!
! !PRIVATE MEMBERS:
  integer :: npfts ! number of plant function types
!
  type emis_eff_t
     real(rk8), pointer :: eff(:) ! emissions efficiency factor
     real(rk8) :: wght            ! molecular weight
     integer :: class_num        ! MEGAN class number
  endtype emis_eff_t
!
  type(emis_eff_t), pointer :: comp_factors_table(:)  ! hash table of MEGAN factors (points to an array of pointers)
  integer, pointer :: hash_table_indices(:)           ! pointer to hash table indices
  integer, parameter :: tbl_hash_sz = 2**16           ! hash table size
!
  character(len=32), allocatable :: comp_names(:)     ! MEGAN compound names
  real(rk8),          allocatable :: comp_molecwghts(:)! MEGAN compound molecular weights
!
! !REVISION HISTORY:
!  28 Oct 2011: Created by Francis Vitt
!
!EOP
!-----------------------------------------------------------------------
contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: megan_factors_get
!
! !INTERFACE:
  subroutine megan_factors_get( comp_name, factors, class_n, molecwght )
!
! !DESCRIPTION:
! Method for getting MEGAN information for a named compound 
!
! !ARGUMENTS:
    character(len=*),intent(in)  :: comp_name      ! MEGAN compound name
    real(rk8),        intent(out) :: factors(npfts) ! vegitation type factors for the compound of intrest
    integer,         intent(out) :: class_n        ! MEGAN class number for the compound of intrest
    real(rk8),        intent(out) :: molecwght      ! molecular weight of the compound of intrest
!
!EOP
!-----------------------------------------------------------------------
! local vars:
    integer :: hashkey, ndx
    character(len=120) :: errmes

    hashkey = gen_hashkey(comp_name)
    ndx = hash_table_indices(hashkey)

    if (ndx<1) then 
       errmes = 'megan_factors_get: '//trim(comp_name)//' compound not found in MEGAN table'
       write(stderr,*) trim(errmes)
       call fatal(__FILE__,__LINE__,errmes)
    endif

    factors(:) = comp_factors_table( ndx )%eff(:)
    class_n    = comp_factors_table( ndx )%class_num
    molecwght  = comp_factors_table( ndx )%wght

  end subroutine megan_factors_get
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: megan_factors_init
!
! !INTERFACE:
  subroutine megan_factors_init( filename )
!
! !DESCRIPTION:
! Initializes the MEGAN factors using data from input file 
!
! !USES:
!
! !ARGUMENTS:
    character(len=*),intent(in) :: filename ! MEGAN factors input file

!EOP
!-----------------------------------------------------------------------
!
    type(clm_filetype) :: ncid             ! netcdf id

    integer :: start(2), count(2)

    integer :: ierr, i, vid
    integer :: n_comps, n_classes, n_pfts
    integer, allocatable :: class_nums(:)

    real(rk8),allocatable :: factors(:)
    real(rk8),allocatable :: comp_factors(:)
    real(rk8),allocatable :: class_factors(:)

    allocate(comp_factors_table(150))
    allocate(hash_table_indices(tbl_hash_sz))


    call clm_openfile(filename,ncid)

    call clm_inqdim(ncid,'Comp_Num',n_comps)
    call clm_inqdim(ncid,'Class_Num',n_classes)
    call clm_inqdim(ncid,'PFT_Num',n_pfts)

    npfts = n_pfts

    allocate( factors(n_pfts) )
    allocate( comp_factors(n_pfts) )
    allocate( class_factors(n_pfts) )

    allocate( comp_names(n_comps) )
    allocate( comp_molecwghts(n_comps) )
    allocate( class_nums(n_comps) )

    call clm_readvar(ncid,'Comp_Name',comp_names)
    call clm_readvar(ncid,'Comp_MW',comp_molecwghts)
    call clm_readvar(ncid,'Class_Num',class_nums)
 
    ! set up hash table where data is stored
    call  bld_hash_table_indices( comp_names )
    do i=1,n_comps
       call clm_readvar(ncid,'Class_EF',comp_factors,nt=i)
       call clm_readvar(ncid,'Comp_EF',class_factors,nt=i)
       factors(:) = comp_factors(:)*class_factors(:)
       call enter_hash_data( trim(comp_names(i)), factors, class_nums(i), comp_molecwghts(i) )
    enddo

    allocate( LDF(n_classes) )
    allocate( Agro(n_classes) )
    allocate( Amat(n_classes) )
    allocate( Anew(n_classes) )
    allocate( Aold(n_classes) )
    allocate( betaT(n_classes) )
    allocate( ct1(n_classes) )
    allocate( ct2(n_classes) )
    allocate( Ceo(n_classes) )

    call clm_readvar(ncid,'LDF',LDF)
    call clm_readvar(ncid,'Agro',Agro)
    call clm_readvar(ncid,'Amat',Amat)
    call clm_readvar(ncid,'Anew',Anew)
    call clm_readvar(ncid,'Aold',Aold)
    call clm_readvar(ncid,'betaT',betaT)
    call clm_readvar(ncid,'ct1',ct1)
    call clm_readvar(ncid,'ct2',ct2)
    call clm_readvar(ncid,'Ceo',Ceo)

    call clm_closefile(ncid)

    deallocate( class_nums, comp_factors,class_factors,factors )

  endsubroutine megan_factors_init
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Private methods...

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine bld_hash_table_indices( names )
    character(len=*),intent(in) :: names(:)

    integer :: n, i, hashkey

    hash_table_indices(:) = 0

    n = size(names)
    do i=1,n
       hashkey = gen_hashkey(names(i))
       hash_table_indices(hashkey) = i
    enddo

  endsubroutine bld_hash_table_indices

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine enter_hash_data( name, data, class_n, molec_wght )
    character(len=*), intent(in) :: name
    real(rk8), intent(in) :: data(:)
    integer,  intent(in) :: class_n
    real(rk8), intent(in) :: molec_wght

    integer :: hashkey, ndx
    integer :: nfactors

    hashkey = gen_hashkey(name)
    nfactors = size(data)

    ndx = hash_table_indices(hashkey)

    allocate (comp_factors_table(ndx)%eff(nfactors))

    comp_factors_table(ndx)%eff(:) = data(:)
    comp_factors_table(ndx)%class_num = class_n
    comp_factors_table(ndx)%wght = molec_wght

  end subroutine enter_hash_data

  !-----------------------------------------------------------------------
  !from cam_history
  !
  ! Purpose: Generate a hash key on the interval [0 .. tbl_hash_sz-1]
  !          given a character string.
  !
  ! Algorithm is a variant of perl's internal hashing function.
  !
  !-----------------------------------------------------------------------
  integer function gen_hashkey(string)

    implicit none
    !
    !  Arguments:
    !
    character(len=*), intent(in) :: string
    !
    !  Local vars
    !
    integer :: hash
    integer :: i

    integer, parameter :: tbl_max_idx = 15  ! 2**N - 1
    integer, parameter :: gen_hash_key_offset = z'000053db'
    integer, dimension(0:tbl_max_idx) :: tbl_gen_hash_key =  (/61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1/)

    hash = gen_hash_key_offset

    if ( len(string) /= 19 ) then
       !
       ! Process arbitrary string length.
       !
       do i = 1, len(string)
          hash = ieor(hash , (ichar(string(i:i)) * tbl_gen_hash_key(iand(i-1,tbl_max_idx))))
       end do
    else
       !
       ! Special case string length = 19
       !
       do i = 1, tbl_max_idx+1
          hash = ieor(hash , ichar(string(i:i))   * tbl_gen_hash_key(i-1)) 
       end do
       do i = tbl_max_idx+2, len(string)
          hash = ieor(hash , ichar(string(i:i))   * tbl_gen_hash_key(i-tbl_max_idx-2)) 
       end do
    end if

    gen_hashkey = iand(hash, tbl_hash_sz-1)

    return

  end function gen_hashkey

end module mod_clm_meganfactors
