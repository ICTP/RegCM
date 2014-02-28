module mod_clm_varctl
  !
  ! Module containing run control variables
  !
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use mod_clm_varpar , only : maxpatch_pft , numpft
  use mod_mppparam
  use mod_runparams

  implicit none
  public :: set_clmvarctl   ! Set variables
  public :: clmvarctl_init  ! Initialize and check values after namelist input

  private

  integer(ik4) , parameter , private :: iundef = -9999999
  integer(ik4) , parameter , private :: rundef = -9999999.D0
  !
  ! Run control variables
  !
  character(len=256) , public :: caseid  = ' ' ! case id
  character(len=256) , public :: ctitle  = ' ' ! case title
  integer(ik4) , public :: nsrest = iundef   ! Type of run
  ! Startup from initial conditions
  integer(ik4) , public , parameter :: nsrStartup  = 0
  ! Continue from restart files
  integer(ik4) , public , parameter :: nsrContinue = 1
  ! true => no valid land points -- do NOT run
  logical , public :: noland = .false.
  ! Hostname of machine running on
  character(len=256) , public :: hostname = ' '
  ! username of user running program
  character(len=256) , public :: username = ' '
  ! description of this source
  character(len=256) , public :: source   = "Community Land Model CLM4.5"
  ! version of program
  character(len=256) , public :: version  = SVN_REV
  ! dataset conventions
  character(len=256) , public :: conventions = "CF-1.0"
  !
  ! Unit Numbers
  !
  ! Output NetCDF files
  !
  ! large file support for output NetCDF files
  logical , public :: outnc_large_files = .true.
  !
  ! Run input files
  !
  ! initial conditions file name
  character(len=256) , public :: finidat    = ' '
  ! surface data file name
  character(len=256) , public :: fsurdat    = ' '
  ! lnd frac file on atm grid
  character(len=256) , public :: fatmlndfrc = ' '
  ! dynamic landuse dataset
  character(len=256) , public :: fpftdyn    = ' '
  ! ASCII data file with PFT physiological constants
  character(len=256) , public :: fpftcon    = ' '
  ! Restart file name
  character(len=256) , public :: nrevsn     = ' '
  ! snow optical properties file name
  character(len=256) , public :: fsnowoptics  = ' '
  ! snow aging parameters file name
  character(len=256) , public :: fsnowaging   = ' '
  !
  ! Irrigate logic
  !
  logical , public :: irrigate = .false. ! do not irrigate by default
  !
  ! Landunit logic
  !
  ! true => separate crop landunit is not created by default
  logical , public :: create_crop_landunit = .false.
  ! true => allocate memory for all possible vegetated pfts on
  ! vegetated landunit if at least one pft has nonzero weight
  logical , public :: allocate_all_vegpfts = .false.
  !
  ! BGC logic and datasets
  !
  ! values of 'prognostic','diagnostic','constant'
  character(len=16) , public :: co2_type = 'constant'
  ! State of the model for the accelerated decomposition (AD) spinup.
  !  0 (default) = normal model; 1 = AD SPINUP
  integer(ik4) , public :: spinup_state = 0
  ! used to override an error check on reading in restart files
  logical , public :: override_bgc_restart_mismatch_dump = .false.
  !
  ! Physics
  !
  !use subgrid fluxes
  integer(ik4) , public :: subgridflag = 1
  ! true => write global average diagnostics to std out
  logical , public :: wrtdia = .false.
  ! atmospheric CO2 molar ratio (by volume) (umol/mol)
  real(rk8) , public :: co2_ppmv = 355.D0
#if (defined LCH4 && defined VERTSOILC)
  ! true => anoxia is applied to heterotrophic respiration
  ! also considered in CH4 model
  logical , public :: anoxia = .true.
#else
  logical , public :: anoxia = .false.
#endif
  !
  ! C isotopes
  !
  ! true => use C-13 model
  logical , public :: use_c13 = .false.
  ! true => use C-14 model
  logical , public :: use_c14 = .false.
  !
  ! upper limit of each class (m)  (set in surfrd)
  real(rk8) , public , allocatable , dimension(:) :: glc_topomax
  ! glc_grid used to determine fglcmask  
  character(len=256) , public :: glc_grid = ' '
  ! glacier mask file name (based on glc_grid)
  character(len=256) , public :: fglcmask = ' '
  !
  ! instance control
  !
  integer(ik4) , public :: inst_index
  character(len=16) , public :: inst_name
  character(len=16) , public :: inst_suffix
  !
  ! Derived variables (run, history and restart file)
  !
  ! directory name for local restart pointer file
  character(len=256) , public :: rpntdir = '.'
  ! file name for local restart pointer file
  character(len=256) , public :: rpntfil = 'rpointer.lnd'
  !
  ! Error growth perturbation limit
  !
  ! perturbation limit when doing error growth test
  real(rk8) , public :: pertlim = 0.0D0
  !
  ! To retrieve namelist
  character(len=256) , public :: NLFilename_in ! Namelist filename
  !
  logical , private :: clmvarctl_isset = .false.

  contains
  !
  ! Set input control variables.
  !
  subroutine set_clmvarctl( caseid_in, ctitle_in, nsrest_in, &
                            version_in, hostname_in, username_in)
    character(len=*) , optional , intent(in) :: caseid_in  ! case id
    character(len=*) , optional , intent(in) :: ctitle_in  ! case title
    ! 0: initial run. 1: restart
    integer(ik4) , optional , intent(in) :: nsrest_in
    ! model version
    character(len=*) , optional , intent(in) :: version_in
    ! hostname running on
    character(len=*) , optional , intent(in) :: hostname_in
    ! username running job
    character(len=*) , optional , intent(in) :: username_in

    character(len=32) :: subname = 'set_clmvarctl'  ! subroutine name

    if ( clmvarctl_isset ) then
      call fatal(__FILE__,__LINE__, &
        subname//' ERROR:: control variables already set -- '// &
        'can not call this subroutine')
    end if
    if ( present(caseid_in       ) ) caseid        = caseid_in
    if ( present(ctitle_in       ) ) ctitle        = ctitle_in
    if ( present(nsrest_in       ) ) nsrest        = nsrest_in
    if ( present(version_in      ) ) version       = version_in
    if ( present(username_in     ) ) username      = username_in
    if ( present(hostname_in     ) ) hostname      = hostname_in
  end subroutine set_clmvarctl
!
! Check that values are correct, and finish setting variables based on
! other variables.
!
  subroutine clmvarctl_init( )
    character(len=32) :: subname = 'clmvarctl_init'  ! subroutine name

    ! landunit generation

    write (0,*) 'START clmvarctl_init'
    write (0,*) maxpatch_pft , numpft+1
    if (maxpatch_pft == numpft+1) then
      allocate_all_vegpfts = .true.
    else
      allocate_all_vegpfts = .false.
#ifdef CROP
      write(stderr,*)'maxpatch_pft = ',maxpatch_pft,&
            ' does NOT equal numpft+1 = ',numpft+1
      call fatal(__FILE__,__LINE__, &
                subname//' ERROR:: Can NOT turn CROP on without all PFTs' )
#endif
    end if
    write (0,*) allocate_all_vegpfts
    write (0,*) 'END clmvarctl_init'

    if ( myid == italk ) then

      ! Consistency settings for co2 type
      if ( co2_type /= 'constant' .and. &
           co2_type /= 'prognostic' .and. &
           co2_type /= 'diagnostic') then
        write(stderr,*)'co2_type = ',co2_type,' is not supported'
        call fatal(__FILE__,__LINE__, &
           subname//' ERROR:: choices are constant, prognostic or diagnostic' )
      end if

      ! Consistency settings for dynamic land use, etc.
      if (fpftdyn /= ' ' .and. create_crop_landunit) &
        call fatal(__FILE__,__LINE__, &
           subname//' ERROR:: dynamic landuse is currently not '// &
           'supported with create_crop_landunit option' )
      if (create_crop_landunit .and. .not.allocate_all_vegpfts) &
        call fatal(__FILE__,__LINE__, &
           subname//' ERROR:: maxpft<numpft+1 is currently not supported '//&
           'with create_crop_landunit option' )
      if (fpftdyn /= ' ') then
#if (defined CNDV)
        call fatal(__FILE__,__LINE__, &
            subname//' ERROR:: dynamic landuse is currently not '//&
            'supported with CNDV option' )
#endif
      end if

      ! Model physics

      if ( (co2_ppmv <= 0.0D0) .or. (co2_ppmv > 3000.0D0) ) &
        call fatal(__FILE__,__LINE__, &
           subname//' ERROR: co2_ppmv is out of a reasonable range' )

      if ( nsrest == nsrStartup ) nrevsn = ' '
      if ( nsrest == nsrContinue ) nrevsn = 'set by restart pointer file file'
      if ( nsrest /= nsrStartup .and. &
           nsrest /= nsrContinue ) &
        call fatal(__FILE__,__LINE__, &
           subname//' ERROR: nsrest NOT set to a valid value' )

#ifndef LCH4
      if ( anoxia ) then
        call fatal(__FILE__,__LINE__, &
           subname//'ERROR:: anoxia is turned on, but this currently '//&
           'requires turning on the CH4 submodel')
      end if
#endif

    endif   ! end of if-masterproc if-block

    clmvarctl_isset = .true.

  end subroutine clmvarctl_init

end module mod_clm_varctl
