module mod_clm_initialize
  !
  ! Performs land model initialization
  !
  use mod_intkinds
  use mod_realkinds
  use mod_runparams
  use mod_date
  use mod_stdio
  use mod_sunorbit
  use mod_mpmessage
  use mod_mppparam
  use mod_dynparam
  use mod_regcm_types
  use mod_clm_nchelper
  use mod_clm_varctl , only : nsrest , nsrStartup , nsrContinue , &
          fsurdat , fatmlndfrc , noland , finidat ,   &
          version , atm_regcm
  use mod_clm_varsur , only : wtxy , vegxy
  use mod_clm_typeinit , only : initClmtype
  use mod_clm_varpar , only : maxpatch , clm_varpar_init
  use mod_clm_varcon , only : clm_varcon_init
  use mod_clm_pftvarcon , only : pftconrd
  use mod_clm_decompinit , only : decompInit_lnd , decompInit_glcp
  use mod_clm_decomp , only : get_proc_bounds
  use mod_clm_domain , only : domain_check , ldomain , domain_init
  use mod_clm_surfrd , only : surfrd_get_grid , surfrd_get_data
  use mod_clm_control , only : control_init , control_print
  use mod_clm_urbaninput , only : UrbanInput
  use mod_clm_atmlnd
  use mod_clm_initgridcells , only : initGridCells
  use mod_clm_filter , only : allocFilters
  use mod_clm_reweight , only : reweightWrapup
#if (defined LCH4)
  use mod_clm_ch4varcon , only : ch4conrd
  use mod_clm_initch4 , only : initch4
#endif
#ifdef CN
  use mod_clm_cnecosystemdyn , only : CNEcosystemDynInit
  use mod_clm_cninitimevar , only : CNiniTimeVar
#endif
  use mod_clm_initslake , only : initSLake
  use mod_clm_mkarbinit , only : mkregcminit
  use mod_clm_pftdyn , only : pftdyn_init , pftdyn_interp
#if (defined CNDV)
  use mod_clm_pftdyn , only : pftwt_init
  use mod_clm_cndvecosystemdynini , only : CNDVEcosystemDynini
#endif
  use mod_clm_staticecosysdyn , only : EcosystemDynini , readAnnualVegetation
  use mod_clm_staticecosysdyn , only : interpMonthlyVeg
  use mod_clm_histflds , only : hist_initFlds
  use mod_clm_histfile , only : hist_htapes_build , htapes_fieldlist
  use mod_clm_restfile , only : restFile_getfile, &
                                 restFile_open, restFile_close, restFile_read
  use mod_clm_accflds , only : initAccFlds , initAccClmtype
  use mod_clm_dust , only : Dustini
  use mod_clm_time_manager, only : get_curr_calday
  use mod_clm_urban , only : UrbanClumpInit
  use mod_clm_urbaninit , only : UrbanInitTimeConst , UrbanInitTimeVar , &
          UrbanInitAero
  use mod_clm_urbaninput , only : UrbanInput
#if (defined LCH4)
#endif
  use mod_clm_drydep , only : n_drydep, drydep_method, DD_XLND
  use mod_clm_initsurfalb , only : initSurfAlb, do_initsurfalb
  use mod_clm_vocemission , only : VOCEmission_init
  use mod_clm_initimeconst , only : iniTimeConst

  implicit none

  private

  save

  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization

  contains
  !
  ! Land model initialization.
  ! o Initializes run control variables via the [clm_inparm] namelist.
  ! o Reads surface data on model grid.
  ! o Defines the multiple plant types and fraction areas for each surface
  !   type.
  ! o Builds the appropriate subgrid <-> grid mapping indices and weights.
  ! o Set up parallel processing.
  ! o Initializes time constant variables.
  ! o Reads restart data for a restart run.
  ! o Reads initial data and initializes the time variant variables for
  !   an initial run.
  ! o Initializes history file output.
  ! o Initializes river routing model.
  ! o Initializes accumulation variables.
  !
  subroutine initialize1
    implicit none
    integer(ik4)  :: ier              ! error status
    integer(ik4)  :: begg , endg      ! beg and ending gridcell indices

    ! -------------------------------------------
    ! Initialize run control variables, timestep
    ! -------------------------------------------

    if ( myid == italk ) then
      write(stdout,*) 'CLM version: 4.5 - RegCM ', trim(version)
    end if

    if ( myid == italk ) then
      write(stdout,*) 'Attempting to initialize the land model .....'
      write(stdout,*) 'Mask given by RegCM model has a total of ', &
              sum(lndcomm%linear_npoint_sg), ' land points'
    endif

    call control_init()
    call control_print()

    call clm_varpar_init()
    call clm_varcon_init()

    ! ------------------
    ! Set decomposition
    ! ------------------

    ! Determine clm decomposition

    call decompInit_lnd

    ! Get grid and land fraction (set ldomain)

    if (myid == italk) then
       write(stdout,*) 'Attempting to read ldomain from ',trim(fatmlndfrc)
    endif
    call surfrd_get_grid(ldomain, fatmlndfrc)
    if (myid == italk) then
       call domain_check(ldomain)
    endif

    ! Initialize urban model input (initialize urbinp data structure)

    call UrbanInput(mode='initialize')

    ! Allocate surface grid dynamic memory (for wtxy and vegxy arrays)
    ! Allocate additional dynamic memory for glacier_mec topo and thickness

    call get_proc_bounds(begg, endg)
    allocate (vegxy(begg:endg,maxpatch), wtxy(begg:endg,maxpatch), stat=ier)
    if (ier /= 0) then
       write(stderr,*)'initialize allocation error'
       call fatal(__FILE__,__LINE__,'clm now stopping')
    endif

    ! Read list of PFTs and their corresponding parameter values
    ! Independent of model resolution, Needs to stay before surfrd_get_data

    call pftconrd()

    ! Read surface dataset and set up vegetation type [vegxy] and
    ! weight [wtxy] arrays for [maxpatch] subgrid patches.

    call surfrd_get_data(ldomain,fsurdat)

    ! Determine decomposition of subgrid scale landunits, columns, pfts

    call decompInit_glcp

#if (defined LCH4)
    ! Set CH4 Model Parameters from namelist.
    call ch4conrd()
    ! Need to do before iniTimeConst so that it knows whether to look for
    ! several optional parameters on surfdata file.
#endif

    ! Allocate memory and initialize values of clmtype data structures

    call initClmtype()

    ! Initialize lnd->atm and atm->lnd data structures

    call init_atm2lnd_type(begg,endg,clm_a2l)
    call init_lnd2atm_type(begg,endg,clm_l2a)

    ! Build hierarchy and topological info for derived types

    call initGridCells()

    ! Deallocate surface grid dynamic memory (for wtxy and vegxy arrays)

    deallocate(vegxy,wtxy)

  end subroutine initialize1
  !
  ! Land model initialization.
  ! o Initializes run control variables via the [clm_inparm] namelist.
  ! o Reads surface data on model grid.
  ! o Defines the multiple plant types and fraction areas for each surface type.
  ! o Builds the appropriate subgrid <-> grid mapping indices and weights.
  ! o Set up parallel processing.
  ! o Initializes time constant variables.
  ! o Reads restart data for a restart run.
  ! o Reads initial data and initializes the time variant variables for an
  !   initial run.
  ! o Initializes history file output.
  ! o Initializes river routing model.
  ! o Initializes accumulation variables.
  !
  subroutine initialize2(rdate)
    implicit none
    character(len=*) , intent(in) :: rdate
    integer(ik4) :: yr      ! current year (0, ...)
    integer(ik4) :: mon     ! current month (1 -> 12)
    integer(ik4) :: day     ! current day (1 -> 31)
    integer(ik4) :: ncsec   ! current time of day [seconds]
    integer(ik4) :: begp , endp   ! beg and ending pft indices
    integer(ik4) :: begc , endc   ! beg and ending column indices
    integer(ik4) :: begl , endl   ! beg and ending landunit indices
    integer(ik4) :: begg , endg   ! beg and ending gridcell indices
    character(len=256) :: fnamer  ! name of netcdf restart file
    real(rk8) :: calday           ! calendar day
    real(rk8) :: caldaym1         ! calendar day for nstep-1
    real(rk8) :: declin           ! solar declination angle in radians
    real(rk8) :: declinm1         ! solar declination angle in radians
    real(rk8) :: eccf             ! earth orbit eccentricity factor

    ! ------------------------------------------------------------------------
    ! Initialize time constant variables
    ! ------------------------------------------------------------------------

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp )

    ! Initialize Ecosystem Dynamics

#if (defined CNDV)
    call CNDVEcosystemDynini()
#elif (!defined CN)
    call EcosystemDynini()
#endif
#if (defined CN) || (defined CNDV)

    ! --------------------------------------------------------------
    ! Initialize CLMSP ecosystem dynamics when drydeposition is used
    ! so that estimates of monthly differences in LAI can be computed
    ! --------------------------------------------------------------

    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
      call EcosystemDynini()
    end if
#endif

    ! Initialize dust emissions model

    call Dustini()

    ! Initialize MEGAN emissions model

    call VOCEmission_init( )

    ! ------------------------------------------------------------------------
    ! Initialize time constant urban variables
    ! ------------------------------------------------------------------------

    call UrbanInitTimeConst()
    call iniTimeConst()

    ! ------------------------------------------------------------------------
    ! Initialize master history list.
    ! ------------------------------------------------------------------------

    call hist_initFlds()
    ! On restart process the history namelist. Later the namelist from the
    ! restart file will be used. But, this allows some basic checking to make
    ! sure you didn't try to change the history namelist on restart.
    if ( nsrest == nsrContinue ) then
      call htapes_fieldlist()
    end if

    ! ------------------------------------------------------------------------
    ! Initialize CN Ecosystem Dynamics
    !   (must be after time-manager initialization)
    ! ------------------------------------------------------------------------
#if (defined CN) || (defined CNDV)
    call CNEcosystemDynInit(begg,endg,begc,endc,begp,endp)
#endif

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! Initialize accumulator fields to be time accumulated for various purposes.
    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed.

    call initAccFlds()

    ! ------------------------------------------------------------------------
    ! Set arbitrary initial conditions for time varying fields
    ! used in coupled carbon-nitrogen code
    ! ------------------------------------------------------------------------

#if (defined CN)
    if ( nsrest == nsrStartup ) then
      call CNiniTimeVar()
    end if
#endif

    ! ------------------------------------------------------------------------
    ! Initialization of dynamic pft weights
    ! ------------------------------------------------------------------------

    ! Determine correct pft weights (interpolate pftdyn dataset if initial run)
    ! Otherwise these are read in for a restart run

#if (defined CNDV)
    call pftwt_init()
#else
#ifdef DYNPFT
    call pftdyn_init()
    call pftdyn_interp( )
#endif
#endif

    ! ------------------------------------------------------------------------
    ! Read restart/initial info
    ! ------------------------------------------------------------------------

    ! No weight related information can be contained in the routines,
    ! "mkregcminit, inicfile and restFile".

    if ( nsrest == nsrContinue ) then
      call restFile_getfile(fnamer, rdate)
      call restFile_read( fnamer )
    else if ( nsrest == nsrStartup ) then
      ! Get initial data from regcm !
      call mkregcminit(adomain)
      call UrbanInitTimeVar( )
    else
      call fatal(__FILE__,__LINE__,'CLM modified to run with RegCM !')
    end if

    ! Attn EK: The calls to initch4 and initSLake combine both setting of
    ! state vars, and constant + flux & diagnostic vars.  This is set up so
    ! that the submodels would be back-compatible with old restart files.
    ! It is intended to work and allow bfb restarts as is, but may not be
    ! consistent style.
    ! See these two routines for structure.  Feel free to modify this but be
    ! careful.
    ! You may want to keep at least the initch4 as is to allow CH4 to be
    ! run even if it wasn't spun up with CH4, as the CH4 sub-model comes to
    ! eq. within a month plus one year for annual mean variables.
    ! Of course if clm_varctl:anoxia is used or NITRIF_DENITRIF is defined,
    ! then CN would no longer be in eq.
    call initSLake(.false.)
#if (defined LCH4)
    call initch4(.false.)
#endif

    ! ------------------------------------------------------------------------
    ! Initialization of model parameterizations that are needed after
    ! restart file is read in
    ! ------------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Initialize history and accumator buffers
    ! ------------------------------------------------------------------------

    ! Initialize active history fields. This is only done if not a restart run.
    ! If a restart run, then this information has already been obtained from
    ! the restart data read above.
    ! Note that routine hist_htapes_build needs time manager information, so
    ! this call must be made after the restart information has been read.

    if ( nsrest == nsrStartup ) then
      call hist_htapes_build()
    end if

    ! Initialize clmtype variables that are obtained from accumulated fields.
    ! This routine is called in an initial run at ktau=0
    ! This routine is also always called for a restart run and must
    ! therefore be called after the restart file is read in

    call initAccClmtype()

    ! --------------------------------------------------------------
    ! Note - everything below this point needs updated weights
    ! --------------------------------------------------------------

    ! Initialize filters

    call allocFilters()
    call reweightWrapup()

    ! Calculate urban "town" roughness length and displacement
    ! height for urban landunits

    call UrbanInitAero()

    ! Initialize urban radiation model - this uses urbinp data structure

    call UrbanClumpInit()

    ! Finalize urban model initialization

    call UrbanInput(mode='finalize')

    !
    ! Even if CN is on, and dry-deposition is active, read CLMSP annual
    ! vegetation to get estimates of monthly LAI
    !
    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
      call readAnnualVegetation()
    end if

    ! End initialization

    if (myid == italk) then
      write(stdout,*) 'Successfully initialized the land model'
      if (nsrest == nsrStartup) then
        write(stdout,*) 'begin initial run at: '
      else
        write(stdout,*) 'begin continuation run at:'
      end if
      call curr_date(idatex, yr, mon, day, ncsec)
      write(stdout,*) '   ktau    = ',ktau
      write(stdout,*) '   year    = ',yr
      write(stdout,*) '   month   = ',mon
      write(stdout,*) '   day     = ',day
      write(stdout,*) '   seconds = ',ncsec
    endif

    if ( nsrest == nsrStartup ) then
      ! Initialize albedos (correct pft filters are needed)
      if (finidat == ' ' .or. do_initsurfalb) then
        calday = get_curr_calday()
        call orb_decl(calday,eccen,mvelpp,lambm0,obliqr,declin,eccf )
        caldaym1 = get_curr_calday(offset=-int(dtsrf))
        call orb_decl(caldaym1,eccen,mvelpp,lambm0,obliqr,declinm1,eccf )
        call initSurfAlb(calday,declin,declinm1)
      else if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
        ! Call interpMonthlyVeg for dry-deposition so that mlaidiff will
        ! be calculated This needs to be done even if CN or CNDV is on!
        call interpMonthlyVeg()
      end if
      ! Determine gridcell averaged properties to send to atm
      call clm_map2gcell(init=.true.)
    end if
  end subroutine initialize2

end module mod_clm_initialize
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
