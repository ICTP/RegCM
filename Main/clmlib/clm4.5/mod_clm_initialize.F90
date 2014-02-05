module mod_clm_initialize
  !
  ! Performs land model initialization
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_clm_varctl , only : nsrest , nsrStartup , nsrContinue , &
          nsrBranch , create_glacier_mec_landunit , fsurdat ,     &
          fatmlndfrc , flndtopo , fglcmask , noland , finidat ,   &
          fpftdyn
  use mod_clm_varsur , only : wtxy , vegxy , topoxy
  use mod_clm_typeinit , only : initClmtype
  use mod_clm_varpar , only : maxpatch , clm_varpar_init
  use mod_clm_varcon , only : clm_varcon_init
  use mod_clm_pftvarcon , only : pftconrd
  use mod_clm_decompInit , only : decompInit_lnd , decompInit_glcp
  use mod_clm_decomp , only : get_proc_bounds
  use mod_clm_domain , only : domain_check , ldomain , domain_init
  use mod_clm_surfrd , only : surfrd_get_globmask , surfrd_get_grid , &
          surfrd_get_topo , surfrd_get_data 
  use mod_clm_control , only : control_init , control_print , nlfilename
  use mod_clm_urbaninput , only : UrbanInput
  use clm_atmlnd
  use mod_clm_initgridcells , only : initGridCells
  use mod_clm_filter , only : allocFilters
  use mod_clm_reweight , only : reweightWrapup
#if (defined LCH4)
  use mod_clm_ch4varcon , only : ch4conrd
  use mod_clm_initch4 , only : initch4
#endif
#ifdef CN
  use mod_clm_cnecosystemdyn , only : CNEcosystemDynInit
#endif
  use mod_clm_initslake , only : initSLake
  use mod_clm_mkarbinit , only : mkarbinit
  use mod_clm_pftdyn , only : pftdyn_init , pftdyn_interp
#if (defined CNDV)
  use mod_clm_pftdyn , only : pftwt_init
  use mod_clm_cndvecosystemdynini , only : CNDVEcosystemDynini
#endif
  use mod_clm_staticecosysdyn , only : EcosystemDynini , readAnnualVegetation
  use mod_clm_staticecosysdyn , only : interpMonthlyVeg

  implicit none

  private

  public :: initialize1  ! Phase one initialization
  public :: initialize2  ! Phase two initialization
  private :: header      ! echo version numbers
  private :: do_restread ! read a restart file

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
  ! o Reads restart data for a restart or branch run.
  ! o Reads initial data and initializes the time variant variables for
  !   an initial run.
  ! o Initializes history file output.
  ! o Initializes river routing model.
  ! o Initializes accumulation variables.
  !
  subroutine initialize1(cl)
    implicit none
    (type (masked_comm) , intent(in) :: cl
    integer(ik4)  :: ier              ! error status
    integer(ik4)  :: i , j , n , k    ! loop indices
    integer(ik4)  :: nl               ! gdc and glo lnd indices
    integer(ik4)  :: ns , ni , nj     ! global grid sizes
    logical  :: isgrid2d              ! true => global grid is regular lat/lon
    integer(ik4)  :: begp , endp      ! beg and ending pft indices
    integer(ik4)  :: begc , endc      ! beg and ending column indices
    integer(ik4)  :: begl , endl      ! beg and ending landunit indices
    integer(ik4)  :: begg , endg      ! beg and ending gridcell indices
    integer(ik4) , pointer  :: amask(:)  ! global land mask
    character(len=32) :: subname = 'initialize1' ! subroutine name

    ! -------------------------------------------
    ! Initialize run control variables, timestep
    ! -------------------------------------------

    call header()

    if ( myid == italk ) then
      write(stdout,*) 'Attempting to initialize the land model .....'
      write(stdout,*)
    endif

    call control_init()
    call clm_varpar_init()
    call clm_varcon_init()

    if ( myid == italk ) call control_print()

    ! ------------------
    ! Set decomposition
    ! ------------------

    ! Determine clm decomposition

    call decompInit_lnd(cl)

    ! Get grid and land fraction (set ldomain)

    if (myid == italk) then
       write(stdout,*) 'Attempting to read ldomain from ',trim(fatmlndfrc)
    endif
    if (create_glacier_mec_landunit) then
       call surfrd_get_grid(ldomain, fatmlndfrc, fglcmask)
    else
       call surfrd_get_grid(ldomain, fatmlndfrc)
    endif
    if (myid == italk) then
       call domain_check(ldomain)
    endif
    ldomain%mask = 1  !!! TODO - is this needed?

    ! Get topo if appropriate (set ldomain%topo)

    if (flndtopo /= " ") then
       if (myid == italk) then
          write(stdout,*) 'Attempting to read atm topo from ',trim(flndtopo)
       endif
       call surfrd_get_topo(ldomain, flndtopo)  
    endif

    ! Initialize urban model input (initialize urbinp data structure)

    call UrbanInput(mode='initialize')

    ! Allocate surface grid dynamic memory (for wtxy and vegxy arrays)
    ! Allocate additional dynamic memory for glacier_mec topo and thickness

    call get_proc_bounds(begg, endg)
    allocate (vegxy(begg:endg,maxpatch), wtxy(begg:endg,maxpatch), stat=ier)   
    if (create_glacier_mec_landunit) then
       allocate (topoxy(begg:endg,maxpatch), stat=ier)
    else
       allocate (topoxy(1,1), stat=ier)
    endif
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

    if (create_glacier_mec_landunit) then
       call decompInit_glcp (ns, ni, nj, ldomain%glcmask)
    else
       call decompInit_glcp (ns, ni, nj)
    endif

#if (defined LCH4)
    ! Set CH4 Model Parameters from namelist.
    call ch4conrd()
    ! Need to do before iniTimeConst so that it knows whether to look for several optional parameters on surfdata file.
#endif

    ! Allocate memory and initialize values of clmtype data structures

    call initClmtype()

    ! Initialize lnd->atm data structure

    call init_atm2lnd_type(begg,endg,clm_a2l)
    call init_lnd2atm_type(begg,endg,clm_l2a)

    ! if (create_glacier_mec_landunit) then
    !    call init_glc2lnd_type(begg, endg, clm_x2s)
    !    call init_lnd2glc_type(begg, endg, clm_s2x)
    ! endif

    ! Build hierarchy and topological info for derived types

    call initGridCells()

    ! Deallocate surface grid dynamic memory (for wtxy and vegxy arrays)

    deallocate(vegxy,wtxy,topoxy)
  end subroutine initialize1
  !
  ! Land model initialization.
  ! o Initializes run control variables via the [clm_inparm] namelist.
  ! o Reads surface data on model grid.
  ! o Defines the multiple plant types and fraction areas for each surface type.
  ! o Builds the appropriate subgrid <-> grid mapping indices and weights.
  ! o Set up parallel processing.
  ! o Initializes time constant variables.
  ! o Reads restart data for a restart or branch run.
  ! o Reads initial data and initializes the time variant variables for an
  !   initial run.
  ! o Initializes history file output.
  ! o Initializes river routing model.
  ! o Initializes accumulation variables.
  !
  subroutine initialize2( )
    use histFldsMod     , only : hist_initFlds
    use histFileMod     , only : hist_htapes_build, htapes_fieldlist
    use restFileMod     , only : restFile_getfile, &
                                 restFile_open, restFile_close, restFile_read 
    use accFldsMod      , only : initAccFlds, initAccClmtype
    use DustMod         , only : Dustini
    use clm_time_manager, only : curr_date, advance_timestep, &
                                 timemgr_init, timemgr_restart_io, timemgr_restart
    use clm_time_manager, only : get_step_size, get_curr_calday
    use fileutils       , only : getfil
    use UrbanMod        , only : UrbanClumpInit
    use UrbanInitMod    , only : UrbanInitTimeConst, UrbanInitTimeVar, UrbanInitAero 
    use UrbanInputMod   , only : UrbanInput
#if (defined LCH4)
#endif
    use clm_glclnd      , only : init_glc2lnd_type, init_lnd2glc_type, &
                                 clm_x2s, clm_s2x
    use seq_drydep_mod  , only : n_drydep, drydep_method, DD_XLND
    use shr_orb_mod        , only : shr_orb_decl
    use initSurfAlbMod     , only : initSurfAlb, do_initsurfalb 
    use clm_varorb         , only : eccen, mvelpp, lambm0, obliqr
    use VOCEmissionMod  , only : VOCEmission_init
    implicit none
    integer(ik4) :: nl , na , nag     ! indices
    integer(ik4) :: i , j , k         ! indices
    integer(ik4) :: yr                ! current year (0, ...)
    integer(ik4) :: mon               ! current month (1 -> 12)
    integer(ik4) :: day               ! current day (1 -> 31)
    integer(ik4) :: ncsec             ! current time of day [seconds]
    integer(ik4) :: begp , endp       ! beg and ending pft indices
    integer(ik4) :: begc , endc       ! beg and ending column indices
    integer(ik4) :: begl , endl       ! beg and ending landunit indices
    integer(ik4) :: begg , endg       ! beg and ending gridcell indices
    character(len=256) :: fnamer      ! name of netcdf restart file 
    character(len=256) :: pnamer      ! full pathname of netcdf restart file
    type(file_desc_t) :: ncid         ! netcdf id
    real(rk8) :: calday               ! calendar day
    real(rk8) :: caldaym1             ! calendar day for nstep-1
    real(rk8) :: declin               ! solar declination angle in radians
    real(rk8) :: declinm1              ! solar declination angle in radians
    real(rk8) :: eccf                  ! earth orbit eccentricity factor
    character(len=32) :: subname = 'initialize2' ! subroutine name
    logical           :: arbinit            ! Make arb init in initSLake

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
    ! Obtain restart file if appropriate
    ! ------------------------------------------------------------------------

    if ( do_restread() ) then
      call restFile_getfile(file=fnamer,path=pnamer)
    end if

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
    if (nsrest == nsrStartup) then
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
    if (fpftdyn /= ' ') then
      call pftdyn_init()
      call pftdyn_interp( )
    end if
#endif

    ! ------------------------------------------------------------------------
    ! Read restart/initial info
    ! ------------------------------------------------------------------------

    ! No weight related information can be contained in the routines,  
    ! "mkarbinit, inicfile and restFile". 

    if ( do_restread() ) then
      if ( myid == italk ) then
        write(stdout,*)'reading restart file ',fnamer
      end if
      call restFile_read( fnamer )

      arbinit = .false.
      call initSLake(arbinit)
#if (defined LCH4)
      arbinit = .false.
      call initch4(arbinit)
#endif
    else if (nsrest == nsrStartup .and. finidat == ' ') then
      call mkarbinit()
      call UrbanInitTimeVar( )

      arbinit = .true.
      call initSLake(arbinit)
#if (defined LCH4)
      arbinit = .true.
      call initch4(arbinit)
#endif
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
    end if

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

    if ( nsrest == nsrStartup .or. nsrest == nsrBranch ) then
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
      call curr_date(yr, mon, day, ncsec)
      write(stdout,*) '   ktau= ',ktau, ' year= ',yr,' month= ',mon,&
            ' day= ',day,' seconds= ',ncsec
      write(stdout,*)
      write(stdout,'(72a1)') ("*",i=1,60)
      write(stdout,*)
    endif

    if ( ktau == 0 .or. nsrest == nsrStartup ) then
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
    ! Initialize sno export state
    !if (create_glacier_mec_landunit) then
    !  call create_clm_s2x(init=.true.)
    !end if
  end subroutine initialize2
  !
  ! Echo and save model version number
  !
  subroutine header()
    use mod_clm_varctl , only : version
    implicit none
    if ( myid == italk ) then
      write(stdout,*) trim(version)
      write(stdout,*)
    end if
  end subroutine header
  !
  ! Determine if restart file will be read
  !
  logical function do_restread( )
    use mod_clm_varctl, only : finidat
    implicit none
    do_restread = ktau > 0
  end function do_restread

end module mod_clm_initialize
