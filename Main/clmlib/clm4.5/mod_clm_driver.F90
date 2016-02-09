module mod_clm_driver
  ! This module provides the main CLM driver physics calling sequence.  Most
  ! computations occurs over gridcells (and associated subgrid
  ! scale entities) assigned to each MPI process.
  !
  ! The main CLM driver physics calling sequence for clm_driver1 is as follows:
  !
  ! + interpMonthlyVeg      interpolate monthly vegetation data [! CN or ! CNDV]
  ! + readMonthlyVegetation read vegetation data for two months [! CN or ! CNDV]
  !
  ! -> dynland_hwcontent   Get initial heat, water content
  !    + pftdyn_interp                                          [pftdyn]
  !    + dynland_hwcontent   Get new heat, water content        [pftdyn]
  !
  ! -> clm_driverInit      save of variables from previous time step
  ! -> Hydrology1          canopy interception and precip on ground
  ! -> FracWet             fraction of wet vegetated surface and dry elai
  ! -> SurfaceRadiation    surface solar radiation
  ! -> UrbanRadiation      surface solar and longwave radiation for Urban land
  ! -> Biogeophysics1      leaf temperature and surface fluxes
  ! -> BareGroundFluxes    surface fluxes for bare soil or snow-covered
  !                        vegetation patches
  ! -> UrbanFluxes         surface fluxes for urban landunits
  ! -> MoninObukIni        first-guess Monin-Obukhov length and wind speed
  ! -> FrictionVelocity    friction velocity and potential temperature and
  !                        humidity profiles
  ! -> CanopyFluxes        leaf temperature and surface fluxes for vegetated
  !                        patches
  ! -> QSat                saturated vapor pressure, specific humidity, &
  !                        derivatives at leaf surface
  ! -> MoninObukIni        first-guess Monin-Obukhov length and wind speed
  ! -> FrictionVelocity    friction velocity and potential temperature and
  !                        humidity profiles
  ! -> Stomata             stomatal resistance and photosynthesis for
  !                        sunlit leaves
  ! -> Stomata             stomatal resistance and photosynthesis for
  !                        shaded leaves
  ! -> QSat                recalculation of saturated vapor pressure,
  !                        specific humidity, & derivatives at leaf surface
  ! + DustEmission         Dust mobilization
  ! + DustDryDep           Dust dry deposition
  ! -> SLakeFluxes         lake surface fluxes
  ! -> SLakeTemperature    lake temperature
  ! + VOCEmission          compute VOC emission                 [VOC]
  ! -> Biogeophysics2      soil/snow & ground temp and update surface fluxes
  ! -> pft2col             Average from PFT level to column level
  ! -> Hydrology2          surface and soil hydrology
  ! -> SLakeHydrology      lake hydrology
  ! -> SnowAge_grain       update snow effective grain size for snow radiative
  !                        transfer
  ! + CNEcosystemDyn       Carbon Nitrogen model ecosystem dynamics: [CN]
  !                        vegetation phenology and soil carbon
  ! + EcosystemDyn         "static" ecosystem dynamics:              [! CN ]
  !                        vegetation phenology and soil carbon
  ! -> BalanceCheck        check for errors in energy and water balances
  ! -> SurfaceAlbedo       albedos for next time step
  ! -> UrbanAlbedo         Urban landunit albedos for next time step
  !
  ! Second phase of the clm main driver, for handling history and restart
  ! file output.
  !
  ! -> updateAccFlds       update accumulated fields
  ! -> hist_update_hbuf    accumulate history fields for time interval
  ! -> htapes_wrapup       write history tapes if appropriate
  ! -> restFile_write      write restart file if appropriate
  ! \end{verbatim}
  !
  ! Optional subroutines are denoted by an plus (+) with the associated
  ! CPP token or variable in brackets at the end of the line.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date
  use mod_mppparam , only : italk
  use mod_runparams
  use mod_clm_type
  use mod_clm_varctl , only : wrtdia
  use mod_clm_decomp , only : get_proc_bounds
  use mod_clm_filter , only : filter
  use mod_clm_reweight , only : reweightWrapup
#if (defined CNDV)
  use mod_clm_cndv , only : dv, histCNDV
  use mod_clm_pftdyn , only : pftwt_interp
#endif
  use mod_clm_pftdyn , only : pftdyn_interp, pftdyn_wbal_init, pftdyn_wbal
#ifdef CN
  use mod_clm_pftdyn , only : pftdyn_cnbal
#endif
  use mod_clm_dynland , only : dynland_hwcontent
  use mod_clm_varcon , only : isturb
  use mod_clm_histfile , only : hist_update_hbuf, hist_htapes_wrapup
  use mod_clm_restfile , only : restFile_write, restFile_filename
  use mod_clm_accflds , only : updateAccFlds
  use mod_clm_driverinit , only : clm_driverInit
  use mod_clm_balancecheck , only : BeginWaterBalance, BalanceCheck
  use mod_clm_surfaceradiation , only : SurfaceRadiation
  use mod_clm_hydrology1 , only : Hydrology1
  use mod_clm_hydrology2 , only : Hydrology2
  use mod_clm_slakefluxes , only : SLakeFluxes
  use mod_clm_slaketemperature , only : SLakeTemperature
  use mod_clm_slakehydrology , only : SLakeHydrology
  use mod_clm_biogeophysics1 , only : Biogeophysics1
  use mod_clm_baregroundfluxes , only : BareGroundFluxes
  use mod_clm_canopyfluxes , only : CanopyFluxes
  use mod_clm_biogeophysics2 , only : Biogeophysics2
  use mod_clm_surfacealbedo , only : SurfaceAlbedo
  use mod_clm_pft2col , only : pft2col
#if (defined CN)
  use mod_clm_cnsetvalue , only : CNZeroFluxes_dwt
  use mod_clm_cnecosystemdyn , only : CNEcosystemDyn
  use mod_clm_cnannualupdate , only : CNAnnualUpdate
  use mod_clm_cnbalancecheck , only : BeginCBalance, BeginNBalance, &
                                      CBalanceCheck, NBalanceCheck
  use mod_clm_cnverticalprofile , only : decomp_vertprofiles
  use mod_clm_cnfire , only : CNFireInterp
#else
  use mod_clm_staticecosysdyn , only : EcosystemDyn
#endif
  use mod_clm_activelayer , only : alt_calc
  use mod_clm_dust , only : DustDryDep, DustEmission
  use mod_clm_vocemission , only : VOCEmission
  use mod_clm_drydep , only : n_drydep, drydep_method, DD_XLND
  use mod_clm_staticecosysdyn , only : interpMonthlyVeg
  use mod_clm_drydepvelocity , only : depvel_compute
#if (defined LCH4)
  use mod_clm_ch4 , only : ch4
#endif
  use mod_clm_urban , only : UrbanAlbedo, UrbanRadiation, UrbanFluxes
  use mod_clm_snicar , only : SnowAge_grain
  use mod_clm_atmlnd , only : clm_map2gcell

  implicit none

  private

  save

  public :: clm_drv            ! clm physics,history, restart writes

  contains
  !
  ! First phase of the clm driver calling the clm physics. An outline of
  ! the calling tree is given in the description of this module.
  !
  subroutine clm_drv(doalb,nextsw_cday,declinp1,declin,rstwr,nlend,nlomon,rdate)
    implicit none
    logical , intent(in) :: doalb     ! true if time for surface albedo calc
    real(rk8) , intent(in) :: nextsw_cday ! calendar day
    real(rk8) , intent(in) :: declinp1 ! declination angle for next time step
    real(rk8) , intent(in) :: declin   ! declination angle for current time step
    logical , intent(in) :: rstwr     ! true => write restart file this step
    logical , intent(in) :: nlend     ! true => end of run on this step
    logical , intent(in) :: nlomon    ! true => end of month on this step
    character(len=*) , intent(in) :: rdate ! restart file time stamp for name
    ! landunit index associated with each column
    integer(ik4) , pointer , dimension(:) :: clandunit
    ! landunit type
    integer(ik4) , pointer , dimension(:) :: itypelun
    integer(ik4)  :: c, l, g   ! indices
    integer(ik4)  :: begg, endg ! beginning and ending gridcell indices
    integer(ik4)  :: begl, endl ! beginning and ending landunit indices
    integer(ik4)  :: begc, endc ! beginning and ending column indices
    integer(ik4)  :: begp, endp ! beginning and ending pft indices
    type(column_type) , pointer :: cptr   ! pointer to column derived subtype
#if (defined CNDV)
    integer(ik4)  :: yr       ! year (0, ...)
    integer(ik4)  :: mon      ! month (1, ..., 12)
    integer(ik4)  :: day      ! day of month (1, ..., 31)
    integer(ik4)  :: ncdate   ! current date
    integer(ik4)  :: nbdate   ! base date (reference date)
    integer(ik4)  :: kyr      ! thousand years, equals 2 at end of first year
    type(rcm_time_and_date) :: nextdate
    type(rcm_time_interval) :: tdif
#endif
    character(len=256) :: filer  ! restart file name
    !FAB
    logical, save  :: lfirstcall
    ! Assign local pointers to derived subtypes components (landunit-level)

    itypelun => clm3%g%l%itype

    ! Assign local pointers to derived subtypes components (column-level)

    clandunit => clm3%g%l%c%landunit

    ! Set pointers into derived type

    cptr => clm3%g%l%c

#ifdef CN
    ! For dry-deposition need to call CLMSP so that mlaidiff is obtained
    if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
      call interpMonthlyVeg()
    endif
#else
    ! ========================================================================
    ! Determine weights for time interpolation of monthly vegetation data.
    ! This also determines whether it is time to read new monthly vegetation
    ! and obtain updated
    !            leaf area index   [mlai1,mlai2]
    !            stem area index   [msai1,msai2]
    !            vegetation top    [mhvt1,mhvt2]
    !            vegetation bottom [mhvb1,mhvb2]
    ! The weights obtained here are used in subroutine ecosystemdyn to obtain
    ! time interpolated values.
    ! ========================================================================
    if ( doalb .or. ( n_drydep > 0 .and. drydep_method == DD_XLND ) ) then
      call interpMonthlyVeg()
    end if
#endif

    ! ========================================================================
    ! Determine boundaries
    ! ========================================================================

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! ========================================================================
    ! change pft weights and compute associated heat & water fluxes
    ! ========================================================================

    ! initialize heat and water content and dynamic balance fields to zero
    do g = begg , endg
      clm3%g%gwf%qflx_liq_dynbal(g) = 0.D0
      clm3%g%gws%gc_liq2(g)         = 0.D0
      clm3%g%gws%gc_liq1(g)         = 0.D0
      clm3%g%gwf%qflx_ice_dynbal(g) = 0.D0
      clm3%g%gws%gc_ice2(g)         = 0.D0
      clm3%g%gws%gc_ice1(g)         = 0.D0
      clm3%g%gef%eflx_dynbal(g)     = 0.D0
      clm3%g%ges%gc_heat2(g)        = 0.D0
      clm3%g%ges%gc_heat1(g)        = 0.D0
    end do

    !--- get initial heat,water content ---
    call dynland_hwcontent(begg,endg,                       &
                           clm3%g%gws%gc_liq1(begg:endg),   &
                           clm3%g%gws%gc_ice1(begg:endg),   &
                           clm3%g%ges%gc_heat1(begg:endg))

    ! ======================================================================
    ! Determine decomp vertical profiles
    !
    ! These routines (alt_calc & decomp_vertprofiles) need to be called
    ! before pftdyn_cnbal, and it appears that they need to be called before
    ! pftdyn_interp and the associated filter updates, too (otherwise we get
    ! a carbon balance error)
    ! =======================================================================

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    call alt_calc(filter%num_soilc,filter%soilc)

#if (defined CN)
    call decomp_vertprofiles(begp,endp,begc,endc, &
                             filter%num_soilc,    &
                             filter%soilc,        &
                             filter%num_soilp,    &
                             filter%soilp)
#endif

#if (!defined CNDV)
#ifdef DYNPFT
    call pftdyn_interp  ! change the pft weights
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! do stuff that needs to be done after changing weights
    ! This call should be made as soon as possible after pftdyn_interp
    call reweightWrapup

    !--- get new heat,water content: (new-old)/dt = flux into lnd model ---
    call dynland_hwcontent(begg,endg,                     &
                           clm3%g%gws%gc_liq2(begg:endg), &
                           clm3%g%gws%gc_ice2(begg:endg), &
                           clm3%g%ges%gc_heat2(begg:endg))
    do g = begg , endg
      clm3%g%gwf%qflx_liq_dynbal(g) = &
              (clm3%g%gws%gc_liq2 (g) - clm3%g%gws%gc_liq1 (g))/dtsrf
      clm3%g%gwf%qflx_ice_dynbal(g) = &
              (clm3%g%gws%gc_ice2 (g) - clm3%g%gws%gc_ice1 (g))/dtsrf
      clm3%g%gef%eflx_dynbal    (g) = &
              (clm3%g%ges%gc_heat2(g) - clm3%g%ges%gc_heat1(g))/dtsrf
    end do
#endif
#endif

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! ======================================================================
    ! Initialize the mass balance checks: water, carbon, and nitrogen
    ! ======================================================================

    call BeginWaterBalance(begc,endc,             &
                           filter%num_nolakec,    &
                           filter%nolakec,        &
                           filter%num_lakec,      &
                           filter%lakec,          &
                           filter%num_hydrologyc, &
                           filter%hydrologyc)

#if (defined CN)
    call BeginCBalance(begc,endc,filter%num_soilc,filter%soilc)
    call BeginNBalance(begc,endc,filter%num_soilc,filter%soilc)
#endif

    ! =======================================================================
    ! Initialize h2ocan_loss to zero
    ! =======================================================================
    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
    call pftdyn_wbal_init(begc,endc)

#if (defined CNDV)
    ! NOTE: Currently CNDV and DYNPFT are incompatible
    call CNZeroFluxes_dwt(begc,endc,begp,endp)
    call pftwt_interp(begp,endp)
    call reweightWrapup
    call pftdyn_wbal(begc,endc,begp,endp)
    call pftdyn_cnbal(begc,endc,begp,endp)
#else

#if (defined CN)
    call CNZeroFluxes_dwt(begc,endc,begp,endp)
#ifdef DYNPFT
    call pftdyn_cnbal(begc,endc,begp,endp)
#endif
#endif
#endif

#if (defined CN)
    ! =======================================================================
    ! Update dynamic N deposition field, on albedo timestep
    ! =======================================================================
    ! PET: switching CN timestep
    call CNFireInterp()
#endif

    ! =======================================================================
    ! Determine proc boundaries
    ! =======================================================================

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! =======================================================================
    ! Initialize variables from previous time step and
    ! Determine canopy interception and precipitation onto ground surface.
    ! Determine the fraction of foliage covered by water and the fraction
    ! of foliage that is dry and transpiring. Initialize snow layer if the
    ! snow accumulation exceeds 10 mm.
    ! =======================================================================

    ! initialize intracellular CO2 (Pa) parameters each timestep for use
    ! in VOCEmission
    clm3%g%l%c%p%pcf%cisun_z(begp:endp,:) = -999.D0
    clm3%g%l%c%p%pcf%cisha_z(begp:endp,:) = -999.D0

    ! initialize declination for current timestep
    do c = begc , endc
      clm3%g%l%c%cps%decl(c) = declin
    end do

    call clm_driverInit(begc,endc,begp,endp, &
                        filter%num_nolakec,  &
                        filter%nolakec)

    ! =======================================================================
    ! Hydrology1
    ! =======================================================================

    call Hydrology1(begc,endc,begp,endp, &
                    filter%num_nolakec,  &
                    filter%nolakec,      &
                    filter%num_nolakep,  &
                    filter%nolakep)

    ! =======================================================================
    ! Surface Radiation
    ! =======================================================================

    ! Surface Radiation for non-urban columns

    call SurfaceRadiation(begp,endp,           &
                          filter%num_nourbanp, &
                          filter%nourbanp)

    ! Surface Radiation for urban columns

    call UrbanRadiation(begl,endl,begp,endp, &
                        filter%num_nourbanl, &
                        filter%nourbanl,     &
                        filter%num_urbanl,   &
                        filter%urbanl,       &
                        filter%num_urbanp,   &
                        filter%urbanp)

    ! =======================================================================
    ! Determine leaf temperature and surface fluxes based on ground
    ! temperature from previous time step.
    ! =======================================================================

    call Biogeophysics1(begg,endg,begc,endc,begp,endp, &
                        filter%num_nolakec,            &
                        filter%nolakec,                &
                        filter%num_nolakep,            &
                        filter%nolakep)

    ! =======================================================================
    ! Determine bare soil or snow-covered vegetation surface temperature and
    ! fluxes. Calculate Ground fluxes (frac_veg_nosno is either 1 or 0)
    ! =======================================================================

    ! BareGroundFluxes for all pfts except lakes and urban landunits

    call BareGroundFluxes(begp,endp,               &
                          filter%num_nolakeurbanp, &
                          filter%nolakeurbanp)

    ! Fluxes for all Urban landunits

    call UrbanFluxes(begp,endp,begl,endl,begc,endc, &
                     filter%num_nourbanl,           &
                     filter%nourbanl,               &
                     filter%num_urbanl,             &
                     filter%urbanl,                 &
                     filter%num_urbanp,             &
                     filter%urbanp)

    ! =======================================================================
    ! Determine non snow-covered vegetation surface temperature and fluxes
    ! Calculate canopy temperature, latent and sensible fluxes from the
    ! canopy, and leaf water change by evapotranspiration
    ! =======================================================================

    call CanopyFluxes(begc,endc,begp,endp,filter%num_nolakep,  &
                      filter%nolakep)

    ! =======================================================================
    ! Determine lake temperature and surface fluxes
    ! =======================================================================

    call SLakeFluxes(begc,endc,begp,endp, &
                     filter%num_lakep,    &
                     filter%lakep)
    call SLakeTemperature(begc,endc,begp,endp, &
                          filter%num_lakec,    &
                          filter%lakec,        &
                          filter%num_lakep,    &
                          filter%lakep)

    ! =======================================================================
    ! DUST and VOC emissions
    ! =======================================================================

    ! Dust mobilization (C. Zender's modified codes)
    call DustEmission(begp,endp,begl,endl, &
                      filter%num_nolakep,  &
                      filter%nolakep)

    ! Dust dry deposition (C. Zender's modified codes)
    call DustDryDep(begp,endp)

    ! VOC emission (A. Guenther's MEGAN (2006) model)
    call VOCEmission(begp,endp,        &
                     filter%num_soilp, &
                     filter%soilp)

    ! =======================================================================
    ! Determine soil/snow temperatures including ground temperature and
    ! update surface fluxes for new ground temperature.
    ! =======================================================================

    call Biogeophysics2(begl,endl,begc,endc,begp,endp, &
                        filter%num_urbanl,             &
                        filter%urbanl,                 &
                        filter%num_nolakec,            &
                        filter%nolakec,                &
                        filter%num_nolakep,            &
                        filter%nolakep)

    ! =======================================================================
    ! Perform averaging from PFT level to column level
    ! =======================================================================

    call pft2col(begc,endc,filter%num_nolakec,filter%nolakec)

    ! =======================================================================
    ! Vertical (column) soil and surface hydrology
    ! =======================================================================

    call Hydrology2(begc,endc,             &
                    filter%num_nolakec,    &
                    filter%nolakec,        &
                    filter%num_hydrologyc, &
                    filter%hydrologyc,     &
                    filter%num_urbanc,     &
                    filter%urbanc,         &
                    filter%num_snowc,      &
                    filter%snowc,          &
                    filter%num_nosnowc,    &
                    filter%nosnowc)

    ! =======================================================================
    ! Lake hydrology
    ! =======================================================================

    call SLakeHydrology(begc,endc,begp,endp, &
                        filter%num_lakec,    &
                        filter%lakec,        &
                        filter%num_lakep,    &
                        filter%lakep)

    ! =======================================================================
    ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
    ! =======================================================================

    do c = begc , endc
      l = clandunit(c)
      if ( itypelun(l) == isturb ) then
        ! Urban landunit use Bonan 1996 (LSM Technical Note)
        cptr%cps%frac_sno(c) = min( cptr%cps%snow_depth(c)/0.05D0, 1.D0)
      end if
    end do

    ! =======================================================================
    ! Snow aging routine based on Flanner and Zender (2006), Linking snowpack
    ! microphysics and albedo evolution, JGR, and Brun (1989), Investigation
    ! of wet-snow metamorphism in respect of liquid-water content,
    ! Ann. Glaciol.
    ! =======================================================================
    ! Note the snow filters here do not include lakes; SnowAge_grain is called
    ! for lakes from SLakeHydrology.

    call SnowAge_grain(begc,endc,          &
                       filter%num_snowc,   &
                       filter%snowc,       &
                       filter%num_nosnowc, &
                       filter%nosnowc)

    ! =======================================================================
    ! Ecosystem dynamics: Uses CN, CNDV, or static parameterizations
    ! =======================================================================

#if (defined CN)
    ! fully prognostic canopy structure and C-N biogeochemistry
    ! - CNDV defined: prognostic biogeography; else prescribed
    ! - crop model:   crop algorithms called from within CNEcosystemDyn
    call CNEcosystemDyn(begc,endc,begp,endp, &
                        filter%num_soilc,    &
                        filter%soilc,        &
                        filter%num_soilp,    &
                        filter%soilp,        &
                        filter%num_pcropp,   &
                        filter%pcropp, doalb)
    call CNAnnualUpdate(begc,endc,begp,endp, &
                        filter%num_soilc,    &
                        filter%soilc,        &
                        filter%num_soilp,    &
                        filter%soilp)
#else
    ! Prescribed biogeography,
    ! prescribed canopy structure, some prognostic carbon fluxes
    call EcosystemDyn(begp,endp,          &
                      filter%num_nolakep, &
                      filter%nolakep, doalb)
#endif

    ! Dry Deposition of chemical tracers (Wesely (1998) parameterizaion)
    call depvel_compute(begp,endp)

    ! =======================================================================
    ! Check the energy and water balance, also carbon and nitrogen balance
    ! =======================================================================

    call BalanceCheck(begp,endp,begc,endc,begl,endl,begg,endg)

#if (defined CN)
!    print*,'FAB', lfirstcall
!    if (ktau < 2) then
    if ( ktau < 2 .or. .not.lfirstcall) then
      if ( myid == italk ) then
        write(stdout,*) &
          '--WARNING-- skipping CN balance check for first timestep'
      end if
    lfirstcall = .true.
    else
      call CBalanceCheck(begc,endc,filter%num_soilc,filter%soilc)
      call NBalanceCheck(begc,endc,filter%num_soilc,filter%soilc)
    end if
#endif

    ! CH4
#if (defined LCH4)
    call ch4(begg,endg,begl,endl,begc,endc,begp,endp, &
             filter%num_soilc, filter%soilc,          &
             filter%num_lakec, filter%lakec,          &
             filter%num_soilp, filter%soilp)
#endif

    ! =======================================================================
    ! Determine albedos for next time step
    ! =======================================================================

    if ( doalb ) then

      ! Albedos for non-urban columns

      call SurfaceAlbedo(begg,endg,begc,endc,begp,endp, &
                         filter%num_nourbanc,           &
                         filter%nourbanc,               &
                         filter%num_nourbanp,           &
                         filter%nourbanp,               &
                         nextsw_cday, declinp1)

      ! Albedos for urban columns

      if ( filter%num_urbanl > 0 ) then
         call UrbanAlbedo(begl,endl,begc,endc,begp,endp, &
                          filter%num_urbanl,             &
                          filter%urbanl,                 &
                          filter%num_urbanc,             &
                          filter%urbanc,                 &
                          filter%num_urbanp,             &
                          filter%urbanp)
      end if
    end if

    ! =======================================================================
    ! Determine gridcell averaged properties to send to atm
    !  (l2as and l2af derived types)
    ! =======================================================================

    call clm_map2gcell( )

    ! =======================================================================
    ! Update accumulators
    ! =======================================================================

    call updateAccFlds()

    ! =======================================================================
    ! Update history buffer
    ! =======================================================================

    call hist_update_hbuf()

    ! =======================================================================
    ! Call dv (dynamic vegetation) at last time step of year
    ! NOTE: monp1, dayp1, and secp1 correspond to ktau+1
    ! =======================================================================

#if (defined CNDV)
    if ( ktau > 0 ) then
      tdif = int(dtsrf + dt/2.0D0)
    else
      tdif = int(dtsrf)
    end if
    nextdate = idatex + tdif
    if ( date_is(nextdate,1,1) .and. &
         time_is(nextdate,0,dtsrf) .and. ktau > 0 ) then
      call split_idate(idatex,yr,mon,day)
      ncdate = yr*10000 + mon*100 + day
      call split_idate(idate0,yr,mon,day)
      nbdate = yr*10000 + mon*100 + day
      kyr = ncdate/10000 - nbdate/10000 + 1
      if ( myid == italk ) then
        write(stdout,*) 'End of year. CNDV called now: ncdate=', &
                       ncdate,' nbdate=',nbdate,' kyr=',kyr,' ktau=', ktau
      end if
      call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
      call dv(begg,endg,begp,endp, &
              filter%num_natvegp,  &
              filter%natvegp, kyr)
    end if
#endif

    ! =======================================================================
    ! Create history and write history tapes if appropriate
    ! =======================================================================

#ifndef _NOIO

    call hist_htapes_wrapup(rstwr,nlend,nlomon)

    ! =======================================================================
    ! Write to CNDV history buffer if appropriate
    ! =======================================================================

#if (defined CNDV)
    if ( date_is(nextdate,1,1) .and. time_is(nextdate,0,dtsrf) .and. &
         ktau > 0 )  then
      call histCNDV()
      if (myid == italk) then
        write(stdout,*) 'Annual CNDV calculations are complete'
      end if
    end if
#endif

    ! =======================================================================
    ! Write restart/initial files if appropriate
    ! =======================================================================

    if ( rstwr ) then
      filer = restFile_filename(rdate=rdate)
      call restFile_write( filer, rdate=rdate )
    end if
#endif

  end subroutine clm_drv

end module mod_clm_driver
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
