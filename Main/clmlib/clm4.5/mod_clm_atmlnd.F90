module mod_clm_atmlnd
  !
  ! Handle atm2lnd, lnd2atm mapping
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_varpar, only : numrad, ndst, nlevsoi
  use mod_clm_varctl, only : use_c13
  use mod_clm_decomp, only : get_proc_bounds

  use mod_clm_drydep, only : n_drydep, drydep_method, DD_XLND
  use mod_clm_megan, only : shr_megan_mechcomps_n

  implicit none

  private

  save

  type atm_domain
    real(rk8), pointer, contiguous, dimension(:) :: xlat => null()
    real(rk8), pointer, contiguous, dimension(:) :: xlon => null()
    real(rk8), pointer, contiguous, dimension(:) :: topo => null()
    real(rk8), pointer, contiguous, dimension(:) :: area => null()
    integer(ik4), pointer, contiguous, dimension(:) :: iveg => null()
    integer(ik4), pointer, contiguous, dimension(:) :: itex => null()
    real(rk8), pointer, contiguous, dimension(:) :: snow => null()
    real(rk8), pointer, contiguous, dimension(:) :: smoist => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: rmoist => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: rts => null()
    real(rk8), pointer, contiguous, dimension(:) :: tgrd => null()
  end type atm_domain

  public :: atm_domain

  !----------------------------------------------------
  ! atmosphere -> land variables structure
  !----------------------------------------------------
  type atm2lnd_type
    !atmospheric temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: forc_t => null()
    !atm wind speed, east direction (m/s)
    real(rk8), pointer, contiguous, dimension(:) :: forc_u => null()
    !atm wind speed, north direction (m/s)
    real(rk8), pointer, contiguous, dimension(:) :: forc_v => null()
    !atmospheric wind speed
    real(rk8), pointer, contiguous, dimension(:) :: forc_wind => null()
    !atmospheric specific humidity (kg/kg)
    real(rk8), pointer, contiguous, dimension(:) :: forc_q => null()
    !atmospheric reference height (m)
    real(rk8), pointer, contiguous, dimension(:) :: forc_hgt => null()
    !obs height of wind [m] (new)
    real(rk8), pointer, contiguous, dimension(:) :: forc_hgt_u => null()
    !obs height of temperature [m] (new)
    real(rk8), pointer, contiguous, dimension(:) :: forc_hgt_t => null()
    !obs height of humidity [m] (new)
    real(rk8), pointer, contiguous, dimension(:) :: forc_hgt_q => null()
    !atmospheric pressure (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: forc_pbot => null()
    !atm potential temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: forc_th => null()
    !density (kg/m**3)
    real(rk8), pointer, contiguous, dimension(:) :: forc_rho => null()
    !atmospheric relative humidity (%)
    real(rk8), pointer, contiguous, dimension(:) :: forc_rh => null()
    !surface pressure (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: forc_psrf => null()
    !CO2 partial pressure (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: forc_pco2 => null()
    !downwrd IR longwave radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: forc_lwrad => null()
    !direct beam radiation (numrad) (vis=forc_sols, nir=forc_soll )
    real(rk8), pointer, contiguous, dimension(:,:) :: forc_solad => null()
    !diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
    real(rk8), pointer, contiguous, dimension(:,:) :: forc_solai => null()
    !incident solar radiation
    real(rk8), pointer, contiguous, dimension(:) :: forc_solar => null()
    !rain rate [mm/s]
    real(rk8), pointer, contiguous, dimension(:) :: forc_rain => null()
    !snow rate [mm/s]
    real(rk8), pointer, contiguous, dimension(:) :: forc_snow => null()
    !nitrogen deposition rate (gN/m2/s)
    real(rk8), pointer, contiguous, dimension(:) :: forc_ndep => null()
    !ALMA rain+snow [mm/s]
    real(rk8), pointer, contiguous, dimension(:) :: rainf => null()
    !C13O2 partial pressure (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: forc_pc13o2 => null()
    !O2 partial pressure (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: forc_po2 => null()
    ! rof flood (mm/s)
    real(rk8), pointer, contiguous, dimension(:) :: forc_flood => null()
    ! rof volr (m3)
    real(rk8), pointer, contiguous, dimension(:) :: volr => null()
    ! aerosol deposition array
    real(rk8), pointer, contiguous, dimension(:,:) :: forc_aer => null()
#ifdef LCH4
    !CH4 partial pressure (Pa)
    real(rk8), pointer, contiguous, dimension(:) :: forc_pch4
#endif
    real(rk8), pointer, contiguous, dimension(:) :: notused => null()
  end type atm2lnd_type

  public :: atm2lnd_type

  !----------------------------------------------------
  ! land -> atmosphere variables structure
  !----------------------------------------------------
  type lnd2atm_type
    !radiative temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_rad => null()
    !vegetation temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_veg => null()
    !2m surface air temperature (Kelvin)
    real(rk8), pointer, contiguous, dimension(:) :: t_ref2m => null()
    !2m surface specific humidity (kg/kg)
    real(rk8), pointer, contiguous, dimension(:) :: q_ref2m => null()
    !10m surface wind speed (m/sec)
    real(rk8), pointer, contiguous, dimension(:) :: u_ref10m => null()
    !snow water (mm H2O)
    real(rk8), pointer, contiguous, dimension(:) :: h2osno => null()
    !(numrad) surface albedo (direct)
    real(rk8), pointer, contiguous, dimension(:,:) :: albd => null()
    !(numrad) surface albedo (diffuse)
    real(rk8), pointer, contiguous, dimension(:,:) :: albi => null()
    !wind stress: e-w (kg/m/s**2)
    real(rk8), pointer, contiguous, dimension(:) :: taux => null()
    !wind stress: n-s (kg/m/s**2)
    real(rk8), pointer, contiguous, dimension(:) :: tauy => null()
    !roughness length over vegetation, momentum
    real(rk8), pointer, contiguous, dimension(:) :: zom => null()
    !roughness length over vegetation, heat
    real(rk8), pointer, contiguous, dimension(:) :: zoh => null()
    !net ground heat flux into ground (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_gnet => null()
    !total latent HF (W/m**2)  [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lh_tot => null()
    !total sensible HF (W/m**2) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: eflx_sh_tot => null()
    !IR (longwave) radiation (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: eflx_lwrad_out => null()
    !qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(rk8), pointer, contiguous, dimension(:) :: qflx_evap_tot => null()
    !solar rad absorbed (total) (W/m**2)
    real(rk8), pointer, contiguous, dimension(:) :: fsa => null()
    !net CO2 flux (kg CO2/m**2/s) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: nee => null()
    !aerodynamical resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: ram1 => null()
    !thermal resistance (s/m)
    real(rk8), pointer, contiguous, dimension(:) :: rah1 => null()
    !bulk Richardson number
    real(rk8), pointer, contiguous, dimension(:) :: br1 => null()
    !friction velocity (m/s) (for dust model)
    real(rk8), pointer, contiguous, dimension(:) :: fv => null()
    !Surface ground emissivity
    real(rk8), pointer, contiguous, dimension(:) :: emg => null()
    !vegrtation emissivity
    real(rk8), pointer, contiguous, dimension(:) :: emv => null()
    !fraction of ground emittimg dust (not vegetated and not snow covered)
    real(rk8), pointer, contiguous, dimension(:) :: vdustfrac => null()
    ! soil water kg/m^2 (water + ice)
    real(rk8), pointer, contiguous, dimension(:,:) :: h2osoi => null()
    ! soil water in first 10 cm
    real(rk8), pointer, contiguous, dimension(:) :: h2o10cm => null()
    ! FAB  soil volumetric water content (m3/m3)
    real(rk8), pointer, contiguous, dimension(:,:) :: h2osoi_vol => null()
    ! soil/snow temperaure profils
    real(rk8), pointer, contiguous, dimension(:,:) :: tsoi => null()
    ! Surface runoff
    real(rk8), pointer, contiguous, dimension(:) :: qflx_surf => null()
    ! Total runoff
    real(rk8), pointer, contiguous, dimension(:) :: qflx_tot => null()
    ! Snow melt
    real(rk8), pointer, contiguous, dimension(:) :: qflx_snow_melt => null()
    ! rof liq forcing
    real(rk8), pointer, contiguous, dimension(:) :: rofliq => null()
    ! rof ice forcing
    real(rk8), pointer, contiguous, dimension(:) :: rofice => null()
    !dust flux (size bins)
    real(rk8), pointer, contiguous, dimension(:,:) :: flxdst => null()
    !dry deposition velocities
    real(rk8), pointer, contiguous, dimension(:,:) :: ddvel => null()
    ! VOC flux (size bins)
    real(rk8), pointer, contiguous, dimension(:,:) :: flxvoc => null()
    ! Total leaf area index at grid level
    real(rk8), pointer, contiguous, dimension(:) :: tlai => null()
#ifdef LCH4
    !net CH4 flux (kg C/m**2/s) [+ to atm]
    real(rk8), pointer, contiguous, dimension(:) :: flux_ch4
#endif
    real(rk8), pointer, contiguous, dimension(:) :: notused => null()
  end type lnd2atm_type

  public :: lnd2atm_type

  type(atm2lnd_type), public, target :: clm_a2l      ! a2l fields on clm grid
  type(lnd2atm_type), public, target :: clm_l2a      ! l2a fields on clm grid

  type(atm_domain), public :: adomain

  public :: init_atm2lnd_type
  public :: init_lnd2atm_type
  public :: clm_map2gcell

  contains
  !
  ! Initialize atmospheric variables required by the land
  !
  subroutine init_atm2lnd_type(ibeg, iend, a2l)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (atm2lnd_type), intent(inout):: a2l
    real(rk8) :: ival   ! initial value

    allocate(a2l%forc_t(ibeg:iend))
    allocate(a2l%forc_u(ibeg:iend))
    allocate(a2l%forc_v(ibeg:iend))
    allocate(a2l%forc_wind(ibeg:iend))
    allocate(a2l%forc_q(ibeg:iend))
    allocate(a2l%forc_rh(ibeg:iend))
    allocate(a2l%forc_hgt(ibeg:iend))
    allocate(a2l%forc_hgt_u(ibeg:iend))
    allocate(a2l%forc_hgt_t(ibeg:iend))
    allocate(a2l%forc_hgt_q(ibeg:iend))
    allocate(a2l%forc_pbot(ibeg:iend))
    allocate(a2l%forc_th(ibeg:iend))
    allocate(a2l%forc_rho(ibeg:iend))
    allocate(a2l%forc_psrf(ibeg:iend))
    allocate(a2l%forc_pco2(ibeg:iend))
    allocate(a2l%forc_lwrad(ibeg:iend))
    allocate(a2l%forc_solad(ibeg:iend,numrad))
    allocate(a2l%forc_solai(ibeg:iend,numrad))
    allocate(a2l%forc_solar(ibeg:iend))
    allocate(a2l%forc_rain(ibeg:iend))
    allocate(a2l%forc_snow(ibeg:iend))
    allocate(a2l%forc_ndep(ibeg:iend))
    allocate(a2l%rainf(ibeg:iend))
    if ( use_c13 ) then
      allocate(a2l%forc_pc13o2(ibeg:iend))
    end if
    allocate(a2l%forc_po2(ibeg:iend))
    allocate(a2l%forc_flood(ibeg:iend))
    allocate(a2l%volr(ibeg:iend))
    allocate(a2l%forc_aer(ibeg:iend,14))
#if (defined LCH4)
    allocate(a2l%forc_pch4(ibeg:iend))
#endif
    allocate(a2l%notused(ibeg:iend))

    ! ival = nan      ! causes core dump in map_maparray, tcx fix
    ival = 0.0_rk8

    a2l%forc_t(ibeg:iend) = ival
    a2l%forc_u(ibeg:iend) = ival
    a2l%forc_v(ibeg:iend) = ival
    a2l%forc_wind(ibeg:iend) = ival
    a2l%forc_q(ibeg:iend) = ival
    a2l%forc_rh(ibeg:iend) = ival
    a2l%forc_hgt(ibeg:iend) = ival
    a2l%forc_hgt_u(ibeg:iend) = ival
    a2l%forc_hgt_t(ibeg:iend) = ival
    a2l%forc_hgt_q(ibeg:iend) = ival
    a2l%forc_pbot(ibeg:iend) = ival
    a2l%forc_th(ibeg:iend) = ival
    a2l%forc_rho(ibeg:iend) = ival
    a2l%forc_psrf(ibeg:iend) = ival
    a2l%forc_pco2(ibeg:iend) = ival
    a2l%forc_lwrad(ibeg:iend) = ival
    a2l%forc_solad(ibeg:iend,1:numrad) = ival
    a2l%forc_solai(ibeg:iend,1:numrad) = ival
    a2l%forc_solar(ibeg:iend) = ival
    a2l%forc_rain(ibeg:iend) = ival
    a2l%forc_snow(ibeg:iend) = ival
    a2l%forc_ndep(ibeg:iend) = ival
    a2l%rainf(ibeg:iend) = nan
    if ( use_c13 ) then
      a2l%forc_pc13o2(ibeg:iend) = ival
    end if
    a2l%forc_po2(ibeg:iend) = ival
    a2l%forc_flood(ibeg:iend) = ival
    a2l%volr(ibeg:iend) = ival
    a2l%forc_aer(ibeg:iend,:) = ival
#if (defined LCH4)
    a2l%forc_pch4(ibeg:iend) = ival
#endif
end subroutine init_atm2lnd_type
  !
  ! Initialize land variables required by the atmosphere
  !
  subroutine init_lnd2atm_type(ibeg,iend,l2a)
    implicit none
    integer(ik4), intent(in) :: ibeg, iend
    type (lnd2atm_type), intent(inout) :: l2a
    real(rk8) :: ival   ! initial value

    allocate(l2a%t_rad(ibeg:iend))
    allocate(l2a%t_veg(ibeg:iend))
    allocate(l2a%t_ref2m(ibeg:iend))
    allocate(l2a%q_ref2m(ibeg:iend))
    allocate(l2a%u_ref10m(ibeg:iend))
    allocate(l2a%h2osno(ibeg:iend))
    allocate(l2a%albd(ibeg:iend,1:numrad))
    allocate(l2a%albi(ibeg:iend,1:numrad))
    allocate(l2a%taux(ibeg:iend))
    allocate(l2a%tauy(ibeg:iend))
    allocate(l2a%zom(ibeg:iend))
    allocate(l2a%zoh(ibeg:iend))
    allocate(l2a%eflx_lwrad_out(ibeg:iend))
    allocate(l2a%eflx_gnet(ibeg:iend))
    allocate(l2a%eflx_sh_tot(ibeg:iend))
    allocate(l2a%eflx_lh_tot(ibeg:iend))
    allocate(l2a%qflx_evap_tot(ibeg:iend))
    allocate(l2a%fsa(ibeg:iend))
    allocate(l2a%nee(ibeg:iend))
    allocate(l2a%ram1(ibeg:iend))
    allocate(l2a%rah1(ibeg:iend))
    allocate(l2a%br1(ibeg:iend))
    allocate(l2a%fv(ibeg:iend))
    allocate(l2a%emg(ibeg:iend))
    allocate(l2a%emv(ibeg:iend))
    allocate(l2a%vdustfrac(ibeg:iend))
    allocate(l2a%h2osoi(ibeg:iend,nlevsoi))
    allocate(l2a%tsoi(ibeg:iend,nlevsoi))
    allocate(l2a%h2osoi_vol(ibeg:iend,nlevsoi))
    allocate(l2a%h2o10cm(ibeg:iend))
    allocate(l2a%qflx_surf(ibeg:iend))
    allocate(l2a%qflx_tot(ibeg:iend))
    allocate(l2a%qflx_snow_melt(ibeg:iend))
    allocate(l2a%rofliq(ibeg:iend))
    allocate(l2a%rofice(ibeg:iend))
    allocate(l2a%tlai(ibeg:iend))
    allocate(l2a%flxdst(ibeg:iend,1:ndst))
#if (defined LCH4)
    allocate(l2a%flux_ch4(ibeg:iend))
#endif
    if ( shr_megan_mechcomps_n > 0 ) then
      allocate(l2a%flxvoc(ibeg:iend,1:shr_megan_mechcomps_n))
    end if
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
      allocate(l2a%ddvel(ibeg:iend,1:n_drydep))
    end if
    allocate(l2a%notused(ibeg:iend))

    ! ival = nan   ! causes core dump in map_maparray, tcx fix
    ival = 0.0_rk8

    l2a%t_rad(ibeg:iend) = ival
    l2a%t_veg(ibeg:iend) = ival
    l2a%t_ref2m(ibeg:iend) = ival
    l2a%q_ref2m(ibeg:iend) = ival
    l2a%u_ref10m(ibeg:iend) = ival
    l2a%h2osno(ibeg:iend) = ival
    l2a%albd(ibeg:iend,1:numrad) = ival
    l2a%albi(ibeg:iend,1:numrad) = ival
    l2a%taux(ibeg:iend) = ival
    l2a%tauy(ibeg:iend) = ival
    l2a%zom(ibeg:iend) = ival
    l2a%zoh(ibeg:iend) = ival
    l2a%eflx_lwrad_out(ibeg:iend) = ival
    l2a%eflx_gnet(ibeg:iend) = ival
    l2a%eflx_sh_tot(ibeg:iend) = ival
    l2a%eflx_lh_tot(ibeg:iend) = ival
    l2a%qflx_evap_tot(ibeg:iend) = ival
    l2a%fsa(ibeg:iend) = ival
    l2a%nee(ibeg:iend) = ival
    l2a%ram1(ibeg:iend) = ival
    l2a%rah1(ibeg:iend) = ival
    l2a%br1(ibeg:iend) = ival
    l2a%fv(ibeg:iend) = ival
    l2a%h2osoi(ibeg:iend,:) = ival
    l2a%emg(ibeg:iend) = ival
    l2a%emv(ibeg:iend) = ival
    l2a%vdustfrac(ibeg:iend) = ival
    l2a%tsoi(ibeg:iend,:) = ival
    l2a%h2o10cm(ibeg:iend) = ival
    l2a%qflx_surf(ibeg:iend) = ival
    l2a%qflx_tot(ibeg:iend) = ival
    l2a%qflx_snow_melt(ibeg:iend) = ival
    l2a%rofliq(ibeg:iend) = ival
    l2a%rofice(ibeg:iend) = ival
    l2a%tlai(ibeg:iend) = ival
    l2a%flxdst(ibeg:iend,1:ndst) = ival
    if ( shr_megan_mechcomps_n > 0 ) then
      l2a%flxvoc(ibeg:iend,1:shr_megan_mechcomps_n) = ival
    end if
#if (defined LCH4)
    l2a%flux_ch4(ibeg:iend) = ival
#endif
    if ( n_drydep > 0 .and. drydep_method == DD_XLND )then
      l2a%ddvel(ibeg:iend, : ) = ival
    end if
  end subroutine init_lnd2atm_type
  !
  ! Compute l2a component of gridcell derived type
  !
  subroutine clm_map2gcell(init)
    use mod_clm_type
    use mod_clm_subgridave
    use mod_clm_varcon , only : sb
    use mod_clm_varpar , only : numrad
#ifdef LCH4
    use mod_clm_ch4varcon  , only : ch4offline
#endif
    implicit none
    save
    ! if true=>only set a subset of arguments
    logical, optional, intent(in) :: init
    ! per-proc beginning and ending pft indices
    integer(ik4) :: begp, endp
    ! per-proc beginning and ending column indices
    integer(ik4) :: begc, endc
    ! per-proc beginning and ending landunit indices
    integer(ik4) :: begl, endl
    ! per-proc gridcell ending gridcell indices
    integer(ik4) :: begg, endg
    integer(ik4) :: g ! indices
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(column_type), pointer :: cptr    ! pointer to column derived subtype
    type(pft_type), pointer :: pptr       ! pointer to pft derived subtype
    real(rk8), parameter :: amC = 12.0_rk8  ! Atomic mass number for Carbon
    real(rk8), parameter :: amO = 16.0_rk8  ! Atomic mass number for Oxygen
    ! Atomic mass number for CO2
    real(rk8), parameter :: amCO2 = amC + 2.0_rk8*amO
    ! The following converts g of C to kg of CO2
    real(rk8), parameter :: convertgC2kgCO2 = 1.0e-3_rk8 * (amCO2/amC)
    real(rk8), allocatable, dimension(:,:) :: tmpc
    real(rk8), allocatable, dimension(:) :: tmp

    ! Set pointers into derived type

    gptr => clm3%g
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Determine processor bounds

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    allocate(tmpc(begc:endc,nlevsoi))
    allocate(tmp(begp:endp))
    !$acc data create(tmpc,tmp)

    ! Compute gridcell averages.

    if ( present(init) ) then
      call c2g(begc,endc,begl,endl,begg,endg,  &
               cptr%cws%h2osno,clm_l2a%h2osno, &
               c2l_scale_type='urbanf',        &
               l2g_scale_type='unity')
      !$acc kernels
      tmpc(:,:) = cptr%cws%h2osoi_liq(:,1:nlevsoi) + &
                  cptr%cws%h2osoi_ice(:,1:nlevsoi)
      !$acc end kernels
      call c2g(begc,endc,begl,endl,begg,endg,nlevsoi, &
               tmpc,clm_l2a%h2osoi, &
               c2l_scale_type='unity', &
               l2g_scale_type='unity')
      call c2g(begc,endc,begl,endl,begg,endg, &
               cptr%cws%h2osoi_liqice_10cm,clm_l2a%h2o10cm, &
               c2l_scale_type='unity', &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,numrad,  &
               pptr%pps%albd,clm_l2a%albd,                      &
               p2c_scale_type='unity',                          &
               c2l_scale_type='urbanf',                         &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,numrad, &
               pptr%pps%albi,clm_l2a%albi,                     &
               p2c_scale_type='unity',                         &
               c2l_scale_type='urbanf',                        &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,        &
               pptr%pef%eflx_lwrad_out,clm_l2a%eflx_lwrad_out, &
               p2c_scale_type='unity',                         &
               c2l_scale_type='urbanf',                        &
               l2g_scale_type='unity')
      do concurrent ( g = begg:endg )
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
      end do
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pps%tlai,clm_l2a%tlai,              &
               p2c_scale_type='unity',                  &
               c2l_scale_type='urbanf',                 &
               l2g_scale_type='unity')
    else
      call c2g(begc,endc,begl,endl,begg,endg,  &
               cptr%cws%h2osno,clm_l2a%h2osno, &
               c2l_scale_type='urbanf',        &
               l2g_scale_type='unity')
      !$acc kernels
      tmpc(:,:) = cptr%cws%h2osoi_liq(:,1:nlevsoi) + &
                  cptr%cws%h2osoi_ice(:,1:nlevsoi)
      !$acc end kernels
      call c2g(begc,endc,begl,endl,begg,endg,nlevsoi, &
               tmpc,clm_l2a%h2osoi, &
               c2l_scale_type='unity', &
               l2g_scale_type='unity')
      !FAB
      !$acc kernels
      tmpc(:,:) = cptr%cws%h2osoi_vol(:,1:nlevsoi)
      !$acc end kernels
      call c2g(begc,endc,begl,endl,begg,endg,nlevsoi, &
               tmpc,clm_l2a%h2osoi_vol, &
               c2l_scale_type='unity', &
               l2g_scale_type='unity')
      !$acc kernels
      tmpc(:,:) = cptr%ces%t_soisno(:,1:nlevsoi)
      !$acc end kernels
      call c2g(begc,endc,begl,endl,begg,endg,nlevsoi, &
               tmpc,clm_l2a%tsoi, &
               c2l_scale_type='unity', &
               l2g_scale_type='unity')
      call c2g(begc,endc,begl,endl,begg,endg, &
               cptr%cws%h2osoi_liqice_10cm,clm_l2a%h2o10cm, &
               c2l_scale_type='unity', &
               l2g_scale_type='unity')
      call c2g(begc,endc,begl,endl,begg,endg, &
               cptr%cps%emg,clm_l2a%emg,      &
               c2l_scale_type='unity',        &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,       &
               pptr%pps%emv,clm_l2a%emv,                      &
               p2c_scale_type='unity',c2l_scale_type='unity', &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,       &
               pptr%pdf%lnd_frc_mbl_dst,clm_l2a%vdustfrac,    &
               p2c_scale_type='unity',c2l_scale_type='unity', &
               l2g_scale_type='unity')
      call c2g(begc,endc,begl,endl,begg,endg,        &
               cptr%cwf%qflx_surf,clm_l2a%qflx_surf, &
               c2l_scale_type='unity',               &
               l2g_scale_type='unity')
      call c2g(begc,endc,begl,endl,begg,endg,         &
               cptr%cwf%qflx_runoff,clm_l2a%qflx_tot, &
               c2l_scale_type='unity',                &
               l2g_scale_type='unity')
      call c2g(begc,endc,begl,endl,begg,endg,                  &
               cptr%cwf%qflx_snow_melt,clm_l2a%qflx_snow_melt, &
               c2l_scale_type='unity',                         &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,numrad, &
               pptr%pps%albd,clm_l2a%albd,                     &
               p2c_scale_type='unity',                         &
               c2l_scale_type='urbanf',                        &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,numrad, &
               pptr%pps%albi,clm_l2a%albi,                     &
               p2c_scale_type='unity',                         &
               c2l_scale_type='urbanf',                        &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pes%t_ref2m,clm_l2a%t_ref2m,        &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pes%t_veg,clm_l2a%t_veg,            &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pes%q_ref2m,clm_l2a%q_ref2m,        &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pps%u10_clm,clm_l2a%u_ref10m,       &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pmf%taux,clm_l2a%taux,              &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pmf%tauy,clm_l2a%tauy,              &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
!FAB: for roughness lenght perform a ln averaging instead of linear
      !$acc kernels
      tmp(:) = pptr%pps%z0mv
      !$acc end kernels
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               tmp,clm_l2a%zom,               &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      !clm_l2a%zom = exp(clm_l2a%zom)
      !tmp(:) = log(pptr%pps%z0hv)
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               tmp,clm_l2a%zoh,               &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      !clm_l2a%zoh = exp(clm_l2a%zoh)
!
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,  &
               pptr%pef%eflx_gnet,clm_l2a%eflx_gnet,     &
               p2c_scale_type='unity',                   &
               c2l_scale_type='urbanf',                  &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,  &
               pptr%pef%eflx_lh_tot,clm_l2a%eflx_lh_tot, &
               p2c_scale_type='unity',                   &
               c2l_scale_type='urbanf',                  &
               l2g_scale_type='unity')
      do concurrent ( g = begg:endg )
        clm_l2a%eflx_sh_tot(g) = gptr%gef%eflx_sh_totg(g)
      end do
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,      &
               pptr%pwf%qflx_evap_tot,clm_l2a%qflx_evap_tot, &
               p2c_scale_type='unity',                       &
               c2l_scale_type='urbanf',                      &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pef%fsa,clm_l2a%fsa,                &
               p2c_scale_type='unity',                  &
               c2l_scale_type='urbanf',                 &
               l2g_scale_type='unity')
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,        &
               pptr%pef%eflx_lwrad_out,clm_l2a%eflx_lwrad_out, &
               p2c_scale_type='unity',                         &
               c2l_scale_type='urbanf',                        &
               l2g_scale_type='unity')
#if (defined CN)
      call c2g(begc,endc,begl,endl,begg,endg, &
               cptr%ccf%nee,clm_l2a%nee,      &
               c2l_scale_type='unity',        &
               l2g_scale_type='unity')
#else
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pcf%fco2,clm_l2a%nee,               &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      ! Note that fco2 in is umolC/m2/sec so units need to be
      ! changed to gC/m2/sec
      do concurrent ( g = begg:endg )
        clm_l2a%nee(g) = clm_l2a%nee(g)*12.011e-6_rk8
      end do
#endif

#if (defined LCH4)
      if ( .not. ch4offline ) then
        ! Adjust flux of CO2 by the net conversion of mineralizing C to CH4
        do concurrent ( g = begg:endg )
          ! nem is in g C/m2/sec nem is calculated in ch4Mod
          ! flux_ch4 is averaged there also.
          clm_l2a%nee(g) = clm_l2a%nee(g) + gptr%gch4%nem(g)
        end do
      end if
#endif
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pps%fv,clm_l2a%fv,                  &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
!FAB for resistance perform conductance linear average
      !$acc kernels
      tmp(:) = 1._rkx/pptr%pps%ram1
      !$acc end kernels
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               tmp,clm_l2a%ram1,              &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      !$acc kernels
      clm_l2a%ram1 = 1._rkx/clm_l2a%ram1
      tmp(:) = 1._rkx/pptr%pps%rah1
      !$acc end kernels
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               tmp,clm_l2a%rah1,              &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      !$acc kernels
      clm_l2a%rah1 = 1._rkx/clm_l2a%rah1
      !$acc end kernels
!
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pps%br1,clm_l2a%br1,                &
               p2c_scale_type='unity',                  &
               c2l_scale_type='unity',                  &
               l2g_scale_type='unity')
      do concurrent ( g = begg:endg )
        clm_l2a%rofliq(g) = gptr%gwf%qflx_runoffg(g)
        clm_l2a%rofice(g) = gptr%gwf%qflx_snwcp_iceg(g)
      end do
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg,ndst, &
               pptr%pdf%flx_mss_vrt_dst,clm_l2a%flxdst,      &
               p2c_scale_type='unity',                       &
               c2l_scale_type='unity',                       &
               l2g_scale_type='unity')
      if ( shr_megan_mechcomps_n > 0 ) then
        call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
                 shr_megan_mechcomps_n,                   &
                 pptr%pvf%vocflx,clm_l2a%flxvoc,          &
                 p2c_scale_type='unity',                  &
                 c2l_scale_type='unity',                  &
                 l2g_scale_type='unity')
      end if
      if ( n_drydep > 0 .and. drydep_method == DD_XLND ) then
        call p2g(begp,endp,begc,endc,begl,endl,begg,endg,n_drydep, &
                 pptr%pdd%drydepvel,clm_l2a%ddvel,                 &
                 p2c_scale_type='unity',                           &
                 c2l_scale_type='unity',                           &
                 l2g_scale_type='unity')
      end if
      ! Convert from gC/m2/s to kgCO2/m2/s
      do concurrent ( g = begg:endg )
        clm_l2a%nee(g) = clm_l2a%nee(g)*convertgC2kgCO2
      end do
      do concurrent ( g = begg:endg )
        clm_l2a%t_rad(g) = sqrt(sqrt(clm_l2a%eflx_lwrad_out(g)/sb))
      end do
      call p2g(begp,endp,begc,endc,begl,endl,begg,endg, &
               pptr%pps%tlai,clm_l2a%tlai,              &
               p2c_scale_type='unity',                  &
               c2l_scale_type='urbanf',                 &
               l2g_scale_type='unity')
    end if
    !$acc end data
    deallocate(tmpc)
    deallocate(tmp)
  end subroutine clm_map2gcell

end module mod_clm_atmlnd
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
