module mod_clm_baregroundfluxes
  !
  ! Compute sensible and latent fluxes and their derivatives with respect
  ! to ground temperature using ground temperatures from previous time step.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_clm_type
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_varpar , only : nlevgrnd
  use mod_clm_varcon , only : cpair , vkc , grav , denice , denh2o , istsoil
  use mod_clm_varcon , only : istcrop , rgas
  use mod_clm_frictionvelocity, only : FrictionVelocity , MoninObukIni
  use mod_clm_qsat , only : QSat
  use mod_clm_varctl , only : use_c13 , use_c14

  implicit none

  private

  save

  public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes

  contains
  !
  ! Compute sensible and latent fluxes and their derivatives with respect
  ! to ground temperature using ground temperatures from previous time step.
  !
  subroutine BareGroundFluxes(lbp, ubp, num_nolakep, filter_nolakep)
    implicit none
    ! pft bounds
    integer(ik4), intent(in) :: lbp, ubp
    ! number of pft non-lake points in pft filter
    integer(ik4), intent(in) :: num_nolakep
    ! pft filter for non-lake points
    integer(ik4), intent(in) :: filter_nolakep(ubp-lbp+1)

    real(rkx), pointer :: t_soisno(:,:)   ! soil temperature (Kelvin)
    integer(ik4) , pointer :: snl(:)      ! number of snow layers
    real(rkx), pointer :: t_h2osfc(:)     ! surface water temperature
    ! sensible heat flux from snow (W/m**2) [+ to atm]
    real(rkx), pointer :: eflx_sh_snow(:)
    ! sensible heat flux from soil (W/m**2) [+ to atm]
    real(rkx), pointer :: eflx_sh_soil(:)
    ! sensible heat flux from soil (W/m**2) [+ to atm]
    real(rkx), pointer :: eflx_sh_h2osfc(:)
    ! specific humidity at snow surface [kg/kg]
    real(rkx), pointer :: qg_snow(:)
    ! specific humidity at soil surface [kg/kg]
    real(rkx), pointer :: qg_soil(:)
    ! specific humidity at h2osfc surface [kg/kg]
    real(rkx), pointer :: qg_h2osfc(:)
    ! evaporation flux from snow (W/m**2) [+ to atm]
    real(rkx), pointer :: qflx_ev_snow(:)
    ! evaporation flux from soil (W/m**2) [+ to atm]
    real(rkx), pointer :: qflx_ev_soil(:)
    ! evaporation flux from h2osfc (W/m**2) [+ to atm]
    real(rkx), pointer :: qflx_ev_h2osfc(:)
    integer(ik4) , pointer :: pcolumn(:)        ! pft's column index
    integer(ik4) , pointer :: pgridcell(:)      ! pft's gridcell index
    integer(ik4) , pointer :: plandunit(:)      ! pft's landunit index
    integer(ik4) , pointer :: ltype(:)          ! landunit type
    ! fraction of vegetation not covered by snow (0 OR 1) [-]
    integer(ik4) , pointer :: frac_veg_nosno(:)
    ! ground surface temperature [K]
    real(rkx), pointer :: t_grnd(:)
    ! intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(rkx), pointer :: thm(:)
    ! specific humidity at ground surface [kg/kg]
    real(rkx), pointer :: qg(:)
    ! virtual potential temperature (kelvin)
    real(rkx), pointer :: thv(:)
    ! temperature derivative of "qg"
    real(rkx), pointer :: dqgdT(:)
    ! latent heat of evaporation (/sublimation) [J/kg]
    real(rkx), pointer :: htvp(:)
    ! coefficient of conective velocity [-]
    real(rkx), pointer :: beta(:)
    ! convective boundary height [m]
    real(rkx), pointer :: zii(:)
    ! atmospheric wind speed in east direction (m/s)
    real(rkx), pointer :: forc_u(:)
    ! atmospheric wind speed in north direction (m/s)
    real(rkx), pointer :: forc_v(:)
    ! atmospheric temperature (Kelvin)
    real(rkx), pointer :: forc_t(:)
    ! atmospheric potential temperature (Kelvin)
    real(rkx), pointer :: forc_th(:)
    ! atmospheric specific humidity (kg/kg)
    real(rkx), pointer :: forc_q(:)
    ! density (kg/m**3)
    real(rkx), pointer :: forc_rho(:)
    ! atmospheric pressure (Pa)
    real(rkx), pointer :: forc_pbot(:)
    ! observational height of wind at pft level [m]
    real(rkx), pointer :: forc_hgt_u_pft(:)
    ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsun(:)
    ! Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsun_wc(:)
    ! RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsun_wj(:)
    ! product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsun_wp(:)
    ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsha(:)
    ! Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsha_wc(:)
    ! RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsha_wj(:)
    ! product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rkx), pointer :: psnsha_wp(:)
    ! roughness length, momentum [m]
    real(rkx), pointer :: z0mg_col(:)
    real(rkx), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(rkx), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(rkx), pointer :: dz(:,:)           ! layer depth (m)
    ! volumetric soil water at saturation (porosity)
    real(rkx), pointer :: watsat(:,:)
    ! fraction of ground covered by snow (0 to 1)
    real(rkx), pointer :: frac_sno(:)
    ! soil wetness relative to field capacity
    real(rkx), pointer :: soilbeta(:)

    real(rkx), pointer :: z0hg_col(:)   ! roughness length, sensible heat [m]
    real(rkx), pointer :: z0qg_col(:)   ! roughness length, latent heat [m]

    ! downward longwave radiation below the canopy [W/m2]
    real(rkx), pointer :: dlrad(:)
    ! upward longwave radiation above the canopy [W/m2]
    real(rkx), pointer :: ulrad(:)
    ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(rkx), pointer :: cgrnds(:)
    ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(rkx), pointer :: cgrndl(:)
    ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(rkx), pointer :: cgrnd(:)
    real(rkx), pointer :: taux(:)    ! wind (shear) stress: e-w (kg/m/s**2)
    real(rkx), pointer :: tauy(:)    ! wind (shear) stress: n-s (kg/m/s**2)
    ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(rkx), pointer :: eflx_sh_grnd(:)
    ! total sensible heat flux (W/m**2) [+ to atm]
    real(rkx), pointer :: eflx_sh_tot(:)
    ! soil evaporation (mm H2O/s) (+ = to atm)
    real(rkx), pointer :: qflx_evap_soi(:)
    ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(rkx), pointer :: qflx_evap_tot(:)
    ! 2 m height surface air temperature (Kelvin)
    real(rkx), pointer :: t_ref2m(:)
    ! 2 m height surface specific humidity (kg/kg)
    real(rkx), pointer :: q_ref2m(:)
    ! Rural 2 m height surface air temperature (Kelvin)
    real(rkx), pointer :: t_ref2m_r(:)
    ! Rural 2 m height surface relative humidity (%)
    real(rkx), pointer :: rh_ref2m_r(:)
    ! 2 m height surface relative humidity (%)
    real(rkx), pointer :: rh_ref2m(:)
    real(rkx), pointer :: t_veg(:)   ! vegetation temperature (Kelvin)
    real(rkx), pointer :: btran(:)   ! transpiration wetness factor (0 to 1)
    real(rkx), pointer :: rssun(:)   ! sunlit stomatal resistance (s/m)
    real(rkx), pointer :: rssha(:)   ! shaded stomatal resistance (s/m)
    real(rkx), pointer :: ram1(:)    ! aerodynamical resistance (s/m)
    real(rkx), pointer :: fpsn(:)    ! photosynthesis (umol CO2 /m**2 /s)
    ! Rubisco-limited photosynthesis (umol CO2 /m**2 /s)
    real(rkx), pointer :: fpsn_wc(:)
    ! RuBP-limited photosynthesis (umol CO2 /m**2 /s)
    real(rkx), pointer :: fpsn_wj(:)
    ! product-limited photosynthesis (umol CO2 /m**2 /s)
    real(rkx), pointer :: fpsn_wp(:)
    ! effective fraction of roots in each soil layer
    real(rkx), pointer :: rootr(:,:)
    ! root resistance by layer (0-1)  (nlevgrnd)
    real(rkx), pointer :: rresis(:,:)
#if (defined LCH4)
    ! tracer conductance for boundary layer [m/s]
    real(rkx), pointer :: grnd_ch4_cond(:)
#endif

    ! maximum number of iterations for surface temperature
    integer(ik4), parameter  :: niters = 3
    integer(ik4)  :: p,c,g,f,j,l        ! indices
    integer(ik4)  :: filterp(ubp-lbp+1) ! pft filter for vegetated pfts
    integer(ik4)  :: fn                 ! number of values in local pft filter
    integer(ik4)  :: fp                 ! lake filter pft index
    integer(ik4)  :: iter               ! iteration index
    ! reference height "minus" zero displacement height [m]
    real(rkx) :: zldis(lbp:ubp)
    real(rkx) :: displa(lbp:ubp)  ! displacement height [m]
    ! dimensionless height used in Monin-Obukhov theory
    real(rkx) :: zeta
    real(rkx) :: wc               ! convective velocity [m/s]
    ! diff of virtual temp. between ref. height and surface
    real(rkx) :: dth(lbp:ubp)
    ! diff of vir. poten. temp. between ref. height and surface
    real(rkx) :: dthv
    ! diff of humidity between ref. height and surface
    real(rkx) :: dqh(lbp:ubp)
    real(rkx) :: obu(lbp:ubp)   ! Monin-Obukhov length (m)
    real(rkx) :: ur(lbp:ubp)    ! wind speed at reference height [m/s]
    real(rkx) :: um(lbp:ubp)    ! wind speed including the stablity effect [m/s]
    real(rkx) :: temp1(lbp:ubp) ! relation for potential temperature profile
    ! relation for potential temperature profile applied at 2-m
    real(rkx) :: temp12m(lbp:ubp)
    real(rkx) :: temp2(lbp:ubp) ! relation for specific humidity profile
    ! relation for specific humidity profile applied at 2-m
    real(rkx) :: temp22m(lbp:ubp)
    real(rkx) :: ustar(lbp:ubp) ! friction velocity [m/s]
    real(rkx) :: tstar          ! temperature scaling parameter
    real(rkx) :: qstar          ! moisture scaling parameter
    ! virtual potential temperature scaling parameter
    real(rkx) :: thvstar
    real(rkx) :: cf           ! heat transfer coefficient from leaves [-]
    real(rkx) :: ram          ! aerodynamical resistance [s/m]
    real(rkx) :: rah          ! thermal resistance [s/m]
    real(rkx) :: raw          ! moisture resistance [s/m]
    real(rkx) :: raih         ! temporary variable [kg/m2/s]
    real(rkx) :: raiw         ! temporary variable [kg/m2/s]
    real(rkx) :: fm(lbp:ubp)  ! needed for BGC only to diagnose 10m wind speed
    real(rkx) :: z0mg_pft(lbp:ubp)
    real(rkx) :: z0hg_pft(lbp:ubp)
    real(rkx) :: z0qg_pft(lbp:ubp)
    real(rkx) :: e_ref2m      ! 2 m height surface saturated vapor pressure [Pa]
    ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(rkx) :: de2mdT
    ! 2 m height surface saturated specific humidity [kg/kg]
    real(rkx) :: qsat_ref2m
    ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    real(rkx) :: dqsat2mdT
    real(rkx) :: www          ! surface soil wetness [-]

    t_soisno       => clm3%g%l%c%ces%t_soisno
    snl            => clm3%g%l%c%cps%snl
    t_h2osfc       => clm3%g%l%c%ces%t_h2osfc
    eflx_sh_snow   => clm3%g%l%c%p%pef%eflx_sh_snow
    eflx_sh_soil   => clm3%g%l%c%p%pef%eflx_sh_soil
    eflx_sh_h2osfc => clm3%g%l%c%p%pef%eflx_sh_h2osfc
    qg_snow        => clm3%g%l%c%cws%qg_snow
    qg_soil        => clm3%g%l%c%cws%qg_soil
    qg_h2osfc      => clm3%g%l%c%cws%qg_h2osfc
    qflx_ev_snow   => clm3%g%l%c%p%pwf%qflx_ev_snow
    qflx_ev_soil   => clm3%g%l%c%p%pwf%qflx_ev_soil
    qflx_ev_h2osfc => clm3%g%l%c%p%pwf%qflx_ev_h2osfc

    ! Assign local pointers to derived type members (gridcell-level)

    forc_u     => clm_a2l%forc_u
    forc_v     => clm_a2l%forc_v

    ! Assign local pointers to derived type members (landunit-level)

    ltype      => clm3%g%l%itype

    ! Assign local pointers to derived type members (column-level)

    forc_th    => clm3%g%l%c%ces%forc_th
    forc_t     => clm3%g%l%c%ces%forc_t
    forc_pbot  => clm3%g%l%c%cps%forc_pbot
    forc_rho   => clm3%g%l%c%cps%forc_rho
    forc_q     => clm3%g%l%c%cws%forc_q
    pcolumn    => clm3%g%l%c%p%column
    pgridcell  => clm3%g%l%c%p%gridcell
    frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
    dlrad  => clm3%g%l%c%p%pef%dlrad
    ulrad  => clm3%g%l%c%p%pef%ulrad
    t_grnd => clm3%g%l%c%ces%t_grnd
    qg     => clm3%g%l%c%cws%qg
    z0mg_col => clm3%g%l%c%cps%z0mg
    z0hg_col => clm3%g%l%c%cps%z0hg
    z0qg_col => clm3%g%l%c%cps%z0qg
    thv    => clm3%g%l%c%ces%thv
    beta   => clm3%g%l%c%cps%beta
    zii    => clm3%g%l%c%cps%zii
    ram1   => clm3%g%l%c%p%pps%ram1
    cgrnds => clm3%g%l%c%p%pef%cgrnds
    cgrndl => clm3%g%l%c%p%pef%cgrndl
    cgrnd  => clm3%g%l%c%p%pef%cgrnd
    dqgdT  => clm3%g%l%c%cws%dqgdT
    htvp   => clm3%g%l%c%cps%htvp
    watsat         => clm3%g%l%c%cps%watsat
    h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
    dz             => clm3%g%l%c%cps%dz
    h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
    frac_sno       => clm3%g%l%c%cps%frac_sno
    soilbeta       => clm3%g%l%c%cws%soilbeta

    ! Assign local pointers to derived type members (pft-level)

    taux => clm3%g%l%c%p%pmf%taux
    tauy => clm3%g%l%c%p%pmf%tauy
    eflx_sh_grnd => clm3%g%l%c%p%pef%eflx_sh_grnd
    eflx_sh_tot => clm3%g%l%c%p%pef%eflx_sh_tot
    qflx_evap_soi => clm3%g%l%c%p%pwf%qflx_evap_soi
    qflx_evap_tot => clm3%g%l%c%p%pwf%qflx_evap_tot
    t_ref2m => clm3%g%l%c%p%pes%t_ref2m
    q_ref2m => clm3%g%l%c%p%pes%q_ref2m
    t_ref2m_r => clm3%g%l%c%p%pes%t_ref2m_r
    rh_ref2m_r => clm3%g%l%c%p%pes%rh_ref2m_r
    plandunit => clm3%g%l%c%p%landunit
    rh_ref2m => clm3%g%l%c%p%pes%rh_ref2m
    t_veg => clm3%g%l%c%p%pes%t_veg
    thm => clm3%g%l%c%p%pes%thm
    btran => clm3%g%l%c%p%pps%btran
    rssun => clm3%g%l%c%p%pps%rssun
    rssha => clm3%g%l%c%p%pps%rssha
    rootr => clm3%g%l%c%p%pps%rootr
    rresis => clm3%g%l%c%p%pps%rresis
    psnsun => clm3%g%l%c%p%pcf%psnsun
    psnsun_wc => clm3%g%l%c%p%pcf%psnsun_wc
    psnsun_wj => clm3%g%l%c%p%pcf%psnsun_wj
    psnsun_wp => clm3%g%l%c%p%pcf%psnsun_wp
    psnsha => clm3%g%l%c%p%pcf%psnsha
    psnsha_wc => clm3%g%l%c%p%pcf%psnsha_wc
    psnsha_wj => clm3%g%l%c%p%pcf%psnsha_wj
    psnsha_wp => clm3%g%l%c%p%pcf%psnsha_wp
    fpsn => clm3%g%l%c%p%pcf%fpsn
    fpsn_wc => clm3%g%l%c%p%pcf%fpsn_wc
    fpsn_wj => clm3%g%l%c%p%pcf%fpsn_wj
    fpsn_wp => clm3%g%l%c%p%pcf%fpsn_wp
    forc_hgt_u_pft => clm3%g%l%c%p%pps%forc_hgt_u_pft
#if (defined LCH4)
    grnd_ch4_cond  => clm3%g%l%c%p%pps%grnd_ch4_cond
#endif

    ! Filter pfts where frac_veg_nosno is zero

    fn = 0
    do fp = 1 , num_nolakep
      p = filter_nolakep(fp)
      if ( frac_veg_nosno(p) == 0 ) then
        fn = fn + 1
        filterp(fn) = p
      end if
    end do

    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step

    do f = 1 , fn
      p = filterp(f)
      c = pcolumn(p)
      g = pgridcell(p)

      ! Initialization variables

      displa(p) = 0._rkx
      dlrad(p)  = 0._rkx
      ulrad(p)  = 0._rkx

      ur(p) = max(1.0_rkx,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
      dth(p) = thm(p)-t_grnd(c)
      dqh(p) = forc_q(c) - qg(c)
      dthv = dth(p)*(1._rkx+0.61_rkx*forc_q(c))+0.61_rkx*forc_th(c)*dqh(p)
      zldis(p) = forc_hgt_u_pft(p)

      ! Copy column roughness to local pft-level arrays

      z0mg_pft(p) = z0mg_col(c)
      z0hg_pft(p) = z0hg_col(c)
      z0qg_pft(p) = z0qg_col(c)

      ! Initialize Monin-Obukhov length and wind speed

      call MoninObukIni(ur(p),thv(c),dthv,zldis(p),z0mg_pft(p),um(p),obu(p))

    end do

    ! Perform stability iteration
    ! Determine friction velocity, and potential temperature and humidity
    ! profiles of the surface boundary layer

    do iter = 1 , niters

      call FrictionVelocity(lbp, ubp, fn, filterp, &
                            displa, z0mg_pft, z0hg_pft, z0qg_pft, &
                            obu, iter, ur, um, ustar, &
                            temp1, temp2, temp12m, temp22m, fm)

      do f = 1 , fn
        p = filterp(f)
        c = pcolumn(p)
        g = pgridcell(p)

        tstar = temp1(p)*dth(p)
        qstar = temp2(p)*dqh(p)
        z0hg_pft(p) = z0mg_pft(p) / &
          exp(0.13_rkx * (ustar(p)*z0mg_pft(p)/1.5e-5_rkx)**0.45_rkx)
        z0qg_pft(p) = z0hg_pft(p)
        thvstar = tstar*(1._rkx+0.61_rkx*forc_q(c)) + 0.61_rkx*forc_th(c)*qstar
        zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

        if ( zeta >= 0._rkx ) then                   !stable
          zeta = min(2._rkx,max(zeta,0.01_rkx))
          um(p) = max(ur(p),0.1_rkx)
        else                                      !unstable
          zeta = max(-100._rkx,min(zeta,-0.01_rkx))
          wc = beta(c)*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333_rkx
          um(p) = sqrt(ur(p)*ur(p) + wc*wc)
        end if
        obu(p) = zldis(p)/zeta
      end do
    end do ! end stability iteration

    do j = 1 , nlevgrnd
      do f = 1 , fn
        p = filterp(f)
        rootr(p,j) = 0._rkx
        rresis(p,j) = 0._rkx
      end do
    end do

    do f = 1 , fn
      p = filterp(f)
      c = pcolumn(p)
      g = pgridcell(p)
      l = plandunit(p)

      ! Determine aerodynamic resistances

      ram  = 1._rkx/(ustar(p)*ustar(p)/um(p))
      rah  = 1._rkx/(temp1(p)*ustar(p))
      raw  = 1._rkx/(temp2(p)*ustar(p))
      raih = forc_rho(c)*cpair/rah
#if (defined LCH4)
      grnd_ch4_cond(p) = 1._rkx/raw
#endif

      ! Soil evaporation resistance
      www = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
      www = min(max(www,0.0_rkx),1._rkx)

      !changed by K.Sakaguchi. Soilbeta is used for evaporation
      if ( dqh(p) > 0._rkx ) then
        !dew  (beta is not applied, just like rsoil used to be)
        raiw = forc_rho(c)/(raw)
      else
        ! Lee and Pielke 1992 beta is applied
        raiw = soilbeta(c)*forc_rho(c)/(raw)
      end if

      ram1(p) = ram  !pass value to global variable

      ! Output to pft-level data structures
      ! Derivative of fluxes with respect to ground temperature

      cgrnds(p) = raih
      cgrndl(p) = raiw*dqgdT(c)
      cgrnd(p) = cgrnds(p) + htvp(c)*cgrndl(p)

      ! Surface fluxes of momentum, sensible and latent heat
      ! using ground temperatures from previous time step

      taux(p) = -forc_rho(c)*forc_u(g)/ram
      tauy(p) = -forc_rho(c)*forc_v(g)/ram
      eflx_sh_grnd(p) = -raih*dth(p)
      eflx_sh_tot(p) = eflx_sh_grnd(p)
      ! compute sensible heat fluxes individually
      eflx_sh_snow(p) = -raih*(thm(p)-t_soisno(c,snl(c)+1))
      eflx_sh_soil(p) = -raih*(thm(p)-t_soisno(c,1))
      eflx_sh_h2osfc(p) = -raih*(thm(p)-t_h2osfc(c))
      qflx_evap_soi(p) = -raiw*dqh(p)
      qflx_evap_tot(p) = qflx_evap_soi(p)
      ! compute latent heat fluxes individually
      qflx_ev_snow(p) = -raiw*(forc_q(c) - qg_snow(c))
      qflx_ev_soil(p) = -raiw*(forc_q(c) - qg_soil(c))
      qflx_ev_h2osfc(p) = -raiw*(forc_q(c) - qg_h2osfc(c))

      ! 2 m height air temperature

      t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._rkx/temp12m(p)-1._rkx/temp1(p))

      ! 2 m height specific humidity

      q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1._rkx/temp22m(p)-1._rkx/temp2(p))

      ! 2 m height relative humidity

      call QSat(t_ref2m(p),forc_pbot(c),e_ref2m,de2mdT,qsat_ref2m,dqsat2mdT)

      rh_ref2m(p) = min(100._rkx, q_ref2m(p) / qsat_ref2m * 100._rkx)

      if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
        rh_ref2m_r(p) = rh_ref2m(p)
        t_ref2m_r(p) = t_ref2m(p)
      end if

      ! Variables needed by history tape

      t_veg(p) = forc_t(c)
      btran(p) = 0._rkx
      cf = forc_pbot(c)/(rgas*0.001_rkx*thm(p))*1.e6_rkx
      rssun(p) = 1._rkx/1.e15_rkx * cf
      rssha(p) = 1._rkx/1.e15_rkx * cf

      ! Add the following to avoid NaN

      psnsun(p) = 0._rkx
      psnsun_wc(p) = 0._rkx
      psnsun_wj(p) = 0._rkx
      psnsun_wp(p) = 0._rkx
      psnsha(p) = 0._rkx
      psnsha_wc(p) = 0._rkx
      psnsha_wj(p) = 0._rkx
      psnsha_wp(p) = 0._rkx
      fpsn(p) = 0._rkx
      fpsn_wc(p) = 0._rkx
      fpsn_wj(p) = 0._rkx
      fpsn_wp(p) = 0._rkx

      ! adding code for isotopes, 8/17/05, PET
      if ( use_c13 ) then
        clm3%g%l%c%p%pps%alphapsnsun(p) = 0._rkx
        clm3%g%l%c%p%pps%alphapsnsha(p) = 0._rkx
        clm3%g%l%c%p%pepv%rc13_canair(p) = 0._rkx
        clm3%g%l%c%p%pepv%rc13_psnsun(p) = 0._rkx
        clm3%g%l%c%p%pepv%rc13_psnsha(p) = 0._rkx
        clm3%g%l%c%p%pc13f%psnsun(p) = 0._rkx
        clm3%g%l%c%p%pc13f%psnsha(p) = 0._rkx
        clm3%g%l%c%p%pc13f%fpsn(p) = 0._rkx
      end if

      if ( use_c14 ) then
        clm3%g%l%c%p%pepv%rc14_atm(p) = 0._rkx
        ! clm3%g%l%c%p%pepv%rc14_canair(p) = 0._rkx
        ! clm3%g%l%c%p%pepv%rc14_psnsun(p) = 0._rkx
        ! clm3%g%l%c%p%pepv%rc14_psnsha(p) = 0._rkx
        clm3%g%l%c%p%pc14f%psnsun(p) = 0._rkx
        clm3%g%l%c%p%pc14f%psnsha(p) = 0._rkx
        clm3%g%l%c%p%pc14f%fpsn(p) = 0._rkx
      end if
    end do
  end subroutine BareGroundFluxes

end module mod_clm_baregroundfluxes
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
