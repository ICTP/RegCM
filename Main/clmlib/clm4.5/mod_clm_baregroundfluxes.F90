module mod_clm_baregroundfluxes

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: BareGroundFluxesMod
!
! !DESCRIPTION:
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
  use mod_realkinds
!
! !PUBLIC TYPES:
   implicit none
   save
!
! !PUBLIC MEMBER FUNCTIONS:
   public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BareGroundFluxes
!
! !INTERFACE:
  subroutine BareGroundFluxes(lbp, ubp, num_nolakep, filter_nolakep)
!
! !DESCRIPTION:
! Compute sensible and latent fluxes and their derivatives with respect
! to ground temperature using ground temperatures from previous time step.
!
! !USES:
    use mod_clm_type
    use mod_clm_atmlnd         , only : clm_a2l
    use mod_clm_varpar         , only : nlevgrnd
    use mod_clm_varcon         , only : cpair, vkc, grav, denice, denh2o, istsoil
    use mod_clm_varcon         , only : istcrop , rgas
    use mod_clm_frictionvelocity, only : FrictionVelocity, MoninObukIni
    use mod_clm_qsat       , only : QSat
    use mod_clm_varctl         , only : use_c13, use_c14

!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp                     ! pft bounds
    integer, intent(in) :: num_nolakep                  ! number of pft non-lake points in pft filter
    integer, intent(in) :: filter_nolakep(ubp-lbp+1)    ! pft filter for non-lake points
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! This routine originally had a long list of parameters, and also a reference to
! the entire clm derived type.  For consistency, only the derived type reference
! is passed (now pointing to the current column and pft), and the other original
! parameters are initialized locally. Using t_grnd instead of tg (tg eliminated
! as redundant).
! 1/23/02, PET: Added pft reference as parameter. All outputs will be written
! to the pft data structures, and averaged to the column level outside of
! this routine.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    real(rk8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)
    integer , pointer :: snl(:)            ! number of snow layers
    real(rk8), pointer :: t_h2osfc(:) 	   ! surface water temperature
    real(rk8), pointer :: eflx_sh_snow(:)   ! sensible heat flux from snow (W/m**2) [+ to atm]
    real(rk8), pointer :: eflx_sh_soil(:)   ! sensible heat flux from soil (W/m**2) [+ to atm]
    real(rk8), pointer :: eflx_sh_h2osfc(:) ! sensible heat flux from soil (W/m**2) [+ to atm]
    real(rk8), pointer :: qg_snow(:)        ! specific humidity at snow surface [kg/kg]
    real(rk8), pointer :: qg_soil(:)        ! specific humidity at soil surface [kg/kg]
    real(rk8), pointer :: qg_h2osfc(:)      ! specific humidity at h2osfc surface [kg/kg]
    real(rk8), pointer :: qflx_ev_snow(:)   ! evaporation flux from snow (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_ev_soil(:)   ! evaporation flux from soil (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_ev_h2osfc(:) ! evaporation flux from h2osfc (W/m**2) [+ to atm]
    integer , pointer :: pcolumn(:)        ! pft's column index
    integer , pointer :: pgridcell(:)      ! pft's gridcell index
    integer , pointer :: plandunit(:)      ! pft's landunit index
    integer , pointer :: ltype(:)          ! landunit type
    integer , pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
    real(rk8), pointer :: t_grnd(:)         ! ground surface temperature [K]
    real(rk8), pointer :: thm(:)            ! intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
    real(rk8), pointer :: qg(:)             ! specific humidity at ground surface [kg/kg]
    real(rk8), pointer :: thv(:)            ! virtual potential temperature (kelvin)
    real(rk8), pointer :: dqgdT(:)          ! temperature derivative of "qg"
    real(rk8), pointer :: htvp(:)           ! latent heat of evaporation (/sublimation) [J/kg]
    real(rk8), pointer :: beta(:)           ! coefficient of conective velocity [-]
    real(rk8), pointer :: zii(:)            ! convective boundary height [m]
    real(rk8), pointer :: forc_u(:)         ! atmospheric wind speed in east direction (m/s)
    real(rk8), pointer :: forc_v(:)         ! atmospheric wind speed in north direction (m/s)
    real(rk8), pointer :: forc_t(:)         ! atmospheric temperature (Kelvin)
    real(rk8), pointer :: forc_th(:)        ! atmospheric potential temperature (Kelvin)
    real(rk8), pointer :: forc_q(:)         ! atmospheric specific humidity (kg/kg)
    real(rk8), pointer :: forc_rho(:)       ! density (kg/m**3)
    real(rk8), pointer :: forc_pbot(:)      ! atmospheric pressure (Pa)
    real(rk8), pointer :: forc_hgt_u_pft(:) ! observational height of wind at pft level [m]
    real(rk8), pointer :: psnsun(:)         ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer :: psnsun_wc(:)      ! Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer :: psnsun_wj(:)      ! RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer :: psnsun_wp(:)      ! product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer :: psnsha(:)         ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer :: psnsha_wc(:)      ! Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer :: psnsha_wj(:)      ! RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer :: psnsha_wp(:)      ! product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
    real(rk8), pointer :: z0mg_col(:)       ! roughness length, momentum [m]
    real(rk8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
    real(rk8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
    real(rk8), pointer :: dz(:,:)           ! layer depth (m)
    real(rk8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity)
    real(rk8), pointer :: frac_sno(:)       ! fraction of ground covered by snow (0 to 1)
    real(rk8), pointer :: soilbeta(:)       ! soil wetness relative to field capacity
!
! local pointers to implicit inout arguments
!
    real(rk8), pointer :: z0hg_col(:)       ! roughness length, sensible heat [m]
    real(rk8), pointer :: z0qg_col(:)       ! roughness length, latent heat [m]
!
! local pointers to implicit out arguments
!
    real(rk8), pointer :: dlrad(:)         ! downward longwave radiation below the canopy [W/m2]
    real(rk8), pointer :: ulrad(:)         ! upward longwave radiation above the canopy [W/m2]
    real(rk8), pointer :: cgrnds(:)        ! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
    real(rk8), pointer :: cgrndl(:)        ! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
    real(rk8), pointer :: cgrnd(:)         ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
    real(rk8), pointer :: taux(:)          ! wind (shear) stress: e-w (kg/m/s**2)
    real(rk8), pointer :: tauy(:)          ! wind (shear) stress: n-s (kg/m/s**2)
    real(rk8), pointer :: eflx_sh_grnd(:)  ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(rk8), pointer :: eflx_sh_tot(:)   ! total sensible heat flux (W/m**2) [+ to atm]
    real(rk8), pointer :: qflx_evap_soi(:) ! soil evaporation (mm H2O/s) (+ = to atm)
    real(rk8), pointer :: qflx_evap_tot(:) ! qflx_evap_soi + qflx_evap_can + qflx_tran_veg
    real(rk8), pointer :: t_ref2m(:)       ! 2 m height surface air temperature (Kelvin)
    real(rk8), pointer :: q_ref2m(:)       ! 2 m height surface specific humidity (kg/kg)
    real(rk8), pointer :: t_ref2m_r(:)     ! Rural 2 m height surface air temperature (Kelvin)
    real(rk8), pointer :: rh_ref2m_r(:)    ! Rural 2 m height surface relative humidity (%)
    real(rk8), pointer :: rh_ref2m(:)      ! 2 m height surface relative humidity (%)
    real(rk8), pointer :: t_veg(:)         ! vegetation temperature (Kelvin)
    real(rk8), pointer :: btran(:)         ! transpiration wetness factor (0 to 1)
    real(rk8), pointer :: rssun(:)         ! sunlit stomatal resistance (s/m)
    real(rk8), pointer :: rssha(:)         ! shaded stomatal resistance (s/m)
    real(rk8), pointer :: ram1(:)          ! aerodynamical resistance (s/m)
    real(rk8), pointer :: fpsn(:)          ! photosynthesis (umol CO2 /m**2 /s)
    real(rk8), pointer :: fpsn_wc(:)       ! Rubisco-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8), pointer :: fpsn_wj(:)       ! RuBP-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8), pointer :: fpsn_wp(:)       ! product-limited photosynthesis (umol CO2 /m**2 /s)
    real(rk8), pointer :: rootr(:,:)       ! effective fraction of roots in each soil layer
    real(rk8), pointer :: rresis(:,:)      ! root resistance by layer (0-1)  (nlevgrnd)	
#if (defined LCH4)
    real(rk8), pointer :: grnd_ch4_cond(:)  ! tracer conductance for boundary layer [m/s]
#endif
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer, parameter  :: niters = 3  ! maximum number of iterations for surface temperature
    integer  :: p,c,g,f,j,l            ! indices
    integer  :: filterp(ubp-lbp+1)     ! pft filter for vegetated pfts
    integer  :: fn                     ! number of values in local pft filter
    integer  :: fp                     ! lake filter pft index
    integer  :: iter                   ! iteration index
    real(rk8) :: zldis(lbp:ubp)         ! reference height "minus" zero displacement height [m]
    real(rk8) :: displa(lbp:ubp)        ! displacement height [m]
    real(rk8) :: zeta                   ! dimensionless height used in Monin-Obukhov theory
    real(rk8) :: wc                     ! convective velocity [m/s]
    real(rk8) :: dth(lbp:ubp)           ! diff of virtual temp. between ref. height and surface
    real(rk8) :: dthv                   ! diff of vir. poten. temp. between ref. height and surface
    real(rk8) :: dqh(lbp:ubp)           ! diff of humidity between ref. height and surface
    real(rk8) :: obu(lbp:ubp)           ! Monin-Obukhov length (m)
    real(rk8) :: ur(lbp:ubp)            ! wind speed at reference height [m/s]
    real(rk8) :: um(lbp:ubp)            ! wind speed including the stablity effect [m/s]
    real(rk8) :: temp1(lbp:ubp)         ! relation for potential temperature profile
    real(rk8) :: temp12m(lbp:ubp)       ! relation for potential temperature profile applied at 2-m
    real(rk8) :: temp2(lbp:ubp)         ! relation for specific humidity profile
    real(rk8) :: temp22m(lbp:ubp)       ! relation for specific humidity profile applied at 2-m
    real(rk8) :: ustar(lbp:ubp)         ! friction velocity [m/s]
    real(rk8) :: tstar                  ! temperature scaling parameter
    real(rk8) :: qstar                  ! moisture scaling parameter
    real(rk8) :: thvstar                ! virtual potential temperature scaling parameter
    real(rk8) :: cf                     ! heat transfer coefficient from leaves [-]
    real(rk8) :: ram                    ! aerodynamical resistance [s/m]
    real(rk8) :: rah                    ! thermal resistance [s/m]
    real(rk8) :: raw                    ! moisture resistance [s/m]
    real(rk8) :: raih                   ! temporary variable [kg/m2/s]
    real(rk8) :: raiw                   ! temporary variable [kg/m2/s]
    real(rk8) :: fm(lbp:ubp)            ! needed for BGC only to diagnose 10m wind speed
    real(rk8) :: z0mg_pft(lbp:ubp)
    real(rk8) :: z0hg_pft(lbp:ubp)
    real(rk8) :: z0qg_pft(lbp:ubp)
    real(rk8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(rk8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(rk8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(rk8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m 
    real(rk8) :: www                    ! surface soil wetness [-]
!------------------------------------------------------------------------------

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
    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       if (frac_veg_nosno(p) == 0) then
          fn = fn + 1
          filterp(fn) = p
       end if
    end do

    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step

    do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)

       ! Initialization variables

       displa(p) = 0.D0
       dlrad(p)  = 0.D0
       ulrad(p)  = 0.D0

       ur(p) = max(1.0D0,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       dth(p) = thm(p)-t_grnd(c)
       dqh(p) = forc_q(c) - qg(c)
       dthv = dth(p)*(1.D0+0.61D0*forc_q(c))+0.61D0*forc_th(c)*dqh(p)
       zldis(p) = forc_hgt_u_pft(p)

       ! Copy column roughness to local pft-level arrays

       z0mg_pft(p) = z0mg_col(c)
       z0hg_pft(p) = z0hg_col(c)
       z0qg_pft(p) = z0qg_col(c)

       ! Initialize Monin-Obukhov length and wind speed

       call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg_pft(p), um(p), obu(p))

    end do

    ! Perform stability iteration
    ! Determine friction velocity, and potential temperature and humidity
    ! profiles of the surface boundary layer

    do iter = 1, niters

       call FrictionVelocity(lbp, ubp, fn, filterp, &
                             displa, z0mg_pft, z0hg_pft, z0qg_pft, &
                             obu, iter, ur, um, ustar, &
                             temp1, temp2, temp12m, temp22m, fm)

       do f = 1, fn
          p = filterp(f)
          c = pcolumn(p)
          g = pgridcell(p)

          tstar = temp1(p)*dth(p)
          qstar = temp2(p)*dqh(p)
          z0hg_pft(p) = z0mg_pft(p)/exp(0.13D0 * (ustar(p)*z0mg_pft(p)/1.5D-5)**0.45D0)
          z0qg_pft(p) = z0hg_pft(p)
          thvstar = tstar*(1.D0+0.61D0*forc_q(c)) + 0.61D0*forc_th(c)*qstar
          zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

          if (zeta >= 0.D0) then                   !stable
             zeta = min(2.D0,max(zeta,0.01D0))
             um(p) = max(ur(p),0.1D0)
          else                                      !unstable
             zeta = max(-100.D0,min(zeta,-0.01D0))
             wc = beta(c)*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333D0
             um(p) = sqrt(ur(p)*ur(p) + wc*wc)
          end if
          obu(p) = zldis(p)/zeta
       end do

    end do ! end stability iteration

     do j = 1, nlevgrnd
       do f = 1, fn
          p = filterp(f)
          rootr(p,j) = 0.D0
          rresis(p,j) = 0.D0
        end do
     end do

    do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)
       l = plandunit(p)

       ! Determine aerodynamic resistances

       ram     = 1.D0/(ustar(p)*ustar(p)/um(p))
       rah     = 1.D0/(temp1(p)*ustar(p))
       raw     = 1.D0/(temp2(p)*ustar(p))
       raih    = forc_rho(c)*cpair/rah
#if (defined LCH4)
       grnd_ch4_cond(p) = 1.D0/raw
#endif

       ! Soil evaporation resistance
       www     = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
       www     = min(max(www,0.0D0),1.D0)

       !changed by K.Sakaguchi. Soilbeta is used for evaporation
       if (dqh(p) .gt. 0.D0) then   !dew  (beta is not applied, just like rsoil used to be)
          raiw    = forc_rho(c)/(raw)
       else
       ! Lee and Pielke 1992 beta is applied
          raiw    = soilbeta(c)*forc_rho(c)/(raw)
       end if

       ram1(p) = ram  !pass value to global variable

       ! Output to pft-level data structures
       ! Derivative of fluxes with respect to ground temperature

       cgrnds(p) = raih
       cgrndl(p) = raiw*dqgdT(c)
       cgrnd(p)  = cgrnds(p) + htvp(c)*cgrndl(p)

       ! Surface fluxes of momentum, sensible and latent heat
       ! using ground temperatures from previous time step

       taux(p)          = -forc_rho(c)*forc_u(g)/ram
       tauy(p)          = -forc_rho(c)*forc_v(g)/ram
       eflx_sh_grnd(p)  = -raih*dth(p)
       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       ! compute sensible heat fluxes individually
       eflx_sh_snow(p)  = -raih*(thm(p)-t_soisno(c,snl(c)+1))
       eflx_sh_soil(p)  = -raih*(thm(p)-t_soisno(c,1))
       eflx_sh_h2osfc(p)  = -raih*(thm(p)-t_h2osfc(c))
       qflx_evap_soi(p) = -raiw*dqh(p)
       qflx_evap_tot(p) = qflx_evap_soi(p)
       ! compute latent heat fluxes individually
       qflx_ev_snow(p) = -raiw*(forc_q(c) - qg_snow(c))
       qflx_ev_soil(p) = -raiw*(forc_q(c) - qg_soil(c))
       qflx_ev_h2osfc(p) = -raiw*(forc_q(c) - qg_h2osfc(c))

       ! 2 m height air temperature

       t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1.D0/temp12m(p) - 1.D0/temp1(p))

       ! 2 m height specific humidity

       q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1.D0/temp22m(p) - 1.D0/temp2(p))

       ! 2 m height relative humidity
                                                                                
       call QSat(t_ref2m(p), forc_pbot(c), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

       rh_ref2m(p) = min(100.D0, q_ref2m(p) / qsat_ref2m * 100.D0)

       if (ltype(l) == istsoil .or. ltype(l) == istcrop) then
         rh_ref2m_r(p) = rh_ref2m(p)
         t_ref2m_r(p) = t_ref2m(p)
       end if

       ! Variables needed by history tape

       t_veg(p) = forc_t(c)
       btran(p) = 0.D0
       cf = forc_pbot(c)/(rgas*0.001D0*thm(p))*1.D06
       rssun(p) = 1.D0/1.D15 * cf
       rssha(p) = 1.D0/1.D15 * cf

       ! Add the following to avoid NaN

       psnsun(p) = 0.D0
       psnsun_wc(p) = 0.D0
       psnsun_wj(p) = 0.D0
       psnsun_wp(p) = 0.D0
       psnsha(p) = 0.D0
       psnsha_wc(p) = 0.D0
       psnsha_wj(p) = 0.D0
       psnsha_wp(p) = 0.D0
       fpsn(p) = 0.D0
       fpsn_wc(p) = 0.D0
       fpsn_wj(p) = 0.D0
       fpsn_wp(p) = 0.D0
       ! adding code for isotopes, 8/17/05, PET

     if ( use_c13 ) then
        clm3%g%l%c%p%pps%alphapsnsun(p) = 0.D0
        clm3%g%l%c%p%pps%alphapsnsha(p) = 0.D0
        clm3%g%l%c%p%pepv%rc13_canair(p) = 0.D0
        clm3%g%l%c%p%pepv%rc13_psnsun(p) = 0.D0
        clm3%g%l%c%p%pepv%rc13_psnsha(p) = 0.D0
        clm3%g%l%c%p%pc13f%psnsun(p) = 0.D0
        clm3%g%l%c%p%pc13f%psnsha(p) = 0.D0
        clm3%g%l%c%p%pc13f%fpsn(p) = 0.D0
     endif

     if ( use_c14 ) then
        clm3%g%l%c%p%pepv%rc14_atm(p) = 0.D0
        ! clm3%g%l%c%p%pepv%rc14_canair(p) = 0.D0
        ! clm3%g%l%c%p%pepv%rc14_psnsun(p) = 0.D0
        ! clm3%g%l%c%p%pepv%rc14_psnsha(p) = 0.D0
        clm3%g%l%c%p%pc14f%psnsun(p) = 0.D0
        clm3%g%l%c%p%pc14f%psnsha(p) = 0.D0
        clm3%g%l%c%p%pc14f%fpsn(p) = 0.D0
     endif

    end do

  end subroutine BareGroundFluxes

end module mod_clm_baregroundfluxes
