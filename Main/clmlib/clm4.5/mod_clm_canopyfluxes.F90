module mod_clm_canopyfluxes
  !
  ! Calculates the leaf temperature and the leaf fluxes,
  ! transpiration, photosynthesis and  updates the dew
  ! accumulation due to evaporation.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use mod_runparams
  use mod_clm_varctl , only: use_c13 , use_c14
  use mod_clm_type
  use mod_clm_atmlnd , only : clm_a2l
  use mod_clm_varpar , only : nlevgrnd, nlevsno, nlevcan
  use mod_clm_varcon , only : sb, cpair, hvap, vkc, grav, denice, &
                    denh2o, tfrz, tlsai_crit, alpha_aero, &
                    isecspday, degpsec , tfrz, c14ratio
  use mod_clm_qsat , only : QSat
  use mod_clm_frictionvelocity , only : FrictionVelocity , MoninObukIni
  use mod_clm_pftvarcon , only : nbrdlf_dcd_tmp_shrub , irrigated
  use mod_clm_pftvarcon , only : nsoybean, nsoybeanirrig, npcropmin , &
    nbrdlf_evr_trp_tree , nbrdlf_dcd_trp_tree
#if (defined CN)
  use mod_clm_cnallocation , only : CNAllocation_Carbon_only
#endif

  implicit none

  private

  save

  public :: CanopyFluxes !Calculates the leaf temperature and leaf fluxes

  ! true => btran is based only on unfrozen soil levels
  logical,  public :: perchroot     = .false.
  ! true  => btran is based on active layer (defined over two years);
  ! false => btran is based on currently unfrozen levels
  logical,  public :: perchroot_alt = .false.

  private :: Photosynthesis ! Leaf stomatal resistance and leaf photosynthesis
  private :: hybrid         ! hybrid solver for ci
  private :: ci_func        ! ci function

  contains
  !
  ! 1. Calculates the leaf temperature:
  ! 2. Calculates the leaf fluxes, transpiration, photosynthesis and
  !    updates the dew accumulation due to evaporation.
  !
  ! Method:
  ! Use the Newton-Raphson iteration to solve for the foliage
  ! temperature that balances the surface energy budget:
  !
  ! f(t_veg) = Net radiation - Sensible - Latent = 0
  ! f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
  !
  ! Note:
  ! (1) In solving for t_veg, t_grnd is given from the previous timestep.
  ! (2) The partial derivatives of aerodynamical resistances, which cannot
  !     be determined analytically, are ignored for d(H)/dT and d(LE)/dT
  ! (3) The weighted stomatal resistance of sunlit and shaded foliage is used
  ! (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
  !                                                          => Ec + Eg = Ea
  ! (5) Energy loss is due to: numerical truncation of energy budget equation
  !     (*); and "ecidif" (see the code) which is dropped into the sensible
  !     heat
  ! (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n)
  !     and del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference
  !     of water flux from the leaf between the iteration step (n+1) and (n)
  !     less than 0.1 W/m2; or the iterative steps over 40.
  !
  subroutine CanopyFluxes(lbc, ubc, lbp, ubp, &
                          num_nolakep, filter_nolakep)
    use mod_clm_varcon , only : csoilc
    implicit none
    integer(ik4), intent(in) :: lbc, ubc ! column bounds
    integer(ik4), intent(in) :: lbp, ubp ! pft bounds
    ! number of column non-lake points in pft filter
    integer(ik4), intent(in) :: num_nolakep
    ! pft filter for non-lake points
    integer(ik4), intent(in) :: filter_nolakep(ubp-lbp+1)

   ! sensible heat flux from snow (W/m**2) [+ to atm]
   real(rkx), pointer :: eflx_sh_snow(:)
   ! sensible heat flux from soil (W/m**2) [+ to atm]
   real(rkx), pointer :: eflx_sh_soil(:)
   ! sensible heat flux from soil (W/m**2) [+ to atm]
   real(rkx), pointer :: eflx_sh_h2osfc(:)
   integer(ik4) , pointer :: snl(:)     ! number of snow layers
   real(rkx), pointer :: t_h2osfc(:)    ! surface water temperature
   real(rkx), pointer :: frac_h2osfc(:) ! fraction of surface water
   real(rkx), pointer :: qg_snow(:) ! specific humidity at snow surface [kg/kg]
   real(rkx), pointer :: qg_soil(:) ! specific humidity at soil surface [kg/kg]
   ! specific humidity at h2osfc surface [kg/kg]
   real(rkx), pointer :: qg_h2osfc(:)
   ! evaporation flux from snow (W/m**2) [+ to atm]
   real(rkx), pointer :: qflx_ev_snow(:)
   ! evaporation flux from soil (W/m**2) [+ to atm]
   real(rkx), pointer :: qflx_ev_soil(:)
   ! evaporation flux from h2osfc (W/m**2) [+ to atm]
   real(rkx), pointer :: qflx_ev_h2osfc(:)
   ! frac of veg not covered by snow (0 OR 1 now) [-]
   integer(ik4) , pointer :: frac_veg_nosno(:)
   integer(ik4) , pointer :: ivt(:)         ! pft vegetation type
   integer(ik4) , pointer :: pcolumn(:)     ! pft's column index
   integer(ik4) , pointer :: plandunit(:)   ! pft's landunit index
   integer(ik4) , pointer :: pgridcell(:)   ! pft's gridcell index
   ! atmospheric potential temperature (Kelvin)
   real(rkx), pointer :: forc_th(:)
   real(rkx), pointer :: t_grnd(:)  ! ground surface temperature [K]
   ! intermediate variable (forc_t+0.0098*forc_hgt_t_pft)
   real(rkx), pointer :: thm(:)
   real(rkx), pointer :: qg(:)   ! specific humidity at ground surface [kg/kg]
   real(rkx), pointer :: thv(:)  ! virtual potential temperature (kelvin)
   ! roughness length over vegetation, momentum [m]
   real(rkx), pointer :: z0mv(:)
   ! roughness length over vegetation, sensible heat [m]
   real(rkx), pointer :: z0hv(:)
   ! roughness length over vegetation, latent heat [m]
   real(rkx), pointer :: z0qv(:)
   ! roughness length of ground, momentum [m]
   real(rkx), pointer :: z0mg(:)
   real(rkx), pointer :: dqgdT(:)  ! temperature derivative of "qg"
   ! latent heat of evaporation (/sublimation) [J/kg]
   real(rkx), pointer :: htvp(:)
   real(rkx), pointer :: emv(:)    ! ground emissivity
   real(rkx), pointer :: emg(:)    ! vegetation emissivity
   real(rkx), pointer :: forc_pbot(:)   ! atmospheric pressure (Pa)
   real(rkx), pointer :: forc_pco2(:)   ! partial pressure co2 (Pa)
   !!! C13
   real(rkx), pointer :: forc_pc13o2(:) ! partial pressure c13o2 (Pa)

   real(rkx), pointer :: forc_po2(:)  ! partial pressure o2 (Pa)
   real(rkx), pointer :: forc_q(:)    ! atmospheric specific humidity (kg/kg)
   ! atmospheric wind speed in east direction (m/s)
   real(rkx), pointer :: forc_u(:)
   ! atmospheric wind speed in north direction (m/s)
   real(rkx), pointer :: forc_v(:)
   !observational height of wind at pft level [m]
   real(rkx), pointer :: forc_hgt_u_pft(:)
   real(rkx), pointer :: forc_rho(:)  ! density (kg/m**3)
   ! downward infrared (longwave) radiation (W/m**2)
   real(rkx), pointer :: forc_lwrad(:)
   real(rkx), pointer :: displa(:)  ! displacement height (m)
   ! one-sided leaf area index with burying by snow
   real(rkx), pointer :: elai(:)
   ! one-sided stem area index with burying by snow
   real(rkx), pointer :: esai(:)
   ! fraction of foliage that is green and dry [-]
   real(rkx), pointer :: fdry(:)
   ! fraction of canopy that is wet (0 to 1)
   real(rkx), pointer :: fwet(:)
   real(rkx), pointer :: laisun(:)      ! sunlit leaf area
   real(rkx), pointer :: laisha(:)      ! shaded leaf area
   ! solar radiation absorbed by vegetation (W/m**2)
   real(rkx), pointer :: sabv(:)
   ! volumetric soil water at saturation (porosity)
   real(rkx), pointer :: watsat(:,:)
   real(rkx), pointer :: watdry(:,:)     ! btran parameter for btran=0
   real(rkx), pointer :: watopt(:,:)     ! btran parameter for btran = 1
   real(rkx), pointer :: h2osoi_ice(:,:) ! ice lens (kg/m2)
   real(rkx), pointer :: h2osoi_liq(:,:) ! liquid water (kg/m2)
   ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] by F. Li and S. Levis
   real(rkx), pointer :: h2osoi_vol(:,:)
   real(rkx), pointer :: dz(:,:)        ! layer depth (m)
   real(rkx), pointer :: t_soisno(:,:)  ! soil temperature (Kelvin)
   real(rkx), pointer :: sucsat(:,:)    ! minimum soil suction (mm)
   real(rkx), pointer :: bsw(:,:)       ! Clapp and Hornberger "b"
   real(rkx), pointer :: rootfr(:,:)    ! fraction of roots in each soil layer
   real(rkx), pointer :: dleaf(:)       ! characteristic leaf dimension (m)
   ! soil water potential at full stomatal opening (mm)
   real(rkx), pointer :: smpso(:)
   ! soil water potential at full stomatal closure (mm)
   real(rkx), pointer :: smpsc(:)
   ! fraction of ground covered by snow (0 to 1)
   real(rkx), pointer :: frac_sno(:)
   real(rkx), pointer :: htop(:)       ! canopy top(m)
   real(rkx), pointer :: snow_depth(:) ! snow height (m)
   ! soil wetness relative to field capacity
   real(rkx), pointer :: soilbeta(:)
   real(rkx), pointer :: lat(:)      ! latitude (radians)
   real(rkx), pointer :: decl(:)     ! declination angle (radians)
   real(rkx), pointer :: max_dayl(:) !maximum daylength for this column (s)
   ! longitude (degrees) (for calculation of local time)
   real(rkx), pointer :: londeg(:)

   ! deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
   real(rkx), pointer :: cgrnds(:)
   ! deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
   real(rkx), pointer :: cgrndl(:)
   real(rkx), pointer :: t_veg(:)   ! vegetation temperature (Kelvin)
   ! 2 m height surface air temperature (Kelvin)
   real(rkx), pointer :: t_ref2m(:)
   ! 2 m height surface specific humidity (kg/kg)
   real(rkx), pointer :: q_ref2m(:)
   ! Rural 2 m height surface air temperature (Kelvin)
   real(rkx), pointer :: t_ref2m_r(:)
   ! 2 m height surface relative humidity (%)
   real(rkx), pointer :: rh_ref2m(:)
   ! Rural 2 m height surface relative humidity (%)
   real(rkx), pointer :: rh_ref2m_r(:)
   real(rkx), pointer :: h2ocan(:)  ! canopy water (mm H2O)

   real(rkx), pointer :: rb1(:)     ! boundary layer resistance (s/m)
   ! deriv. of soil energy flux wrt to soil temp [w/m2/k]
   real(rkx), pointer :: cgrnd(:)
   ! downward longwave radiation below the canopy [W/m2]
   real(rkx), pointer :: dlrad(:)
   ! upward longwave radiation above the canopy [W/m2]
   real(rkx), pointer :: ulrad(:)
   real(rkx), pointer :: ram1(:)    ! aerodynamical resistance (s/m)
   real(rkx), pointer :: btran(:)   ! transpiration wetness factor (0 to 1)
   real(rkx), pointer :: btran2(:)  !F. Li and S. Levis
   real(rkx), pointer :: rssun(:)   ! sunlit stomatal resistance (s/m)
   real(rkx), pointer :: rssha(:)   ! shaded stomatal resistance (s/m)
   real(rkx), pointer :: rhal(:)
   real(rkx), pointer :: vpdal(:)
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
   !!! C13
   ! sunlit leaf photosynthesis (umol 13CO2 /m**2/ s)
   real(rkx), pointer :: c13_psnsun(:)
   ! shaded leaf photosynthesis (umol 13CO2 /m**2/ s)
   real(rkx), pointer :: c13_psnsha(:)
   real(rkx), pointer :: rc13_canair(:) !C13O2/C12O2 in canopy air
   real(rkx), pointer :: rc13_psnsun(:) !C13O2/C12O2 in sunlit canopy psn flux
   real(rkx), pointer :: rc13_psnsha(:) !C13O2/C12O2 in shaded canopy psn flux
   !fractionation factor in sunlit canopy psn flux
   real(rkx), pointer :: alphapsnsun(:)
   !fractionation factor in shaded canopy psn flux
   real(rkx), pointer :: alphapsnsha(:)
   !!! C14
   ! sunlit leaf photosynthesis (umol 14CO2 /m**2/ s)
   real(rkx), pointer :: c14_psnsun(:)
   ! shaded leaf photosynthesis (umol 14CO2 /m**2/ s)
   real(rkx), pointer :: c14_psnsha(:)
   real(rkx), pointer :: rc14_atm(:)  ! C14O2/C12O2 in atmosphere

   ! vegetation transpiration (mm H2O/s) (+ = to atm)
   real(rkx), pointer :: qflx_tran_veg(:)
   ! change in t_veg, last iteration (Kelvin)
   real(rkx), pointer :: dt_veg(:)
   ! vegetation evaporation (mm H2O/s) (+ = to atm)
   real(rkx), pointer :: qflx_evap_veg(:)
   ! sensible heat flux from leaves (W/m**2) [+ to atm]
   real(rkx), pointer :: eflx_sh_veg(:)
   real(rkx), pointer :: taux(:)      ! wind (shear) stress: e-w (kg/m/s**2)
   real(rkx), pointer :: tauy(:)      ! wind (shear) stress: n-s (kg/m/s**2)
   ! sensible heat flux from ground (W/m**2) [+ to atm]
   real(rkx), pointer :: eflx_sh_grnd(:)
   ! soil evaporation (mm H2O/s) (+ = to atm)
   real(rkx), pointer :: qflx_evap_soi(:)
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
   !KO
   ! fractional humidity of canopy air [dimensionless]
   real(rkx), pointer :: rhaf(:)
   !KO

#if (defined LCH4)
   ! tracer conductance for boundary layer [m/s]
   real(rkx), pointer :: grnd_ch4_cond(:)
   !tracer conductance for canopy [m/s]
   real(rkx), pointer :: canopy_cond(:)
#endif
   ! current irrigation rate [mm/s]
   real(rkx), pointer :: irrig_rate(:)
   ! # of time steps for which we still need to irrigate today
   integer(ik4), pointer  :: n_irrig_steps_left(:)
   ! maximum annual depth of thaw
   integer(ik4), pointer  :: altmax_indx(:)
   ! prior year maximum annual depth of thaw
   integer(ik4), pointer  :: altmax_lastyear_indx(:)

   real(rkx), parameter :: btran0 = 0.0_rkx  ! initial value
   ! convective boundary layer height [m]
   real(rkx), parameter :: zii = 1000.0_rkx
   ! coefficient of conective velocity [-]
   real(rkx), parameter :: beta = 1.0_rkx
   ! maxchange in  leaf temperature [K]
   real(rkx), parameter :: delmax = 1.0_rkx
   ! max limit for energy flux convergence [w/m2]
   real(rkx), parameter :: dlemin = 0.1_rkx
   ! max limit for temperature convergence [K]
   real(rkx), parameter :: dtmin = 0.01_rkx
   ! maximum number of iteration [-]
   integer(ik4) , parameter :: itmax = 40
   ! minimum number of iteration [-]
   integer(ik4) , parameter :: itmin = 2
   ! Minimum LAI for irrigation
   real(rkx), parameter :: irrig_min_lai = 0.0_rkx
   ! Irrigate when btran falls below 0.999999 rather than 1 to
   ! allow for round-off error
   real(rkx), parameter :: irrig_btran_thresh = 0.999999_rkx
   ! Time of day to check whether we need irrigation, seconds (0 = midnight).
   ! We start applying the irrigation in the time step FOLLOWING this time,
   ! since we won't begin irrigating until the next call to Hydrology1
   integer(ik4) , parameter :: irrig_start_time = isecspday/4
   ! Desired amount of time to irrigate per day (sec). Actual time may
   ! differ if this is not a multiple of dtsrf. Irrigation won't work properly
   ! if dtsrf > secsperday
   integer(ik4) , parameter :: irrig_length = isecspday/6
   ! Determines target soil moisture level for irrigation. If h2osoi_liq_so
   ! is the soil moisture level at which stomata are fully open and
   ! h2osoi_liq_sat is the soil moisture level at saturation (eff_porosity),
   ! then the target soil moisture level is
   !     (h2osoi_liq_so + irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)).
   ! A value of 0 means that the target soil moisture level is h2osoi_liq_so.
   ! A value of 1 means that the target soil moisture level is h2osoi_liq_sat
   real(rkx), parameter :: irrig_factor = 0.7_rkx

   !added by K.Sakaguchi for litter resistance
   ! placeholder for (dry) plant litter area index (m2/m2)
   real(rkx), parameter :: lai_dl = 0.5_rkx
   ! placeholder for (dry) litter layer thickness (m)
   real(rkx), parameter :: z_dl = 0.05_rkx
   !added by K.Sakaguchi for stability formulation
   ! free parameter for stable formulation
   ! (currently = 0.5, "gamma" in Sakaguchi&Zeng,2008)
   real(rkx), parameter :: ria  = 0.5_rkx
   ! reference height "minus" zero displacement height [m]
   real(rkx) :: zldis(lbp:ubp)
   ! dimensionless height used in Monin-Obukhov theory
   real(rkx) :: zeta
   real(rkx) :: wc             ! convective velocity [m/s]
   ! diff of virtual temp. between ref. height and surface
   real(rkx) :: dth(lbp:ubp)
   ! diff of vir. poten. temp. between ref. height and surface
   real(rkx) :: dthv(lbp:ubp)
   ! diff of humidity between ref. height and surface
   real(rkx) :: dqh(lbp:ubp)
   real(rkx) :: obu(lbp:ubp) ! Monin-Obukhov length (m)
   real(rkx) :: um(lbp:ubp)  ! wind speed including the stablity effect [m/s]
   real(rkx) :: ur(lbp:ubp)  ! wind speed at reference height [m/s]
   real(rkx) :: uaf(lbp:ubp) ! velocity of air within foliage [m/s]
   real(rkx) :: temp1(lbp:ubp)   ! relation for potential temperature profile
   ! relation for potential temperature profile applied at 2-m
   real(rkx) :: temp12m(lbp:ubp)
   real(rkx) :: temp2(lbp:ubp)   ! relation for specific humidity profile
   ! relation for specific humidity profile applied at 2-m
   real(rkx) :: temp22m(lbp:ubp)
   real(rkx) :: ustar(lbp:ubp)   ! friction velocity [m/s]
   real(rkx) :: tstar         ! temperature scaling parameter
   real(rkx) :: qstar         ! moisture scaling parameter
   real(rkx) :: thvstar       ! virtual potential temperature scaling parameter
   real(rkx) :: taf(lbp:ubp)  ! air temperature within canopy space [K]
   real(rkx) :: qaf(lbp:ubp)  ! humidity of canopy air [kg/kg]
   ! fraction of potential evaporation from leaf [-]
   real(rkx) :: rpp
   ! fraction of potential evaporation through transp [-]
   real(rkx) :: rppdry
   ! heat transfer coefficient from leaves [-]
   real(rkx) :: cf
   ! leaf boundary layer resistance [s/m]
   real(rkx) :: rb(lbp:ubp)
   real(rkx) :: rah(lbp:ubp,2)  ! thermal resistance [s/m]
   real(rkx) :: raw(lbp:ubp,2)  ! moisture resistance [s/m]
   real(rkx) :: wta           ! heat conductance for air [m/s]
   real(rkx) :: wtg(lbp:ubp)  ! heat conductance for ground [m/s]
   real(rkx) :: wtl           ! heat conductance for leaf [m/s]
   real(rkx) :: wta0(lbp:ubp) ! normalized heat conductance for air [-]
   real(rkx) :: wtl0(lbp:ubp) ! normalized heat conductance for leaf [-]
   real(rkx) :: wtg0          ! normalized heat conductance for ground [-]
   ! normalized heat conductance for air and leaf [-]
   real(rkx) :: wtal(lbp:ubp)
   real(rkx) :: wtga          ! normalized heat cond. for air and ground  [-]
   real(rkx) :: wtaq          ! latent heat conductance for air [m/s]
   real(rkx) :: wtlq          ! latent heat conductance for leaf [m/s]
   real(rkx) :: wtgq(lbp:ubp) ! latent heat conductance for ground [m/s]
   ! normalized latent heat conductance for air [-]
   real(rkx) :: wtaq0(lbp:ubp)
   ! normalized latent heat conductance for leaf [-]
   real(rkx) :: wtlq0(lbp:ubp)
   real(rkx) :: wtgq0          ! normalized heat conductance for ground [-]
   ! normalized latent heat cond. for air and leaf [-]
   real(rkx) :: wtalq(lbp:ubp)
   ! normalized latent heat cond. for air and ground [-]
   real(rkx) :: wtgaq
   real(rkx) :: el(lbp:ubp)      ! vapor pressure on leaf surface [pa]
   real(rkx) :: deldT            ! derivative of "el" on "t_veg" [pa/K]
   real(rkx) :: qsatl(lbp:ubp)   ! leaf specific humidity [kg/kg]
   real(rkx) :: qsatldT(lbp:ubp) ! derivative of "qsatl" on "t_veg"
   ! 2 m height surface saturated vapor pressure [Pa]
   real(rkx) :: e_ref2m
   ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
   real(rkx) :: de2mdT
   ! 2 m height surface saturated specific humidity [kg/kg]
   real(rkx) :: qsat_ref2m
   ! derivative of 2 m height surface saturated specific humidity on t_ref2m
   real(rkx) :: dqsat2mdT
   ! atmos. radiation temporay set
   real(rkx) :: air(lbp:ubp),bir(lbp:ubp),cir(lbp:ubp)
   real(rkx) :: dc1,dc2          ! derivative of energy flux [W/m2/K]
   real(rkx) :: delt             ! temporary
   real(rkx) :: delq(lbp:ubp)    ! temporary
   ! absolute change in leaf temp in current iteration [K]
   real(rkx) :: del(lbp:ubp)
   ! change in leaf temperature in previous iteration [K]
   real(rkx) :: del2(lbp:ubp)
   ! change in latent heat flux from leaf [K]
   real(rkx) :: dele(lbp:ubp)
   ! change in leaf temperature in current iteration [K]
   real(rkx) :: dels
   ! maximum leaf temp. change in two consecutive iter [K]
   real(rkx) :: det(lbp:ubp)
   ! latent heat flux from leaf (previous iter) [mm/s]
   real(rkx) :: efeb(lbp:ubp)
   ! latent heat flux from leaf (previous iter) [mm/s]
   real(rkx) :: efeold
   ! potential latent energy flux [kg/m2/s]
   real(rkx) :: efpot
   real(rkx) :: efe(lbp:ubp)    ! water flux from leaf [mm/s]
   real(rkx) :: efsh            ! sensible heat from leaf [mm/s]
   real(rkx) :: obuold(lbp:ubp) ! monin-obukhov length from previous iteration
   real(rkx) :: tlbef(lbp:ubp)  ! leaf temperature from previous iteration [K]
   real(rkx) :: ecidif          ! excess energies [W/m2]
   real(rkx) :: err(lbp:ubp)    ! balance error
   real(rkx) :: erre            ! balance error
   real(rkx) :: co2(lbp:ubp)    ! atmospheric co2 partial pressure (pa)
   !!! C13
   real(rkx) :: c13o2(lbp:ubp)  ! atmospheric c13o2 partial pressure (pa)

   real(rkx) :: o2(lbp:ubp)     ! atmospheric o2 partial pressure (pa)
   real(rkx) :: svpts(lbp:ubp)  ! saturation vapor pressure at t_veg (pa)
   real(rkx) :: eah(lbp:ubp)    ! canopy air vapor pressure (pa)
   real(rkx) :: s_node          ! vol_liq/eff_porosity
   real(rkx) :: smp_node        ! matrix potential
   real(rkx) :: smp_node_lf     ! F. Li and S. Levis
   real(rkx) :: vol_ice         ! partial volume of ice lens in layer
   real(rkx) :: eff_porosity    ! effective porosity in layer
   real(rkx) :: vol_liq         ! partial volume of liquid water in layer
   integer(ik4)  :: itlef       ! counter for leaf temperature iteration [-]
   ! number of times stability changes sign
   integer(ik4)  :: nmozsgn(lbp:ubp)
   real(rkx) :: w               ! exp(-LSAI)
   ! interpolated csoilc for less than dense canopies
   real(rkx) :: csoilcn
   real(rkx) :: fm(lbp:ubp)   ! needed for BGC only to diagnose 10m wind speed
   ! sensible heat resistance for air, grnd and leaf [-]
   real(rkx) :: wtshi
   ! latent heat resistance for air, grnd and leaf [-]
   real(rkx) :: wtsqi
   integer(ik4)  :: j          ! soil/snow level index
   integer(ik4)  :: p          ! pft index
   integer(ik4)  :: c          ! column index
   integer(ik4)  :: l          ! landunit index
   integer(ik4)  :: g          ! gridcell index
   integer(ik4)  :: fp         ! lake filter pft index
   integer(ik4)  :: fn         ! number of values in pft filter
   integer(ik4)  :: fnorig     ! number of values in pft filter copy
   integer(ik4)  :: fnold      ! temporary copy of pft count
   integer(ik4)  :: f          ! filter index
   integer(ik4)  :: filterp(ubp-lbp+1)    ! temporary filter
   integer(ik4)  :: fporig(ubp-lbp+1)     ! temporary filter
   real(rkx) :: displa_loc(lbp:ubp)   ! temporary copy
   real(rkx) :: z0mv_loc(lbp:ubp)     ! temporary copy
   real(rkx) :: z0hv_loc(lbp:ubp)     ! temporary copy
   real(rkx) :: z0qv_loc(lbp:ubp)     ! temporary copy
   logical  :: found           ! error flag for canopy above forcing hgt
   integer(ik4)  :: index      ! pft index for error
   real(rkx) :: egvf           ! effective green vegetation fraction
   real(rkx) :: lt     ! elai+esai
   real(rkx) :: ri     ! stability parameter for under canopy air (unitless)
   ! turbulent transfer coefficient over bare soil (unitless)
   real(rkx) :: csoilb
   ! modified transfer coefficient under dense canopy (unitless)
   real(rkx) :: ricsoilc
   ! critical snow depth to cover plant litter (m)
   real(rkx) :: snow_depth_c
   ! dry litter layer resistance for water vapor  (s/m)
   real(rkx) :: rdl
   real(rkx) :: elai_dl       ! exposed (dry) plant litter area index
   real(rkx) :: fsno_dl       ! effective snow cover over plant litter
   real(rkx) :: dayl          ! daylength (s)
   real(rkx) :: temp          ! temporary, for daylength calculation
   real(rkx) :: fpeav         ! temporary, for avoid fpe
   ! scalar (0-1) for daylength effect on Vcmax
   real(rkx) :: dayl_factor(lbp:ubp)
   ! Rootfraction defined for unfrozen layers only.
   ! If no unfrozen layers, put all in the top layer.
   real(rkx) :: rootfr_unf(lbp:ubp,1:nlevgrnd)
   real(rkx) :: rootsum(lbp:ubp)
   real(rkx) :: delt_snow
   real(rkx) :: delt_soil
   real(rkx) :: delt_h2osfc
   real(rkx) :: lw_grnd
   real(rkx) :: delq_snow
   real(rkx) :: delq_soil
   real(rkx) :: delq_h2osfc
   ! time UTC at start of time step (seconds after solar midnight)
   integer(ik4)  :: time
   ! local time at start of time step (seconds after solar midnight)
   integer(ik4)  :: local_time
   ! number of time steps per day in which we irrigate
   integer(ik4)  :: irrig_nsteps_per_day
   ! where do we need to check soil moisture to see if we need to irrigate?
   logical  :: check_for_irrig(lbp:ubp)
   ! set to true if we have encountered a frozen soil layer
   logical  :: frozen_soil(lbc:ubc)
   ! partial volume of liquid water in layer for which smp_node = smpso
   real(rkx) :: vol_liq_so
   ! liquid water corresponding to vol_liq_so for this layer [kg/m2]
   real(rkx) :: h2osoi_liq_so
   ! liquid water corresponding to eff_porosity for this layer [kg/m2]
   real(rkx) :: h2osoi_liq_sat
   ! difference between desired soil moisture level for this layer and
   ! current soil moisture level [kg/m2]
   real(rkx) :: deficit

   ! Assign local pointers to derived type members (gridcell-level)

   eflx_sh_snow   => clm3%g%l%c%p%pef%eflx_sh_snow
   eflx_sh_soil   => clm3%g%l%c%p%pef%eflx_sh_soil
   eflx_sh_h2osfc => clm3%g%l%c%p%pef%eflx_sh_h2osfc
   snl            => clm3%g%l%c%cps%snl
   t_h2osfc       => clm3%g%l%c%ces%t_h2osfc
   frac_h2osfc    => clm3%g%l%c%cps%frac_h2osfc
   qg_snow        => clm3%g%l%c%cws%qg_snow
   qg_soil        => clm3%g%l%c%cws%qg_soil
   qg_h2osfc      => clm3%g%l%c%cws%qg_h2osfc
   qflx_ev_snow   => clm3%g%l%c%p%pwf%qflx_ev_snow
   qflx_ev_soil   => clm3%g%l%c%p%pwf%qflx_ev_soil
   qflx_ev_h2osfc => clm3%g%l%c%p%pwf%qflx_ev_h2osfc

   forc_lwrad     => clm_a2l%forc_lwrad
   forc_pco2      => clm_a2l%forc_pco2
   if ( use_c13 ) then
     forc_pc13o2    => clm_a2l%forc_pc13o2
   end if
   forc_po2       => clm_a2l%forc_po2
   forc_q         => clm_a2l%forc_q
   forc_pbot      => clm_a2l%forc_pbot
   forc_u         => clm_a2l%forc_u
   forc_v         => clm_a2l%forc_v
   forc_th        => clm_a2l%forc_th
   forc_rho       => clm_a2l%forc_rho
   lat            => clm3%g%lat
   londeg         => clm3%g%londeg

   ! Assign local pointers to derived type members (column-level)

   t_soisno       => clm3%g%l%c%ces%t_soisno
   watsat         => clm3%g%l%c%cps%watsat
   watdry         => clm3%g%l%c%cps%watdry
   watopt         => clm3%g%l%c%cps%watopt
   h2osoi_ice     => clm3%g%l%c%cws%h2osoi_ice
   h2osoi_vol     => clm3%g%l%c%cws%h2osoi_vol
   dz             => clm3%g%l%c%cps%dz
   h2osoi_liq     => clm3%g%l%c%cws%h2osoi_liq
   sucsat         => clm3%g%l%c%cps%sucsat
   bsw            => clm3%g%l%c%cps%bsw
   emg            => clm3%g%l%c%cps%emg
   t_grnd         => clm3%g%l%c%ces%t_grnd
   qg             => clm3%g%l%c%cws%qg
   thv            => clm3%g%l%c%ces%thv
   dqgdT          => clm3%g%l%c%cws%dqgdT
   htvp           => clm3%g%l%c%cps%htvp
   z0mg           => clm3%g%l%c%cps%z0mg
   frac_sno       => clm3%g%l%c%cps%frac_sno_eff
   snow_depth         => clm3%g%l%c%cps%snow_depth
   soilbeta       => clm3%g%l%c%cws%soilbeta
   decl           => clm3%g%l%c%cps%decl
   max_dayl       => clm3%g%l%c%cps%max_dayl

   ! Assign local pointers to derived type members (pft-level)

   rb1            => clm3%g%l%c%p%pps%rb1
   ivt            => clm3%g%l%c%p%itype
   pcolumn        => clm3%g%l%c%p%column
   plandunit      => clm3%g%l%c%p%landunit
   pgridcell      => clm3%g%l%c%p%gridcell
   frac_veg_nosno => clm3%g%l%c%p%pps%frac_veg_nosno
   btran          => clm3%g%l%c%p%pps%btran
   btran2          => clm3%g%l%c%p%pps%btran2
   rootfr         => clm3%g%l%c%p%pps%rootfr
   rootr          => clm3%g%l%c%p%pps%rootr
   rresis         => clm3%g%l%c%p%pps%rresis
   emv            => clm3%g%l%c%p%pps%emv
   t_veg          => clm3%g%l%c%p%pes%t_veg
   displa         => clm3%g%l%c%p%pps%displa
   z0mv           => clm3%g%l%c%p%pps%z0mv
   z0hv           => clm3%g%l%c%p%pps%z0hv
   z0qv           => clm3%g%l%c%p%pps%z0qv
   ram1           => clm3%g%l%c%p%pps%ram1
   htop           => clm3%g%l%c%p%pps%htop
   rssun          => clm3%g%l%c%p%pps%rssun
   rssha          => clm3%g%l%c%p%pps%rssha
   rhal           => clm3%g%l%c%p%pps%rhal
   vpdal          => clm3%g%l%c%p%pps%vpdal
   psnsun         => clm3%g%l%c%p%pcf%psnsun
   psnsun_wc      => clm3%g%l%c%p%pcf%psnsun_wc
   psnsun_wj      => clm3%g%l%c%p%pcf%psnsun_wj
   psnsun_wp      => clm3%g%l%c%p%pcf%psnsun_wp
   psnsha         => clm3%g%l%c%p%pcf%psnsha
   psnsha_wc      => clm3%g%l%c%p%pcf%psnsha_wc
   psnsha_wj      => clm3%g%l%c%p%pcf%psnsha_wj
   psnsha_wp      => clm3%g%l%c%p%pcf%psnsha_wp
   if ( use_c13 ) then
     c13_psnsun     => clm3%g%l%c%p%pc13f%psnsun
     c13_psnsha     => clm3%g%l%c%p%pc13f%psnsha
     rc13_canair    => clm3%g%l%c%p%pepv%rc13_canair
     rc13_psnsun    => clm3%g%l%c%p%pepv%rc13_psnsun
     rc13_psnsha    => clm3%g%l%c%p%pepv%rc13_psnsha
     alphapsnsun    => clm3%g%l%c%p%pps%alphapsnsun
     alphapsnsha    => clm3%g%l%c%p%pps%alphapsnsha
   end if
   if ( use_c14 ) then
     c14_psnsun     => clm3%g%l%c%p%pc14f%psnsun
     c14_psnsha     => clm3%g%l%c%p%pc14f%psnsha
     rc14_atm       => clm3%g%l%c%p%pepv%rc14_atm
   end if
   elai           => clm3%g%l%c%p%pps%elai
   esai           => clm3%g%l%c%p%pps%esai
   fdry           => clm3%g%l%c%p%pps%fdry
   laisun         => clm3%g%l%c%p%pps%laisun
   laisha         => clm3%g%l%c%p%pps%laisha
   qflx_tran_veg  => clm3%g%l%c%p%pwf%qflx_tran_veg
   fwet           => clm3%g%l%c%p%pps%fwet
   h2ocan         => clm3%g%l%c%p%pws%h2ocan
   dt_veg         => clm3%g%l%c%p%pps%dt_veg
   sabv           => clm3%g%l%c%p%pef%sabv
   qflx_evap_veg  => clm3%g%l%c%p%pwf%qflx_evap_veg
   eflx_sh_veg    => clm3%g%l%c%p%pef%eflx_sh_veg
   taux           => clm3%g%l%c%p%pmf%taux
   tauy           => clm3%g%l%c%p%pmf%tauy
   eflx_sh_grnd   => clm3%g%l%c%p%pef%eflx_sh_grnd
   qflx_evap_soi  => clm3%g%l%c%p%pwf%qflx_evap_soi
   t_ref2m        => clm3%g%l%c%p%pes%t_ref2m
   q_ref2m        => clm3%g%l%c%p%pes%q_ref2m
   t_ref2m_r      => clm3%g%l%c%p%pes%t_ref2m_r
   rh_ref2m_r     => clm3%g%l%c%p%pes%rh_ref2m_r
   rh_ref2m       => clm3%g%l%c%p%pes%rh_ref2m
   dlrad          => clm3%g%l%c%p%pef%dlrad
   ulrad          => clm3%g%l%c%p%pef%ulrad
   cgrnds         => clm3%g%l%c%p%pef%cgrnds
   cgrndl         => clm3%g%l%c%p%pef%cgrndl
   cgrnd          => clm3%g%l%c%p%pef%cgrnd
   fpsn           => clm3%g%l%c%p%pcf%fpsn
   fpsn_wc        => clm3%g%l%c%p%pcf%fpsn_wc
   fpsn_wj        => clm3%g%l%c%p%pcf%fpsn_wj
   fpsn_wp        => clm3%g%l%c%p%pcf%fpsn_wp
   forc_hgt_u_pft => clm3%g%l%c%p%pps%forc_hgt_u_pft
   thm            => clm3%g%l%c%p%pes%thm
!KO
   rhaf           => clm3%g%l%c%p%pps%rhaf
!KO
#if (defined LCH4)
   grnd_ch4_cond  => clm3%g%l%c%p%pps%grnd_ch4_cond
   canopy_cond    => clm3%g%l%c%p%pps%canopy_cond
#endif
   irrig_rate           => clm3%g%l%c%cps%irrig_rate
   n_irrig_steps_left   => clm3%g%l%c%cps%n_irrig_steps_left
   altmax_lastyear_indx => clm3%g%l%c%cps%altmax_lastyear_indx
   altmax_indx          => clm3%g%l%c%cps%altmax_indx

   ! Assign local pointers to derived type members (ecophysiological)

   dleaf          => pftcon%dleaf
   smpso          => pftcon%smpso
   smpsc          => pftcon%smpsc

   ! Determine step size

   irrig_nsteps_per_day = int((dble(irrig_length) + &
                                   (dtsrf-1.0_rkx))/dtsrf)  ! round up

   ! Filter pfts where frac_veg_nosno is non-zero

   fn = 0
   do fp = 1 , num_nolakep
     p = filter_nolakep(fp)
     if (frac_veg_nosno(p) /= 0) then
       fn = fn + 1
       filterp(fn) = p
     end if
   end do

   ! Initialize

   do f = 1, fn
     p = filterp(f)
     del(p)    = 0._rkx  ! change in leaf temperature from previous iteration
     efeb(p)   = 0._rkx  ! latent head flux from leaf for previous iteration
     wtlq0(p)  = 0._rkx
     wtalq(p)  = 0._rkx
     wtgq(p)   = 0._rkx
     wtaq0(p)  = 0._rkx
     obuold(p) = 0._rkx
     btran(p)  = btran0
     btran2(p)  = btran0
   end do

   ! calculate daylength control for Vcmax
   do f = 1, fn
     p = filterp(f)
     c = pcolumn(p)
     g = pgridcell(p)
     ! calculate daylength
     temp = -(sin(lat(g))*sin(decl(c)))/(cos(lat(g)) * cos(decl(c)))
     temp = min(1._rkx,max(-1._rkx,temp))
     dayl = 2.0_rkx * 13750.9871_rkx * acos(temp)
     ! calculate dayl_factor as the ratio of (current:max dayl)^2
     ! set a minimum of 0.01 (1%) for the dayl_factor
     dayl_factor(p) = min(1._rkx,max(0.01_rkx, &
             (dayl*dayl)/(max_dayl(c)*max_dayl(c))))
   end do

   rb1(lbp:ubp) = 0._rkx

   ! Define rootfraction for unfrozen soil only
   if (perchroot .or. perchroot_alt) then
     if (perchroot_alt) then
       ! use total active layer
       ! (defined ass max thaw depth for current and prior year)
       do j = 1,nlevgrnd
         do f = 1, fn
           p = filterp(f)
           c = pcolumn(p)
           if ( j <= max(altmax_lastyear_indx(c), altmax_indx(c), 1) ) then
             rootfr_unf(p,j) = rootfr(p,j)
           else
             rootfr_unf(p,j) = 0._rkx
           end if
         end do
       end do
     else
       ! use instantaneous temperature
       do j = 1,nlevgrnd
         do f = 1, fn
           p = filterp(f)
           c = pcolumn(p)
           if (t_soisno(c,j) >= tfrz) then
             rootfr_unf(p,j) = rootfr(p,j)
           else
             rootfr_unf(p,j) = 0._rkx
           end if
         end do
       end do
     end if ! perchroot_alt

     ! sum unfrozen roots
     do j = 1,nlevgrnd
       do f = 1, fn
         p = filterp(f)
         c = pcolumn(p)
         if (j == 1) rootsum(p) = 0._rkx
         rootsum(p) = rootsum(p) + rootfr_unf(p,j)
       end do
     end do

     ! normalize rootfr to total unfrozen depth
     do j = 1,nlevgrnd
       do f = 1, fn
         p = filterp(f)
         c = pcolumn(p)
         if (rootsum(p) > 0._rkx) then
           rootfr_unf(p,j) = rootfr_unf(p,j) / rootsum(p)
         end if
       end do
     end do
   end if ! perchroot

   ! Effective porosity of soil, partial volume of ice and liquid
   ! (needed for btran) and root resistance factors

   do j = 1 , nlevgrnd
     do f = 1 , fn
       p = filterp(f)
       c = pcolumn(p)
       l = plandunit(p)

       ! Root resistance factors

       vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
       eff_porosity = watsat(c,j)-vol_ice
       vol_liq = min(eff_porosity, h2osoi_liq(c,j)/(dz(c,j)*denh2o))
       if ( vol_liq <= 0._rkx .or. t_soisno(c,j) <= tfrz-2._rkx ) then
         rootr(p,j) = 0._rkx
       else
         s_node = max(vol_liq/eff_porosity,0.01_rkx)
         smp_node = max(smpsc(ivt(p)), -sucsat(c,j)*s_node**(-bsw(c,j)))

         rresis(p,j) = min( (eff_porosity/watsat(c,j))* &
                        (smp_node - smpsc(ivt(p))) / &
                        (smpso(ivt(p)) - smpsc(ivt(p))), 1._rkx)
         if (.not. (perchroot .or. perchroot_alt) ) then
           rootr(p,j) = rootfr(p,j)*rresis(p,j)
         else
           rootr(p,j) = rootfr_unf(p,j)*rresis(p,j)
         end if
         btran(p)    = btran(p) + rootr(p,j)
         smp_node_lf = max(smpsc(ivt(p)), &
                 -sucsat(c,j)*(h2osoi_vol(c,j)/watsat(c,j))**(-bsw(c,j)))
         btran2(p)   = btran2(p) + &
                 rootfr(p,j)*min((smp_node_lf - smpsc(ivt(p))) / &
                 (smpso(ivt(p)) - smpsc(ivt(p))), 1._rkx)
       endif
     end do
   end do

   ! Normalize root resistances to get layer contribution to ET

   do j = 1 , nlevgrnd
     do f = 1, fn
       p = filterp(f)
       if ( btran(p) > btran0 ) then
         rootr(p,j) = rootr(p,j)/btran(p)
       else
         rootr(p,j) = 0._rkx
       end if
     end do
   end do

   ! Determine if irrigation is needed (over irrigated soil columns)

   ! First, determine in what grid cells we need to bother 'measuring'
   ! soil water, to see if we need irrigation
   ! Also set n_irrig_steps_left for these grid cells
   ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
   ! in this case, we'll irrigate by 0 for the given number of time steps
   ! get time as of beginning of time step
   time = idatex%second_of_day
   do f = 1, fn
     p = filterp(f)
     c = pcolumn(p)
     g = pgridcell(p)
     if (irrigated(ivt(p)) == 1._rkx .and. &
         elai(p) > irrig_min_lai .and.   &
         btran(p) < irrig_btran_thresh) then
       ! see if it's the right time of day to start irrigating:
       local_time = modulo(time + nint(londeg(g)/degpsec), isecspday)
       if (modulo(local_time - irrig_start_time, isecspday) < dtsrf) then
         ! it's time to start irrigating
         check_for_irrig(p)    = .true.
         n_irrig_steps_left(c) = irrig_nsteps_per_day
         irrig_rate(c)         = 0._rkx  ! reset; we'll add to this later
       else
         check_for_irrig(p)    = .false.
       end if
     else  ! non-irrig pft or elai<=irrig_min_lai or btran>irrig_btran_thresh
       check_for_irrig(p)       = .false.
     end if
   end do

   ! Now 'measure' soil water for the grid cells identified above and see
   ! if the soil is dry enough to warrant irrigation
   frozen_soil(:) = .false.
   do j = 1,nlevgrnd
     do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       if ( check_for_irrig(p) .and. .not. frozen_soil(c) ) then
         ! if level L was frozen, then we don't look at any levels below L
         if ( t_soisno(c,j) <= tfrz ) then
           frozen_soil(c) = .true.
         else if ( rootfr(p,j) > 0.0_rkx ) then
           ! determine soil water deficit in this layer:
           ! Calculate vol_liq_so - i.e., vol_liq at which smp_node = smpso
           ! - by inverting the above equations for the root resistance factors
           ! this duplicates the above equation for vol_ice
           vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
           ! this duplicates the above equation for eff_porosity
           eff_porosity = watsat(c,j)-vol_ice
           vol_liq_so   = eff_porosity * &
                   (-smpso(ivt(p))/sucsat(c,j))**(-1.0_rkx/bsw(c,j))
           ! Translate vol_liq_so and eff_porosity into h2osoi_liq_so and
           ! h2osoi_liq_sat and calculate deficit
           h2osoi_liq_so  = vol_liq_so * denh2o * dz(c,j)
           h2osoi_liq_sat = eff_porosity * denh2o * dz(c,j)
           deficit = max((h2osoi_liq_so + &
                   irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)) - &
                   h2osoi_liq(c,j), 0._rkx)
           ! Add deficit to irrig_rate, converting units from mm to mm/sec
           irrig_rate(c) = irrig_rate(c)+deficit/(dtsrf*irrig_nsteps_per_day)

         end if  ! else if (rootfr(p,j) > 0)
       end if    ! if (check_for_irrig(p) .and. .not. frozen_soil(c))
     end do      ! do f
   end do        ! do j

   ! Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
   do f = 1, fn
     p = filterp(f)
     c = pcolumn(p)

     lt = min(elai(p)+esai(p), tlsai_crit)
     egvf = (1._rkx - alpha_aero * exp(-lt)) / &
            (1._rkx - alpha_aero * exp(-tlsai_crit))
     displa(p) = egvf * displa(p)
     z0mv(p)   = exp(egvf * log(z0mv(p)) + (1._rkx - egvf) * log(z0mg(c)))
     z0hv(p)   = z0mv(p)
     z0qv(p)   = z0mv(p)
   end do

   found = .false.
   do f = 1, fn
     p = filterp(f)
     c = pcolumn(p)
     g = pgridcell(p)

     ! Net absorbed longwave radiation by canopy and ground
     ! =air+bir*t_veg**4+cir*t_grnd(c)**4

     air(p) =   emv(p) * (1._rkx+(1._rkx-emv(p))*(1._rkx-emg(c))) * forc_lwrad(g)
     bir(p) = - (2._rkx-emv(p)*(1._rkx-emg(c))) * emv(p) * sb
     cir(p) =   emv(p)*emg(c)*sb

     ! Saturated vapor pressure, specific humidity, and their derivatives
     ! at the leaf surface

     call QSat (t_veg(p), forc_pbot(g), el(p), deldT, qsatl(p), qsatldT(p))

     ! Determine atmospheric co2 and o2

     co2(p) = forc_pco2(g)
     o2(p)  = forc_po2(g)

     if ( use_c13 ) then
       c13o2(p) = forc_pc13o2(g)
     end if

     ! Initialize flux profile

     nmozsgn(p) = 0

     taf(p) = (t_grnd(c) + thm(p)) / 2._rkx
     qaf(p) = (forc_q(g) + qg(c)) / 2._rkx

     ur(p) = max(1.0_rkx,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
     dth(p) = thm(p) - taf(p)
     dqh(p) = forc_q(g) - qaf(p)
     delq(p) = qg(c) - qaf(p)
     dthv(p) = dth(p)*(1._rkx+0.61_rkx*forc_q(g))+0.61_rkx*forc_th(g)*dqh(p)
     zldis(p) = forc_hgt_u_pft(p) - displa(p)

     ! Check to see if the forcing height is below the canopy height
     if (zldis(p) < 0._rkx) then
       found = .true.
       index = p
       write(stderr,*) 'At pft index ', p
       write(stderr,*) 'Canopy height : ', displa(p)
       write(stderr,*) 'Forcing model : ', forc_hgt_u_pft(p)
     end if
   end do

   if (found) then
     write(stderr,*) &
           'Error: Forcing height is below canopy height for pft index ',index
     call fatal(__FILE__,__LINE__,'clm now stopping')
   end if

   do f = 1, fn
     p = filterp(f)
     c = pcolumn(p)
     ! Initialize Monin-Obukhov length and wind speed
     call MoninObukIni(ur(p),thv(c),dthv(p),zldis(p),z0mv(p),um(p),obu(p))
   end do

   ! Set counter for leaf temperature iteration (itlef)

   itlef = 0
   fnorig = fn
   fporig(1:fn) = filterp(1:fn)

   ! Make copies so that array sections are not passed in function
   ! calls to friction velocity

   do f = 1, fn
     p = filterp(f)
     displa_loc(p) = displa(p)
     z0mv_loc(p) = z0mv(p)
     z0hv_loc(p) = z0hv(p)
     z0qv_loc(p) = z0qv(p)
   end do

   ! Begin stability iteration

   ITERATION : do while (itlef <= itmax .and. fn > 0)

     ! Determine friction velocity, and potential temperature and humidity
     ! profiles of the surface boundary layer

     call FrictionVelocity (lbp, ubp, fn, filterp, &
                            displa_loc, z0mv_loc, z0hv_loc, z0qv_loc, &
                            obu, itlef+1, ur, um, ustar, &
                            temp1, temp2, temp12m, temp22m, fm)

     do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)

       tlbef(p) = t_veg(p)
       del2(p) = del(p)

       ! Determine aerodynamic resistances

       ram1(p)  = 1._rkx/(ustar(p)*ustar(p)/um(p))
       rah(p,1) = 1._rkx/(temp1(p)*ustar(p))
       raw(p,1) = 1._rkx/(temp2(p)*ustar(p))

       ! Bulk boundary layer resistance of leaves

       uaf(p) = um(p)*sqrt( 1._rkx/(ram1(p)*um(p)) )
       cf  = 0.01_rkx/(sqrt(uaf(p))*sqrt(dleaf(ivt(p))))
       rb(p)  = 1._rkx/(cf*uaf(p))
       rb1(p) = rb(p)

       ! Parameterization for variation of csoilc with canopy density from
       ! X. Zeng, University of Arizona

       w = exp(-(elai(p)+esai(p)))

       ! changed by K.Sakaguchi from here
       ! transfer coefficient over bare soil is changed to a local variable
       ! just for readability of the code (from line 680)
       csoilb = (vkc/(0.13_rkx*(z0mg(c)*uaf(p)/1.5e-5_rkx)**0.45_rkx))

       !compute the stability parameter for ricsoilc
       !  ("S" in Sakaguchi&Zeng,2008)

       ri = ( grav*htop(p)*(taf(p) - t_grnd(c)) )/(taf(p) * uaf(p)**2.00_rkx)

       !! modify csoilc value (0.004) if the under-canopy is in stable condition

       if ( (taf(p) - t_grnd(c) ) > 0._rkx) then
         ! decrease the value of csoilc by dividing it with
         ! (1+gamma*min(S, 10.0))
         ! ria ("gmanna" in Sakaguchi&Zeng, 2008) is a constant (=0.5)
         ricsoilc = csoilc / (1.00_rkx + ria*min( ri, 10.0_rkx) )
         csoilcn = csoilb*w + ricsoilc*(1._rkx-w)
       else
         csoilcn = csoilb*w + csoilc*(1._rkx-w)
       end if

       !! Sakaguchi changes for stability formulation ends here

       rah(p,2) = 1._rkx/(csoilcn*uaf(p))
       raw(p,2) = rah(p,2)
#if (defined LCH4)
       grnd_ch4_cond(p) = 1._rkx/(raw(p,1)+raw(p,2))
#endif

       ! Stomatal resistances for sunlit and shaded fractions of canopy.
       ! Done each iteration to account for differences in eah, tv.

       svpts(p) = el(p)                         ! pa
       eah(p) = forc_pbot(g) * qaf(p) / 0.622_rkx ! pa
!KO
       rhaf(p) = eah(p)/svpts(p)
!KO
       rhal(p) = rhaf(p)
       vpdal(p) = svpts(p) - eah(p)
     end do

     call Photosynthesis (fn, filterp, lbp, ubp, svpts, eah, &
             o2, co2, rb, dayl_factor, phase='sun')
     if ( use_c13 ) then
       call Fractionation (lbp, ubp, fn, filterp, phase='sun')
     endif
     call Photosynthesis (fn, filterp, lbp, ubp, svpts, eah, &
             o2, co2, rb, dayl_factor, phase='sha')
     if ( use_c13 ) then
       call Fractionation (lbp, ubp, fn, filterp, phase='sha')
     end if
     do f = 1, fn
       p = filterp(f)
       c = pcolumn(p)
       g = pgridcell(p)

       ! Sensible heat conductance for air, leaf and ground
       ! Moved the original subroutine in-line...

       wta    = 1._rkx/rah(p,1)             ! air
       wtl    = (elai(p)+esai(p))/rb(p)   ! leaf
       wtg(p) = 1._rkx/rah(p,2)             ! ground
       wtshi  = 1._rkx/(wta+wtl+wtg(p))

       wtl0(p) = wtl*wtshi         ! leaf
       wtg0    = wtg(p)*wtshi      ! ground
       wta0(p) = wta*wtshi         ! air

       wtga    = wta0(p)+wtg0      ! ground + air
       wtal(p) = wta0(p)+wtl0(p)   ! air + leaf

       ! Fraction of potential evaporation from leaf

       if ( fdry(p) > 0.0_rkx ) then
         rppdry  = fdry(p)*rb(p)*(laisun(p)/(rb(p)+rssun(p)) + &
                                  laisha(p)/(rb(p)+rssha(p)))/elai(p)
       else
         rppdry = 0.0_rkx
       end if

#if (defined LCH4)
       ! Calculate canopy conductance for methane / oxygen
       ! (e.g. stomatal conductance & leaf bdy cond)
       canopy_cond(p) = (laisun(p)/(rb(p)+rssun(p)) + &
               laisha(p)/(rb(p)+rssha(p)))/max(elai(p), 0.01_rkx)
#endif
       efpot = forc_rho(g)*wtl*(qsatl(p)-qaf(p))

       if ( efpot > 0._rkx ) then
         if ( btran(p) > btran0 ) then
           qflx_tran_veg(p) = efpot * rppdry
           rpp = rppdry + fwet(p)
         else
           ! No transpiration if btran below 1.e-10
           rpp = fwet(p)
           qflx_tran_veg(p) = 0._rkx
         end if
         ! Check total evapotranspiration from leaves
         rpp = min(rpp, (qflx_tran_veg(p)+h2ocan(p)/dtsrf)/efpot)
       else
         !No transpiration if potential evaporation less than zero
         rpp = 1._rkx
         qflx_tran_veg(p) = 0._rkx
       end if

       ! Update conductances for changes in rpp
       ! Latent heat conductances for ground and leaf.
       ! Air has same conductance for both sensible and latent heat.
       ! Moved the original subroutine in-line...

       wtaq = frac_veg_nosno(p)/raw(p,1)                        ! air
       wtlq = frac_veg_nosno(p)*(elai(p)+esai(p))/rb(p) * rpp   ! leaf

       !Litter layer resistance. Added by K.Sakaguchi
       ! critical depth for 100% litter burial by snow (=litter thickness)
       snow_depth_c = z_dl
       ! effective snow cover for (dry)plant litter
       fsno_dl = snow_depth(c)/snow_depth_c
       ! exposed (dry)litter area index
       elai_dl = lai_dl*(1._rkx - min(fsno_dl,1._rkx))
       ! dry litter layer resistance
       rdl = ( 1._rkx - exp(-elai_dl) ) / ( 0.004_rkx*uaf(p))

       ! add litter resistance and Lee and Pielke 1992 beta
       if (delq(p) < 0._rkx) then
         !dew. Do not apply beta for negative flux (follow old rsoil)
         wtgq(p) = frac_veg_nosno(p)/(raw(p,2)+rdl)
       else
         wtgq(p) = soilbeta(c)*frac_veg_nosno(p)/(raw(p,2)+rdl)
       end if

       wtsqi   = 1._rkx/(wtaq+wtlq+wtgq(p))

       wtgq0    = wtgq(p)*wtsqi      ! ground
       wtlq0(p) = wtlq*wtsqi         ! leaf
       wtaq0(p) = wtaq*wtsqi         ! air

       wtgaq    = wtaq0(p)+wtgq0     ! air + ground
       wtalq(p) = wtaq0(p)+wtlq0(p)  ! air + leaf

       dc1 = forc_rho(g)*cpair*wtl
       dc2 = hvap*forc_rho(g)*wtlq

       efsh   = dc1*(wtga*t_veg(p)-wtg0*t_grnd(c)-wta0(p)*thm(p))
       efe(p) = dc2*(wtgaq*qsatl(p)-wtgq0*qg(c)-wtaq0(p)*forc_q(g))

       ! Evaporation flux from foliage

       erre = 0._rkx

       if ( abs(efe(p)) < 1.e-20_rkx ) efe(p)  = 0.0_rkx
       if ( abs(efeb(p)) < 1.e-20_rkx ) efeb(p) = 0.0_rkx
       if (efe(p)*efeb(p) < 0._rkx) then
         efeold = efe(p)
         efe(p) = 0.1_rkx*efeold
         erre = efe(p) - efeold
       end if
       ! fractionate ground emitted longwave
       lw_grnd = (frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
              +(1._rkx-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
              +frac_h2osfc(c)*t_h2osfc(c)**4)

       dt_veg(p) = (sabv(p) + air(p) + bir(p)*t_veg(p)**4 + &
              cir(p)*lw_grnd - efsh - efe(p)) / &
              (- 4._rkx*bir(p)*t_veg(p)**3 +dc1*wtga +dc2*wtgaq*qsatldT(p))
       t_veg(p) = tlbef(p) + dt_veg(p)
       dels = dt_veg(p)
       del(p)  = abs(dels)
       err(p) = 0._rkx
       if (del(p) > delmax) then
         dt_veg(p) = delmax*dels/del(p)
         t_veg(p) = tlbef(p) + dt_veg(p)
         err(p) = sabv(p) + air(p) + bir(p)*tlbef(p)**3*(tlbef(p) + &
                 4._rkx*dt_veg(p)) + cir(p)*lw_grnd - &
                 (efsh + dc1*wtga*dt_veg(p)) - (efe(p) + &
                 dc2*wtgaq*qsatldT(p)*dt_veg(p))
       end if

       ! Fluxes from leaves to canopy space
       ! "efe" was limited as its sign changes frequently.  This limit may
       ! result in an imbalance in "hvap*qflx_evap_veg" and
       ! "efe + dc2*wtgaq*qsatdt_veg"

       efpot = forc_rho(g)*wtl*(wtgaq*(qsatl(p)+qsatldT(p)*dt_veg(p)) &
            -wtgq0*qg(c)-wtaq0(p)*forc_q(g))
       qflx_evap_veg(p) = rpp*efpot

       ! Calculation of evaporative potentials (efpot) and
       ! interception losses; flux in kg m**-2 s-1.  ecidif
       ! holds the excess energy if all intercepted water is evaporated
       ! during the timestep.  This energy is later added to the
       ! sensible heat flux.

       ecidif = 0._rkx
       if (efpot > 0._rkx .and. btran(p) > btran0) then
         qflx_tran_veg(p) = efpot*rppdry
       else
         qflx_tran_veg(p) = 0._rkx
       end if
       if ( h2ocan(p) > 1.0e-20_rkx ) then
         fpeav = h2ocan(p)/dtsrf
       else
         fpeav = 0.0_rkx
       end if
       ecidif = max(0._rkx, qflx_evap_veg(p)-qflx_tran_veg(p)-fpeav)
       qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+fpeav)

       ! The energy loss due to above two limits is added to
       ! the sensible heat flux.

       eflx_sh_veg(p) = efsh + dc1*wtga*dt_veg(p) + &
               err(p) + erre + hvap*ecidif

       ! Re-calculate saturated vapor pressure, specific humidity, and their
       ! derivatives at the leaf surface

       call QSat(t_veg(p), forc_pbot(g), el(p), deldT, qsatl(p), qsatldT(p))

       ! Update vegetation/ground surface temperature, canopy air
       ! temperature, canopy vapor pressure, aerodynamic temperature, and
       ! Monin-Obukhov stability parameter for next iteration.

       taf(p) = wtg0*t_grnd(c) + wta0(p)*thm(p) + wtl0(p)*t_veg(p)
       qaf(p) = wtlq0(p)*qsatl(p) + wtgq0*qg(c) + forc_q(g)*wtaq0(p)

       ! Update Monin-Obukhov length and wind speed including the
       ! stability effect

       dth(p) = thm(p)-taf(p)
       dqh(p) = forc_q(g)-qaf(p)
       delq(p) = wtalq(p)*qg(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(g)

       tstar = temp1(p)*dth(p)
       qstar = temp2(p)*dqh(p)

       thvstar = tstar*(1._rkx+0.61_rkx*forc_q(g)) + 0.61_rkx*forc_th(g)*qstar
       zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

       if (zeta >= 0._rkx) then     !stable
         zeta = min(2._rkx,max(zeta,0.01_rkx))
         um(p) = max(ur(p),0.1_rkx)
       else                     !unstable
         zeta = max(-100._rkx,min(zeta,-0.01_rkx))
         wc = beta*(-grav*ustar(p)*thvstar*zii/thv(c))**0.333_rkx
         um(p) = sqrt(ur(p)*ur(p)+wc*wc)
       end if
       obu(p) = zldis(p)/zeta

       if (obuold(p)*obu(p) < 0._rkx) nmozsgn(p) = nmozsgn(p)+1
       if (nmozsgn(p) >= 4) obu(p) = zldis(p)/(-0.01_rkx)
       obuold(p) = obu(p)

     end do   ! end of filtered pft loop

     ! Test for convergence

     itlef = itlef+1
     if (itlef > itmin) then
       do f = 1, fn
         p = filterp(f)
         dele(p) = abs(efe(p)-efeb(p))
         efeb(p) = efe(p)
         det(p)  = max(del(p),del2(p))
       end do
       fnold = fn
       fn = 0
       do f = 1, fnold
         p = filterp(f)
         if (.not. (det(p) < dtmin .and. dele(p) < dlemin)) then
           fn = fn + 1
           filterp(fn) = p
         end if
       end do
     end if
   end do ITERATION     ! End stability iteration

   fn = fnorig
   filterp(1:fn) = fporig(1:fn)

   do f = 1, fn
     p = filterp(f)
     c = pcolumn(p)
     g = pgridcell(p)

     ! Energy balance check in canopy

     lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
           +(1._rkx-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
           +frac_h2osfc(c)*t_h2osfc(c)**4)

     err(p) = sabv(p) + air(p) + bir(p)*tlbef(p)**3*(tlbef(p) + &
             4._rkx*dt_veg(p)) + cir(p)*lw_grnd - eflx_sh_veg(p) - &
             hvap*qflx_evap_veg(p)
!         + cir(p)*t_grnd(c)**4 - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)

     ! Fluxes from ground to canopy space

     delt    = wtal(p)*t_grnd(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
     taux(p) = -forc_rho(g)*forc_u(g)/ram1(p)
     tauy(p) = -forc_rho(g)*forc_v(g)/ram1(p)
     eflx_sh_grnd(p) = cpair*forc_rho(g)*wtg(p)*delt
     ! compute individual sensible heat fluxes
     delt_snow = wtal(p)*t_soisno(c,snl(c)+1)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
     eflx_sh_snow(p) = cpair*forc_rho(g)*wtg(p)*delt_snow

     delt_soil  = wtal(p)*t_soisno(c,1)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
     eflx_sh_soil(p) = cpair*forc_rho(g)*wtg(p)*delt_soil

     delt_h2osfc  = wtal(p)*t_h2osfc(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
     eflx_sh_h2osfc(p) = cpair*forc_rho(g)*wtg(p)*delt_h2osfc
     qflx_evap_soi(p) = forc_rho(g)*wtgq(p)*delq(p)

     ! compute individual latent heat fluxes
     delq_snow = wtalq(p)*qg_snow(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(g)
     qflx_ev_snow(p) = forc_rho(g)*wtgq(p)*delq_snow

     delq_soil = wtalq(p)*qg_soil(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(g)
     qflx_ev_soil(p) = forc_rho(g)*wtgq(p)*delq_soil

     delq_h2osfc = wtalq(p)*qg_h2osfc(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(g)
     qflx_ev_h2osfc(p) = forc_rho(g)*wtgq(p)*delq_h2osfc

     ! 2 m height air temperature

     t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._rkx/temp12m(p) - 1._rkx/temp1(p))
     t_ref2m_r(p) = t_ref2m(p)

     ! 2 m height specific humidity

     q_ref2m(p) = forc_q(g) + temp2(p)*dqh(p)*(1._rkx/temp22m(p) - 1._rkx/temp2(p))

     ! 2 m height relative humidity

     call QSat(t_ref2m(p), forc_pbot(g), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
     rh_ref2m(p) = min(100._rkx, q_ref2m(p) / qsat_ref2m * 100._rkx)
     rh_ref2m_r(p) = rh_ref2m(p)

     ! Downward longwave radiation below the canopy

     dlrad(p) = (1._rkx-emv(p))*emg(c)*forc_lwrad(g) + &
         emv(p)*emg(c)*sb*tlbef(p)**3*(tlbef(p) + 4._rkx*dt_veg(p))

     ! Upward longwave radiation above the canopy

     ulrad(p) = ((1._rkx-emg(c))*(1._rkx-emv(p))*(1._rkx-emv(p))*forc_lwrad(g) &
         + emv(p)*(1._rkx+(1._rkx-emg(c))*(1._rkx-emv(p)))* &
         sb*tlbef(p)**3*(tlbef(p) + &
         4._rkx*dt_veg(p)) + emg(c)*(1._rkx-emv(p))*sb*lw_grnd)

     ! Derivative of soil energy flux with respect to soil temperature

     cgrnds(p) = cgrnds(p) + cpair*forc_rho(g)*wtg(p)*wtal(p)
     cgrndl(p) = cgrndl(p) + forc_rho(g)*wtgq(p)*wtalq(p)*dqgdT(c)
     cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

     ! Update dew accumulation (kg/m2)

     h2ocan(p) = max(0._rkx,h2ocan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtsrf)

     ! total photosynthesis

     fpsn(p)    = psnsun(p)   *laisun(p) + psnsha(p)   *laisha(p)
     fpsn_wc(p) = psnsun_wc(p)*laisun(p) + psnsha_wc(p)*laisha(p)
     fpsn_wj(p) = psnsun_wj(p)*laisun(p) + psnsha_wj(p)*laisha(p)
     fpsn_wp(p) = psnsun_wp(p)*laisun(p) + psnsha_wp(p)*laisha(p)

#if (defined CN)
     if ( use_c13 ) then
       rc13_canair(p) = c13o2(p)/(co2(p)-c13o2(p))
       rc13_psnsun(p) = rc13_canair(p)/alphapsnsun(p)
       rc13_psnsha(p) = rc13_canair(p)/alphapsnsha(p)
       c13_psnsun(p) = psnsun(p) * (rc13_psnsun(p)/(1._rkx+rc13_psnsun(p)))
       c13_psnsha(p) = psnsha(p) * (rc13_psnsha(p)/(1._rkx+rc13_psnsha(p)))

       ! use fixed c13 ratio with del13C of -25 to test the
       ! overall c13 structure
       ! c13_psnsun(p) = 0.01095627 * psnsun(p)
       ! c13_psnsha(p) = 0.01095627 * psnsha(p)
     endif

     if ( use_c14 ) then
       c14_psnsun(p) = rc14_atm(p) * psnsun(p)
       c14_psnsha(p) = rc14_atm(p) * psnsha(p)
     endif
#endif
   end do

   ! Filter out pfts which have small energy balance errors; report others

   fnold = fn
   fn = 0
   do f = 1, fnold
     p = filterp(f)
     if (abs(err(p)) > 0.1_rkx) then
       fn = fn + 1
       filterp(fn) = p
     end if
   end do

   do f = 1, fn
     p = filterp(f)
     write(stderr,*) 'energy balance in canopy ',p,', err=',err(p)
   end do

   contains
    !
    ! C13 fractionation during photosynthesis is calculated here
    ! after the nitrogen limitation is taken into account in the
    ! CNAllocation module.
    !
    subroutine Fractionation(lbp, ubp, fn, filterp, phase)
      implicit none
      integer(ik4), intent(in) :: lbp, ubp    ! pft bounds
      integer(ik4), intent(in) :: fn          ! size of pft filter
      integer(ik4), intent(in) :: filterp(fn) ! pft filter
      character(len=*), intent(in) :: phase   ! 'sun' or 'sha'

      integer(ik4) , pointer :: pgridcell(:)  ! pft's gridcell index
      integer(ik4) , pointer :: ivt(:)        ! pft vegetation type
      ! photosynthetic pathway: 0. = c4, 1. = c3
      real(rkx), pointer :: c3psn(:)
      ! number of canopy layers, above snow for radiative transfer
      integer(ik4) , pointer :: nrad(:)
      ! par absorbed per unit lai for canopy layer (w/m**2)
      real(rkx), pointer :: par_z(:,:)
      ! partial pressure co2 (Pa)
      real(rkx), pointer :: forc_pco2(:)
      ! fractional reduction in GPP due to N limitation (DIM)
      real(rkx), pointer :: downreg(:)
      real(rkx), pointer :: alphapsn(:)
      real(rkx), pointer :: forc_pbot(:) ! atmospheric pressure (Pa)
      ! net leaf photosynthesis (umol CO2/m**2/s)
      real(rkx), pointer :: an(:,:)
      ! leaf stomatal conductance (umol H2O/m**2/s)
      real(rkx), pointer :: gs_mol(:,:)
      ! leaf boundary layer conductance (umol H2O/m**2/s)
      real(rkx), pointer :: gb_mol(:)

      integer(ik4) :: f,p,g,iv  ! indices
      real(rkx) :: co2(lbp:ubp) ! atmospheric co2 partial pressure (pa)
      real(rkx) :: ci

      pgridcell      => clm3%g%l%c%p%gridcell
      nrad           => clm3%g%l%c%p%pps%nrad
      forc_pbot      => clm_a2l%forc_pbot
      forc_pco2      => clm_a2l%forc_pco2
      c3psn          => pftcon%c3psn
      ivt            => clm3%g%l%c%p%itype
      downreg        => clm3%g%l%c%p%pepv%downreg
      an             => clm3%g%l%c%p%ppsyns%an
      gb_mol         => clm3%g%l%c%p%ppsyns%gb_mol
      gs_mol         => clm3%g%l%c%p%ppsyns%gs_mol

      if (phase == 'sun') then
        par_z       => clm3%g%l%c%p%pef%parsun_z
        alphapsn    => clm3%g%l%c%p%pps%alphapsnsun
      else if (phase == 'sha') then
        par_z       => clm3%g%l%c%p%pef%parsha_z
        alphapsn    => clm3%g%l%c%p%pps%alphapsnsha
      end if

      do f = 1, fn
        p = filterp(f)
        g= pgridcell(p)
        co2(p) = forc_pco2(g)
        do iv = 1,nrad(p)
          if ( par_z(p,iv) <= 0._rkx ) then           ! night time
            alphapsn(p) = 1._rkx
          else                                      ! day time
            ci = co2(p) - ((an(p,iv) * (1._rkx-downreg(p)) )* &
                   forc_pbot(g) * (1.4_rkx*gs_mol(p,iv)+1.6_rkx*gb_mol(p)) / &
                   (gb_mol(p)*gs_mol(p,iv)))
            alphapsn(p) = 1._rkx + (((c3psn(ivt(p)) * (4.4_rkx + &
                 (22.6_rkx*(ci/co2(p))))) + ((1._rkx - c3psn(ivt(p))) * &
                 4.4_rkx))/1000._rkx)
          end if
        end do
      end do
    end subroutine Fractionation

  end subroutine CanopyFluxes
  !
  ! Leaf photosynthesis and stomatal conductance calculation as described by
  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
  ! a multi-layer canopy
  !
  subroutine Photosynthesis(fn, filterp, lbp, ubp, esat_tv, eair, oair, cair, &
                            rb, dayl_factor, phase)
    use mod_clm_varcon , only : rgas , tfrz
    implicit none
    integer(ik4) , intent(in)    :: fn            ! size of pft filter
    integer(ik4) , intent(in)    :: filterp(fn)   ! pft filter
    integer(ik4) , intent(in)    :: lbp, ubp      ! pft bounds
    ! saturation vapor pressure at t_veg (Pa)
    real(rkx), intent(in)    :: esat_tv(lbp:ubp)
    ! vapor pressure of canopy air (Pa)
    real(rkx), intent(in)    :: eair(lbp:ubp)
    ! Atmospheric O2 partial pressure (Pa)
    real(rkx), intent(in)    :: oair(lbp:ubp)
    ! Atmospheric CO2 partial pressure (Pa)
    real(rkx), intent(in)    :: cair(lbp:ubp)
    ! boundary layer resistance (s/m)
    real(rkx), intent(inout) :: rb(lbp:ubp)
    ! scalar (0-1) for daylength
    real(rkx), intent(in)    :: dayl_factor(lbp:ubp)
    character(len=*), intent(in) :: phase  ! 'sun' or 'sha'

    integer(ik4) , pointer :: pgridcell(:)! pft's gridcell index
    integer(ik4) , pointer :: ivt(:)      ! pft vegetation type
    real(rkx), pointer :: forc_pbot(:)! atmospheric pressure (Pa)
    real(rkx), pointer :: t_veg(:)    ! vegetation temperature (Kelvin)
    real(rkx), pointer :: btran(:)    ! transpiration wetness factor (0 to 1)
    ! air temperature at agcm reference height (kelvin)
    real(rkx), pointer :: tgcm(:)
    ! photosynthetic pathway: 0. = c4, 1. = c3
    real(rkx), pointer :: c3psn(:)
    ! specific leaf area at top of canopy, projected area basis [m^2/gC]
    real(rkx), pointer :: slatop(:)
    ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
    real(rkx), pointer :: flnr(:)
    real(rkx), pointer :: fnitr(:)  ! foliage nitrogen limitation factor (-)
    real(rkx), pointer :: leafcn(:) ! leaf C:N (gC/gN)
    ! number of canopy layers, above snow for radiative transfer
    integer(ik4) , pointer :: nrad(:)
    real(rkx), pointer :: tlai_z(:,:) ! total leaf area index for canopy layer
    ! leaf area index for canopy layer, sunlit or shaded
    real(rkx), pointer :: lai_z(:,:)
    ! par absorbed per unit lai for canopy layer (w/m**2)
    real(rkx), pointer :: par_z(:,:)
    real(rkx), pointer :: vcmaxcint(:) ! leaf to canopy scaling coefficient
    !KO
    ! 10-day running mean of the 2 m temperature (K)
    real(rkx), pointer :: t10(:)
    !KO

    !!! C13
    real(rkx), pointer :: alphapsn(:) ! 13C fractionation factor for PSN ()

    ! canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(rkx), pointer :: psn_z(:,:)
    ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rkx), pointer :: lmr_z(:,:)
    ! canopy layer: leaf stomatal resistance (s/m)
    real(rkx), pointer :: rs_z(:,:)
    real(rkx), pointer :: ci_z(:,:)   ! intracellular leaf CO2 (Pa)
    ! foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(rkx), pointer :: psn(:)
    ! Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(rkx), pointer :: psn_wc(:)
    ! RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(rkx), pointer :: psn_wj(:)
    ! product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(rkx), pointer :: psn_wp(:)
    ! leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rkx), pointer :: lmr(:)
    ! leaf stomatal resistance (s/m)
    real(rkx), pointer :: rs(:)
    !KO
    ! fractional humidity at leaf surface (dimensionless)
    real(rkx), pointer :: rh_leaf(:)
    !KO

    ! Leaf photosynthesis parameters

    ! maximum rate of carboxylation (umol co2/m**2/s)
    real(rkx),pointer :: vcmax_z(:,:)
    ! maximum electron transport rate (umol electrons/m**2/s)
    real(rkx) :: jmax_z(lbp:ubp,nlevcan)
    ! triose phosphate utilization rate (umol CO2/m**2/s)
    real(rkx),pointer :: tpu_z(:,:)
    ! initial slope of CO2 response curve (C4 plants)
    real(rkx),pointer :: kp_z(:,:)

    logical,pointer  :: c3flag(:) ! true if C3 and false if C4
    real(rkx) :: lnc(lbp:ubp)     ! leaf N concentration (gN leaf/m^2)
    real(rkx),pointer :: kc(:)    ! Michaelis-Menten constant for CO2 (Pa)
    real(rkx),pointer :: ko(:)    ! Michaelis-Menten constant for O2 (Pa)
    real(rkx),pointer :: cp(:)    ! CO2 compensation point (Pa)
    ! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
    real(rkx) :: bbbopt(lbp:ubp)
    ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(rkx),pointer :: bbb(:)
    ! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    real(rkx) :: mbbopt(lbp:ubp)
    ! Ball-Berry slope of conductance-photosynthesis relationship
    real(rkx),pointer :: mbb(:)
    ! leaf nitrogen decay coefficient
    real(rkx) :: kn(lbp:ubp)
    ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(rkx) :: vcmax25top
    ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(rkx) :: jmax25top
    ! canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(rkx) :: tpu25top
    ! canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(rkx) :: lmr25top
    ! canopy top: initial slope of CO2 response curve (C4 plants) at 25C
    real(rkx) :: kp25top

    ! leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(rkx) :: vcmax25
    ! leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(rkx) :: jmax25
    ! leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(rkx) :: tpu25
    ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(rkx) :: lmr25
    ! leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
    real(rkx) :: kp25
    real(rkx) :: kc25   ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(rkx) :: ko25   ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(rkx) :: cp25   ! CO2 compensation point at 25C (Pa)

    real(rkx) :: vcmaxha  ! activation energy for vcmax (J/mol)
    real(rkx) :: jmaxha   ! activation energy for jmax (J/mol)
    real(rkx) :: tpuha    ! activation energy for tpu (J/mol)
    real(rkx) :: lmrha    ! activation energy for lmr (J/mol)
    real(rkx) :: kcha     ! activation energy for kc (J/mol)
    real(rkx) :: koha     ! activation energy for ko (J/mol)
    real(rkx) :: cpha     ! activation energy for cp (J/mol)

    real(rkx) :: vcmaxhd  ! deactivation energy for vcmax (J/mol)
    real(rkx) :: jmaxhd   ! deactivation energy for jmax (J/mol)
    real(rkx) :: tpuhd    ! deactivation energy for tpu (J/mol)
    real(rkx) :: lmrhd    ! deactivation energy for lmr (J/mol)

    real(rkx) :: vcmaxse  ! entropy term for vcmax (J/mol/K)
    real(rkx) :: jmaxse   ! entropy term for jmax (J/mol/K)
    real(rkx) :: tpuse    ! entropy term for tpu (J/mol/K)
    real(rkx) :: lmrse    ! entropy term for lmr (J/mol/K)

    ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(rkx) :: vcmaxc
    ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(rkx) :: jmaxc
    ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(rkx) :: tpuc
    ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(rkx) :: lmrc

    ! quantum efficiency, used only for C4 (mol CO2 / mol photons)
    real(rkx),pointer :: qe(:)
    ! fraction of light absorbed by non-photosynthetic pigments
    real(rkx) :: fnps
    ! empirical curvature parameter for electron transport rate
    real(rkx) :: theta_psii

    ! empirical curvature parameter for ac, aj photosynthesis co-limitation
    real(rkx),pointer :: theta_cj(:)
    ! empirical curvature parameter for ap photosynthesis co-limitation
    real(rkx) :: theta_ip

    integer(ik4)  :: f,p,g,iv      ! indices
    real(rkx) :: cf       ! s m**2/umol -> s/m
    real(rkx) :: rsmax0   ! maximum stomatal resistance [s/m]
    real(rkx) :: gb       ! leaf boundary layer conductance (m/s)
    real(rkx) :: cs       ! CO2 partial pressure at leaf surface (Pa)
    real(rkx) :: gs       ! leaf stomatal conductance (m/s)
    real(rkx) :: hs       ! fractional humidity at leaf surface (dimensionless)
    real(rkx) :: sco      ! relative specificity of rubisco
    ! scaling factor for high temperature inhibition (25 C = 1.0)
    real(rkx) :: gs_mol_err  ! gs_mol for error check
    real(rkx) :: je          ! electron transport rate (umol electrons/m**2/s)
    real(rkx) :: qabs        ! PAR absorbed by PS II (umol photons/m**2/s)
    real(rkx) :: aquad,bquad,cquad ! terms for quadratic equations
    real(rkx) :: r1,r2         ! roots of quadratic equation
    real(rkx) :: ceair         ! vapor pressure of air, constrained (Pa)
    real(rkx) :: fnr           ! (gRubisco/gN in Rubisco)
    real(rkx) :: act25         ! (umol/mgRubisco/min) Rubisco activity at 25 C
    real(rkx) :: nscaler       ! leaf nitrogen scaling coefficient

    ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: ac(:,:)
    ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: aj(:,:)
    ! product-limited (C3) or CO2-limited (C4) gross
    ! photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: ap(:,:)
    ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: ag(:,:)
    ! net leaf photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: an(:,:)
    ! leaf stomatal conductance (umol H2O/m**2/s)
    real(rkx),pointer :: gs_mol(:,:)
    ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(rkx),pointer :: gb_mol(:)

    ! Rubisco-limited contribution to psn_z (umol CO2/m**2/s)
    real(rkx) :: psn_wc_z(lbp:ubp,nlevcan)
    ! RuBP-limited contribution to psn_z (umol CO2/m**2/s)
    real(rkx) :: psn_wj_z(lbp:ubp,nlevcan)
    ! product-limited contribution to psn_z (umol CO2/m**2/s)
    real(rkx) :: psn_wp_z(lbp:ubp,nlevcan)

    real(rkx) :: psncan       ! canopy sum of psn_z
    real(rkx) :: psncan_wc    ! canopy sum of psn_wc_z
    real(rkx) :: psncan_wj    ! canopy sum of psn_wj_z
    real(rkx) :: psncan_wp    ! canopy sum of psn_wp_z
    real(rkx) :: lmrcan       ! canopy sum of lmr_z
    real(rkx) :: gscan        ! canopy sum of leaf conductance
    real(rkx) :: laican       ! canopy sum of lai_z
    real(rkx) :: rh_can

    ! Assign local pointers to derived type members (gridcell-level)
    forc_pbot => clm_a2l%forc_pbot

    ! Assign local pointers to derived type members (pft-level)
    pgridcell => clm3%g%l%c%p%gridcell
    ivt       => clm3%g%l%c%p%itype
    t_veg     => clm3%g%l%c%p%pes%t_veg
    btran     => clm3%g%l%c%p%pps%btran
    tgcm      => clm3%g%l%c%p%pes%thm
    nrad      => clm3%g%l%c%p%pps%nrad
    tlai_z    => clm3%g%l%c%p%pps%tlai_z
    if (phase == 'sun') then
      lai_z  => clm3%g%l%c%p%pps%laisun_z
      par_z  => clm3%g%l%c%p%pef%parsun_z
      psn_z  => clm3%g%l%c%p%pcf%psnsun_z
      ci_z   => clm3%g%l%c%p%pcf%cisun_z
      lmr_z  => clm3%g%l%c%p%pcf%lmrsun_z
      rs_z   => clm3%g%l%c%p%pps%rssun_z
      psn    => clm3%g%l%c%p%pcf%psnsun
      psn_wc => clm3%g%l%c%p%pcf%psnsun_wc
      psn_wj => clm3%g%l%c%p%pcf%psnsun_wj
      psn_wp => clm3%g%l%c%p%pcf%psnsun_wp
      lmr    => clm3%g%l%c%p%pcf%lmrsun
      rs     => clm3%g%l%c%p%pps%rssun
      vcmaxcint => clm3%g%l%c%p%pps%vcmaxcintsun
      if ( use_c13 ) then
        alphapsn  => clm3%g%l%c%p%pps%alphapsnsun
      end if
    else if (phase == 'sha') then
      lai_z  => clm3%g%l%c%p%pps%laisha_z
      par_z  => clm3%g%l%c%p%pef%parsha_z
      psn_z  => clm3%g%l%c%p%pcf%psnsha_z
      ci_z   => clm3%g%l%c%p%pcf%cisha_z
      lmr_z  => clm3%g%l%c%p%pcf%lmrsha_z
      rs_z   => clm3%g%l%c%p%pps%rssha_z
      psn    => clm3%g%l%c%p%pcf%psnsha
      psn_wc => clm3%g%l%c%p%pcf%psnsha_wc
      psn_wj => clm3%g%l%c%p%pcf%psnsha_wj
      psn_wp => clm3%g%l%c%p%pcf%psnsha_wp
      lmr    => clm3%g%l%c%p%pcf%lmrsha
      rs     => clm3%g%l%c%p%pps%rssha
      vcmaxcint => clm3%g%l%c%p%pps%vcmaxcintsha
      if ( use_c13 ) then
        alphapsn  => clm3%g%l%c%p%pps%alphapsnsha
      end if
    end if
    !KO
    t10       => clm3%g%l%c%p%pes%t10
    rh_leaf   => clm3%g%l%c%p%pps%rh_leaf
    !KO

    ! Assign local pointers to pft constants
    c3psn     => pftcon%c3psn
    leafcn    => pftcon%leafcn
    flnr      => pftcon%flnr
    fnitr     => pftcon%fnitr
    slatop    => pftcon%slatop

    c3flag => clm3%g%l%c%p%ppsyns%c3flag
    ac     => clm3%g%l%c%p%ppsyns%ac
    aj     => clm3%g%l%c%p%ppsyns%aj
    ap     => clm3%g%l%c%p%ppsyns%ap
    ag     => clm3%g%l%c%p%ppsyns%ag
    an     => clm3%g%l%c%p%ppsyns%an
    gb_mol => clm3%g%l%c%p%ppsyns%gb_mol
    gs_mol => clm3%g%l%c%p%ppsyns%gs_mol
    vcmax_z=> clm3%g%l%c%p%ppsyns%vcmax_z
    cp     => clm3%g%l%c%p%ppsyns%cp
    kc     => clm3%g%l%c%p%ppsyns%kc
    ko     => clm3%g%l%c%p%ppsyns%ko
    qe     => clm3%g%l%c%p%ppsyns%qe
    tpu_z  => clm3%g%l%c%p%ppsyns%tpu_z
    kp_z   => clm3%g%l%c%p%ppsyns%kp_z
    theta_cj=>clm3%g%l%c%p%ppsyns%theta_cj
    bbb     =>clm3%g%l%c%p%ppsyns%bbb
    mbb    => clm3%g%l%c%p%ppsyns%mbb

    !==========================================================
    ! Photosynthesis and stomatal conductance parameters, from:
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
    !==========================================================

    ! vcmax25 parameters, from CN

    fnr = 7.16_rkx
    act25 = 3.6_rkx   !umol/mgRubisco/min
    ! Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
    act25 = act25 * 1000.0_rkx / 60.0_rkx

    ! Activation energy, from:
    ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
    ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
    ! except TPU from:
    !     Harley et al (1992) Plant, Cell and Environment 15:271-282

    kcha    = 79430._rkx
    koha    = 36380._rkx
    cpha    = 37830._rkx
    !KO   vcmaxha = 65330._rkx
    !KO   jmaxha  = 43540._rkx
    !KO   tpuha   = 53100._rkx
    !KO
    vcmaxha = 72000._rkx
    jmaxha  = 50000._rkx
    tpuha   = 72000._rkx
    !KO
    lmrha   = 46390._rkx

    ! High temperature deactivation, from:
    ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
    ! The factor "c" scales the deactivation to a value of 1.0 at 25C

    !KO   vcmaxhd = 149250._rkx
    !KO   jmaxhd  = 152040._rkx
    !KO   tpuhd   = 150650._rkx
    !KO
    vcmaxhd = 200000._rkx
    jmaxhd  = 200000._rkx
    tpuhd   = 200000._rkx
    !KO
    lmrhd   = 150650._rkx

    !KO   vcmaxse = 485._rkx
    !KO   jmaxse  = 495._rkx
    !KO   tpuse   = 490._rkx
    lmrse   = 490._rkx

    !KO   vcmaxc = fth25 (vcmaxhd, vcmaxse)
    !KO   jmaxc  = fth25 (jmaxhd, jmaxse)
    !KO   tpuc   = fth25 (tpuhd, tpuse)
    lmrc   = fth25 (lmrhd, lmrse)

    ! Miscellaneous parameters, from Bonan et al (2011)
    ! JGR, 116, doi:10.1029/2010JG001593

    fnps = 0.15_rkx
    theta_psii = 0.7_rkx
    theta_ip = 0.95_rkx

    do f = 1 , fn
      p = filterp(f)
      g = pgridcell(p)

      ! Modification for shrubs proposed by X.D.Z
      ! Why does he prefer this line here instead of in subr.
      ! CanopyFluxes? (slevis)
      ! Equivalent modification for soy following AgroIBIS
#if (defined CNDV)
      if (ivt(p) == nbrdlf_dcd_tmp_shrub) then
        btran(p) = min(1._rkx, btran(p) * 3.33_rkx)
      end if
#endif
      if (ivt(p) == nsoybean .or. ivt(p) == nsoybeanirrig) then
        btran(p) = min(1._rkx, btran(p) * 1.25_rkx)
      end if

      ! C3 or C4 photosynthesis logical variable

      if (nint(c3psn(ivt(p))) == 1) then
        c3flag(p) = .true.
      else if (nint(c3psn(ivt(p))) == 0) then
        c3flag(p) = .false.
      end if

      ! C3 and C4 dependent parameters

      if (c3flag(p)) then
        qe(p) = 0._rkx
        theta_cj(p) = 0.98_rkx
        if ( ivt(p) == nbrdlf_evr_trp_tree ) then
          bbbopt(p) = 80000._rkx
        else if ( ivt(p) == nbrdlf_dcd_trp_tree ) then
          bbbopt(p) = 1000._rkx
        else
          bbbopt(p) = 10000._rkx
        end if
        mbbopt(p) = 9._rkx
      else
        qe(p) = 0.05_rkx
        theta_cj(p) = 0.80_rkx
        bbbopt(p) = 40000._rkx
        mbbopt(p) = 4._rkx
      end if

      ! Soil water stress applied to Ball-Berry parameters

      bbb(p) = max (bbbopt(p)*btran(p), 1._rkx)
      mbb(p) = mbbopt(p)

      ! kc, ko, cp, from: Bernacchi et al (2001)
      ! Plant, Cell and Environment 24:253-259
      !
      !       kc25 = 404.9 umol/mol
      !       ko25 = 278.4 mmol/mol
      !       cp25 = 42.75 umol/mol
      !
      ! Derive sco from cp and O2 using present-day O2 (0.209 mol/mol)
      ! and re-calculate cp to account for variation in O2 using
      ! cp = 0.5 O2 / sco
      !

      kc25 = (404.9_rkx / 1.e6_rkx) * forc_pbot(g)
      ko25 = (278.4_rkx / 1.e3_rkx) * forc_pbot(g)
      sco  = 0.5_rkx * 0.209_rkx / (42.75_rkx / 1.e6_rkx)
      cp25 = 0.5_rkx * oair(p) / sco

      kc(p) = kc25 * ft(t_veg(p), kcha)
      ko(p) = ko25 * ft(t_veg(p), koha)
      cp(p) = cp25 * ft(t_veg(p), cpha)

    end do

    ! Multi-layer parameters scaled by leaf nitrogen profile.
    ! Loop through each canopy layer to calculate nitrogen profile using
    ! cumulative lai at the midpoint of the layer

    do f = 1, fn
      p = filterp(f)

      ! Leaf nitrogen concentration at the top of the canopy
      ! (g N leaf / m**2 leaf)

      lnc(p) = 1._rkx / (slatop(ivt(p)) * leafcn(ivt(p)))

      ! vcmax25 at canopy top, as in CN but using lnc at top of the canopy

      vcmax25top = lnc(p) * flnr(ivt(p)) * fnr * act25 * dayl_factor(p)
#ifndef CN
      vcmax25top = vcmax25top * fnitr(ivt(p))
#else
      if ( CNAllocation_Carbon_only() ) vcmax25top = vcmax25top * fnitr(ivt(p))
#endif

      ! Parameters derived from vcmax25top.
      ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      ! used jmax25 = 1.97 vcmax25, from
      ! Wullschleger (1993) Journal of Experimental Botany 44:907-920.

      !KO  jmax25top = 1.97_rkx * vcmax25top
      !KO
      jmax25top = (2.59_rkx - &
        0.035_rkx*min(max((t10(p)-tfrz),11._rkx),35._rkx)) * vcmax25top
      !KO
      tpu25top = 0.167_rkx * vcmax25top
      kp25top = 20000._rkx * vcmax25top

      ! Nitrogen scaling factor.
      ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
      ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al
      ! (2010) Biogeosciences, 7, 1833-1859
      ! Remove daylength factor from vcmax25 so that kn is based on
      ! maximum vcmax25
      ! But not used as defined here if using sun/shade big leaf code. Instead,
      ! will use canopy integrated scaling factors from SurfaceAlbedo.

      if (dayl_factor(p) == 0._rkx) then
        kn(p) =  0._rkx
      else
        kn(p) = exp(0.00963_rkx * vcmax25top/dayl_factor(p) - 2.43_rkx)
      end if

#if (defined CN)
      ! Leaf maintenance respiration to match the base rate used in CN
      ! but with the new temperature functions for C3 and C4 plants.
      !
      ! Base rate for maintenance respiration is from:
      ! M. Ryan, 1991. Effects of climate change on plant respiration.
      ! Ecological Applications, 1(2), 157-167.
      ! Original expression is br = 0.0106 molC/(molN h)
      ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
      !
      ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
      !
      ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
      ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s
      !
      ! Then scale this value at the top of the canopy for canopy depth

      lmr25top = 2.525e-6_rkx * (1.5_rkx ** ((25._rkx - 20._rkx)/10._rkx))
      lmr25top = lmr25top * lnc(p) / 12.e-6_rkx

#else
      ! Leaf maintenance respiration in proportion to vcmax25top

      if (c3flag(p)) then
        lmr25top = vcmax25top * 0.015_rkx
      else
        lmr25top = vcmax25top * 0.025_rkx
      end if

#endif

      ! Loop through canopy layers (above snow). Respiration needs to be
      ! calculated every timestep. Others are calculated only if daytime

      laican = 0._rkx
      do iv = 1, nrad(p)

        ! Cumulative lai at middle of layer

        if (iv == 1) then
          laican = 0.5_rkx * tlai_z(p,iv)
        else
          laican = laican + 0.5_rkx * (tlai_z(p,iv-1)+tlai_z(p,iv))
        end if

        ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
        ! profile. If sun/shade big leaf code, use canopy integrated factor.

        if (nlevcan == 1) then
          nscaler = vcmaxcint(p)
        else if (nlevcan > 1) then
          nscaler = exp(-kn(p) * laican)
        end if

        ! Maintenance respiration

        lmr25 = lmr25top * nscaler
        if (c3flag(p)) then
          lmr_z(p,iv) = lmr25 * ft(t_veg(p), lmrha) * &
                                fth(t_veg(p), lmrhd, lmrse, lmrc)
        else
          lmr_z(p,iv) = lmr25 * 2._rkx**((t_veg(p)-(tfrz+25._rkx))/10._rkx)
          lmr_z(p,iv) = lmr_z(p,iv) / &
            (1._rkx + exp( 1.3_rkx*(t_veg(p)-(tfrz+55._rkx)) ))
        end if

        if (par_z(p,iv) <= 0._rkx) then           ! night time

          vcmax_z(p,iv) = 0._rkx
          jmax_z(p,iv) = 0._rkx
          tpu_z(p,iv) = 0._rkx
          kp_z(p,iv) = 0._rkx

          if ( use_c13 ) then
            alphapsn(p) = 1._rkx
          end if

        else                                     ! day time

          vcmax25 = vcmax25top * nscaler
          jmax25 = jmax25top * nscaler
          tpu25 = tpu25top * nscaler
          kp25 = kp25top * nscaler

          ! Adjust for temperature
          !KO
          vcmaxse = 668.39_rkx - 1.07_rkx * min(max((t10(p)-tfrz),11._rkx),35._rkx)
          jmaxse  = 659.70_rkx - 0.75_rkx * min(max((t10(p)-tfrz),11._rkx),35._rkx)
          tpuse = vcmaxse
          vcmaxc = fth25 (vcmaxhd, vcmaxse)
          jmaxc  = fth25 (jmaxhd, jmaxse)
          tpuc   = fth25 (tpuhd, tpuse)
          !KO
          vcmax_z(p,iv) = vcmax25 * ft(t_veg(p), vcmaxha) * &
            fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
          jmax_z(p,iv) = jmax25 * ft(t_veg(p), jmaxha) * &
            fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
          tpu_z(p,iv) = tpu25 * ft(t_veg(p), tpuha) * &
            fth(t_veg(p), tpuhd, tpuse, tpuc)

          if (.not. c3flag(p)) then
            vcmax_z(p,iv) = vcmax25 * 2._rkx**((t_veg(p)-(tfrz+25._rkx))/10._rkx)
            vcmax_z(p,iv) = vcmax_z(p,iv) / &
              (1._rkx + exp( 0.2_rkx*((tfrz+15._rkx)-t_veg(p)) ))
            vcmax_z(p,iv) = vcmax_z(p,iv) / &
              (1._rkx + exp( 0.3_rkx*(t_veg(p)-(tfrz+40._rkx)) ))
          end if

          kp_z(p,iv) = kp25 * 2._rkx**((t_veg(p)-(tfrz+25._rkx))/10._rkx)

        end if

        ! Adjust for soil water

        vcmax_z(p,iv) = vcmax_z(p,iv) * btran(p)
        lmr_z(p,iv) = lmr_z(p,iv) * btran(p)

      end do       ! canopy layer loop
    end do          ! pft loop

    !====================================================
    ! Leaf-level photosynthesis and stomatal conductance
    !====================================================

    rsmax0 = 2.e4_rkx

    do f = 1, fn
      p = filterp(f)
      g = pgridcell(p)

      ! Leaf boundary layer conductance, umol/m**2/s

      cf = forc_pbot(g)/(rgas*1.e-3_rkx*tgcm(p))*1.e6_rkx
      gb = 1._rkx / rb(p)
      gb_mol(p) = gb * cf

      ! Loop through canopy layers (above snow).
      ! Only do calculations if daytime

      do iv = 1, nrad(p)

        if (par_z(p,iv) <= 0._rkx) then           ! night time

          ac(p,iv) = 0._rkx
          aj(p,iv) = 0._rkx
          ap(p,iv) = 0._rkx
          ag(p,iv) = 0._rkx
          an(p,iv) = ag(p,iv) - lmr_z(p,iv)
          psn_z(p,iv) = 0._rkx
          psn_wc_z(p,iv) = 0._rkx
          psn_wj_z(p,iv) = 0._rkx
          psn_wp_z(p,iv) = 0._rkx
          rs_z(p,iv) = min(rsmax0, 1._rkx/bbb(p) * cf)
          ci_z(p,iv) = 0._rkx
          !KO
          rh_leaf(p) = 0._rkx
          !KO

        else                                     ! day time

          !now the constraint is no longer needed, Jinyun Tang
          ceair = min( eair(p),  esat_tv(p) )
          rh_can = ceair / esat_tv(p)

          ! Electron transport rate for C3 plants. Convert par from W/m2 to
          ! umol photons/m**2/s using the factor 4.6

          qabs = 0.5_rkx * (1._rkx - fnps) * par_z(p,iv) * 4.6_rkx
          aquad = theta_psii
          bquad = -(qabs + jmax_z(p,iv))
          cquad = qabs * jmax_z(p,iv)
          call quadratic (aquad, bquad, cquad, r1, r2)
          je = min(r1,r2)

          ! Iterative loop for ci beginning with initial guess

          if (c3flag(p)) then
            ci_z(p,iv) = 0.7_rkx * cair(p)
          else
            ci_z(p,iv) = 0.4_rkx * cair(p)
          end if

          !find ci and stomatal conductance
          call hybrid(ci_z(p,iv), p, iv, g, gb_mol(p), je, cair(p), oair(p), &
                lmr_z(p,iv), par_z(p,iv), rh_can, gs_mol(p,iv))

          ! End of ci iteration.  Check for an < 0, in which case gs_mol = bbb

          if ( an(p,iv) < 0._rkx ) gs_mol(p,iv) = bbb(p)

          ! Final estimates for cs and ci
          ! (needed for early exit of ci iteration when an < 0)

          cs = cair(p) - 1.4_rkx/gb_mol(p) * an(p,iv) * forc_pbot(g)
          cs = max(cs, 1.e-6_rkx)
          ci_z(p,iv) = cair(p) - an(p,iv) * forc_pbot(g) * &
            (1.4_rkx*gs_mol(p,iv)+1.6_rkx*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv))

          ! Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)

          gs = gs_mol(p,iv) / cf
          rs_z(p,iv) = min(1._rkx/gs, rsmax0)

          ! Photosynthesis. Save rate-limiting photosynthesis

          psn_z(p,iv) = ag(p,iv)

          psn_wc_z(p,iv) = 0._rkx
          psn_wj_z(p,iv) = 0._rkx
          psn_wp_z(p,iv) = 0._rkx
          if (ac(p,iv) <= aj(p,iv) .and. ac(p,iv) <= ap(p,iv)) then
            psn_wc_z(p,iv) =  psn_z(p,iv)
          else if (aj(p,iv) < ac(p,iv) .and. aj(p,iv) <= ap(p,iv)) then
            psn_wj_z(p,iv) =  psn_z(p,iv)
          else if (ap(p,iv) < ac(p,iv) .and. ap(p,iv) < aj(p,iv)) then
            psn_wp_z(p,iv) =  psn_z(p,iv)
          end if

          ! Make sure iterative solution is correct

          if (gs_mol(p,iv) < 0._rkx) then
            write (stderr,*) 'Negative stomatal conductance:'
            write (stderr,*) gs_mol(p,iv)
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if

          ! Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b

          hs = (gb_mol(p)*ceair + gs_mol(p,iv)*esat_tv(p)) / &
            ((gb_mol(p)+gs_mol(p,iv))*esat_tv(p))
          !KO
          rh_leaf(p) = hs
          !KO
          gs_mol_err = mbb(p)*max(an(p,iv), 0._rkx)*hs/cs*forc_pbot(g) + bbb(p)

          if (abs(gs_mol(p,iv)-gs_mol_err) > 1.e-1_rkx) then
            write (stderr,*) 'Ball-Berry error check - stomatal conductance :'
            write (stderr,*) gs_mol(p,iv), gs_mol_err
          end if

        end if    ! night or day if
      end do       ! canopy layer loop
    end do          ! pft loop

    !================================================
    ! Canopy photosynthesis and stomatal conductance
    !================================================

    ! Sum canopy layer fluxes and then derive effective leaf-level fluxes (per
    ! unit leaf area), which are used in other parts of the model. Here, laican
    ! sums to either laisun or laisha.

    do f = 1, fn
      p = filterp(f)
      psncan = 0._rkx
      psncan_wc = 0._rkx
      psncan_wj = 0._rkx
      psncan_wp = 0._rkx
      lmrcan = 0._rkx
      gscan = 0._rkx
      laican = 0._rkx
      do iv = 1, nrad(p)
        psncan = psncan + psn_z(p,iv) * lai_z(p,iv)
        psncan_wc = psncan_wc + psn_wc_z(p,iv) * lai_z(p,iv)
        psncan_wj = psncan_wj + psn_wj_z(p,iv) * lai_z(p,iv)
        psncan_wp = psncan_wp + psn_wp_z(p,iv) * lai_z(p,iv)
        lmrcan = lmrcan + lmr_z(p,iv) * lai_z(p,iv)
        gscan = gscan + lai_z(p,iv) / (rb(p)+rs_z(p,iv))
        laican = laican + lai_z(p,iv)
      end do
      if (laican > 0._rkx) then
        psn(p) = psncan / laican
        psn_wc(p) = psncan_wc / laican
        psn_wj(p) = psncan_wj / laican
        psn_wp(p) = psncan_wp / laican
        lmr(p) = lmrcan / laican
        rs(p) = laican / gscan - rb(p)
      else
        psn(p) =  0._rkx
        psn_wc(p) =  0._rkx
        psn_wj(p) =  0._rkx
        psn_wp(p) =  0._rkx
        lmr(p) = 0._rkx
        rs(p) = 0._rkx
      end if
    end do

    contains
    !
    ! Temperature and soil water response functions
    !
    ! photosynthesis temperature response
    !
    pure real(rkx) function ft(tl,ha)
      use mod_clm_varcon , only : rgas , tfrz
      implicit none
      ! leaf temperature in photosynthesis temperature function (K)
      real(rkx), intent(in) :: tl
      ! activation energy in photosynthesis temperature function (J/mol)
      real(rkx), intent(in) :: ha
      ft = exp( ha / (rgas*1.e-3_rkx*(tfrz+25._rkx)) * (1._rkx - (tfrz+25._rkx)/tl) )
    end function ft
    !
    ! photosynthesis temperature inhibition
    !
    pure real(rkx) function fth(tl,hd,se,cc)
      use mod_clm_varcon , only : rgas
      implicit none
      ! leaf temperature in photosynthesis temperature function (K)
      real(rkx), intent(in) :: tl
      ! deactivation energy in photosynthesis temperature function (J/mol)
      real(rkx), intent(in) :: hd
      ! entropy term in photosynthesis temperature function (J/mol/K)
      real(rkx), intent(in) :: se
      ! scaling factor for high temperature inhibition (25 C = 1.0)
      real(rkx), intent(in) :: cc
      fth = cc / ( 1._rkx + exp( (-hd+se*tl) / (rgas*1.e-3_rkx*tl) ) )
    end function fth
    !
    ! scaling factor for photosynthesis temperature inhibition
    !
    pure real(rkx) function fth25(hd,se)
      use mod_clm_varcon , only : rgas , tfrz
      implicit none
      ! deactivation energy in photosynthesis temperature function (J/mol)
      real(rkx), intent(in) :: hd
      ! entropy term in photosynthesis temperature function (J/mol/K)
      real(rkx), intent(in) :: se
      fth25 = 1._rkx + exp( (-hd+se*(tfrz+25._rkx)) / (rgas*1.e-3_rkx*(tfrz+25._rkx)))
    end function fth25

  end subroutine Photosynthesis
  !
  ! Evaluate the function
  ! f(ci) = ci - (ca - (1.37rb+1.65rs))*patm*an
  !
  ! remark:  I am attempting to maintain the original code structure, also
  ! considering one may be interested to output relevant variables for the
  ! photosynthesis model, I have decided to add these relevant variables to
  ! the clmtype structure.
  !
  subroutine ci_func(ci,fval,p,iv,g,gb_mol,je,cair,oair,lmr_z,par_z, &
                     rh_can,gs_mol)
    implicit none
    real(rkx), intent(in) :: ci     ! intracellular leaf CO2 (Pa)
    ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rkx), intent(in) :: lmr_z
    ! par absorbed per unit lai for canopy layer (w/m**2)
    real(rkx), intent(in) :: par_z
    ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(rkx), intent(in) :: gb_mol
    ! electron transport rate (umol electrons/m**2/s)
    real(rkx), intent(in) :: je
    ! Atmospheric CO2 partial pressure (Pa)
    real(rkx), intent(in) :: cair
    ! Atmospheric O2 partial pressure (Pa)
    real(rkx), intent(in) :: oair
    ! canopy air realtive humidity
    real(rkx), intent(in) :: rh_can
    ! pft, vegetation type and column indexes
    integer(ik4),  intent(in) :: p, iv, g
    real(rkx), intent(out) :: fval   !return function of the value f(ci)
    ! leaf stomatal conductance (umol H2O/m**2/s)
    real(rkx), intent(out) :: gs_mol
    logical, pointer :: c3flag(:)    ! true if C3 and false if C4
    ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: ac(:,:)
    ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: aj(:,:)
    ! product-limited (C3) or CO2-limited (C4) gross
    ! photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: ap(:,:)
    ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: ag(:,:)
    ! net leaf photosynthesis (umol CO2/m**2/s)
    real(rkx),pointer :: an(:,:)
    ! maximum rate of carboxylation (umol co2/m**2/s)
    real(rkx),pointer :: vcmax_z(:,:)
    real(rkx),pointer :: cp(:)  ! CO2 compensation point (Pa)
    real(rkx),pointer :: kc(:)  ! Michaelis-Menten constant for CO2 (Pa)
    real(rkx),pointer :: ko(:)  ! Michaelis-Menten constant for O2 (Pa)
    ! quantum efficiency, used only for C4 (mol CO2 / mol photons)
    real(rkx),pointer :: qe(:)
    ! triose phosphate utilization rate (umol CO2/m**2/s)
    real(rkx),pointer :: tpu_z(:,:)
    ! initial slope of CO2 response curve (C4 plants)
    real(rkx),pointer :: kp_z(:,:)
    ! empirical curvature parameter for ac, aj photosynthesis co-limitation
    real(rkx),pointer :: theta_cj(:)
    real(rkx),pointer :: forc_pbot(:) ! atmospheric pressure (Pa)
    ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
    real(rkx),pointer :: bbb(:)
    ! Ball-Berry slope of conductance-photosynthesis relationship
    real(rkx),pointer :: mbb(:)

    real(rkx) :: ai ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
    real(rkx) :: cs ! CO2 partial pressure at leaf surface (Pa)

    real(rkx) :: aquad, bquad, cquad  ! terms for quadratic equations
    real(rkx) :: r1, r2  ! roots of quadratic equation
    ! fraction of light absorbed by non-photosynthetic pigments
    real(rkx) :: fnps
    ! empirical curvature parameter for electron transport rate
    real(rkx) :: theta_psii
    ! empirical curvature parameter for ap photosynthesis co-limitation
    real(rkx) :: theta_ip

    !grid level
    forc_pbot => clm_a2l%forc_pbot

    !pft level
    c3flag => clm3%g%l%c%p%ppsyns%c3flag
    ac     => clm3%g%l%c%p%ppsyns%ac
    aj     => clm3%g%l%c%p%ppsyns%aj
    ap     => clm3%g%l%c%p%ppsyns%ap
    ag     => clm3%g%l%c%p%ppsyns%ag
    an     => clm3%g%l%c%p%ppsyns%an
    vcmax_z=> clm3%g%l%c%p%ppsyns%vcmax_z
    cp     => clm3%g%l%c%p%ppsyns%cp
    kc     => clm3%g%l%c%p%ppsyns%kc
    ko     => clm3%g%l%c%p%ppsyns%ko
    qe     => clm3%g%l%c%p%ppsyns%qe
    tpu_z  => clm3%g%l%c%p%ppsyns%tpu_z
    kp_z   => clm3%g%l%c%p%ppsyns%kp_z
    theta_cj=>clm3%g%l%c%p%ppsyns%theta_cj
    bbb     =>clm3%g%l%c%p%ppsyns%bbb
    mbb    => clm3%g%l%c%p%ppsyns%mbb

    ! Miscellaneous parameters, from
    ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
    fnps = 0.15_rkx
    theta_psii = 0.7_rkx
    theta_ip = 0.95_rkx

    if ( c3flag(p) ) then

      ! C3: Rubisco-limited photosynthesis
      ac(p,iv) = vcmax_z(p,iv) * &
        max(ci-cp(p), 0._rkx) / (ci+kc(p)*(1._rkx+oair/ko(p)))

      ! C3: RuBP-limited photosynthesis
      aj(p,iv) = je * max(ci-cp(p), 0._rkx) / (4._rkx*ci+8._rkx*cp(p))

      ! C3: Product-limited photosynthesis
      ap(p,iv) = 3._rkx * tpu_z(p,iv)

    else

      ! C4: Rubisco-limited photosynthesis
      ac(p,iv) = vcmax_z(p,iv)

      ! C4: RuBP-limited photosynthesis
      aj(p,iv) = qe(p) * par_z * 4.6_rkx

      ! C4: PEP carboxylase-limited (CO2-limited)
      ap(p,iv) = kp_z(p,iv) * max(ci, 0._rkx) / forc_pbot(g)

    end if

    ! Gross photosynthesis. First co-limit ac and aj. Then co-limit ap

    aquad = theta_cj(p)
    bquad = -(ac(p,iv) + aj(p,iv))
    cquad = ac(p,iv) * aj(p,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ai = min(r1,r2)

    aquad = theta_ip
    bquad = -(ai + ap(p,iv))
    cquad = ai * ap(p,iv)
    call quadratic (aquad, bquad, cquad, r1, r2)
    ag(p,iv) = min(r1,r2)

    ! Net photosynthesis. Exit iteration if an < 0

    an(p,iv) = ag(p,iv) - lmr_z
    if (an(p,iv) < 0._rkx) then
      fval = 0._rkx
      return
    endif
    ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
    ! With an <= 0, then gs_mol = bbb

    cs = cair - 1.4_rkx/gb_mol * an(p,iv) * forc_pbot(g)
    cs = max(cs,1.e-6_rkx)
    aquad = cs
    bquad = cs*(gb_mol - bbb(p)) - mbb(p)*an(p,iv)*forc_pbot(g)
    cquad = -gb_mol*(cs*bbb(p) + mbb(p)*an(p,iv)*forc_pbot(g)*rh_can)
    call quadratic (aquad, bquad, cquad, r1, r2)
    gs_mol = max(r1,r2)

    ! Derive new estimate for ci

    fval = ci - cair + an(p,iv) * forc_pbot(g) * &
      (1.4_rkx*gs_mol+1.6_rkx*gb_mol) / (gb_mol*gs_mol)

  end subroutine ci_func
  !
  !=========================================================================!
  !-------------- Solve quadratic equation for its two roots ---------------!
  !=========================================================================!
  ! Solution from Press et al (1986) Numerical Recipes: The Art of Scientific
  ! Computing (Cambridge University Press, Cambridge), pp. 145.
  !
  subroutine quadratic (a, b, c, r1, r2)
    implicit none
    real(rkx), intent(in)  :: a,b,c  ! Terms for quadratic equation
    real(rkx), intent(out) :: r1,r2  ! Roots of quadratic equation
    real(rkx) :: q                   ! Temporary term for quadratic solution

    if ( abs(a) < 1.e-20_rkx ) then
      write (stderr,*) 'Quadratic solution error: a = ',a
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    if ( b >= 0._rkx ) then
      q = -0.5_rkx * (b + sqrt(b*b - 4._rkx*a*c))
    else
      q = -0.5_rkx * (b - sqrt(b*b - 4._rkx*a*c))
    end if

    r1 = q / a
    if ( q /= 0._rkx ) then
      r2 = c / q
    else
      r2 = 1.e36_rkx
    end if
  end subroutine quadratic
  !
  ! use a hybrid solver to find the root of equation
  ! f(x) = x - h(x),
  ! we want to find x, s.t. f(x) = 0.
  ! the hybrid approach combines the strength of the newton secant
  ! approach (find the solution domain) and the bisection approach
  ! implemented with the Brent's method to guarrantee convergence.
  !
  subroutine hybrid(x0, p, iv, g, gb_mol, je, cair, oair, lmr_z, par_z,&
                    rh_can, gs_mol)
    implicit none
    !initial guess and final value of the solution
    real(rkx), intent(inout) :: x0
    ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(rkx), intent(in) :: lmr_z
    ! par absorbed per unit lai for canopy layer (w/m**2)
    real(rkx), intent(in) :: par_z
    ! canopy air relative humidity
    real(rkx), intent(in) :: rh_can
    ! leaf boundary layer conductance (umol H2O/m**2/s)
    real(rkx), intent(in) :: gb_mol
    ! electron transport rate (umol electrons/m**2/s)
    real(rkx), intent(in) :: je
    ! Atmospheric CO2 partial pressure (Pa)
    real(rkx), intent(in) :: cair
    ! Atmospheric O2 partial pressure (Pa)
    real(rkx), intent(in) :: oair
    ! pft, c3/c4, and column index
    integer(ik4),  intent(in) :: p, iv, g
    ! leaf stomatal conductance (umol H2O/m**2/s)
    real(rkx), intent(out) :: gs_mol
    !number of iterations used, for record only
    integer(ik4) :: iter

    real(rkx) :: x1, f0, f1
    real(rkx) :: x, dx
    real(rkx) , parameter :: eps = 1.e-2_rkx !relative accuracy
    real(rkx) , parameter :: eps1= 1.e-4_rkx
    integer(ik4),  parameter :: itmax = 40 !maximum number of iterations
    real(rkx) :: tol,minx,minf

    call ci_func(x0, f0, p, iv, g, gb_mol, je, cair, &
                 oair, lmr_z, par_z, rh_can, gs_mol)
    if ( abs(f0) < 1.e-20_rkx ) return

    minx = x0
    minf = f0
    x1 = x0 * 0.99_rkx
    call ci_func(x1, f1, p, iv, g, gb_mol, je, cair, &
                 oair, lmr_z, par_z, rh_can, gs_mol)

    if ( abs(f1) < 1.e-20_rkx ) then
      x0 = x1
      return
    end if
    if ( f1 < minf ) then
      minx = x1
      minf = f1
    end if

    !first use the secant approach, then use the brent approach as a backup
    iter = 0
    do
      iter = iter + 1
      dx = - f1 * (x1-x0)/(f1-f0)
      x = x1 + dx
      tol = abs(x) * eps
      if ( abs(dx) < tol ) then
        x0 = x
        exit
      end if
      x0 = x1
      f0 = f1
      x1 = x
      call ci_func(x1, f1, p, iv, g, gb_mol, je, cair, oair, &
                   lmr_z, par_z, rh_can, gs_mol)
      if ( f1 < minf ) then
        minx = x1
        minf = f1
      end if
      if ( abs(f1) <= eps1 ) then
        x0 = x1
        exit
      end if
      !
      ! if a root zone is found, use the brent method for a
      ! robust backup strategy
      !
      if ( f1 * f0 < 0._rkx ) then
        call brent(x, x0, x1, f0, f1, tol, p, iv, g, gb_mol, je, cair, oair, &
                   lmr_z, par_z, rh_can, gs_mol)
        x0 = x
        exit
      end if
      if ( iter > itmax ) then
        ! in case of failing to converge within itmax iterations
        ! stop at the minimum function
        ! this happens because of some other issues besides the
        ! stomatal conductance calculation and it happens usually
        ! in very dry places and more likely with c4 plants.
        call ci_func(minx, f1, p, iv, g, gb_mol, je, cair, oair, &
                     lmr_z, par_z, rh_can, gs_mol)
        exit
      end if
    end do

    contains
    !
    ! Use Brent's method to find the root of a single variable function
    ! ci_func, which is known to exist between x1 and x2.
    ! The found root will be updated until its accuracy is tol.
    !
    subroutine brent(x,x1,x2,f1,f2,tol,ip,iv,ig,gb_mol,je,cair,oair, &
                     lmr_z,par_z,rh_can,gs_mol)
      implicit none
      ! indepedent variable of the single value function ci_func(x)
      real(rkx), intent(out) :: x
      ! minimum and maximum of the variable domain to search for the
      ! solution ci_func(x1) = f1, ci_func(x2)=f2
      real(rkx), intent(in) :: x1, x2, f1, f2
      real(rkx), intent(in) :: tol    !the error tolerance

      ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
      real(rkx), intent(in) :: lmr_z
      ! par absorbed per unit lai for canopy layer (w/m**2)
      real(rkx), intent(in) :: par_z
      ! leaf boundary layer conductance (umol H2O/m**2/s)
      real(rkx), intent(in) :: gb_mol
      ! electron transport rate (umol electrons/m**2/s)
      real(rkx), intent(in) :: je
      real(rkx), intent(in) :: cair   ! Atmospheric CO2 partial pressure (Pa)
      real(rkx), intent(in) :: oair   ! Atmospheric O2 partial pressure (Pa)
      real(rkx), intent(in) :: rh_can ! inside canopy relative humidity
      ! pft, c3/c4, and column index
      integer(ik4),  intent(in) :: ip, iv, ig
      ! leaf stomatal conductance (umol H2O/m**2/s)
      real(rkx), intent(out) :: gs_mol

      integer(ik4), parameter :: itmax = 20  !maximum number of iterations
      real(rkx), parameter :: eps = 1.e-2_rkx    !relative error tolerance

      integer(ik4) :: iter
      real(rkx) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

      a = x1
      b = x2
      fa = f1
      fb = f2
      if ( (fa > 0._rkx .and. fb > 0._rkx) .or. &
           (fa < 0._rkx .and. fb < 0._rkx) ) then
        write(stderr,*) 'root must be bracketed for brent'
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
      c = b
      fc = fb
      iter = 0
      do
        if ( iter == itmax ) exit
        iter = iter + 1
        if ( (fb > 0._rkx .and. fc > 0._rkx) .or. &
             (fb < 0._rkx .and. fc < 0._rkx) ) then
          c = a   !Rename a, b, c and adjust bounding interval d.
          fc = fa
          d = b-a
          e = d
        end if
        if ( abs(fc) < abs(fb) ) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
        end if
        tol1 = 2._rkx*eps*abs(b)+0.5_rkx*tol  !Convergence check.
        xm = 0.5_rkx*(c-b)
        if ( abs(xm) <= tol1 .or. fb == 0.0_rkx ) then
          x = b
          return
        end if
        if ( abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then
          s = fb/fa !Attempt inverse quadratic interpolation.
          if ( a == c ) then
            p = 2._rkx*xm*s
            q = 1._rkx-s
          else
            q = fa/fc
            r = fb/fc
            p = s*(2._rkx*xm*q*(q-r)-(b-a)*(r-1._rkx))
            q = (q-1._rkx)*(r-1._rkx)*(s-1._rkx)
          end if
          if ( p > 0._rkx ) q = -q !Check whether in bounds.
          p = abs(p)
          if ( 2._rkx*p < min(3._rkx*xm*q-abs(tol1*q),abs(e*q)) ) then
            e = d !Accept interpolation.
            d = p/q
          else
            d = xm  !Interpolation failed, use bisection.
            e = d
          end if
        else !Bounds decreasing too slowly, use bisection.
          d = xm
          e = d
        end if
        a = b !Move last best guess to a.
        fa = fb
        if ( abs(d) > tol1 ) then !Evaluate new trial root.
          b = b+d
        else
          b = b+sign(tol1,xm)
        end if
        call ci_func(b, fb, ip, iv, ig, gb_mol, je, cair, &
                     oair, lmr_z, par_z, rh_can, gs_mol)
        if ( fb == 0._rkx ) exit
      end do
      if ( iter == itmax ) then
        write(stderr,*) 'brent exceeding maximum iterations', b, fb
      end if
      x = b
    end subroutine brent

  end subroutine hybrid

end module mod_clm_canopyfluxes
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
