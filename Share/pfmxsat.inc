  ! Computes saturation mixing ratio over a mixed-phase condensate
  ! and the effective latent heat for use in the Clausius-Clapeyron
  ! denominator of the temperature adjustment iteration.
  !
  ! Phase blending follows the quadratic interpolation of Chaboureau
  ! and Bechtold (2002, JAS) and is consistent with the formulation
  ! used in Bretherton et al. (2004, MWR) for the UW TCM:
  !
  !   alpha(T) = [(T - T_ice)/(T_liq - T_ice)]^2   for T_ice <= T <= T_liq
  !   alpha(T) = 1                                   for T >= T_liq   (pure liquid)
  !   alpha(T) = 0                                   for T <= T_ice   (pure ice)
  !
  ! where T_liq = 273.15 K (0°C) and T_ice = 233.15 K (-40°C).
  !
  !   es_eff  = alpha * es_liq  + (1-alpha) * es_ice   [Pa]
  !   L_eff   = alpha * Lv      + (1-alpha) * Ls        [J/kg]
  !   ws      = ep2 * es_eff / (p - es_eff)            [kg/kg]
  !
  ! Saturation pressures use the Flatau et al. (1992, JAM) 8th-order
  ! polynomial fits, accurate to < 0.03% from -75°C to +100°C.
  !
  ! ARGUMENTS
  !   t   - temperature [K]
  !   p   - total pressure [Pa]
  !
  ! RETURNS (as a subroutine so both ws and L_eff can be returned)
  !   ws    - saturation mixing ratio [kg/kg]
  !   Leff  - effective latent heat   [J/kg]
  !
  ! USAGE IN THE CLAUSIUS-CLAPEYRON ITERATION:
  !
  !   call pfmxsat(temps, preshl(k), rvls, Leff)
  !   cpoxleff = cp(k) / Leff
  !   deltat = ((templ - temps)*cpoxleff + qwx(k) - rvls) /
  !            (cpoxleff + ep2*Leff*rvls/(rgas*temps*temps))
  !
  ! REFERENCES
  !   Flatau, P.J., Walko, R.L. and Cotton, W.R. (1992), Polynomial
  !     Fits to Saturation Vapor Pressure. J. Appl. Meteor., 31, 1507-1513.
  !   Chaboureau, J-P. and Bechtold, P. (2002), A Simple Cloud Parameterization
  !     Derived from Cloud Resolving Model Data. J. Atmos. Sci., 59, 2362-2372.
  !   Bretherton, C.S., McCaa, J.R. and Grenier, H. (2004), A New
  !     Parameterization for Shallow Cumulus Convection. Mon. Wea. Rev., 132.
  !
!DIR$ ATTRIBUTES FORCEINLINE :: pfmxsat
  pure elemental subroutine pfmxsat(t, p, ws, Leff)
!$acc routine seq
    use mod_constants, only : ep2, tzero
    implicit none
    real(rkx), intent(in)  :: t    ! Temperature [K]
    real(rkx), intent(in)  :: p    ! Pressure [Pa]
    real(rkx), intent(out) :: ws   ! Saturation mixing ratio [kg/kg]
    real(rkx), intent(out) :: Leff ! Effective latent heat [J/kg]

    ! Flatau et al. (1992) liquid-water coefficients (0°C to 100°C)
    real(rk8), parameter :: a0 =  0.611213476e+03_rk8
    real(rk8), parameter :: a1 =  0.444007856e+02_rk8
    real(rk8), parameter :: a2 =  0.143064234e+01_rk8
    real(rk8), parameter :: a3 =  0.264461437e-01_rk8
    real(rk8), parameter :: a4 =  0.305903558e-03_rk8
    real(rk8), parameter :: a5 =  0.196237241e-05_rk8
    real(rk8), parameter :: a6 =  0.892344772e-08_rk8
    real(rk8), parameter :: a7 = -0.373208410e-10_rk8
    real(rk8), parameter :: a8 =  0.209339997e-13_rk8

    ! Flatau et al. (1992) ice coefficients (-75°C to 0°C)
    real(rk8), parameter :: c0 =  0.611123516e+03_rk8
    real(rk8), parameter :: c1 =  0.503109514e+02_rk8
    real(rk8), parameter :: c2 =  0.188369801e+01_rk8
    real(rk8), parameter :: c3 =  0.420547422e-01_rk8
    real(rk8), parameter :: c4 =  0.614396778e-03_rk8
    real(rk8), parameter :: c5 =  0.602780717e-05_rk8
    real(rk8), parameter :: c6 =  0.387940929e-07_rk8
    real(rk8), parameter :: c7 =  0.149436277e-09_rk8
    real(rk8), parameter :: c8 =  0.262655803e-12_rk8

    ! Phase-blending temperature limits
    real(rkx), parameter :: T_liq  = 273.15_rkx   ! 0°C: pure liquid above
    real(rkx), parameter :: T_ice  = 233.15_rkx   ! -40°C: pure ice below
    real(rkx), parameter :: dTmx   = T_liq - T_ice ! = 40 K
    real(rkx), parameter :: rdTmx  = d_one / dTmx

    ! Reference latent heats at 0°C [J/kg]
    real(rkx), parameter :: Lv0 = 2.501e6_rkx   ! vaporisation
    real(rkx), parameter :: Ls0 = 2.835e6_rkx   ! sublimation

    real(rk8)  :: td, es_liq, es_ice, es_eff
    real(rkx)  :: alpha, es_rkx

    td = real(t,rk8) - real(tzero,rk8)
    td = min(max(td, -75.0_rk8), 100.0_rk8)

    ! Saturation pressure over liquid water [Pa]
    es_liq = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
           + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    es_liq = min(es_liq, 0.15_rk8*real(p,rk8))

    ! Saturation pressure over ice [Pa]
    es_ice = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
           + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    es_ice = min(es_ice, 0.15_rk8*real(p,rk8))

    ! Quadratic liquid fraction: 1 at T_liq, 0 at T_ice
    if ( t >= T_liq ) then
      alpha = d_one
    else if ( t <= T_ice ) then
      alpha = d_zero
    else
      alpha = ((t - T_ice) * rdTmx) ** 2
    end if

    ! Blended effective saturation pressure [Pa]
    es_eff = alpha * es_liq + (d_one - real(alpha,rk8)) * es_ice

    ! Exact saturation mixing ratio [kg/kg]
    es_rkx = real(es_eff, rkx)
    ws = ep2 * es_rkx / (p - es_rkx)

    ! Effective latent heat [J/kg], blended consistently with es_eff
    ! so that d(ws)/dT = ep2*p*Leff*ws / (Rd*(p-es)*T^2) holds.
    Leff = alpha * Lv0 + (d_one - alpha) * Ls0

  end subroutine pfmxsat

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
