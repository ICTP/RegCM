!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_pbl_myj_sflayer
  !
  ! Janjic Eta surface layer for use with the MYJ PBL scheme.
  !
  ! The scheme computes Monin-Obukhov exchange coefficients (akms, akhs)
  ! and updates the viscous-sublayer reference values (uz0, vz0, thz0, qz0)
  ! used as lower boundary conditions by the MYJ diffusion routines.
  !
  ! Surface stress and heat/moisture fluxes (uvdrag, hfx, qfx) from RegCM's
  ! existing land/ocean surface scheme provide the anchor for ustar and the
  ! Obukhov length, so the MYJ surface layer is consistent with — rather than
  ! replacing — the BATS/CLM/COARE computation.  Exchange coefficients akms
  ! and akhs are then derived self-consistently with that ustar.
  !
  ! REFERENCES:
  !   Janjic (1994), Mon. Wea. Rev., 122, 927-945.
  !   Janjic (2002), NCEP Office Note 437.
  !   Mellor & Yamada (1982), Rev. Geophys. Space Phys., 20, 851-875.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_regcm_types
  use mod_runparams, only : rcmtimer, iqv

  implicit none

  private

  public :: myj_sfclayer

  !--------------------------------------------------------------------
  ! Constants shared with the Janjic viscous-sublayer parameterisation.
  ! Numerical values are identical to those in mod_pbl_myj.F90 so that
  ! the two modules are consistent.
  !--------------------------------------------------------------------

  ! Minimum / maximum Obukhov length to avoid singularities
  ! THIS CAN BE TUNED.
  real(rkx), parameter :: obmin =  1.0_rkx    ! [m]  minimum |L| (stable)
  real(rkx), parameter :: obmax = -1.0_rkx    ! [m]  maximum L  (unstable)
  real(rkx), parameter :: minak =  0.1_rkx    ! 0.01 if over diffusive

  ! Von Karman constant (must match mod_constants value)
  ! vonkar is imported from mod_constants via mod_dynparam.

  ! Neutral drag-coefficient lower bound
  real(rkx), parameter :: cdmin = 1.0e-6_rkx  ! [m/s]  lower bound on akms

  ! Turbulent Prandtl number at neutral stability (used as Ch/Cm ratio)
  real(rkx), parameter :: prt0 = 1.0_rkx

  ! Coefficients for Businger-Dyer-Pandolfo stability functions
  ! (Janjic 1994, eqs. 22-27; same as WRF sf_myjsfc).
  real(rkx), parameter :: alpha4 = 4.7_rkx     ! stable   psi_m/h coefficient
  real(rkx), parameter :: alpha5 = 16.0_rkx    ! unstable psi_m coefficient
  real(rkx), parameter :: alpha6 = 16.0_rkx    ! unstable psi_h coefficient

  ! Viscous sublayer constants (identical to mod_pbl_myj.F90)
  real(rkx), parameter :: visc  = 1.5e-5_rkx  ! kinematic viscosity [m^2/s]
  real(rkx), parameter :: tvisc = 2.1e-5_rkx  ! thermal diffusivity [m^2/s]
  real(rkx), parameter :: qvisc = 2.1e-5_rkx  ! moisture diffusivity [m^2/s]
  real(rkx), parameter :: glkbs = 30.0_rkx ! Kader-Yaglom constant (scalar)
  real(rkx), parameter :: glkbr = 10.0_rkx ! Kader-Yaglom constant (momentum)
  real(rkx), parameter :: sqpr  = 0.84_rkx ! sqrt(Pr)
  real(rkx), parameter :: sqsc  = 0.84_rkx ! sqrt(Sc)
  real(rkx), parameter :: small = 0.35_rkx ! = cziv / glkbs
  real(rkx), parameter :: ustc  = 0.7_rkx  ! viscous→wave transition [m/s]
  real(rkx), parameter :: ustr  = 0.225_rkx   ! smooth→transitional   [m/s]
  real(rkx), parameter :: ustfc = 0.018_rkx/egrav  ! Charnock coefficient
  real(rkx), parameter :: seafc = 0.98_rkx    ! sea-surface activity

  ! Derived viscous-sublayer constants
  real(rkx), parameter :: rvisc  = d_one/visc
  real(rkx), parameter :: rtvisc = d_one/tvisc
  real(rkx), parameter :: rqvisc = d_one/qvisc
  real(rkx), parameter :: cziv   = small*glkbs
  real(rkx), parameter :: grrs   = glkbr/glkbs
  real(rkx), parameter :: zqrzt  = sqsc/sqpr
  real(rkx), parameter :: fzu1   = cziv*visc
  real(rkx), parameter :: fzt1   = rvisc*tvisc*sqpr
  real(rkx), parameter :: fzt2   = cziv*grrs*tvisc*sqpr
  real(rkx), parameter :: fzq1   = rtvisc*qvisc*zqrzt
  real(rkx), parameter :: fzq2   = fzq1

  contains

  ! ===========================================================================
  subroutine myj_sfclayer(m2p, ustar2d, akms2d, akhs2d)
  ! ===========================================================================
  !
  ! PURPOSE
  !   Compute the MYJ surface-layer exchange coefficients akms and akhs
  !   [m s^-1], update the viscous-sublayer reference values uz0, vz0,
  !   thz0, qz0, and return the friction velocity ustar.
  !
  ! INPUTS  (via m2p struct)
  !   uvdrag  - surface drag [kg m^-2 s^-1]  = rho * Cd * |U|
  !   uxatm   - lowest model-level u [m/s]
  !   vxatm   - lowest model-level v [m/s]
  !   tatm    - lowest model-level temperature [K]
  !   patm    - lowest model-level pressure [Pa]
  !   patmf   - interface pressures [Pa]  (patmf(:,:,kzp1) = p_sfc)
  !   qxatm   - moisture mixing ratios [kg/kg]
  !   hfx     - upward sensible heat flux [W m^-2]  (positive upward)
  !   qfx     - upward moisture flux [kg m^-2 s^-1] (positive upward)
  !   tg      - surface (skin) temperature [K]
  !   q2m     - near-surface specific humidity [kg/kg] (land)
  !   zq      - half-level heights AGL + orography [m]
  !   ht      - surface orography height [m]
  !   rhox2d  - surface air density [kg m^-3]
  !   ldmsk   - land mask (0 = ocean, 1 = land)
  !   uz0,vz0,thz0,qz0  - viscous-sublayer values (updated in place)
  !
  ! OUTPUTS
  !   ustar2d - friction velocity [m/s]
  !   akms2d  - momentum exchange coefficient [m/s]
  !   akhs2d  - heat/moisture exchange coefficient [m/s]
  !   m2p%uz0, vz0, thz0, qz0 updated in place.
  !
  ! ALGORITHM
  !   1.  Compute ustar from the surface stress already calculated by the
  !       land/ocean scheme:
  !           tau = uvdrag * |U_kz|     (uvdrag = rho * Cd * |U|)
  !           ustar = sqrt( tau / rho ) = sqrt( uvdrag * |U| / rho )
  !
  !   2.  Compute the Obukhov length L from ustar and hfx:
  !           L = -rho * cp * ustar^3 * theta_v /
  !               (vonkar * g * hfx_tv)
  !       where hfx_tv = hfx + rho*cp*ep1*theta_v*qfx/rho is the virtual
  !       heat flux.
  !
  !   3.  Compute the stability parameter zeta = za / L, where za is the
  !       height of the lowest model level AGL.
  !
  !   4.  Evaluate Businger-Dyer-Pandolfo integrated stability functions
  !       psi_m(zeta) and psi_h(zeta).
  !
  !   5.  Derive akms and akhs:
  !           akms = vonkar * ustar / (log(za/z0) - psi_m)
  !           akhs = vonkar * ustar / (log(za/zt) - psi_h)
  !       where z0 = max(ustfc*ustar^2, z0_min) (Charnock + smooth floor)
  !       and zt ≈ z0 (Janjic assumes scalar and momentum roughness equal
  !       at the resolved-scale interface; viscous sublayer handles the
  !       actual scalar-momentum separation).
  !
  !   6.  Update viscous-sublayer values uz0, vz0, thz0, qz0 following
  !       the three-regime Janjic parameterisation (smooth, transitional,
  !       rough/wavy) exactly as in the original mod_pbl_myj.F90 code,
  !       now driven by the self-consistent ustar from step 1.
  !
  implicit none

  type(mod_2_pbl), intent(inout) :: m2p
  real(rkx), dimension(jci1:jci2,ici1:ici2), intent(out) :: ustar2d
  real(rkx), dimension(jci1:jci2,ici1:ici2), intent(out) :: akms2d
  real(rkx), dimension(jci1:jci2,ici1:ici2), intent(out) :: akhs2d

  integer(ik4) :: i, j
  real(rkx) :: uatm, vatm, wspd, tau_over_rho, ustar, ust2
  real(rkx) :: hfxv, thv, theta_v, lob, zeta, za
  real(rkx) :: z0, logza_z0, psi_m, psi_h, x, x2, x4
  real(rkx) :: akms, akhs
  real(rkx) :: psfc, rexnsfc, tg, thsk, qsfc
  real(rkx) :: tha, qha, ratiomx
  real(rkx) :: zu, zt, zq, wght, wghtt, wghtq
  real(rkx) :: exner_kz, elocp_loc

  ! elocp = eliwv/cpd  (cannot use the private parameter from mod_pbl_myj)
  elocp_loc = eliwv / cpd

  do concurrent ( j = jci1:jci2, i = ici1:ici2 )

    !----------------------------------------------------------------
    ! Step 1.  Friction velocity from surface stress
    !
    ! uvdrag [kg m^-2 s^-1] = rho * Cd * |U|
    ! tau    [Pa]            = rho * Cd * |U|^2  = uvdrag * |U|
    ! ustar  [m/s]           = sqrt(tau / rho)   = sqrt(uvdrag * |U| / rho)
    !----------------------------------------------------------------
    uatm = m2p%uxatm(j,i,kz)
    vatm = m2p%vxatm(j,i,kz)
    wspd = max(sqrt(uatm*uatm + vatm*vatm), 0.01_rkx)

    ! tau/rho = uvdrag * wspd / rho
    tau_over_rho = max(m2p%uvdrag(j,i) * wspd / m2p%rhox2d(j,i), 1.0e-8_rkx)
    ustar  = sqrt(tau_over_rho)   ! [m/s]
    ust2   = ustar * ustar         ! [m^2/s^2] = Cd * wspd^2

    ustar2d(j,i) = ustar

    !----------------------------------------------------------------
    ! Step 2.  Obukhov length
    !
    ! L = -rho * cp * ustar^3 * theta_v / (vonkar * g * hfx_v)
    !
    ! where hfx_v is the virtual (buoyancy) heat flux [W m^-2]:
    !   hfx_v = hfx + rho * cp * ep1 * theta_v * (qfx / rho)
    !         = hfx + cp * ep1 * theta_v * qfx
    !
    ! theta_v at the lowest model level.
    !----------------------------------------------------------------
    exner_kz = (m2p%patm(j,i,kz) / p00) ** rovcp
    ratiomx  = m2p%qxatm(j,i,kz,iqv)
    qha      = ratiomx / (d_one + ratiomx)  ! specific humidity at kz
    theta_v  = (m2p%tatm(j,i,kz) / exner_kz) * (d_one + ep1 * qha)

    ! Virtual heat flux [W m^-2]
    hfxv = m2p%hfx(j,i) + cpd * ep1 * theta_v * m2p%qfx(j,i)

    ! Obukhov length [m]; protect against zero flux (neutral limit)
    if ( abs(hfxv) > 1.0e-4_rkx ) then
      lob = -(m2p%rhox2d(j,i) * cpd * ustar**3 * theta_v) / &
            (vonkar * egrav * hfxv)
      ! Apply realistic bounds to avoid singularity
      if ( lob >= d_zero ) then
        lob = max(lob, obmin)
      else
        lob = min(lob, obmax)
      end if
    else
      lob = 1.0e6_rkx   ! effectively neutral
    end if

    !----------------------------------------------------------------
    ! Step 3.  Stability parameter zeta = za / L
    !
    ! za is the height of the lowest half-level above the surface.
    ! zq(:,:,kzp1) is the surface interface height (AGL + orography),
    ! zq(:,:,kz)   is the lowest half-level height.
    ! za = zq(kz) - zq(kzp1)  = height above local surface [m].
    !----------------------------------------------------------------
    za = max(m2p%zq(j,i,kz) - m2p%zq(j,i,kzp1), 1.0_rkx)
    zeta = za / lob

    !----------------------------------------------------------------
    ! Step 4.  Roughness length and integrated stability functions
    !
    ! Charnock formula for z0 (Janjic 1994, eq. 8):
    !   z0 = max( ustfc * ustar^2 / g , z0_min )
    ! ustfc = 0.018/g  =>  ustfc*ustar^2 = 0.018*ustar^2/g
    ! z0_min = 1.59e-5 m (smooth viscous sublayer limit)
    !
    ! Businger-Dyer-Pandolfo functions (Janjic 1994, eqs. 22-27):
    !   Stable   (zeta > 0): psi_m = psi_h = -alpha4 * zeta
    !   Unstable (zeta < 0): psi_m = 2*ln((1+x)/2) + ln((1+x^2)/2)
    !                                 - 2*atan(x) + pi/2
    !                         psi_h = 2*ln((1+x^2)/2)
    !                         x = (1 - alpha5*zeta)^(1/4)  [for psi_m]
    !                         x = (1 - alpha6*zeta)^(1/2)  [for psi_h]
    !----------------------------------------------------------------
    z0 = max(ustfc * ust2, 1.59e-5_rkx)

    ! log(za/z0) — the neutral part of the profile integral
    logza_z0 = log(za / z0)

    if ( zeta >= d_zero ) then
      ! Stable branch: linear correction (Businger et al. 1971)
      psi_m = -alpha4 * min(zeta, 1.0_rkx)
      psi_h = psi_m
    else
      ! Unstable branch: Businger-Dyer-Pandolfo
      ! Momentum
      x4    = max(d_one - alpha5 * zeta, d_one)   ! x^4
      x     = x4 ** 0.25_rkx                       ! x = (1-16*zeta)^0.25
      psi_m = d_two * log((d_one + x) * d_half) + &
              log((d_one + x * x) * d_half) - &
              d_two * atan(x) + d_half * mathpi
      ! Heat/moisture
      x2    = max(d_one - alpha6 * zeta, d_one)   ! x^2
      x2    = sqrt(x2)                              ! x = (1-16*zeta)^0.5
      psi_h = d_two * log((d_one + x2) * d_half)
    end if

    !----------------------------------------------------------------
    ! Step 5.  Exchange coefficients
    !
    ! Denominator = log(za/z0) - psi_m (or psi_h)
    ! bounded below to prevent blow-up in extreme unstable conditions.
    !
    ! akms = vonkar * ustar / denom_m   [m/s]
    ! akhs = vonkar * ustar / denom_h   [m/s]
    !
    ! Physical interpretation:
    !   tau   = rho * akms * (U_kz - uz0)     [Pa]
    !   H     = rho * cp * akhs * (theta_kz - thz0)  [W m^-2]
    !   E     = rho * akhs * (q_kz - qz0)     [kg m^-2 s^-1]
    !----------------------------------------------------------------
    akms = vonkar * ustar / max(logza_z0 - psi_m, minak)
    akhs = vonkar * ustar / max(logza_z0 - psi_h, minak)

    ! Apply a physically motivated lower bound: akms >= Cd_min * wspd
    akms = max(akms, cdmin * wspd)
    akhs = max(akhs, cdmin * wspd)

    akms2d(j,i) = akms
    akhs2d(j,i) = akhs

    !----------------------------------------------------------------
    ! Step 6.  Viscous sublayer values uz0, vz0, thz0, qz0
    !
    ! The three-regime Janjic parameterisation (identical structure to
    ! the original mod_pbl_myj.F90 lines 355-396), but now driven by
    ! the self-consistent ustar from step 1 and akhs from step 5.
    !
    ! Regime boundaries:
    !   ustar <  ustr : smooth/viscous sublayer
    !   ustr  <= ustar < ustc : transitional
    !   ustar >= ustc : rough/wavy
    !
    ! Surface values for the tridiagonal lower BC:
    !   thz0 = potential temperature at z=z0  [K]
    !   qz0  = specific humidity    at z=z0  [kg/kg]
    !   uz0  = u-velocity           at z=z0  [m/s]
    !   vz0  = v-velocity           at z=z0  [m/s]
    !----------------------------------------------------------------
    psfc     = m2p%patmf(j,i,kzp1)
    rexnsfc  = (p00 / psfc) ** rovcp   ! 1/Exner at surface
    tg       = m2p%tg(j,i)
    thsk     = tg * rexnsfc            ! surface potential temperature [K]

    if ( m2p%ldmsk(j,i) == 0 ) then
      qsfc = seafc * pfqsat(tg, psfc)  ! ocean: 98% of saturation
    else
      qsfc = m2p%q2m(j,i)              ! land:  2-m specific humidity
    end if

    ! Lowest model-level thermodynamic values
    tha = m2p%tatm(j,i,kz) / exner_kz   ! potential temperature at kz
    ! qha already computed above

    if ( ustar < ustr ) then
      !------------------------------------------------------------------
      ! Smooth/viscous regime
      ! Roughness lengths for momentum, heat and moisture
      ! in the viscous sublayer (Janjic 1994, eqs. 12-14):
      !   zu = fzu1 * sqrt(sqrt(z0 * ustar / nu)) / ustar
      !   zt = fzt1 * zu
      !   zq = fzq1 * zt
      ! These are viscous sublayer thicknesses.
      !
      ! The weighting factor wght blends the sublayer value toward the
      ! interior model value as akms*zu/nu increases.
      !------------------------------------------------------------------
      zu = fzu1 * sqrt(sqrt(z0 * ustar * rvisc)) / ustar
      wght  = akms * zu * rvisc
      wght  = wght / (d_one + wght)

      ! uz0, vz0: viscous sublayer velocity (blended with previous step)
      m2p%uz0(j,i) = d_half * (uatm * wght + m2p%uz0(j,i))
      m2p%vz0(j,i) = d_half * (vatm * wght + m2p%vz0(j,i))

      zt = fzt1 * zu
      zq = fzq1 * zt
      wghtt = akhs * zt * rtvisc
      wghtq = akhs * zq * rqvisc

      if ( rcmtimer%lcount < 1 ) then
        ! First call: no blending
        m2p%thz0(j,i) = (wghtt * tha + thsk) / (wghtt + d_one)
        m2p%qz0(j,i)  = (wghtq * qha + qsfc) / (wghtq + d_one)
      else
        m2p%thz0(j,i) = d_half * ((wghtt * tha + thsk) / &
                                   (wghtt + d_one) + m2p%thz0(j,i))
        m2p%qz0(j,i)  = d_half * ((wghtq * qha + qsfc) / &
                                   (wghtq + d_one) + m2p%qz0(j,i))
      end if

    else if ( ustar >= ustr .and. ustar < ustc ) then
      !------------------------------------------------------------------
      ! Transitional regime
      ! uz0 = vz0 = 0 (no-slip lost; wave motion sets surface velocity).
      ! Scalar roughness follows Kader-Yaglom parameterisation.
      !------------------------------------------------------------------
      m2p%uz0(j,i) = d_zero
      m2p%vz0(j,i) = d_zero

      zt = fzt2 * sqrt(sqrt(z0 * ustar * rvisc)) / ustar
      zq = fzq2 * zt
      wghtt = akhs * zt * rtvisc
      wghtq = akhs * zq * rqvisc

      if ( rcmtimer%lcount < 1 ) then
        m2p%thz0(j,i) = (wghtt * tha + thsk) / (wghtt + d_one)
        m2p%qz0(j,i)  = (wghtq * qha + qsfc) / (wghtq + d_one)
      else
        m2p%thz0(j,i) = d_half * ((wghtt * tha + thsk) / &
                                   (wghtt + d_one) + m2p%thz0(j,i))
        m2p%qz0(j,i)  = d_half * ((wghtq * qha + qsfc) / &
                                   (wghtq + d_one) + m2p%qz0(j,i))
      end if

    else
      !------------------------------------------------------------------
      ! Rough/wavy regime (ustar >= ustc)
      ! Surface values set directly to their bulk values.
      !------------------------------------------------------------------
      m2p%uz0(j,i)  = d_zero
      m2p%vz0(j,i)  = d_zero
      m2p%thz0(j,i) = thsk
      m2p%qz0(j,i)  = qsfc
    end if

  end do

  contains

  ! pfqsat is needed for qsfc over ocean
#include <pfqsat.inc>

  end subroutine myj_sfclayer

end module mod_pbl_myj_sflayer

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
