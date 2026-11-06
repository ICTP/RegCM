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
  ! Surface-layer quantities for the MYJ and Shin-Hong PBL schemes.
  !
  ! The surface stress and the heat and moisture fluxes (uvdrag, hfx, qfx)
  ! are the ones of the RegCM land and ocean surface schemes, which the PBL
  ! schemes apply as they are. From them this module derives the friction
  ! velocity and the Monin-Obukhov stability parameter.
  !
  ! REFERENCES:
  !   Janjic (1994), Mon. Wea. Rev., 122, 927-945.
  !   Janjic (2002), NCEP Office Note 437.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_regcm_types
  use mod_runparams, only : iqv

  implicit none

  private

  public :: myj_sfclayer, myj_sfclayer_zol

  contains

  subroutine myj_sfclayer(m2p, ustar2d)
  !
  ! PURPOSE
  !   The friction velocity from the surface stress of the surface scheme:
  !     uvdrag [kg m^-2 s^-1] = rho * Cd * |U|
  !     tau    [Pa]           = uvdrag * |U|
  !     ustar  [m/s]          = sqrt(tau / rho)
  !
  implicit none

  type(mod_2_pbl), intent(in) :: m2p
  real(rkx), dimension(jci1:jci2,ici1:ici2), intent(out) :: ustar2d

  integer(ik4) :: i, j
  real(rkx) :: uatm, vatm, wspd, tau_over_rho

  do concurrent ( j = jci1:jci2, i = ici1:ici2 )
    uatm = m2p%uxatm(j,i,kz)
    vatm = m2p%vxatm(j,i,kz)
    wspd = max(sqrt(uatm*uatm + vatm*vatm), 0.01_rkx)
    tau_over_rho = max(m2p%uvdrag(j,i) * wspd / m2p%rhox2d(j,i), 1.0e-8_rkx)
    ustar2d(j,i) = sqrt(tau_over_rho)
  end do

  end subroutine myj_sfclayer

  subroutine myj_sfclayer_zol(m2p, ustar2d, zol2d)
  !
  ! PURPOSE
  !   Surface-layer diagnostic for the Shin-Hong scheme: the friction
  !   velocity and the Obukhov stability parameter zeta = za/L.
  !
  ! OUTPUTS
  !   ustar2d  - friction velocity [m/s]
  !   zol2d    - stability parameter za/L  (> 0 stable, < 0 unstable)
  !
  ! METHOD
  !   ustar is derived from the surface stress (uvdrag), as in
  !   myj_sfclayer.
  !
  !   The Obukhov length is:
  !     L = -rho * cp * ustar^3 * theta_v / (vonkar * g * H_v)
  !   where H_v = hfx + cp * ep1 * theta_v * qfx  (virtual heat flux [W/m^2]).
  !   Hence:
  !     za/L = -vonkar * g * za * H_v / (rho * cp * ustar^3 * theta_v)
  !          = -vonkar * g * za * hfxv_kin / (ustar^3 * theta_v)
  !   where hfxv_kin = H_v / (rho * cp)  [K m/s] is the kinematic virtual flux.
  !
  !   Result is bounded to the physically safe interval [-100, 100].
  !
  implicit none

  type(mod_2_pbl), intent(in)  :: m2p
  real(rkx), dimension(jci1:jci2,ici1:ici2), intent(out) :: ustar2d
  real(rkx), dimension(jci1:jci2,ici1:ici2), intent(out) :: zol2d

  integer(ik4) :: i, j
  real(rkx) :: uatm, vatm, wspd, tau_over_rho, ustar, ust3
  real(rkx) :: hfxv_kin, theta_v, exner_kz, ratiomx, qha, za

  do concurrent ( j = jci1:jci2, i = ici1:ici2 )

    !------------------------------------------------------------
    ! Friction velocity from surface stress  (as in myj_sfclayer)
    !------------------------------------------------------------
    uatm = m2p%uxatm(j,i,kz)
    vatm = m2p%vxatm(j,i,kz)
    wspd = max(sqrt(uatm*uatm + vatm*vatm), 0.01_rkx)
    tau_over_rho = max(m2p%uvdrag(j,i) * wspd / m2p%rhox2d(j,i), 1.0e-8_rkx)
    ustar = sqrt(tau_over_rho)
    ustar2d(j,i) = ustar
    ust3  = ustar**3

    !------------------------------------------------------------
    ! Virtual potential temperature at the lowest model level
    !------------------------------------------------------------
    exner_kz = (m2p%patm(j,i,kz) / p00) ** rovcp
    ratiomx  = m2p%qxatm(j,i,kz,iqv)
    qha      = ratiomx / (d_one + ratiomx)
    theta_v  = (m2p%tatm(j,i,kz) / exner_kz) * (d_one + ep1 * qha)

    !------------------------------------------------------------
    ! Kinematic virtual heat flux [K m/s]
    !   hfxv_kin = (hfx + cp * ep1 * theta_v * qfx) / (rho * cp)
    !------------------------------------------------------------
    hfxv_kin = (m2p%hfx(j,i) + cpd * ep1 * theta_v * m2p%qfx(j,i)) / &
               (m2p%rhox2d(j,i) * cpd)

    !------------------------------------------------------------
    ! Stability parameter zeta = za/L
    !   za is the height of the lowest model level above the surface
    !------------------------------------------------------------
    za = max(m2p%za(j,i,kz), 1.0_rkx)

    if ( abs(hfxv_kin) > 1.0e-6_rkx ) then
      ! zol = za/L = -vonkar*g*za*hfxv / (ustar^3 * theta_v)
      zol2d(j,i) = -vonkar * egrav * za * hfxv_kin / (ust3 * theta_v)
      zol2d(j,i) = max(min(zol2d(j,i), 100.0_rkx), -100.0_rkx)
    else
      zol2d(j,i) = 0.0_rkx   ! neutral
    end if

  end do

  end subroutine myj_sfclayer_zol

end module mod_pbl_myj_sflayer

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
