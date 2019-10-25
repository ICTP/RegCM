  ! Potential Evapotranspiration
  ! Shuttleworth, J., Putting the vap' into evaporation (1993)
  ! Return value in kg m-2 s-1
  pure elemental real(rkx) function evpt(p,t,u,e,es,m,sw,lw)
    implicit none
    real(rkx) , intent(in) :: p   ! Pressure in Pa
    real(rkx) , intent(in) :: t   ! Temperature in K
    real(rkx) , intent(in) :: u   ! Wind speed in m/s
    real(rkx) , intent(in) :: e   ! Vapor pressure in Pa
    real(rkx) , intent(in) :: es  ! Saturation vapor pressure in Pa
    real(rkx) , intent(in) :: m   ! Saturation vapor pressure derivative Pa K-1
    real(rkx) , intent(in) :: sw  ! Surface Solar net Shortwave
    real(rkx) , intent(in) :: lw  ! Surface net Longwave

    real(rkx) :: mkpa , pkpa , xgamma , rn , des , lath , lambd

    lath = wlh(t)
    lambd = lath * 1.0e-6
    ! Port to kPa
    pkpa = p * 1.0e-3
    mkpa = m * 1.0e-3
    des = max(((es-e) * 1.0e-3),d_zero)
    rn = max(sw-lw,0.0_rk8) * 1.e-6_rk8      ! Mj s-1
    ! Units kg m-2 day-1
    rn = (rn/lath)*86400.0_rkx
    ! Compute psychrometric constant (kPa K-1)
    xgamma = (0.0016286_rkx*pkpa)/lambd
    ! Potential evapotranspiration in kg m-2 day-1
    evpt = ( (xgamma * 6.43_rkx * (1.0_rkx + 0.536_rkx*u) * des)/lambd + &
             mkpa*rn ) / (mkpa+xgamma)
    evpt = evpt / 86400.0_rkx
  end function evpt

  ! ASCE-EWRI for Reference Evapotranspiration for grass
  pure real(rkx) function evpt_fao(p,t,u,e,es,m,sw,lw)
    implicit none
    real(rkx) , intent(in) :: p   ! Pressure in Pa
    real(rkx) , intent(in) :: t   ! 2m Temperature in K
    real(rkx) , intent(in) :: u   ! 10 m Wind speed in m/s
    real(rkx) , intent(in) :: e   ! Vapor pressure in Pa
    real(rkx) , intent(in) :: es  ! Saturation vapor pressure in Pa
    real(rkx) , intent(in) :: m   ! Saturation vapor pressure derivative Pa K-1
    real(rkx) , intent(in) :: sw  ! Surface Solar net Shortwave
    real(rkx) , intent(in) :: lw  ! Surface net Longwave

    real(rk8) :: pkpa , rn , u2 , delta , xgamma , ut , dt , pt , tt
    real(rk8) :: ekpa , eskpa , de , etr , etw

    rn = (sw-lw) * 1.e-6_rk8           ! Mj s-1
    pkpa = p * 0.001_rk8               ! kPa
    ekpa = e * 0.001_rk8               ! kPa
    eskpa = es * 0.001_rk8             ! kPa
    de = eskpa - ekpa                  ! kPa
    u2 = max(u * 0.7342_rk8,0.5_rk8)   ! 4.87/ln(67.8*h-5.42) m s-1
    rn = 0.408_rk8 * rn                ! mm s-1
    delta = m * 0.001_rk8              ! kPa K-1
    xgamma = 0.000665_rk8 * pkpa       ! Psychrometric constant kPa K-1
    if ( rn > 0.0_rkx ) then
      ut = (1.0_rk8 + 0.24_rk8 * u2)   ! Modificator for grass
      rn = rn * 0.9_rk8
    else
      ut = (1.0_rk8 + 0.96_rk8 * u2)   ! Modificator for grass
      rn = rn * 0.5_rk8
    end if
    dt = delta / (delta + xgamma * ut)
    pt = xgamma / (delta + xgamma * ut)
    tt = 0.010417_rk8 * ( u2/t )       ! 900/86400
    etw = pt * tt * de
    etr = dt * rn
    evpt_fao = (etw+etr)               ! mm s-1
  end function evpt_fao

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
