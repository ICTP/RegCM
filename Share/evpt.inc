  ! Potential Evapotranspiration
  ! Shuttleworth, J., Putting the vap' into evaporation (1993)
  ! Return value in kg m-2 s-1
  pure elemental real(rkx) function evpt(p,t,u,e,es,m,sw,lw)
!$acc routine seq
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

  ! ASCE-EWRI for Reference Evapotranspiration
  ! https://ascelibrary.org/doi/book/10.1061/9780784408056
  pure real(rkx) function evpt_fao(p,t2,u10,e,es,m,sw,lw)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: p   ! Pressure in Pa
    real(rkx) , intent(in) :: t2  ! 2m Temperature in K
    real(rkx) , intent(in) :: u10 ! 10m Wind speed in m/s
    real(rkx) , intent(in) :: e   ! Vapor pressure in Pa
    real(rkx) , intent(in) :: es  ! Saturation vapor pressure in Pa
    real(rkx) , intent(in) :: m   ! Saturation vapor pressure derivative Pa K-1
    real(rkx) , intent(in) :: sw  ! Surface Solar net Shortwave
    real(rkx) , intent(in) :: lw  ! Surface net Longwave

    real(rk8) :: pkpa , rn , u2 , delta , xgamma , ut , dt , pt , tt
    real(rk8) :: de , etr , etw

    rn = (sw-lw) * 0.00360_rk8         ! MJ m-2 h-1
    pkpa = p * 0.001_rk8               ! kPa
    de = (es - e) * 0.001_rk8          ! Saturation deficit kPa
    delta = m * 0.001_rk8              ! kPa K-1
    xgamma = 0.000665_rk8 * pkpa       ! Psychrometric constant kPa K-1
    u2 = max(u10 * 0.7342_rk8,0.1_rk8) ! 4.87/ln(67.8*h-5.42) m s-1
    evpt_fao = ((0.408_rk8*delta*rn+xgamma*(37.0_rk8/t2)*u2*de) / &
                (delta + xgamma*(1.0_rkx+0.34_rk8*u2)))/3600.0_rkx
  end function evpt_fao

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
