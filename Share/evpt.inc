  ! Potential Evapotranspiration
  ! Shuttleworth, J., Putting the vap' into evaporation (1993)
  ! Return value in kg m-2 s-1
  pure real(rkx) function evpt(p,t,u,e,es,m,d)
    implicit none
    real(rkx) , intent(in) :: p   ! Pressure in Pa
    real(rkx) , intent(in) :: t   ! Pressure in K
    real(rkx) , intent(in) :: u   ! Wind speed in m/s
    real(rkx) , intent(in) :: e   ! Vapor pressure in Pa
    real(rkx) , intent(in) :: es  ! Saturation vapor pressure in Pa
    real(rkx) , intent(in) :: m   ! Saturation vapor pressure derivative Pa K-1
    real(rk8) , intent(in) :: d   ! Daily solar irradiance W m-2 s-1

    ! Latent heat of vaporization
    ! (MJ kg-1 , FAO-56 estimation for Penman-Monteith Eq.
    real(rkx) , parameter :: lambd = 2.45_rkx
    real(rkx) :: mkpa , pkpa , xgamma , rn , des

    ! Port to kPa
    pkpa = p * 1.0e-3
    mkpa = m * 1.0e-3
    des = max((es-e) * 1.0e-3,d_zero)
    ! Units kg m-2 day-1
    rn = d/wlh(t)*86400.0_rkx
    ! Compute psychrometric constant (kPa K-1)
    xgamma = (0.0016286_rkx*pkpa)/lambd
    ! Potential evapotranspiration in kg m-2 day-1
    evpt = ( (xgamma*6.43_rkx*(1.0_rkx+0.536_rkx*u)*des)/lambd + &
             mkpa*rn ) / (mkpa+xgamma)
    evpt = evpt / 86400_rkx
  end function evpt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
