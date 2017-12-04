  ! Potential Evapotranspiration
  ! Shuttleworth, J., Putting the vap' into evaporation (1993)
  ! Return value in kg m-2 s-1
  pure real(rkx) function evpt(t,p,e,es,m,d,l,u)
    implicit none
    real(rkx) , intent(in) :: t   ! Temperature in K
    real(rkx) , intent(in) :: p   ! Pressure in Pa
    real(rkx) , intent(in) :: e   ! Vapor pressure in Pa
    real(rkx) , intent(in) :: es  ! Saturation vapor pressure in Pa
    real(rkx) , intent(in) :: m   ! Saturation vapor pressure derivative Pa K-1
    real(rkx) , intent(in) :: d   ! Declination angle (rad)
    real(rkx) , intent(in) :: l   ! Latitude (rad)
    real(rkx) , intent(in) :: u   ! Wind speed in m/s

    ! Latent heat of vaporization
    ! (MJ kg-1 , FAO-56 estimation for Penman-Monteith Eq.
    real(rkx) , parameter :: lambd = 2.45_rkx
    real(rkx) :: mkpa , pkpa , gamma , rn , hdl , temp

    ! Port to kPa
    pkpa = p * 1.0e-3
    mkpa = m * 1.0e-3
    ! Half day lenght
    temp = -(sin(l)*sin(d))/(cos(l)*cos(d))
    temp = min(1.0_rkx,max(-1.0_rkx,temp))
    hdl = acos(temp)
    ! Daily insolation (Bristow and Campbell 1984)
    rn = (86400.0_rkx*1360.0_rkx*(hdl*sin(l)*sin(d) + &
                                  cos(l)*cos(d)*sin(hdl)))/3.1416_rkx
    ! Port to MJ m-2 day-1
    rn = 1.0e-6_rkx * rn
    ! Compute psychrometric constant (kPa K-1)
    gamma = (0.0016286_rkx*pkpa)/lambd
    ! Potential evapotranspiration in kg m-2 day-1
    evpt = (mkpa*rn + gamma*6.43_rkx*(1.0_rkx+0.536_rkx*u)*(es-e)) / &
            (lambd*(mkpa+gamma))
    evpt = evpt / 86400_rkx
  end function evpt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
