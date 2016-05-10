  !
  ! Computes derivative in temperature of saturation pressure
  !
  pure real(rkx) function pfesdt(t) result(esdt)
    implicit none
    real(rkx), intent(in)  :: t     ! Temperature (K)

    real(rkx) :: td , t_limit
    !
    ! For derivative:water vapor
    !
    real(rkx), parameter :: b0 =  0.444017302_rkx
    real(rkx), parameter :: b1 =  0.286064092e-01_rkx
    real(rkx), parameter :: b2 =  0.794683137e-03_rkx
    real(rkx), parameter :: b3 =  0.121211669e-04_rkx
    real(rkx), parameter :: b4 =  0.103354611e-06_rkx
    real(rkx), parameter :: b5 =  0.404125005e-09_rkx
    real(rkx), parameter :: b6 = -0.788037859e-12_rkx
    real(rkx), parameter :: b7 = -0.114596802e-13_rkx
    real(rkx), parameter :: b8 =  0.381294516e-16_rkx
    !
    ! For derivative:ice
    !
    real(rkx), parameter :: d0 =  0.503277922_rkx
    real(rkx), parameter :: d1 =  0.377289173e-01_rkx
    real(rkx), parameter :: d2 =  0.126801703e-02_rkx
    real(rkx), parameter :: d3 =  0.249468427e-04_rkx
    real(rkx), parameter :: d4 =  0.313703411e-06_rkx
    real(rkx), parameter :: d5 =  0.257180651e-08_rkx
    real(rkx), parameter :: d6 =  0.133268878e-10_rkx
    real(rkx), parameter :: d7 =  0.394116744e-13_rkx
    real(rkx), parameter :: d8 =  0.498070196e-16_rkx

    t_limit = t - tzero
    if ( t_limit > 100.0_rkx ) t_limit = 100.0_rkx
    if ( t_limit < -75.0_rkx ) t_limit = -75.0_rkx
    td = t_limit
    if ( td >= 0.0_rkx ) then
      esdt = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
           + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
      esdt = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
           + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    end if
    esdt = esdt * 100.0_rkx ! pa/K
  end function pfesdt

  pure real(rkx) function pfqsdt(t,p,e,dedt) result(qsdt)
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) , intent(in) , optional :: dedt ! derivative of e in dt (Pa/K)
    real(rkx) :: es , esdt , vp , vp1 , vp2
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    if ( present(dedt) ) then
      esdt = dedt
    else
      esdt = pfesdt(t)
    end if
    vp  = 1.0_rkx / (p - 0.378_rkx*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    qsdt = esdt * vp2 * p ! 1 / K
  end function pfqsdt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
