  !
  ! Computes derivative in temperature of saturation pressure
  !
  pure real(rk8) function pfesdt(t) result(esdt)
    implicit none
    real(rk8), intent(in)  :: t     ! Temperature (K)

    real(rk8) :: td , t_limit
    !
    ! For derivative:water vapor
    !
    real(rk8), parameter :: b0 =  0.444017302D0
    real(rk8), parameter :: b1 =  0.286064092D-01
    real(rk8), parameter :: b2 =  0.794683137D-03
    real(rk8), parameter :: b3 =  0.121211669D-04
    real(rk8), parameter :: b4 =  0.103354611D-06
    real(rk8), parameter :: b5 =  0.404125005D-09
    real(rk8), parameter :: b6 = -0.788037859D-12
    real(rk8), parameter :: b7 = -0.114596802D-13
    real(rk8), parameter :: b8 =  0.381294516D-16
    !
    ! For derivative:ice
    !
    real(rk8), parameter :: d0 =  0.503277922D0
    real(rk8), parameter :: d1 =  0.377289173D-01
    real(rk8), parameter :: d2 =  0.126801703D-02
    real(rk8), parameter :: d3 =  0.249468427D-04
    real(rk8), parameter :: d4 =  0.313703411D-06
    real(rk8), parameter :: d5 =  0.257180651D-08
    real(rk8), parameter :: d6 =  0.133268878D-10
    real(rk8), parameter :: d7 =  0.394116744D-13
    real(rk8), parameter :: d8 =  0.498070196D-16

    t_limit = t - tzero
    if ( t_limit > 100.0D0 ) t_limit = 100.0D0
    if ( t_limit < -75.0D0 ) t_limit = -75.0D0
    td = t_limit
    if ( td >= 0.0D0 ) then
      esdt = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
           + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
      esdt = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
           + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    end if
    esdt = esdt * 100.D0 ! pa/K
  end function pfesdt

  pure real(rk8) function pfqsdt(t,p,e,dedt) result(qsdt)
    implicit none
    real(rk8) , intent(in) :: t             ! Temperature (K)
    real(rk8) , intent(in) :: p             ! Pressure (Pa)
    real(rk8) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rk8) , intent(in) , optional :: dedt ! derivative of e in dt (Pa/K)
    real(rk8) :: es , esdt , vp , vp1 , vp2
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
    vp  = 1.0D0 / (p - 0.378D0*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    qsdt = esdt * vp2 * p ! 1 / K
  end function pfqsdt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
