  !
  ! Computes derivative in temperature of saturation pressure
  !
!DIR$ ATTRIBUTES FORCEINLINE :: pfesdt
  pure elemental real(rkx) function pfesdt(t) result(esdt)
    implicit none
    real(rkx), intent(in)  :: t     ! Temperature (K)

    real(rk8) :: td , t_limit , esatdt
    !
    ! For derivative:water vapor
    !
    real(rk8), parameter :: b0 =  0.444017302_rk8
    real(rk8), parameter :: b1 =  0.286064092e-01_rk8
    real(rk8), parameter :: b2 =  0.794683137e-03_rk8
    real(rk8), parameter :: b3 =  0.121211669e-04_rk8
    real(rk8), parameter :: b4 =  0.103354611e-06_rk8
    real(rk8), parameter :: b5 =  0.404125005e-09_rk8
    real(rk8), parameter :: b6 = -0.788037859e-12_rk8
    real(rk8), parameter :: b7 = -0.114596802e-13_rk8
    real(rk8), parameter :: b8 =  0.381294516e-16_rk8
    !
    ! For derivative:ice
    !
    real(rk8), parameter :: d0 =  0.503277922_rk8
    real(rk8), parameter :: d1 =  0.377289173e-01_rk8
    real(rk8), parameter :: d2 =  0.126801703e-02_rk8
    real(rk8), parameter :: d3 =  0.249468427e-04_rk8
    real(rk8), parameter :: d4 =  0.313703411e-06_rk8
    real(rk8), parameter :: d5 =  0.257180651e-08_rk8
    real(rk8), parameter :: d6 =  0.133268878e-10_rk8
    real(rk8), parameter :: d7 =  0.394116744e-13_rk8
    real(rk8), parameter :: d8 =  0.498070196e-16_rk8

    t_limit = t - tzero
    if ( t_limit > 100.0_rk8 ) t_limit = 100.0_rk8
    if ( t_limit < -75.0_rk8 ) t_limit = -75.0_rk8
    td = t_limit
    if ( td >= 0.0_rk8 ) then
      esatdt = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
           + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
      esatdt = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
           + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    end if
    esdt = real(esatdt * 100.0_rk8,rkx) ! pa/K
  end function pfesdt

!DIR$ ATTRIBUTES FORCEINLINE :: pfesdt_water
  pure elemental real(rkx) function pfesdt_water(t) result(esdt)
    implicit none
    real(rkx), intent(in)  :: t     ! Temperature (K)

    real(rk8) :: td , esatdt
    !
    ! For derivative:water vapor
    !
    real(rk8), parameter :: b0 =  0.444017302_rk8
    real(rk8), parameter :: b1 =  0.286064092e-01_rk8
    real(rk8), parameter :: b2 =  0.794683137e-03_rk8
    real(rk8), parameter :: b3 =  0.121211669e-04_rk8
    real(rk8), parameter :: b4 =  0.103354611e-06_rk8
    real(rk8), parameter :: b5 =  0.404125005e-09_rk8
    real(rk8), parameter :: b6 = -0.788037859e-12_rk8
    real(rk8), parameter :: b7 = -0.114596802e-13_rk8
    real(rk8), parameter :: b8 =  0.381294516e-16_rk8

    td = min(t-tzero,100.0_rk8)
    esatdt = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
           + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    esdt = real(esatdt * 100.0_rk8,rkx) ! pa/K
  end function pfesdt_water

!DIR$ ATTRIBUTES FORCEINLINE :: pfesdt_ice
  pure elemental real(rkx) function pfesdt_ice(t) result(esdt)
    implicit none
    real(rkx), intent(in)  :: t     ! Temperature (K)

    real(rk8) :: td , esatdt
    !
    ! For derivative:ice
    !
    real(rk8), parameter :: d0 =  0.503277922_rk8
    real(rk8), parameter :: d1 =  0.377289173e-01_rk8
    real(rk8), parameter :: d2 =  0.126801703e-02_rk8
    real(rk8), parameter :: d3 =  0.249468427e-04_rk8
    real(rk8), parameter :: d4 =  0.313703411e-06_rk8
    real(rk8), parameter :: d5 =  0.257180651e-08_rk8
    real(rk8), parameter :: d6 =  0.133268878e-10_rk8
    real(rk8), parameter :: d7 =  0.394116744e-13_rk8
    real(rk8), parameter :: d8 =  0.498070196e-16_rk8

    td = max(t-tzero,-75.0_rk8)
    esatdt = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
           + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    esdt = real(esatdt * 100.0_rk8,rkx) ! pa/K
  end function pfesdt_ice

