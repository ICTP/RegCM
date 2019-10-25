  ! Computes saturation pressurre
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
!DIR$ ATTRIBUTES FORCEINLINE :: esatliq
  pure elemental real(rkx) function esatliq(t) result(es)
    implicit none
    real(rkx) , intent(in) :: t     ! Temperature (K)

    real(rk8) :: td , t_limit , esat
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8) , parameter :: a0 =  6.11213476_rk8
    real(rk8) , parameter :: a1 =  0.444007856_rk8
    real(rk8) , parameter :: a2 =  0.143064234e-01_rk8
    real(rk8) , parameter :: a3 =  0.264461437e-03_rk8
    real(rk8) , parameter :: a4 =  0.305903558e-05_rk8
    real(rk8) , parameter :: a5 =  0.196237241e-07_rk8
    real(rk8) , parameter :: a6 =  0.892344772e-10_rk8
    real(rk8) , parameter :: a7 = -0.373208410e-12_rk8
    real(rk8) , parameter :: a8 =  0.209339997e-15_rk8

    t_limit = t - tzero
    if ( t_limit > 100.0_rk8 ) t_limit = 100.0_rk8
    td = t_limit
    esat = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    es = real(esat * 100.0_rk8,rkx) ! pa
  end function esatliq

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
