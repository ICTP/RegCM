  ! Computes saturation pressurre
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
!DIR$ ATTRIBUTES FORCEINLINE :: esatice
  pure elelmental real(rkx) function esatice(t) result(es)
    implicit none
    real(rkx) , intent(in) :: t     ! Temperature (K)

    real(rk8) :: td , t_limit , esat
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8) , parameter :: c0 =  6.11123516_rk8
    real(rk8) , parameter :: c1 =  0.503109514_rk8
    real(rk8) , parameter :: c2 =  0.188369801e-01_rk8
    real(rk8) , parameter :: c3 =  0.420547422e-03_rk8
    real(rk8) , parameter :: c4 =  0.614396778e-05_rk8
    real(rk8) , parameter :: c5 =  0.602780717e-07_rk8
    real(rk8) , parameter :: c6 =  0.387940929e-09_rk8
    real(rk8) , parameter :: c7 =  0.149436277e-11_rk8
    real(rk8) , parameter :: c8 =  0.262655803e-14_rk8

    t_limit = t - tzero
    if ( t_limit < -75.0_rk8 ) t_limit = -75.0_rk8
    td = t_limit
    esat = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    es = real(esat * 100.0_rk8,rkx) ! pa
  end function esatice

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
