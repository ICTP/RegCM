  !
  ! Computes moist lapse rate
  !
!DIR$ ATTRIBUTES FORCEINLINE :: mlrate
  pure elemental real(rkx) function mlrate(t,q) result(lr)
!$acc routine seq
    use mod_constants , only : egrav, cpd, rdry, rwat
    implicit none
    ! Temperature (K), q Mixing Ratio kg kg-1
    real(rkx) , intent(in) :: t, q
    real(rk8) :: tc , lheat, wlh, wlhsq
    tc = t - tzero
    if ( tc > 0.0_rk8 ) then
      lheat = 2500.79_rk8 - 2.36418_rk8*tc + &
            1.58927e-3_rk8*tc*tc - 6.14342e-5_rk8*tc*tc*tc
    else
      lheat = 2834.1_rk8 - 0.29_rk8*tc - 0.004_rk8*tc*tc
    end if
    wlh = lheat*1.0e3_rk8
    wlhsq = wlh*wlh
    lr = egrav * (d_one + (wlh*q)/(rdry*t))/(cpd+(wlhsq*q)/(rwat*t*t))
  end function mlrate

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
