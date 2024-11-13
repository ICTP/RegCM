  !
  ! Computes moist lapse rate
  !
!DIR$ ATTRIBUTES FORCEINLINE :: mlrate
  pure elemental real(rkx) function mlrate(t,q) result(lr)
!$acc routine seq
    use mod_constants , only : egrav , wlhv , cpd , rdry , rwat
    implicit none
    ! Temperature (K) , Pressure (Pa) , q Mixing Ratio kg kg-1
    real(rkx) , intent(in) :: t , p , q
    real(rkx) , parameter :: wlhvsq = wlhv*wlhv

    lr = egrav * (d_one + (wlhv*q)/(rdry*t))/(cpd+(wlhvsq*q)/(rwat*t*t))
  end function mlrate

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
