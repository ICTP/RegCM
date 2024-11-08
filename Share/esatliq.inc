  ! Computes saturation pressure over flat surface of liquid water
  ! Reference: Maarten H. P. Ambaum
  !    Accurate, simple equation for saturated vapour pressure over
  !    water and ice , QJRMS, 2020, https://doi.org/10.1002/qj.3899
  !
!DIR$ ATTRIBUTES FORCEINLINE :: esatliq
  pure elemental real(rkx) function esatliq(t) result(es)
    use mod_constants , only : rwat , wlhv
    implicit none
    real(rkx) , intent(in) :: t     ! Temperature (K)

    real(rkx) , parameter :: es0 = 611.655_rkx ! Pa
    real(rkx) , parameter :: cplmcpv = 2180.0_rkx ! J kg-1 K-1
    real(rkx) , parameter :: arg = cplmcpv/rwat
    real(rkx) , parameter :: t0 = 273.16_rkx
    real(rkx) , parameter :: lv0 = 2.501e6_rkx
    real(rkx) , parameter :: al0 = lv0/(rwat*t0)
    real(rkx) :: lv

    lv = lv0-(cplmcpv)*(t-t0)
    es = es0 * (t0/t)**arg * exp(al0-lv/(rwat*t))
  end function esatliq

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
