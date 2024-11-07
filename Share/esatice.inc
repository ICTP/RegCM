  ! Computes saturation pressure over flat surface of liquid water
  ! Reference: Maarten H. P. Ambaum
  !    Accurate, simple equation for saturated vapour pressure over
  !    water and ice , QJRMS, 2020, https://doi.org/10.1002/qj.3899
  !
!DIR$ ATTRIBUTES FORCEINLINE :: esatice
  pure elemental real(rkx) function esatice(t) result(es)
    use mod_constants , only : rwat
    implicit none
    real(rkx) , intent(in) :: t     ! Temperature (K)

    real(rkx) , parameter :: ei0 = 611.655_rkx ! Pa
    real(rkx) , parameter :: cpimcpv = 212.0_rkx ! J kg-1 K-1
    real(rkx) , parameter :: arg = cpimcpv/rwat
    real(rkx) , parameter :: t0 = 273.16_rkx
    real(rkx) , parameter :: li0 = 2.835e6_rkx ! J kg-1
    real(rkx) , parameter :: ai0 = li0/(rwat*t0)
    real(rkx) :: li
    li = li0-(cpimcpv)*(t-t0)
    es = ei0 * (t0/t)**arg * exp(ai0-li/(rwat*t))
  end function esatice

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
