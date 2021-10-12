!DIR$ ATTRIBUTES FORCEINLINE :: lapse_rate
  pure elemental real(rkx) function lapse_rate(jday,lat) result(gammac)
    use mod_constants , only : mathpi
    implicit none
    real(rkx) , intent(in) :: jday , lat
    real(rk8) :: coeday
    if ( lat > 30.0_rkx ) then
      coeday = 0.5_rk8*(1.0_rk8-cos(jday*2.0_rk8*mathpi/365.0_rk8))   ! N.H.
    else if ( lat < -30.0_rkx ) then
      coeday = 0.5_rk8*(1.0_rk8+cos(jday*2.0_rk8*mathpi/365.0_rk8))   ! S.H.
    else
      coeday = 0.5_rk8                        ! tropics
    endif
    gammac = 4.5e-3_rk8 + 2.8e-3_rk8*coeday
  end function lapse_rate

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
