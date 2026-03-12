  ! Temperature dependency for cloud water content
  ! in g/m3 (Lemus et al., 1997)
  ! NOTE : THIS IS IN-CLOUD VARIABLE.
!DIR$ ATTRIBUTES FORCEINLINE :: clwfromt
  pure elemental real(rkx) function clwfromt(t) result(clw)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t
    real(rk8) :: tcel, xclw
    tcel = t - tzero
    if ( tcel < -50.0_rk8 ) then
      xclw = 0.001_rk8
    else
      xclw = 0.127_rk8 + 6.78e-03_rk8 * tcel +    &
                        1.29e-04_rk8 * tcel**2 + &
                        8.68e-07_rk8 * tcel**3
      ! Diferrent formulas
      !xclw = exp(-0.48911_rk8 + 0.01244_rk8 * tcel)
      !xclw = max(3.0e-4_rk8,0.18_rk8 * min(1.0_rk8,(tcel+53.15_rk8)/60.0_rk8))
    end if
    clw = real(max(min(xclw,0.1_rk8),0.001_rk8),rkx)
  end function clwfromt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
