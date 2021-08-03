  ! Temperature dependency for cloud water content
  ! in g/m3 (Lemus et al., 1997)
  ! NOTE : THIS IS IN-CLOUD VARIABLE.
!DIR$ ATTRIBUTES FORCEINLINE :: clwfromt
  pure elemental real(rkx) function clwfromt(t) result(clw)
    implicit none
    real(rkx) , intent(in) :: t
    real(rkx) :: tcel
    tcel = t - tzero
    if ( tcel < -50.0_rkx ) then
      clw = 0.001_rkx
    else
      clw = 0.127_rkx + 6.78e-03_rkx * tcel +    &
                        1.29e-04_rkx * tcel**2 + &
                        8.68e-07_rkx * tcel**3
    end if
    ! Diferrent formulas
    !clw = exp(-0.48911_rkx + 0.01244_rkx * tcel)
    !clw = max(3.0e-4_rkx,0.18_rkx * min(1.0_rkx,(tcel+53.15_rkx)/60.0_rkx))
  end function clwfromt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
