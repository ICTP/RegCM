  ! Temperature dependency for cloud water content
  ! in g/m3 (Lemus et al., 1997)
  ! NOTE : THIS IS IN-CLOUD VARIABLE.
  pure real(rkx) function clwfromt(t) result(clw)
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
  end function clwfromt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
