  ! Heat capacity for moist air at constant pressure
  pure real(rk8) function cpmf(qv)
    implicit none
    real(rk8) , intent(in) :: qv
    cpmf = cpd * (d_one + 0.856D0 * qv)
    !cpmf = cpd * (d_one + 0.84D0 * qv)
    !cpmf = cpd + 1.820D0 * qv
    !cpmf = cpd * (d_one + 0.80D0 * qv)
  end function cpmf

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
