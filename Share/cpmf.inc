  ! Heat capacity for moist air at constant pressure
  pure real(rkx) function cpmf(qv)
    implicit none
    real(rkx) , intent(in) :: qv
    !cpmf = cpd * (d_one + 0.856_rkx * qv)
    !cpmf = cpd * (d_one + 0.84_rkx * qv)
    !cpmf = cpd + 1.820_rkx * qv
    cpmf = cpd * (d_one + 0.80_rkx * qv)
  end function cpmf

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
