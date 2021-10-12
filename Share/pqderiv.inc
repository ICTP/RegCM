!DIR$ ATTRIBUTES FORCEINLINE :: pfqsdt
  pure real(rkx) function pfqsdt(t,p,e,dedt) result(qsdt)
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) , intent(in) , optional :: dedt ! derivative of e in dt (Pa/K)
    real(rkx) :: es , esdt , vp , vp1 , vp2
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    if ( present(dedt) ) then
      esdt = dedt
    else
      esdt = pfesdt(t)
    end if
    vp  = 1.0_rkx / (p - 0.378_rkx*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    qsdt = esdt * vp2 * p ! 1 / K
  end function pfqsdt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
