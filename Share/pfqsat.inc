  pure real(rk8) function pfqsat(t,p,e) result(qs)
    implicit none
    real(rk8) , intent(in) :: t             ! Temperature (K)
    real(rk8) , intent(in) :: p             ! Pressure (Pa)
    real(rk8) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rk8) :: es , vp , vp1 , vp2
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    ! Bolton 1980
    vp  = 1.0D0 / (p - 0.378D0*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    qs = max(es * vp1, minqq)  ! kg/kg
  end function pfqsat
