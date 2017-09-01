!DIR$ ATTRIBUTES FORCEINLINE :: pfqsat
  pure real(rkx) function pfqsat(t,p,e) result(qs) ! In kg/kg
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) :: es , ws
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    ws = ep2 * es / (p - es)
    qs = ws / ( 1.0_rkx + ws )
  end function pfqsat
