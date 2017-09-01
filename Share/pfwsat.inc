!DIR$ ATTRIBUTES FORCEINLINE :: pfwsat
  pure real(rkx) function pfwsat(t,p,e) result(ws) ! In kg/kg
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) :: es
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    ws = ep2 * (es / (p - es))
  end function pfwsat
