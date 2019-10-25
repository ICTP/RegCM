!DIR$ ATTRIBUTES FORCEINLINE :: pfwsat
  pure elemental real(rkx) function pfwsat(t,p,e) result(ws) ! In kg/kg
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) :: es ! , qs
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    ws = ep2 * (es / (p - es))
    ! Bolton 1980
    ! qs = ep2 * es / (p - (0.378_rkx * es))
    ! ws = qs * ( d_one - qs)
  end function pfwsat
