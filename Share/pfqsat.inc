!DIR$ ATTRIBUTES FORCEINLINE :: pfqsat
  pure elemental real(rkx) function pfqsat(t,p,e) result(qs) ! In kg/kg
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
    ws = ep2 * (es / (p - es))
    qs = ws / ( d_one + ws )
    ! Bolton 1980
    ! qs = ep2 * es / (p - (0.378_rkx * es))
  end function pfqsat
