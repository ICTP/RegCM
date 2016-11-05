  pure real(rkx) function pfwsat(t,p,e) result(ws)
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) :: es , vp , vp1 , vp2 , qs
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t)
    end if
    ! Bolton 1980
    vp  = 1.0_rkx / (p - 0.378_rkx*es)
    vp1 = ep2 * vp
    vp2 = vp1 * vp
    qs = es * vp1  ! kg/kg
    ws = max(qs/(d_one - qs),minqq)
  end function pfwsat
