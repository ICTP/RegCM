  ! Computes saturation pressurre
  ! Reference: Tetens, O., 1930:
  ! Uber einige meteorologische Begriffe. Zeitschrift fur Geophysik, Vol. 6:297.
  !
!DIR$ ATTRIBUTES FORCEINLINE :: pfesat
  pure elemental real(rkx) function pfesat(t) result(es)
    implicit none
    real(rkx) , intent(in) :: t     ! Temperature (K)

    if ( t > wattp ) then
      es = c1es*exp(c3les*((t-wattp)/(t-c4les)))
    else
      es = c1es*exp(c3ies*((t-wattp)/(t-c4ies)))
    end if
  end function pfesat

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
