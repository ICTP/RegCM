  !
  ! Computes derivative in temperature of saturation pressure
  ! Tetens formula
  !
!DIR$ ATTRIBUTES FORCEINLINE :: pfesdt
  pure elemental real(rkx) function pfesdt(t) result(esdt)
    implicit none
    real(rkx), intent(in)  :: t     ! Temperature (K)

    real(rkx) :: tc , es

    tc = t - tzero
    es = pfesat(t)

    if ( tc > 0.0_rkx ) then
      esdt = es * c3les/((tc-c4les)**2)
    else
      esdt = es * c3ies/((tc-c4ies)**2)
    end if
  end function pfesdt

!DIR$ ATTRIBUTES FORCEINLINE :: pfqsdt
  pure real(rkx) function pfqsdt(t,p) result(qsdt)
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) :: es , qs

    es = pfesat(t)
    qs = ep2 * es / (p - (0.378_rkx * es))
    if ( t < tzero ) then
      qsdt = qs * (c3ies * (tzero-c4ies) * (d_one/(t-c4ies))**d_two)
    else
      qsdt = qs * (c3les * (tzero-c4les) * (d_one/(t-c4les))**d_two)
    end if
  end function pfqsdt

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
