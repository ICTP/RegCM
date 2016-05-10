module mod_clm_qsat
  !
  ! Computes saturation mixing ratio and the change in saturation
  !
  implicit none

  private

  save

  public :: qsat

  contains
  !
  ! Computes saturation mixing ratio and the change in saturation
  ! mixing ratio with respect to temperature.
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
  subroutine qsat(t,p,es,esdt,qs,qsdt)
    use mod_realkinds
    use mod_constants
    implicit none
    real(rkx), intent(in)  :: t     ! Temperature (K)
    real(rkx), intent(in)  :: p     ! surface atmospheric pressure (pa)
    real(rkx), intent(out) :: es    ! vapor pressure (pa)
    real(rkx), intent(out) :: esdt  ! d(es)/d(T)
    real(rkx), intent(out) :: qs    ! humidity (kg/kg)
    real(rkx), intent(out) :: qsdt  ! d(qs)/d(T)

    real(rkx) :: t_limit
    real(rkx) :: td , vp , vp1 , vp2
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rkx), parameter :: a0 =  6.11213476_rkx
    real(rkx), parameter :: a1 =  0.444007856_rkx
    real(rkx), parameter :: a2 =  0.143064234e-01_rkx
    real(rkx), parameter :: a3 =  0.264461437e-03_rkx
    real(rkx), parameter :: a4 =  0.305903558e-05_rkx
    real(rkx), parameter :: a5 =  0.196237241e-07_rkx
    real(rkx), parameter :: a6 =  0.892344772e-10_rkx
    real(rkx), parameter :: a7 = -0.373208410e-12_rkx
    real(rkx), parameter :: a8 =  0.209339997e-15_rkx
    !
    ! For derivative:water vapor
    !
    real(rkx), parameter :: b0 =  0.444017302_rkx
    real(rkx), parameter :: b1 =  0.286064092e-01_rkx
    real(rkx), parameter :: b2 =  0.794683137e-03_rkx
    real(rkx), parameter :: b3 =  0.121211669e-04_rkx
    real(rkx), parameter :: b4 =  0.103354611e-06_rkx
    real(rkx), parameter :: b5 =  0.404125005e-09_rkx
    real(rkx), parameter :: b6 = -0.788037859e-12_rkx
    real(rkx), parameter :: b7 = -0.114596802e-13_rkx
    real(rkx), parameter :: b8 =  0.381294516e-16_rkx
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rkx), parameter :: c0 =  6.11123516_rkx
    real(rkx), parameter :: c1 =  0.503109514_rkx
    real(rkx), parameter :: c2 =  0.188369801e-01_rkx
    real(rkx), parameter :: c3 =  0.420547422e-03_rkx
    real(rkx), parameter :: c4 =  0.614396778e-05_rkx
    real(rkx), parameter :: c5 =  0.602780717e-07_rkx
    real(rkx), parameter :: c6 =  0.387940929e-09_rkx
    real(rkx), parameter :: c7 =  0.149436277e-11_rkx
    real(rkx), parameter :: c8 =  0.262655803e-14_rkx
    !
    ! For derivative:ice
    !
    real(rkx), parameter :: d0 =  0.503277922_rkx
    real(rkx), parameter :: d1 =  0.377289173e-01_rkx
    real(rkx), parameter :: d2 =  0.126801703e-02_rkx
    real(rkx), parameter :: d3 =  0.249468427e-04_rkx
    real(rkx), parameter :: d4 =  0.313703411e-06_rkx
    real(rkx), parameter :: d5 =  0.257180651e-08_rkx
    real(rkx), parameter :: d6 =  0.133268878e-10_rkx
    real(rkx), parameter :: d7 =  0.394116744e-13_rkx
    real(rkx), parameter :: d8 =  0.498070196e-16_rkx

    t_limit = t - tzero
    if ( t_limit > 100.0_rkx ) t_limit = 100.0_rkx
    if ( t_limit < -75.0_rkx ) t_limit = -75.0_rkx
    td = t_limit

    if ( td >= 0.0_rkx ) then
      es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
           + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
      esdt = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
           + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
      es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
           + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
      esdt = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
           + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    end if

    es   = es   * 100._rkx ! pa
    esdt = esdt * 100._rkx ! pa/K

    vp  = 1.0_rkx / (p - 0.378_rkx*es)
    vp1 = 0.622_rkx * vp
    vp2 = vp1 * vp

    qs   = es   * vp1     ! kg/kg
    qsdt = esdt * vp2 * p ! 1 / K

  end subroutine qsat

end module mod_clm_qsat
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
