module mod_clm_qsat
  !
  ! Computes saturation specific humidity and the change in saturation
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
    real(rk8), intent(in)  :: t     ! Temperature (K)
    real(rk8), intent(in)  :: p     ! surface atmospheric pressure (pa)
    real(rk8), intent(out) :: es    ! vapor pressure (pa)
    real(rk8), intent(out) :: esdt  ! d(es)/d(T)
    real(rk8), intent(out) :: qs    ! humidity (kg/kg)
    real(rk8), intent(out) :: qsdt  ! d(qs)/d(T)

    real(rk8) :: t_limit
    real(rk8) :: td , vp , vp1 , vp2
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8), parameter :: a0 =  6.11213476_rk8
    real(rk8), parameter :: a1 =  0.444007856_rk8
    real(rk8), parameter :: a2 =  0.143064234e-01_rk8
    real(rk8), parameter :: a3 =  0.264461437e-03_rk8
    real(rk8), parameter :: a4 =  0.305903558e-05_rk8
    real(rk8), parameter :: a5 =  0.196237241e-07_rk8
    real(rk8), parameter :: a6 =  0.892344772e-10_rk8
    real(rk8), parameter :: a7 = -0.373208410e-12_rk8
    real(rk8), parameter :: a8 =  0.209339997e-15_rk8
    !
    ! For derivative:water vapor
    !
    real(rk8), parameter :: b0 =  0.444017302_rk8
    real(rk8), parameter :: b1 =  0.286064092e-01_rk8
    real(rk8), parameter :: b2 =  0.794683137e-03_rk8
    real(rk8), parameter :: b3 =  0.121211669e-04_rk8
    real(rk8), parameter :: b4 =  0.103354611e-06_rk8
    real(rk8), parameter :: b5 =  0.404125005e-09_rk8
    real(rk8), parameter :: b6 = -0.788037859e-12_rk8
    real(rk8), parameter :: b7 = -0.114596802e-13_rk8
    real(rk8), parameter :: b8 =  0.381294516e-16_rk8
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8), parameter :: c0 =  6.11123516_rk8
    real(rk8), parameter :: c1 =  0.503109514_rk8
    real(rk8), parameter :: c2 =  0.188369801e-01_rk8
    real(rk8), parameter :: c3 =  0.420547422e-03_rk8
    real(rk8), parameter :: c4 =  0.614396778e-05_rk8
    real(rk8), parameter :: c5 =  0.602780717e-07_rk8
    real(rk8), parameter :: c6 =  0.387940929e-09_rk8
    real(rk8), parameter :: c7 =  0.149436277e-11_rk8
    real(rk8), parameter :: c8 =  0.262655803e-14_rk8
    !
    ! For derivative:ice
    !
    real(rk8), parameter :: d0 =  0.503277922_rk8
    real(rk8), parameter :: d1 =  0.377289173e-01_rk8
    real(rk8), parameter :: d2 =  0.126801703e-02_rk8
    real(rk8), parameter :: d3 =  0.249468427e-04_rk8
    real(rk8), parameter :: d4 =  0.313703411e-06_rk8
    real(rk8), parameter :: d5 =  0.257180651e-08_rk8
    real(rk8), parameter :: d6 =  0.133268878e-10_rk8
    real(rk8), parameter :: d7 =  0.394116744e-13_rk8
    real(rk8), parameter :: d8 =  0.498070196e-16_rk8

    t_limit = t - tzero
    if ( t_limit > 100.0_rk8 ) t_limit = 100.0_rk8
    if ( t_limit < -75.0_rk8 ) t_limit = -75.0_rk8
    td = t_limit

    if ( td >= 0.0_rk8 ) then
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

    es   = es   * 100._rk8 ! pa
    esdt = esdt * 100._rk8 ! pa/K

    vp  = 1.0_rk8 / (p - 0.378_rk8*es)
    vp1 = 0.622_rk8 * vp
    vp2 = vp1 * vp

    qs   = es   * vp1     ! kg/kg
    qsdt = esdt * vp2 * p ! 1 / K

  end subroutine qsat

end module mod_clm_qsat
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
