  ! Computes saturation pressurre
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
  pure real(rk8) function pfesat(t) result(es)
    implicit none
    real(rk8) , intent(in) :: t     ! Temperature (K)

    real(rk8) :: td , t_limit
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8) , parameter :: a0 =  6.11213476D0
    real(rk8) , parameter :: a1 =  0.444007856D0
    real(rk8) , parameter :: a2 =  0.143064234D-01
    real(rk8) , parameter :: a3 =  0.264461437D-03
    real(rk8) , parameter :: a4 =  0.305903558D-05
    real(rk8) , parameter :: a5 =  0.196237241D-07
    real(rk8) , parameter :: a6 =  0.892344772D-10
    real(rk8) , parameter :: a7 = -0.373208410D-12
    real(rk8) , parameter :: a8 =  0.209339997D-15
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8) , parameter :: c0 =  6.11123516D0
    real(rk8) , parameter :: c1 =  0.503109514D0
    real(rk8) , parameter :: c2 =  0.188369801D-01
    real(rk8) , parameter :: c3 =  0.420547422D-03
    real(rk8) , parameter :: c4 =  0.614396778D-05
    real(rk8) , parameter :: c5 =  0.602780717D-07
    real(rk8) , parameter :: c6 =  0.387940929D-09
    real(rk8) , parameter :: c7 =  0.149436277D-11
    real(rk8) , parameter :: c8 =  0.262655803D-14

    t_limit = t - tzero
    if ( t_limit > 100.0D0 ) t_limit = 100.0D0
    if ( t_limit < -75.0D0 ) t_limit = -75.0D0
    td = t_limit
    if ( td >= 0.0D0 ) then
      es = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    else
      es = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    end if
    es   = es * d_100 ! pa
  end function pfesat

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
