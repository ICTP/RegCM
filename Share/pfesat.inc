  ! Computes saturation pressurre
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
!DIR$ ATTRIBUTES FORCEINLINE :: pfesat
  pure elemental real(rkx) function pfesat(t,p) result(es)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t , p ! Temperature (K) , Pressure (Pa)

    real(rk8) :: td , t_limit , esat
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8) , parameter :: a0 =  0.611213476e+03_rk8
    real(rk8) , parameter :: a1 =  0.444007856e+02_rk8
    real(rk8) , parameter :: a2 =  0.143064234e+01_rk8
    real(rk8) , parameter :: a3 =  0.264461437e-01_rk8
    real(rk8) , parameter :: a4 =  0.305903558e-03_rk8
    real(rk8) , parameter :: a5 =  0.196237241e-05_rk8
    real(rk8) , parameter :: a6 =  0.892344772e-08_rk8
    real(rk8) , parameter :: a7 = -0.373208410e-10_rk8
    real(rk8) , parameter :: a8 =  0.209339997e-13_rk8
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8) , parameter :: c0 =  0.611123516e+03_rk8
    real(rk8) , parameter :: c1 =  0.503109514e+02_rk8
    real(rk8) , parameter :: c2 =  0.188369801e+01_rk8
    real(rk8) , parameter :: c3 =  0.420547422e-01_rk8
    real(rk8) , parameter :: c4 =  0.614396778e-03_rk8
    real(rk8) , parameter :: c5 =  0.602780717e-05_rk8
    real(rk8) , parameter :: c6 =  0.387940929e-07_rk8
    real(rk8) , parameter :: c7 =  0.149436277e-09_rk8
    real(rk8) , parameter :: c8 =  0.262655803e-12_rk8

    t_limit = t - tzero
    if ( t_limit > 100.0_rk8 ) t_limit = 100.0_rk8
    if ( t_limit < -75.0_rk8 ) t_limit = -75.0_rk8
    td = t_limit
    if ( td >= 0.0_rk8 ) then
      esat = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    else
      esat = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    end if
    es = real(min(esat,0.15_rk8*p),rkx) ! pa
  end function pfesat

!DIR$ ATTRIBUTES FORCEINLINE :: pfesat_water
  pure elemental real(rkx) function pfesat_water(t,p) result(es)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t , p ! Temperature (K) , Pressure (Pa)

    real(rk8) :: td , esat
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8) , parameter :: a0 =  0.611213476e+03_rk8
    real(rk8) , parameter :: a1 =  0.444007856e+02_rk8
    real(rk8) , parameter :: a2 =  0.143064234e+01_rk8
    real(rk8) , parameter :: a3 =  0.264461437e-01_rk8
    real(rk8) , parameter :: a4 =  0.305903558e-03_rk8
    real(rk8) , parameter :: a5 =  0.196237241e-05_rk8
    real(rk8) , parameter :: a6 =  0.892344772e-08_rk8
    real(rk8) , parameter :: a7 = -0.373208410e-10_rk8
    real(rk8) , parameter :: a8 =  0.209339997e-13_rk8

    td = min(t - tzero,100.0_rk8)
    esat = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    es = real(min(esat,0.15_rk8*p),rkx) ! pa
  end function pfesat_water

!DIR$ ATTRIBUTES FORCEINLINE :: pfesat_ice
  pure elemental real(rkx) function pfesat_ice(t,p) result(es)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t , p ! Temperature (K) , Pressure (Pa)

    real(rk8) :: td , esat
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8) , parameter :: c0 =  0.611123516e+03_rk8
    real(rk8) , parameter :: c1 =  0.503109514e+02_rk8
    real(rk8) , parameter :: c2 =  0.188369801e+01_rk8
    real(rk8) , parameter :: c3 =  0.420547422e-01_rk8
    real(rk8) , parameter :: c4 =  0.614396778e-03_rk8
    real(rk8) , parameter :: c5 =  0.602780717e-05_rk8
    real(rk8) , parameter :: c6 =  0.387940929e-07_rk8
    real(rk8) , parameter :: c7 =  0.149436277e-09_rk8
    real(rk8) , parameter :: c8 =  0.262655803e-12_rk8

    td = max(t - tzero,-75.0_rk8)
    esat = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    es = real(min(esat,0.15_rk8*p),rkx) ! pa
  end function pfesat_ice

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
