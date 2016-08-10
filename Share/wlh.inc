  ! Latent heat for condensation/sublimation of water in J/kg
  pure real(rkx) function wlh(t)
    implicit none
    real(rkx) , intent(in) :: t
    real(rk8) :: tc , lheat
    tc = t - tzero
    if ( tc > d_zero ) then
      lheat = 2500.79_rk8 - 2.36418_rk8*tc + &
            1.58927e-3_rk8*tc*tc - 6.14342e-5_rk8*tc*tc*tc
    else
      lheat = 2834.1_rk8 - 0.29_rk8*tc - 0.004_rk8*tc*tc
    end if
    wlh = real(lheat*1.0e3_rk8,rkx)
  end function wlh

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
