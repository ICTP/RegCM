  ! Latent heat for condensation/sublimation of water in J/kg
  pure real(rkx) function wlh(t)
    implicit none
    real(rkx) , intent(in) :: t
    real(rkx) :: tc
    tc = t - tzero
    if ( tc > d_zero ) then
      wlh = 2500.79_rkx - 2.36418_rkx*tc + &
            1.58927e-3_rkx*tc*tc - 6.14342e-5_rkx*tc*tc*tc
    else
      wlh = 2834.1_rkx - 0.29_rkx*tc - 0.004_rkx*tc*tc
    end if
    wlh = wlh*1.0e3_rkx
  end function wlh

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
