  ! Latent heat for condensation/sublimation of water in J/kg
  pure real(rk8) function wlh(t)
    implicit none
    real(rk8) , intent(in) :: t
    real(rk8) :: tc
    tc = t - tzero
    if ( tc > d_zero ) then
      wlh = 2500.79D0 - 2.36418D0*tc + 1.58927D-3*tc*tc - 6.14342D-5*tc*tc*tc
    else
      wlh = 2834.1D0 - 0.29D0*tc - 0.004D0*tc*tc
    end if
    wlh = wlh*1.0D3
  end function wlh

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
