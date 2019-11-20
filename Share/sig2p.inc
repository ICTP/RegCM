!DIR$ ATTRIBUTES FORCEINLINE :: sig2p
  pure elemental real(rkx) function sig2p(ps,sigma,ptop) result(p)
    implicit none
    real(rkx) , intent(in) :: ps , sigma , ptop
    !!!!!!!!!! Assume input is in Pa !!!!!!!!!!
    p = (sigma * ( ps - ptop ) + ptop) ! Pressure in Pa
  end function sig2p

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
