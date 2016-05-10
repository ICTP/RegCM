  pure real(rkx) function sig2p(ps,sigma,ptop) result(p)
    implicit none
    real(rkx) , intent(in) :: ps , sigma , ptop
    !!!!!!!!!! Assume input is in cbar !!!!!!!!!!
    p = (sigma * ( ps - ptop ) + ptop) * d_1000 ! Pressure in Pa
    !!!!!!!!!! Assume input is in cbar !!!!!!!!!!
  end function sig2p

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
