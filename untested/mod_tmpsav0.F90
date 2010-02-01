      module tmpsav0
      use regcm_param
      implicit none
#ifdef MPP1
      real(8) , dimension(ix,kx*4+2,jxp) :: sav0
      real(8) , dimension(ix,kx,jxp) :: sav0s
      real(8) , dimension(ix,kx*4+2,mjx) :: sav_0
      real(8) , dimension(ix,kx,mjx) :: sav_0s
#endif
      end module tmpsav0
