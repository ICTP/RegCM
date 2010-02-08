      module tmpsav
      use regcm_param
      implicit none
#ifdef MPP1
      real(8) , dimension(ix,kx*4+2,jxp) :: sav0
      real(8) , dimension(ix,kx+nnsg+5,jxp) :: sav0a
      real(8) , dimension(ix,kx+1,jxp) :: sav0b
      real(8) , dimension(ix,kx*2,jxp) :: sav0c
      real(8) , dimension(ix,kx,jxp) :: sav0s
      real(8) , dimension(ix,kx*4+2,mjx) :: sav_0
      real(8) , dimension(ix,kx+nnsg+5,mjx) :: sav_0a
      real(8) , dimension(ix,kx+1,mjx) :: sav_0b
      real(8) , dimension(ix,kx*2,mjx) :: sav_0c
      real(8) , dimension(ix,kx,mjx) :: sav_0s
      real(8) , dimension(ilx,kx*4+(kx+1)*(kx+2),jxp) :: sav1
      real(8) , dimension(ilx,kx*4+(kx+1)*(kx+2),mjx) :: sav_1
      real(8) , dimension(ilx,nnsg*4+4,jxp) :: sav2
      real(8) , dimension(ilx,nnsg*5+1,jxp) :: sav2a
      real(8) , dimension(ilx,nnsg*4+4,mjx) :: sav_2
      real(8) , dimension(ilx,nnsg*5+1,mjx) :: sav_2a
      real(8) , dimension(ix,ntr*(kx*4+1),jxp) :: sav4
      real(8) , dimension(ilx,7,jxp) :: sav4a
      real(8) , dimension(ix,ntr*(kx*4+1),mjx) :: sav_4
      real(8) , dimension(ilx,7,mjx) :: sav_4a
#endif
      end module tmpsav
