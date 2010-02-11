      module mod_tmpsav
      use mod_regcm_param
      implicit none

#ifdef MPP1
      real(8) , dimension(ix,kx*4+2,jxp) :: sav0
      real(8) , dimension(ix,kx+nnsg+5,jxp) :: sav0a
      real(8) , dimension(ix,kxp1,jxp) :: sav0b
      real(8) , dimension(ix,kx*2,jxp) :: sav0c
      real(8) , dimension(ix,kx,jxp) :: sav0s
      real(8) , dimension(ix,kx*4+2,jx) :: sav_0
      real(8) , dimension(ix,kx+nnsg+5,jx) :: sav_0a
      real(8) , dimension(ix,kxp1,jx) :: sav_0b
      real(8) , dimension(ix,kx*2,jx) :: sav_0c
      real(8) , dimension(ix,kx,jx) :: sav_0s
      real(8) , dimension(ix,nsplit*2,jxp) :: sav0d
      real(8) , dimension(ix,nsplit*2,jx) :: sav_0d
      real(8) , dimension(ixm1,kx*4+(kxp1)*(kxp2),jxp) :: sav1
      real(8) , dimension(ixm1,kx*4+(kxp1)*(kxp2),jx) :: sav_1
      real(8) , dimension(ixm1,nnsg*4+4,jxp) :: sav2
      real(8) , dimension(ixm1,nnsg*5+1,jxp) :: sav2a
      real(8) , dimension(ixm1,nnsg*4+4,jx) :: sav_2
      real(8) , dimension(ixm1,nnsg*5+1,jx) :: sav_2a
      real(8) , dimension(ix,ntr*(kx*4+1),jxp) :: sav4
      real(8) , dimension(ixm1,7,jxp) :: sav4a
      real(8) , dimension(ix,ntr*(kx*4+1),jx) :: sav_4
      real(8) , dimension(ixm1,7,jx) :: sav_4a
      real(8) , dimension(kx,8,jxp) :: sav6
      real(8) , dimension(kx,8,jx) :: sav_6
#endif

      end module mod_tmpsav
