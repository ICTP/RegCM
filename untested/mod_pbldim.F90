      module pbldim

      use regcm_param

      implicit none
!
! COMMON /PBLDIM/
!
#ifdef MPP1
      real(8) , dimension(ix,kx,jxp) :: dzq , thvx , thx3d , za
      real(8) , dimension(ix,jxp) :: rhox2d
      real(8) , dimension(ix,kxp1) :: zq
#else
      real(8) , dimension(ix,kx,jx) :: dzq , thvx , thx3d , za
      real(8) , dimension(ix,jx) :: rhox2d
      real(8) , dimension(ix,kxp1) :: zq
#endif
      end module pbldim
