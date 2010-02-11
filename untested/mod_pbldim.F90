      module mod_pbldim

      use mod_regcm_param

      implicit none
!
! COMMON /PBLDIM/
!
      real(8) , dimension(ix,kxp1) :: zq
#ifdef MPP1
      real(8) , dimension(ix,kx,jxp) :: dzq , thvx , thx3d , za
      real(8) , dimension(ix,jxp) :: rhox2d
#else
      real(8) , dimension(ix,kx,jx) :: dzq , thvx , thx3d , za
      real(8) , dimension(ix,jx) :: rhox2d
#endif
      end module mod_pbldim
