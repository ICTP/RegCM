      module rad

      use regcm_param

      implicit none
!
! COMMON /RADIATION/
!
#ifdef MPP1
      real(8) , dimension(ilx,kx) :: cldfra , cldlwc
      real(8) , dimension(ilx,kx,jxp) :: heatrt
      real(8) , dimension(ilx,kxp1,jxp) :: o3prof
#else
      real(8) , dimension(ilx,kx) :: cldfra , cldlwc
      real(8) , dimension(ilx,kx,jlx) :: heatrt
      real(8) , dimension(ilx,kxp1,jlx) :: o3prof
#endif
      end module rad
