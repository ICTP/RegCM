      module mod_rad

      use mod_regcm_param

      implicit none

#ifdef MPP1
      real(8) , dimension(ixm1,kx) :: cldfra , cldlwc
      real(8) , dimension(ixm1,kx,jxp) :: heatrt
      real(8) , dimension(ixm1,kxp1,jxp) :: o3prof
#else
      real(8) , dimension(ixm1,kx) :: cldfra , cldlwc
      real(8) , dimension(ixm1,kx,jxm1) :: heatrt
      real(8) , dimension(ixm1,kxp1,jxm1) :: o3prof
#endif

      end module mod_rad
