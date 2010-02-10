      module mod_rad

      use mod_regcm_param

      implicit none
!
! COMMON /RADIATION/
!
#ifdef MPP1
      real(8) , dimension(ixm1,kx) :: cldfra , cldlwc
      real(8) , dimension(ixm1,kx,jxp) :: heatrt
      real(8) , dimension(ixm1,kxp1,jxp) :: o3prof
#else
      real(8) , dimension(ixm1,kx) :: cldfra , cldlwc
      real(8) , dimension(ixm1,kx,jxm1) :: heatrt
      real(8) , dimension(ixm1,kxp1,jxm1) :: o3prof
#endif

#ifdef MPP1
!
! COMMON /RAD1IO/
!
      real(8) , dimension(ixm1,kx,jxm1) :: heatrt_io
      real(8) , dimension(ixm1,kxp1,jxm1) :: o3prof_io
#endif

      end module mod_rad
