      module mod_rad

      use mod_regcm_param

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

#ifdef MPP1
!
! COMMON /RAD1IO/
!
      real(8) , dimension(ilx,kx,mjx-1) :: heatrt_io
      real(8) , dimension(ilx,kxp1,mjx-1) :: o3prof_io
#endif

      end module mod_rad
