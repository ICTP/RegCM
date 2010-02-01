      module outrad
      use regcm_param
      implicit none

      integer , parameter :: nrad2d = 21 , nrad3d = 5
#ifdef MPP1
      real(4) , dimension(jxp,ix-2,nrad2d) :: frad2d
      real(4) , dimension(jxp,ix-2,kx,nrad3d) :: frad3d
#else
      real(4) , dimension(jx-2,ix-2,nrad2d) :: frad2d
      real(4) , dimension(jx-2,ix-2,kx,nrad3d) :: frad3d
#endif

      end module outrad
