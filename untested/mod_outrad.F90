      module mod_outrad
      use mod_regcm_param
      implicit none

      integer , parameter :: nrad2d = 21 , nrad3d = 5
      integer :: nrcrad
#ifdef MPP1
      real(4) , dimension(jxp,ixm2,nrad2d) :: frad2d
      real(4) , dimension(jxp,ixm2,kx,nrad3d) :: frad3d
      real(4) , dimension(jxm2,ixm2,nrad2d) :: frad2d_io
      real(4) , dimension(jxm2,ixm2,kx,nrad3d) :: frad3d_io
#else
      real(4) , dimension(jxm2,ixm2,nrad2d) :: frad2d
      real(4) , dimension(jxm2,ixm2,kx,nrad3d) :: frad3d
#endif

      end module mod_outrad
