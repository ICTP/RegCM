      module bxq

      use regcm_param

      implicit none
#ifdef MPP1
!
! COMMON /bxq_aaa/
!
      real(8) , dimension(ix,jxp,nsplit) :: ddsum
      real(8) , dimension(ix,jxp,nsplit,3) :: deld
      real(8) , dimension(ix,0:jxp,nsplit,3) :: delh
      real(8) , dimension(ix,0:jxp,nsplit) :: dhsum
      real(8) , dimension(ix,jxp) :: psdot
      real(8) , dimension(ix,jxp,3) :: work
!
! COMMON /bxq_tmp/
!
      real(8) , dimension(ix,jxp+1) :: uu , vv
!
! COMMON /bxq_tmq/
!
      real(8) , dimension(ix,kx,jxp+1) :: uuu , vvv
#else
!
! COMMON /bxq_aaa/
!
      real(8) , dimension(ix,jx,nsplit) :: ddsum , dhsum
      real(8) , dimension(ix,jx,nsplit,3) :: deld , delh
      real(8) , dimension(ix,jx) :: psdot
      real(8) , dimension(ix,jx,3) :: work
!
! COMMON /bxq_tmp/
!
      real(8) , dimension(ix,jx) :: uu , vv
!
! COMMON /bxq_tmq/
!
      real(8) , dimension(ix,kx,jx) :: uuu , vvv
#endif

      end module bxq
