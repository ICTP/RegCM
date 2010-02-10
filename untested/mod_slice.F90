      module mod_slice

      use mod_regcm_param

      implicit none
!
! COMMON /CBSLICE/
!
#ifdef MPP1
      real(8) , dimension(ix,kx,-1:jxp+2,ntr) :: chib3d
      real(8) , dimension(ix,kx,jxp) :: pb3d , qsb3d , rhb3d , rhob3d , &
                                      & ubx3d , vbx3d
      real(8) , dimension(ix,kx,-1:jxp+2) :: qcb3d , qvb3d , tb3d ,     &
           & ubd3d , vbd3d
#else
      real(8) , dimension(ix,kx,jx,ntr) :: chib3d
      real(8) , dimension(ix,kx,jx) :: pb3d , qcb3d , qsb3d , qvb3d ,   &
                                     & rhb3d , rhob3d , tb3d , ubd3d ,  &
                                     & ubx3d , vbd3d , vbx3d
#endif

      end module mod_slice
