      module mainchem

      use regcm_param

      implicit none
!
! COMMON /TRAC/
!
#ifdef MPP1
      real(8) , dimension(ix,jxp,12,ntr) :: chemsrc
      real(8) , dimension(ix,kx,-1:jxp+2,ntr) :: chia , chib
      real(8) , dimension(ix,jxp,ntr) :: srclp2
#else
      real(8) , dimension(ix,jx,12,ntr) :: chemsrc
      real(8) , dimension(ix,kx,jx,ntr) :: chia , chib
      real(8) , dimension(ix,jx,ntr) :: srclp2
#endif
!
! COMMON /TRACFLUX/
!
#ifdef MPP1
      real(8) , dimension(ix,jxp,ntr) :: ddsfc , dtrace , wdcvc ,       &
                                       & wdlsc , wxaq , wxsg
#else
      real(8) , dimension(ix,jx,ntr) :: ddsfc , dtrace , wdcvc , wdlsc ,&
                                      & wxaq , wxsg
#endif

#ifdef MPP1
!
! COMMON /IO/
!
      real(8) , dimension(ix,mjx,12,ntr) :: chemsrc_io
      real(8) , dimension(ix,mjx,ntr) :: ddsfc_io , dtrace_io ,         &
           & wdcvc_io , wdlsc_io
      real(8) , dimension(ilx,jxbb) :: pptc_io , pptnc_io , prca2d_io , &
                                     & prnca2d_io
#endif

      integer :: nrcchem
#ifdef MPP1
      real(4) , dimension(mjx-2,ix-2) :: fchem
      real(8) , dimension(ix,12,ntr,jxp) :: src0
      real(8) , dimension(ix,12,ntr,mjx) :: src_0
#else
      real(4) , dimension(jx-2,ix-2) :: fchem
#endif

      end module mainchem
