      module mod_mainchem

      use mod_regcm_param

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

      real(4) , dimension(jxm2,ixm2) :: fchem

#ifdef MPP1
      real(8) , dimension(ix,12,ntr,jxp) :: src0
      real(8) , dimension(ix,12,ntr,jx) :: src_0
#endif

      end module mod_mainchem
