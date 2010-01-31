      module dust

      use regcm_param

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: nsoil = 152 , nats = 12 , mode = 3 ,       &
                           & jsoilm = 1 , jfs = 0 , ust = 1
!
! COMMON /DUST/
!
#ifdef MPP1
      real(8) , dimension(ix,nats,jxp) :: clay2row2 , sand2row2 ,       &
           & silt2row2
      real(8) , dimension(ix,jxp) :: clayrow2 , dustsotex , sandrow2
      real(8) , dimension(ix,jxp,nsoil) :: srel2d
#else
      real(8) , dimension(ix,nats,jx) :: clay2row2 , sand2row2 ,        &
           & silt2row2
      real(8) , dimension(ix,jx) :: clayrow2 , dustsotex , sandrow2
      real(8) , dimension(ix,jx,nsoil) :: srel2d
#endif
      real(8) , dimension(nsoil) :: dp

      end module dust
