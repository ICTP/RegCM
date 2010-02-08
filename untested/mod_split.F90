      module split

      use regcm_param

      implicit none
!
! COMMON /CDIAG/
!
      real(8) , dimension(kx,kx) :: a , a1 , a2 , a3 , a4 , d1 , d2 ,   &
                                  & e1 , e2 , e3 , g1 , g2 , g3 , s1 ,  &
                                  & s2 , w1 , w2 , x1
      integer , dimension(kx) :: iw2
      real(8) , dimension(kxp1) :: tbarf , thetaf
      real(8) , dimension(kx) :: thetah , tweigh
      real(8) , dimension(kxp1,kx) :: w3
!
! COMMON /CVERT/
!
      real(8) :: alpha1 , alpha2 , pd , ps , pt , r
      real(8) , dimension(kx) :: cpfac , dsigma , hbar , hweigh , tbarh
      real(8) , dimension(kx,kxp1) :: hydroc , varpa1
      real(8) , dimension(kx,kx) :: hydror , hydros , tau , zmatx ,     &
                                  & zmatxr
      real(8) , dimension(kxp1) :: sigmah
      real(8) , dimension(kxp1,kxp1) :: varpa2
!
! COMMON /DPASS/
!
      real(8) , dimension(kx,nsplit) :: am
      real(8) , dimension(nsplit) :: an
#ifdef MPP1
      real(8) , dimension(ix,0:jxp+1,nsplit) :: dstor , hstor
#else
      real(8) , dimension(ix,jx,nsplit) :: dstor , hstor
#endif
      integer , dimension(nsplit) :: m

#ifdef MPP1
!
! COMMON /SPLITIO/
!
      real(8) , dimension(ix,mjx,nsplit) :: dstor_io , hstor_io
#endif

      end module split
