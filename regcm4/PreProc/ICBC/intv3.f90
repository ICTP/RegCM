      subroutine intv3(fsccm,fccm,psrccm,sccm,ptop,ni,nj,kccm)
      use mod_constants , only : rgas , gti , lrate
      implicit none
!
! PARAMETER definitions
!
      real , parameter :: rgas2 = rgas/2.
      real , parameter :: b1 = -gti/lrate
!
! Dummy arguments
!
      integer :: kccm , ni , nj
      real :: ptop
      real , dimension(ni,nj,kccm) :: fccm
      real , dimension(ni,nj) :: fsccm , psrccm
      real , dimension(kccm) :: sccm
      intent (in) fccm , kccm , ni , nj , psrccm , ptop , sccm
      intent (out) fsccm
!
! Local variables
!
      real :: a1 , rc , rc1 , sc
      integer :: i , j , k , k1 , k1p
!
!**   INTV3 IS FOR VERTICAL INTERPOLATION OF TSCCM.  THE INTERPOLATION
!     IS LINEAR IN LOG P.  WHERE EXTRAPOLATION UPWARD IS NECESSARY,
!     THE T FIELD IS CONSIDERED TO HAVE 0 VERTICAL DERIVATIVE.
!     WHERE EXTRAPOLATION DOWNWARD IS NECESSARY, THE T FIELD IS
!     CONSIDERED TO HAVE A LAPSE RATE OF RLAPSE (K/M), AND THE
!     THICKNESS IS DETERMINED HYDROSTATICALLY FROM THE MEAN OF THE
!     TWO EXTREME TEMPERATUES IN THE LAYER.
!
      do i = 1 , ni
        do j = 1 , nj
          sc = (psrccm(i,j)+ptop)/100.
          do k = 1 , kccm - 1
            if ( sc<=sccm(k+1) .and. sc>=sccm(k) ) k1 = k
          end do
 
!
!         CONDITION FOR SC .GT. SCCM(KCCM) FOLLOWS
!
          if ( sc>sccm(kccm) ) then
            a1 = rgas2*log(sc/sccm(kccm))
            fsccm(i,j) = fccm(i,j,kccm+1-kccm)*(b1-a1)/(b1+a1)
            cycle
          end if
!
!         CONDITION FOR SC .LT. SCCM(KCCM) FOLLOWS
!
          k1p = k1 + 1
          rc = log(sc/sccm(k1))/log(sccm(k1)/sccm(k1p))
          rc1 = rc + 1.
          fsccm(i,j) = rc1*fccm(i,j,kccm+1-k1) - rc*fccm(i,j,kccm+1-k1p)
!
        end do
      end do
 
      end subroutine intv3
